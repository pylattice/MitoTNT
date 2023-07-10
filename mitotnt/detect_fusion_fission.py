import numpy as np
import pandas as pd
import igraph as ig
from tqdm.notebook import trange

def list_to_str(list1):
    string = ""
    if len(list1) == 0:
        return string

    for i in range(len(list1)-1):
        string += str(list1[i])
        string += " "
    string += str(list1[-1])

    return string

# return filtered and paired two lists of same size
def remove_untracked_entries(frag_list1, frag_list2):
    new_list1, new_list2 = [], []
    for i in range(np.min([len(frag_list1), len(frag_list2)])):
        if np.isnan(frag_list1[i]) or np.isnan(frag_list2[i]):
            continue
        else:
            new_list1.append(frag_list1[i])
            new_list2.append(frag_list2[i])
    return new_list1, new_list2


def overlap_score(list1, list2, min_frames):
    if len(list1) != len(list2):
        raise Exception('Two lists should be the same size!')

    num_frames = len(list1)
    if num_frames >= min_frames:
        score = np.sum(np.array(list1)==np.array(list2))
        return 1 - score / num_frames

    else: # less than min_frames
        return np.nan

def find_site_nodes(graph, event_edges, max_edge_gap=5):
    all_site_nodes = []
    merged_edges = []

    for event_edge in event_edges:
        if event_edge not in merged_edges:
            a, b = graph.es[event_edge].source, graph.es[event_edge].target
            neighs_a = graph.neighborhood(vertices=a, order=max_edge_gap)
            neighs_b = graph.neighborhood(vertices=b, order=max_edge_gap)
            all_neighs = np.unique(neighs_a + neighs_b)

            # look at local site and collect event edges
            # local graph need to use attribute for consistent indexing
            local_graph = graph.induced_subgraph(all_neighs)
            site_edges = [e for e in local_graph.es if e['index'] in event_edges]

            # get node id from a list of edges
            site_nodes = []
            for e in site_edges:
                site_nodes += [local_graph.vs[e.source]['frame_node_id'], local_graph.vs[e.target]['frame_node_id']]
            all_site_nodes.append(np.unique(site_nodes).tolist())

            # store event edges that are already merged
            merged_edges += [e for e in local_graph.es['index'] if e in event_edges]

    return all_site_nodes

def detect(input_dir, output_dir, analy_remodeling_dir,
           start_frame, end_frame,
           stride_size=1, half_win_size=4, min_tracked_frames=2):

    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']
    total_number_of_frames = len(full_graph_all_frames)

    if start_frame < half_win_size:
        raise Exception('start_frame must be >= half_win_size but start_frame is '+str(start_frame)+' and half_win_size is '+str(half_win_size))
    if end_frame > total_number_of_frames - half_win_size:
        raise Exception('end_frame must be <= total_number_of_frames - half_win_size but end_frame is '+str(end_frame)+' and total_number_of_frames - half_win_size is '+str(total_number_of_frames-half_win_size))

    tracks = pd.read_csv(output_dir+'final_node_tracks.csv')
    tracks.loc[:,'frame_id'] = tracks.loc[:,'frame_id'].astype(int)
    tracks.loc[:,'frame_frag_id'] = tracks.loc[:,'frame_frag_id'].astype(int)

    all_fragments = []
    for frame in range(0, total_number_of_frames):
        full_graph = full_graph_all_frames[frame]
        frags = full_graph.components()
        all_fragments.append(frags)

    event_list = []
    node_list = []

    print('Performing fusion and fission event detection ...')
    for current_frame in trange(start_frame, end_frame, stride_size):

        frame_tracks = tracks[tracks['frame_id']==current_frame]
        frame_node_id = frame_tracks['frame_node_id'].tolist()
        unique_nodes = frame_tracks['unique_node_id'].tolist()
        frame_nodes = frame_tracks['frame_node_id'].tolist()
        frame_to_unique = {frame_nodes[i]:unique_nodes[i] for i in range(len(frame_nodes))}
        unique_to_frame = {unique_nodes[i]:frame_nodes[i] for i in range(len(unique_nodes))}

        current_fragments = all_fragments[current_frame]
        for frag_id in range(len(current_fragments)):
            all_frag_frame_nodes = current_fragments[frag_id]

            # frame node id that are both tracked and belong to this fragment
            frag_frame_nodes = [n for n in frame_node_id if n in all_frag_frame_nodes]

            # find the unique id for tracked nodes on this fragment
            frag_unique_nodes = [frame_to_unique[f] for f in frag_frame_nodes]

            # for fission events (fusion is reverse)
            forward_track_node_list, backward_track_node_list = [], []
            forward_track_frag_list, backward_track_frag_list = [], []
            for track_id in frag_unique_nodes:
                node_track = tracks[tracks['unique_node_id']==track_id]
                track_frames = node_track['frame_id'].to_numpy()

                # get the frag list for the next half_win_size frames
                track_node_list, track_frag_list = [], []
                for frame in range(current_frame+1, current_frame+half_win_size+1, 1):
                    if frame in track_frames:
                        node_frame_track = node_track[track_frames==frame]
                        node_id = node_frame_track['frame_node_id'].tolist()[0]
                        track_node_list.append(node_id)
                        frame_frag_id = node_frame_track['frame_frag_id'].tolist()[0]
                        track_frag_list.append(frame_frag_id)
                    else:
                        track_node_list.append(np.nan)
                        track_frag_list.append(np.nan)
                forward_track_node_list.append(track_node_list)
                forward_track_frag_list.append(track_frag_list)

                # get the frag list for the past half_win_size frames
                track_node_list, track_frag_list = [], []
                for frame in range(current_frame-1, current_frame-half_win_size-1, -1):
                    if frame in track_frames:
                        node_frame_track = node_track[track_frames==frame]
                        node_id = node_frame_track['frame_node_id'].tolist()[0]
                        track_node_list.append(node_id)
                        frame_frag_id = node_frame_track['frame_frag_id'].tolist()[0]
                        track_frag_list.append(frame_frag_id)
                    else:
                        track_node_list.append(np.nan)
                        track_frag_list.append(np.nan)
                backward_track_node_list.append(track_node_list)
                backward_track_frag_list.append(track_frag_list)

            # make graph
            frag_track_graph = ig.Graph()
            frag_track_graph.add_vertices(len(frag_unique_nodes))
            frag_track_graph.vs['frame_node_id'] = frag_frame_nodes
            frag_track_graph.vs['unique_node_id'] = frag_unique_nodes
            unique_to_index = {frag_unique_nodes[i]:i for i in range(len(frag_unique_nodes))}

            for node in frag_unique_nodes:
                this_track = frame_tracks[frame_tracks['unique_node_id']==node]
                neighs_str = this_track['connected_unique_node_id'].tolist()[0]

                neighs = []
                if pd.isna(neighs_str):
                    continue
                if len(neighs_str) > 0:
                    temp = neighs_str.split(' ')
                    neighs = [int(n) for n in temp]

                for neigh in neighs:
                    index_node, index_neigh = unique_to_index[node], unique_to_index[neigh]
                    frag_track_graph.add_edge(index_node, index_neigh)

                    this_edge = frag_track_graph.es[-1]
                    this_edge['forward_node_list'] = (forward_track_node_list[index_node], forward_track_node_list[index_neigh])
                    this_edge['forward_frag_list'] = (forward_track_frag_list[index_node], forward_track_frag_list[index_neigh])
                    frag_list1, frag_list2 = remove_untracked_entries(forward_track_frag_list[index_node], forward_track_frag_list[index_neigh])
                    this_edge['forward_jaccard'] = overlap_score(frag_list1, frag_list2, min_frames=min_tracked_frames)

                    this_edge = frag_track_graph.es[-1]
                    this_edge['backward_node_list'] = (backward_track_node_list[index_node], backward_track_node_list[index_neigh])
                    this_edge['backward_frag_list'] = (backward_track_frag_list[index_node], backward_track_frag_list[index_neigh])
                    frag_list1, frag_list2 = remove_untracked_entries(backward_track_frag_list[index_node], backward_track_frag_list[index_neigh])
                    this_edge['backward_jaccard'] = overlap_score(frag_list1, frag_list2, min_frames=min_tracked_frames)

            # since edges are added twice
            frag_track_graph.simplify(combine_edges='max')
            for i, e in enumerate(frag_track_graph.es):
                e['index'] = i

            # now extract unique edges
            fission_edges, fusion_edges = [], []
            for e in frag_track_graph.es:
                if e['forward_jaccard'] == 1.0 and e['backward_jaccard'] == 0.0:
                    fission_edges.append(e['index'])
                if e['forward_jaccard'] == 0.0 and e['backward_jaccard'] == 1.0:
                    fusion_edges.append(e['index'])

            # cluster nodes found on valid edges that are close by, as frame node ids
            fission_clusters = find_site_nodes(frag_track_graph, fission_edges)
            fusion_clusters = find_site_nodes(frag_track_graph, fusion_edges)

            # for fission, also find node_id and frag_id for frame after
            fission_next_frame_id, fission_clusters_next_frame, fission_frags_next_frame = [], [], []
            for cluster in fission_clusters:
                cluster_next_frame_id, cluster_next_frame_node_id, cluster_next_frame_frag_id = [], [], []
                for node in cluster:
                    unique_id = frame_to_unique[node]
                    node_track = tracks[tracks['unique_node_id']==unique_id]
                    track_frames = node_track['frame_id'].to_numpy()

                    next_frame = track_frames[int(np.argwhere(track_frames==current_frame)) + 1]
                    next_track = node_track[track_frames==next_frame]
                    next_frame_node_id = next_track['frame_node_id'].tolist()[0]
                    next_frame_frag_id = next_track['frame_frag_id'].tolist()[0]

                    cluster_next_frame_id.append(next_frame)
                    cluster_next_frame_node_id.append(next_frame_node_id)
                    cluster_next_frame_frag_id.append(next_frame_frag_id)

                fission_next_frame_id.append(cluster_next_frame_id)
                fission_clusters_next_frame.append(cluster_next_frame_node_id)
                fission_frags_next_frame.append(cluster_next_frame_frag_id)

            # for fusion, also find node_id and frag_id for frame before
            fusion_last_frame_id, fusion_clusters_last_frame, fusion_frags_last_frame = [], [], []
            for cluster in fusion_clusters:
                cluster_last_frame_id, cluster_last_frame_node_id, cluster_last_frame_frag_id = [], [], []
                for node in cluster:
                    unique_id = frame_to_unique[node]
                    node_track = tracks[tracks['unique_node_id']==unique_id]
                    track_frames = node_track['frame_id'].to_numpy()

                    last_frame = track_frames[int(np.argwhere(track_frames==current_frame)) - 1]
                    last_track = node_track[track_frames==last_frame]
                    last_frame_node_id = last_track['frame_node_id'].tolist()[0]
                    last_frame_frag_id = last_track['frame_frag_id'].tolist()[0]

                    cluster_last_frame_id.append(last_frame)
                    cluster_last_frame_node_id.append(last_frame_node_id)
                    cluster_last_frame_frag_id.append(last_frame_frag_id)

                fusion_last_frame_id.append(cluster_last_frame_id)
                fusion_clusters_last_frame.append(cluster_last_frame_node_id)
                fusion_frags_last_frame.append(cluster_last_frame_frag_id)

            # FISSION
            for i in range(len(fission_clusters)):
                event_list.append({'type':'fission',
                                   'frame_id': current_frame,
                                   'frame_id_before':list_to_str([current_frame]*len(fission_next_frame_id[i])),
                                   'frame_id_after':list_to_str(fission_next_frame_id[i]),
                                   'node_id_before':list_to_str(fission_clusters[i]),
                                   'node_id_after':list_to_str(fission_clusters_next_frame[i]),
                                   'frag_id_before':list_to_str([frag_id]*len(fission_frags_next_frame[i])),
                                   'frag_id_after':list_to_str(fission_frags_next_frame[i]),
                                   'unique_node_id':list_to_str([frame_to_unique[n] for n in fission_clusters[i]])})
            # FUSION
            for i in range(len(fusion_clusters)):
                event_list.append({'type':'fusion',
                                   'frame_id': current_frame-1,
                                   'frame_id_before':list_to_str(fusion_last_frame_id[i]),
                                   'frame_id_after':list_to_str([current_frame]*len(fusion_last_frame_id[i])),
                                   'node_id_before':list_to_str(fusion_clusters_last_frame[i]),
                                   'node_id_after':list_to_str(fusion_clusters[i]),
                                   'frag_id_before':list_to_str(fusion_frags_last_frame[i]),
                                   'frag_id_after':list_to_str([frag_id]*len(fusion_frags_last_frame[i])),
                                   'unique_node_id':list_to_str([frame_to_unique[n] for n in fusion_clusters[i]])})

            # Add reaction nodes
            for i in range(len(fission_clusters)):
                for j in range(len(fission_clusters[i])):
                    node_list.append({'type':'fission', 'frame_id':current_frame,
                                      'frame_node_id':fission_clusters[i][j],
                                      'unique_node_id':frame_to_unique[fission_clusters[i][j]],
                                      'frag_id':frag_id})

            for i in range(len(fusion_clusters)):
                for j in range(len(fusion_clusters_last_frame[i])):
                    node_list.append({'type':'fusion', 'frame_id':fusion_last_frame_id[i][j],
                                      'frame_node_id':fusion_clusters_last_frame[i][j],
                                      'unique_node_id':frame_to_unique[fusion_clusters[i][j]],
                                      'frag_id':fusion_frags_last_frame[i][j]})

    event_list = pd.DataFrame.from_dict(event_list)
    event_list.sort_values('frame_id', inplace=True)
    event_list.to_csv(analy_remodeling_dir+'remodeling_events.csv', index=False)

#     node_list = pd.DataFrame.from_dict(node_list)
#     node_list.sort_values('frame_id', inplace=True)
#     node_list.to_csv(analy_remodeling_dir+'remodeling_nodes.csv', index=False)

    print('The data is saved at', analy_remodeling_dir+'remodeling_events.csv')
#     print('The data is saved at', analy_remodeling_dir+'remodeling_nodes.csv')