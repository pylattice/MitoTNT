import time, os
import numpy as np
import pandas as pd
import igraph as ig
from tqdm.notebook import trange
from scipy.optimize import linear_sum_assignment as lap_solver
from fastdist import fastdist
import warnings
warnings.filterwarnings('ignore')

def dissimilarity_score(sel_edge_len_m, sel_edge_len_n):

    num_edge_diff = len(sel_edge_len_n) - len(sel_edge_len_m)
    if num_edge_diff > 0:
        sel_edge_len_m += [0] * abs(num_edge_diff)
    else:
        sel_edge_len_n += [0] * abs(num_edge_diff)

    num_edges = max(len(sel_edge_len_m), len(sel_edge_len_n))
    cost_mat = np.zeros((num_edges, num_edges))

    # fill cost matrix
    for i in range(num_edges):
        for j in range(num_edges):
            cost_mat[i,j] = abs(sel_edge_len_m[i]-sel_edge_len_n[j]) / max(sel_edge_len_m[i], sel_edge_len_n[j]) # should never have two zeros

    assigned = lap_solver(cost_mat)[1]
    min_sum_cost = 0
    for i in range(num_edges):
        min_sum_cost += cost_mat[i,assigned[i]]

    return min_sum_cost


def local_graph_comparison_score(depth, node_i, node_j, contracted_graphs_m, contracted_graphs_n):

    frag_m, frag_n = contracted_graphs_m[node_i].copy(), contracted_graphs_n[node_j].copy()
    root_m, root_n = frag_m.vs['index'].index(node_i), frag_n.vs['index'].index(node_j)

    node_mapping = {root_m:root_n}

    # iterate each level
    for n_level in range(depth):

        visited_nodes_m = node_mapping.keys()
        visited_nodes_n = node_mapping.values()

        if n_level > 0:
            last_level_m = frag_m.neighborhood(vertices=root_m, order=n_level-1, mindist=n_level-1)
            last_level_n = frag_n.neighborhood(vertices=root_n, order=n_level-1, mindist=n_level-1)
        else:
            last_level_m = []
            last_level_n = []

        # abort once the graph is fully mapped
        current_level_m = frag_m.neighborhood(vertices=root_m, order=n_level, mindist=n_level)
        if len(current_level_m) == 0:
            break

        for node_m in current_level_m:
            node_n = node_mapping[node_m] # use mapping determined from last level

            neighbors_m = frag_m.neighbors(node_m)
            neighbors_n = frag_n.neighbors(node_n)

            # replace each cycle edge with two pseudo-edges of same lengths and add two pseudo-nodes
            for neigh in neighbors_m:
                if neigh in visited_nodes_m and neigh not in last_level_m:
                    dist = frag_m.es[frag_m.get_eid(node_m, neigh)]['distance']
                    frag_m.delete_edges(frag_m.get_eid(node_m, neigh))
                    frag_m.add_vertices(2)
                    frag_m.add_edges([[node_m, frag_m.vs[-2].index]], {'distance':dist})
                    frag_m.add_edges([[neigh, frag_m.vs[-1].index]], {'distance':dist})

            for neigh in neighbors_n:
                if neigh in visited_nodes_n and neigh not in last_level_n:
                    dist = frag_n.es[frag_n.get_eid(node_n, neigh)]['distance']
                    frag_n.delete_edges(frag_n.get_eid(node_n, neigh))
                    frag_n.add_vertices(2)
                    frag_n.add_edges([[node_n, frag_n.vs[-2].index]], {'distance':dist})
                    frag_n.add_edges([[neigh, frag_n.vs[-1].index]], {'distance':dist})


            # add pseudo-nodes and pseudo-edges of 0 at this level
            num_node_diff = len(neighbors_m) - len(neighbors_n)
            if num_node_diff >= 0:
                for it in range(num_node_diff):
                    frag_n.add_vertices(1)
                    frag_n.add_edges([[node_n, frag_n.vs[-1].index]], {'distance':0})
            else:
                for it in range(-1*num_node_diff):
                    frag_m.add_vertices(1)
                    frag_m.add_edges([[node_m, frag_m.vs[-1].index]], {'distance':0})


            # update neighbor list to include pseudo-nodes
            neighbors_m = frag_m.neighbors(node_m)
            neighbors_n = frag_n.neighbors(node_n)
            # remember to exclude parents
            for neigh in neighbors_m:
                if neigh in last_level_m:
                    neighbors_m.remove(neigh)
            for neigh in neighbors_n:
                if neigh in last_level_n:
                    neighbors_n.remove(neigh)


            # map index of cost matrix to real node ids
            index_mapping_m = {i:neighbors_m[i] for i in range(len(neighbors_m))}
            index_mapping_n = {i:neighbors_n[i] for i in range(len(neighbors_n))}

            num_node = max(len(neighbors_m), len(neighbors_n))
            cost_mat = np.zeros((num_node, num_node))

            # fill cost matrix with dissimilarity scores of each node m and n
            for i in range(num_node):
                for j in range(num_node):
                    sel_edge_len_m = frag_m.es[frag_m.incident(index_mapping_m[i])]['distance']
                    sel_edge_len_n = frag_n.es[frag_n.incident(index_mapping_n[j])]['distance']
                    cost_mat[i,j] = dissimilarity_score(sel_edge_len_m, sel_edge_len_n)

            nc = lap_solver(cost_mat)[1]

            # get node correspondence between local graphs m and n
            for a in range(len(nc)):
                node_mapping[index_mapping_m[a]] = index_mapping_n[nc[a]]

    # CALCULATE ADJACENCY MATRIX
    # n-to-m mapping
    reverse_node_mapping = {}
    for a in node_mapping.keys():
        b = node_mapping[a]
        reverse_node_mapping[b] = a

    # fill adjacency matrix with distances between nodes with ordering given by node correspondence
    total_nodes = len(node_mapping.keys())
    weighted_adj_mat_m = np.zeros([total_nodes, total_nodes])
    weighted_adj_mat_n = np.zeros([total_nodes, total_nodes])

    index_mapping_m = {sorted(node_mapping.keys())[i]:i for i in range(len(node_mapping.keys()))}
    index_mapping_n = {sorted(node_mapping.values())[i]:i for i in range(len(node_mapping.values()))}

    for i in node_mapping.keys():
        for j in frag_m.neighbors(i):
            if j in node_mapping.keys():
                weighted_adj_mat_m[index_mapping_m[i],index_mapping_m[j]] = \
                frag_m.es[frag_m.get_eid(i,j)]['distance']

    for i in node_mapping.keys():
        mapped_i = node_mapping[i]
        for mapped_j in frag_n.neighbors(mapped_i):
            if mapped_j in node_mapping.values():
                j = reverse_node_mapping[mapped_j]
                weighted_adj_mat_n[index_mapping_m[i],index_mapping_m[j]] = \
                frag_n.es[frag_n.get_eid(mapped_i,mapped_j)]['distance']

    # calculate euclidean distance between the two weighted adjacency matrices
    weight_diff = weighted_adj_mat_m - weighted_adj_mat_n

    score = np.sum((weight_diff.ravel())**2)

    return score

def list_to_str(list1):

    string = ""
    if len(list1) == 0:
        return string

    for i in range(len(list1)-1):
        string += str(list1[i])
        string += " "
    string += str(list1[-1])

    return string

# =============================================================================
# FRAME-TO-FRAME TRACKING
# =============================================================================

def frametoframe_tracking(input_dir, output_dir, start_frame, end_frame, frame_interval, tracking_interval=1,
                          cutoff_num_neighbor=10, cutoff_speed=None,
                          graph_matching_depth=2, dist_exponent=1, top_exponent=1):

    print('Data loading ...\n')
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)

    # store the data for all frames for easy access
    full_graph_all_frames = inputs['full_graphs']
    classic_graph_per_node_all_frames = inputs['classic_graphs_per_node']
    segment_node_all_frames = inputs['segment_nodes']

    linked_nodes, terminated_nodes, initiated_nodes = [], [], []
    terminated_tracks, ongoing_tracks = [], []

    # declare useful data holders
    previous_costs = []
    for frame in trange(start_frame, end_frame, tracking_interval, desc='Tracking in progress'):

        start = time.time()
        print('\nStart tracking frame {} and {} ...'.format(frame, frame+tracking_interval))

        ### Load data ###

        # load full graph
        full_graph_m = full_graph_all_frames[frame]
        full_graph_n = full_graph_all_frames[frame+tracking_interval]
        cc_m, cc_n = full_graph_m.components(), full_graph_n.components()

        # get number of nodes and coordinates
        number_m, number_n = len(full_graph_m.vs), len(full_graph_n.vs)
        coords_m, coords_n = full_graph_m.vs['coordinate'], full_graph_n.vs['coordinate']

        # get properties
        intensity_m, intensity_n = full_graph_m.vs['intensity'], full_graph_n.vs['intensity']
        width_m, width_n = full_graph_m.vs['width'], full_graph_n.vs['width']

        # load contracted graphs
        classic_graphs_m, classic_graphs_n = classic_graph_per_node_all_frames[frame], classic_graph_per_node_all_frames[frame+tracking_interval]

        # load nodes for each segment
        segment_nodes_m, segment_nodes_n = segment_node_all_frames[frame], segment_node_all_frames[frame+tracking_interval]

        # store branching nodes and ignore them for segments
        branching_nodes_m, branching_nodes_n = [], []
        for i in range(number_m):
            if full_graph_m.vs[i].degree() > 2:
                branching_nodes_m.append(i)
        for j in range(number_n):
            if full_graph_n.vs[j].degree() > 2:
                branching_nodes_n.append(j)

        # know which node belongs to which segment
        node_to_segment_m = {}
        for segment_id, segment in enumerate(segment_nodes_m): # segment consists of of segment nodes
            for b in segment:
                if b in branching_nodes_m:
                    node_to_segment_m[b] = np.nan
                else:
                    node_to_segment_m[b] = segment_id
        node_to_segment_n = {}
        for segment_id, segment in enumerate(segment_nodes_n): # segment consists of of segment nodes
            for b in segment:
                if b in branching_nodes_n:
                    node_to_segment_n[b] = np.nan
                else:
                    node_to_segment_n[b] = segment_id
        ### Finish data loading ###


        ### Calculate distance cost matrix ###
        cost_start = time.time()

        coords_m_mat = np.array(coords_m); coords_n_mat = np.array(coords_n)
        dist_cost_mat = fastdist.matrix_to_matrix_distance(coords_m_mat, coords_n_mat, fastdist.euclidean, "euclidean")

        min_dists = []

        for i in range(number_m):

            row = dist_cost_mat[i,:]

            neighbor_cutoff = sorted(row)[cutoff_num_neighbor]

            row[row > neighbor_cutoff] = np.nan
            dist_cost_mat[i,:] = row

            min_dists.append(np.nanmin(row))

        if cutoff_speed is None:
            speed_cutoff = np.mean(min_dists) + 3 * np.std(min_dists)
        else:
            speed_cutoff = cutoff_speed * frame_interval

        dist_cost_mat[dist_cost_mat > speed_cutoff] = np.nan

        valid_node_pairs = np.argwhere(~np.isnan(dist_cost_mat))
        cost_end = time.time()
        print('Distance cost matrix takes {:.2f} s'.format(cost_end - cost_start))
        ### Distance matrix complete ###

        ### Calculate topology cost matrix ###
        cost_start = time.time()

        topology_cost_mat = np.empty([number_m, number_n])
        topology_cost_mat[:] = np.nan

        for i, j in valid_node_pairs:
            topology_cost_mat[i,j] = local_graph_comparison_score(graph_matching_depth, i, j, classic_graphs_m, classic_graphs_n)

        cost_end = time.time()
        print('Topology cost matrix takes {:.2f} s'.format(cost_end - cost_start))
        ### Topology matrix complete ###


        ### Build final cost matrix ###
        # add pseudo-counts for zero scores
        dist_cost_mat += 0.01 * np.nanmax(dist_cost_mat.ravel())
        topology_cost_mat += 0.01 * np.nanmax(topology_cost_mat.ravel())

        # construct linking cost matrix based on three matrics and relative scaling
        cost_m_n = dist_cost_mat**dist_exponent * topology_cost_mat**top_exponent

        # construct termination cost matrix
        cost_m_m = np.empty([number_m, number_m])
        cost_m_m[:] = np.nan
        for i in range(number_m):
            row = cost_m_n[i,:]
            if np.isnan(row).all():
                cost_m_m[i,i] = 0 # must be assigned to itself since all other nodes exceed max radius
            else:
                if len(previous_costs) == 0:
                    min_cost = np.nanmin(row)
                    cost_m_m[i,i] = 2 * min_cost
                else:
                    cost_m_m[i,i] = np.percentile(previous_costs, 98)

        # construct initiation cost matrix
        cost_n_n = np.empty([number_n, number_n])
        cost_n_n[:] = np.nan
        for j in range(number_n):
            column = cost_m_n[:,j]
            if np.isnan(column).all():
                cost_n_n[j,j] = 0 # must be assigned to itself since all other nodes exceed max radius
            else:
                if len(previous_costs) == 0:
                    min_cost = np.nanmin(column)
                    cost_n_n[j,j] = 2 * min_cost
                else:
                    cost_n_n[j,j] = np.percentile(previous_costs, 98)

        # construct auxiliary block
        cost_n_m = cost_m_n.T.copy()
        cost_n_m[:] = np.nanmin(cost_m_n) # this matrix is needed for LAP solver but not used for tracking

        # assemble into one matrix
        left_block = np.concatenate((cost_m_n, cost_n_n), axis=0)
        right_block = np.concatenate((cost_m_m, cost_n_m), axis=0)
        cost_matrix = np.concatenate((left_block, right_block), axis=1)
        cost_matrix[np.isnan(cost_matrix)] = np.inf # blocking values
        ### Final matrix is complete ###

        ### Solve LAP ###
        assignment = lap_solver(cost_matrix)[1]

        ### Filter out unlikely links ###
        # find linked nodes at frame m, n
        assigned_m = assignment[:number_m]

        linked_m, linked_n = [], []
        for i in range(len(assigned_m)):
            if assigned_m[i] < number_n:
                linked_m.append(i) # first being index for frame t and second for frame t+tracking_interval
                linked_n.append(assigned_m[i])

        filtered_nodes = []
        for full_segment in segment_nodes_m: # segment consists of of segment nodes

            # find only linked nodes in the segment
            segment_nodes = np.array([node for node in full_segment if node in linked_m and node not in branching_nodes_m])

            if len(segment_nodes) > 0:
                # find the assigned nodes at frame n
                assigned_nodes = np.array([assignment[n] for n in segment_nodes])

                # find the seg id of the assigned nodes
                assigned_segment_ids = []
                for n in assigned_nodes:
                    if n in node_to_segment_n.keys():
                        assigned_segment_ids.append(node_to_segment_n[n])
                    else: # if the node doesn't belong to any segment (lone node)
                        assigned_segment_ids.append(np.nan)
                assigned_segment_ids = np.array(assigned_segment_ids)

                # seg id to node mapping
                seg_to_nodes = {}
                for index, segment_id in enumerate(assigned_segment_ids):
                    if segment_id not in seg_to_nodes.keys():
                        seg_to_nodes[segment_id] = [segment_nodes[index]]
                    else:
                        seg_to_nodes[segment_id] += [segment_nodes[index]]

                # node count of each segment
                counts = {}
                for segment_id in assigned_segment_ids:
                    if segment_id not in counts.keys():
                        counts[segment_id] = 1
                    else:
                        counts[segment_id] += 1

                max_segment_id = max(counts, key=counts.get)

                if counts[max_segment_id] == 1:
                    filtered_nodes += list(segment_nodes)

                else:
                    # filter by displacement vector norm
                    majority_nodes = segment_nodes[assigned_segment_ids == max_segment_id]
                    mean_majority_dist = np.mean([dist_cost_mat[i,j] for i,j in list(zip(majority_nodes, [assignment[n] for n in majority_nodes]))])
                    other_nodes = [n for n in segment_nodes if n not in majority_nodes]
                    for node in other_nodes:
                        if dist_cost_mat[node, assignment[node]] > 3 * mean_majority_dist: # cutoff is here
                            filtered_nodes.append(node)

                    # remove outliers pointing to other segments without partners
                    unremoved_nodes = [n for n in segment_nodes if n not in filtered_nodes]
                    for index, node in enumerate(unremoved_nodes):
                        seg_id = assigned_segment_ids[index]
                        if np.isnan(seg_id):
                            continue
                        elif seg_id != max_segment_id:
                            neighs = full_graph_m.neighbors(node)
                            count = np.sum([neigh in seg_to_nodes[seg_id] for neigh in neighs])
                            if count == 0: # if no neigh pointing to the same segment
                                filtered_nodes.append(node)

        # print('Number of filtered nodes is', len(filtered_nodes))

        # update assignment after filtering
        assignment_filtered = assignment.copy()

        for index_m in filtered_nodes:
            index_n = linked_n[linked_m.index(index_m)]
            # set the linked node to initated
            assignment_filtered[number_m + index_n] = index_n
            # set the node to terminated
            assignment_filtered[index_m] = number_n + index_m
        ### Unlikely links are filtered ###

        ### Report assignments ###
        # interpretate assignment
        assigned_m, assigned_n = assignment_filtered[:number_m], assignment_filtered[number_m:]

        linked, terminated, initiated = [], [], []
        for i in range(len(assigned_m)):
            if assigned_m[i] < number_n:
                linked.append([i, assigned_m[i]]) # first being index for frame t and second for frame t+tracking_interval
            else:
                terminated.append(i)

        for i in range(len(assigned_n)):
            if assigned_n[i] < number_n:
                initiated.append(i)

        linked = np.array(linked); terminated = np.array(terminated); initiated = np.array(initiated)
        linked_nodes.append(linked); terminated_nodes.append(terminated); initiated_nodes.append(initiated)

        assignment_dist = [dist_cost_mat[a,b] for (a,b) in linked]
        previous_costs += list(cost_matrix[np.arange(number_m), assignment[:number_m]])
        avg_cost = cost_matrix[np.arange(number_m),assignment[:number_m]].sum() / number_m
        end = time.time()

        # output stats
        print('Number of nodes at frame {}, {} are {}, {}'.format(frame, frame+tracking_interval, number_m, number_n))
        print('Number of nodes linked, terminated at frame {}: {}, {}'.format(frame, len(linked), len(terminated)))
        print('Number of nodes initiated at frame {}: {}'.format(frame+tracking_interval, len(initiated)))
        print('Mean, std of distances for tracked nodes: {:2f}, {:2f}'.format(np.mean(assignment_dist), np.std(assignment_dist)))
        print('Average cost: {:.1f}'.format(avg_cost))
        print('Tracking is complete and takes {:.2f} s\n'.format(end - start))
        ### Assignments reported ###

        ### Update tracks ###
        nodes_m, nodes_n = linked[:,0].tolist(), linked[:,1].tolist()
        tracks_to_remove = []

        if frame == start_frame:
            for i in range(len(nodes_m)):
                # initiate with first two frames
                ongoing_tracks.append([[frame, frame+tracking_interval],
                                       [nodes_m[i], nodes_n[i]],
                                       [node_to_segment_m[nodes_m[i]], node_to_segment_n[nodes_n[i]]],
                                       [cc_m.membership[nodes_m[i]], cc_n.membership[nodes_n[i]]],
                                       [coords_m[nodes_m[i]], coords_n[nodes_n[i]]],
                                       [intensity_m[nodes_m[i]], intensity_n[nodes_n[i]]],
                                       [width_m[nodes_m[i]], width_n[nodes_n[i]]],])

        else:
            for idx, track in enumerate(ongoing_tracks):
                # if terminated, remove track; else append linked node
                if track[1][-1] in terminated:
                    terminated_tracks.append(track)
                    tracks_to_remove.append(idx)
                else:
                    linked_index = nodes_m.index(track[1][-1])
                    linked_node = nodes_n[linked_index]
                    track[0].append(frame+tracking_interval)
                    track[1].append(linked_node)
                    track[2].append(node_to_segment_n[linked_node])
                    track[3].append(cc_n.membership[linked_node])
                    track[4].append(coords_n[linked_node])
                    track[5].append(intensity_n[linked_node])
                    track[6].append(width_n[linked_node])

        # delete terminated tracks from ongoing tracks
        ongoing_tracks = [t for i,t in enumerate(ongoing_tracks) if i not in tracks_to_remove]

        for init_node in initiated:
            ongoing_tracks.append([[frame+tracking_interval],
                                   [init_node],
                                   [node_to_segment_n[init_node]],
                                   [cc_n.membership[init_node]],
                                   [coords_n[init_node]],
                                   [intensity_n[init_node]],
                                   [width_n[init_node]]])

    ### Save frame-to-frame tracking ###
    path = output_dir+'frametoframe_tracking_outputs.npz'
    data = {}
    data['linked_nodes'] = np.array(linked_nodes, dtype=object)
    data['terminated_nodes'] = np.array(terminated_nodes, dtype=object)
    data['initiated_nodes'] = np.array(initiated_nodes, dtype=object)
    np.savez(path, **data)

    # each element in all_tracks is a track with 1) frame numbers; 2) node indices; 3) segment ids of the node, 4) frag ids of the node; 4) node coords; 5) node intensities; 6) node widths
    terminated_tracks = np.array(terminated_tracks, dtype=object)
    ongoing_tracks = np.array(ongoing_tracks, dtype=object)
    all_tracks = np.concatenate([terminated_tracks, ongoing_tracks])
    np.save(output_dir+'all_tracks.npy', np.array(all_tracks, dtype=object))
    ### tracking saved ###


# =============================================================================
# GAP CLOSING
# =============================================================================

def gap_closing(input_dir, output_dir, start_frame, end_frame, tracking_interval=1,
                graph_matching_depth=2, dist_exponent=1, top_exponent=1,
                min_track_size=5, max_gap_size=3, memory_efficient_gap_closing=False):

    # reload all tracks
    all_tracks = np.load(output_dir+'/all_tracks.npy', allow_pickle=True)

    # filter out very short tracks
    very_short_tracks = []
    for i in range(len(all_tracks)):
        if len(all_tracks[i][0]) < min_track_size:
            very_short_tracks.append(i)

    all_tracks = np.delete(all_tracks, very_short_tracks, axis=0)
    all_tracks = sorted(all_tracks, key=lambda track : track[0][0]) # sort by start frame
    num_tracks = len(all_tracks)

    print('Gap closing initiates ...\n')

    # get track disps
    all_track_disps = []
    for t in all_tracks:
        track_coords = t[4] # use index for the coordinates
        all_track_disps.append([np.linalg.norm(track_coords[t+1]-track_coords[t]) for t in range(len(track_coords)-1)])

    ### load data again
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']
    classic_graph_per_node_all_frames = inputs['classic_graphs_per_node']
    ### done loading

    all_track_assignments = {}
    partition_start = 0
    iter_num = 1
    while partition_start < num_tracks:
        num_frame = end_frame-start_frame

        if not memory_efficient_gap_closing:
            partition_size = num_tracks
            overlap_size = 0
        else:
            partition_size = int(num_tracks / num_frame) * 50
            overlap_size =  int(num_tracks / num_frame) * 10

        partition_end = partition_start + partition_size
        if partition_end > num_tracks:
            partition_end = num_tracks
        print('Block '+str(iter_num)+' index:', partition_start, 'to', partition_end)
        print('Computing cost terms for block '+str(iter_num))

        track_cost_m_n = np.empty([partition_size, partition_size])
        track_cost_m_n[:] = np.nan

        for i in trange(partition_start, partition_end):

            track_frames_m, track_nodes_m, track_coords_m = all_tracks[i][0], all_tracks[i][1], all_tracks[i][4]
            end_frame_m, end_node_m, end_coord_m = track_frames_m[-1], track_nodes_m[-1], track_coords_m[-1]

            classic_graphs_m = classic_graph_per_node_all_frames[end_frame_m]
            disps_m = all_track_disps[i]

            for j in range(i+1, partition_end): # no need to check index less than i because the start frame is sorted

                track_frames_n, track_nodes_n, track_coords_n = all_tracks[j][0], all_tracks[j][1], all_tracks[j][4]
                start_frame_n, start_node_n, start_coord_n = track_frames_n[0], track_nodes_n[0], track_coords_n[0]

                gap_size = start_frame_n - end_frame_m - 1

                # check only those within max gap size
                if 1 <= gap_size <= max_gap_size: # if gap == 1 it should have been linked before - so skip it

                    # load data
                    classic_graphs_n = classic_graph_per_node_all_frames[start_frame_n]

                    # compute distance cutoff based on the two tracks
                    disps_n = all_track_disps[j]
                    comb_disps = disps_m + disps_n

                    dist_cutoff = (gap_size + 1) * (3 * np.std(comb_disps)) ** 2

                    # compute node-to-node distance
                    dist = end_coord_m - start_coord_n
                    dist_cost = np.sum(dist**2)

                    # filter by distance cutoff
                    if dist_cost > dist_cutoff:
                        continue

                    # compute topology cost
                    topology_cost = local_graph_comparison_score(graph_matching_depth, end_node_m, start_node_n, classic_graphs_m, classic_graphs_n)

                    # assign g.c. cost
                    track_cost_m_n[i-partition_start,j-partition_start] = dist_cost**dist_exponent * topology_cost**top_exponent

        # construct termination cost matrix
        track_cost_m_m = np.empty([partition_size, partition_size])
        track_cost_m_m[:] = np.nan

        for i in range(partition_size):
            row = track_cost_m_n[i,:]

            if np.isnan(row).all():
                track_cost_m_m[i,i] = 0 # must be assigned to itself since all other nodes exceed max radius
            else:
                min_cost = np.nanmin(row)
                track_cost_m_m[i,i] = 2 * min_cost

        # assemble into one matrix
        track_cost_matrix = np.concatenate((track_cost_m_n, track_cost_m_m), axis=1)
        track_cost_matrix[np.isnan(track_cost_matrix)] = np.inf

        # evaluate memory usage
        print('Cost matrix memory usage: {:.1f} MB\n'.format(track_cost_matrix.nbytes / 1024**2))

        # solve LAP and store linking results
        assignment = lap_solver(track_cost_matrix)[1]

        linked = []
        for i in range(len(assignment)):
            if assignment[i] < partition_size:
                linked.append([i, assignment[i]]) # first being index for frame t and second for frame t+tracking_interval

        for pair in linked:
            all_track_assignments[partition_start + pair[0]] = partition_start + pair[1] # offset by start index of the partition

        # go to next partition
        partition_start = partition_start + partition_size - overlap_size
        iter_num += 1

    # convert dictionary to array and overwrite assignment by next partition's assignment for the overlapped region
    linked_tracks = np.zeros([len(all_track_assignments.keys()), 2], dtype=int)
    for i, a in enumerate(all_track_assignments.keys()):
        linked_tracks[i,0] = a
        linked_tracks[i,1] = all_track_assignments[a]


    # combine tracks for gap closing
    print('Start combining closed tracks ...')
    if linked_tracks.shape[0] > 0:
        linked_tracks_for_update = linked_tracks.copy() # used to record appended tracks
        tracks_of_track = [] # list of linked tracks
        all_linked_tracks = [] # record tracks that are closed

        # recursive function for finding linked tracks
        def find_all_linked_tracks(tracks, track_id):
            for i in range(0, len(linked_tracks_for_update)):
                if linked_tracks_for_update[i,0] == track_id: # find the first track
                    tracks.append(linked_tracks_for_update[i,1])
                    linked_tracks_for_update[i,0] = -1 # note that the track is already appended
                    find_all_linked_tracks(tracks, linked_tracks_for_update[i,1]) # go to the next track and find linked_tracks track

        # for each track find the series of linked tracks
        for index in trange(len(linked_tracks)):
            track_id = linked_tracks[index,0]
            if track_id in linked_tracks_for_update[:,0]:
                tracks = [track_id]
                find_all_linked_tracks(tracks, track_id)
                tracks_of_track.append(tracks)
                all_linked_tracks += tracks

        # concatenate data for closed tracks
        all_closed_tracks = []
        for tot in tracks_of_track:
            all_closed_tracks.append([sum([all_tracks[t][0] for t in tot], []),
                                      sum([all_tracks[t][1] for t in tot], []),
                                      sum([all_tracks[t][2] for t in tot], []),
                                      sum([all_tracks[t][3] for t in tot], []),
                                      sum([all_tracks[t][4] for t in tot], []),
                                      sum([all_tracks[t][5] for t in tot], []),
                                      sum([all_tracks[t][6] for t in tot], [])])

        # add unclosed tracks back
        for t in range(num_tracks):
            if t not in all_linked_tracks:
                all_closed_tracks.append(all_tracks[t].tolist())

        # sort tracks
        sort_by_length = sorted(all_closed_tracks, key=lambda track : len(track[0]), reverse=True) # first sort by size of track
        sort_by_start = sorted(sort_by_length, key=lambda track : track[0][0]) # then sort by start frame
        all_closed_tracks = sort_by_start

    else:
        all_closed_tracks = all_tracks

    np.save(output_dir+'all_closed_tracks.npy', np.array(all_closed_tracks, dtype=object))

    print('Number of tracks and average track length before gap closing: {}, {:.2f}\n'.format(
          len(all_tracks), np.mean([len(track[0]) for track in all_tracks])))

    print('Number of tracks and average track length after gap closing: {}, {:.2f}\n'.format(
          len(all_closed_tracks), np.mean([len(track[0]) for track in all_closed_tracks])))


    # Save tracks in the form of one node per row ###
    print('Saving final node trajectory file ... This might take a while.\n')
    tracks = pd.DataFrame(columns={'frame_id', 'unique_node_id', 'frame_node_id', 'frame_seg_id', 'frame_frag_id', 'x', 'y', 'z','intensity','width'})
    tracks = tracks[['frame_id', 'unique_node_id', 'frame_node_id', 'frame_seg_id', 'frame_frag_id', 'x', 'y', 'z','intensity','width']] # reorder the columns

    df_index = 0
    for track_id, track in enumerate(all_closed_tracks):
        track_frames, track_nodes, track_segs, track_frags, track_coords, track_ints, track_widths = track[0], track[1], track[2], track[3], track[4], track[5], track[6]

        for f in range(len(track_frames)):
            x, y, z = track_coords[f][0], track_coords[f][1], track_coords[f][2]
            tracks.loc[df_index] = {'frame_id': track_frames[f], 'unique_node_id': track_id,
                            'frame_node_id': track_nodes[f], 'frame_seg_id': track_segs[f], 'frame_frag_id': track_frags[f],
                            'x': x, 'y': y, 'z': z, 'intensity': track_ints[f], 'width': track_widths[f]}
            df_index += 1

    tracks.sort_values(['unique_node_id', 'frame_id'], inplace=True, ignore_index=True)

    # add connected nodes using unique indexing #
    new_tracks = pd.DataFrame()
    for frame in range(start_frame, end_frame, tracking_interval):

        full_graph = full_graph_all_frames[frame]
        tracks_frame = tracks[tracks['frame_id'] == frame]

        tracked_frame_nodes = tracks_frame['frame_node_id'].astype('int').tolist()
        unique_nodes = tracks_frame['unique_node_id'].astype('int').tolist()
        frame_to_unique = {tracked_frame_nodes[i]:unique_nodes[i] for i in range(len(tracks_frame))}

        def find_connected_unique_nodes(this_node, visited_nodes):

            neighs = full_graph.neighbors(this_node)
            for visited_node in visited_nodes:
                if visited_node in neighs:
                    neighs.remove(visited_node)

            visited_nodes.append(this_node)

            for neigh in neighs:
                # if the frame node is tracked for this frame, add to list
                if neigh in tracked_frame_nodes:
                    connected_unique_nodes.append(frame_to_unique[neigh])
                else:
                    find_connected_unique_nodes(neigh, visited_nodes)
            return

        all_connected_unique_nodes = []
        for node in tracked_frame_nodes:
            connected_unique_nodes = []
            find_connected_unique_nodes(node, [node])
            all_connected_unique_nodes.append(connected_unique_nodes)

        tracks_frame.insert(2, 'connected_unique_node_id', [list_to_str(a) for a in all_connected_unique_nodes])
        new_tracks = pd.concat([new_tracks, tracks_frame]) # accumulate tracks from each frame

    # reorder the columns
    new_tracks = new_tracks[['frame_id', 'unique_node_id', 'frame_node_id', 'frame_seg_id', 'frame_frag_id', 'connected_unique_node_id', 'x', 'y', 'z', 'intensity', 'width']]
    new_tracks.to_csv(output_dir+'final_node_tracks.csv', index=False)

    print('File saved. Tracking is complete.')