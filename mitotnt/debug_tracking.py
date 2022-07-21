from network_tracking import *
import time, os
import numpy as np
import pandas as pd
import igraph as ig
from scipy.optimize import linear_sum_assignment as lap_solver
from fastdist import fastdist
# D:\Python Scripts\Mito\MitoTNT\mitotnt

#%%
work_dir = 'D:/Python Scripts/Mito/MitoTNT/test_data/'
data_dir = work_dir+'mitograph/'
input_dir = work_dir+'tracking_input/'
if not os.path.isdir(input_dir):
    os.mkdir(input_dir)
output_dir = work_dir+'tracking_output/'
if not os.path.isdir(output_dir):
    os.mkdir(output_dir)

start_frame = 10
# end_frame = len([frame for frame in os.listdir(data_dir) if 'FRAME' in frame])
end_frame = 20

tracking_interval = 1
frame_interval = 3.253

distance_cutoff_mode = 2
cutoff_num_neighbor = 8
cutoff_speed = 3

graph_matching_depth = 2
dist_exponent, top_exponent = 1, 1

#%%
print('Data loading ...')
inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)

# store the data for all frames for easy access
full_graph_all_frames = inputs['full_graphs']
classic_graph_per_node_all_frames = inputs['classic_graphs_per_node']
segment_node_all_frames = inputs['segment_nodes']

#%%
linked_nodes, terminated_nodes, initiated_nodes = [], [], []
terminated_tracks, ongoing_tracks = [], []

# declare useful data holders
previous_costs = []
all_cutoff_dict = {}
for frame in range(start_frame, end_frame-tracking_interval, tracking_interval):

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

    all_cutoff = []
    ongoing_nodes = [t[-1] for t in ongoing_tracks]
    for i in range(number_m):

        dist_cutoff = np.inf
        row = dist_cost_mat[i,:]

        if cutoff_speed is not None:
            dist_cutoff = cutoff_speed * frame_interval # fixed distance cutoff of default 1 um/s * frame interval (s)

        elif cutoff_num_neighbor is not None:
            dist_cutoff = sorted(row)[cutoff_num_neighbor]

        else:
            raise Exception('Either cutoff speed or cutoff number of neighbors needs to be provided')

        row[row > dist_cutoff] = np.nan
        dist_cost_mat[i,:] = row
        all_cutoff.append(dist_cutoff)

    all_cutoff_dict[frame] = all_cutoff
    print('Mean and std of the cutoffs: {}, {}'.format(np.mean(all_cutoff), np.std(all_cutoff)))

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

            # first filter by displacement vector norm
            majority_nodes = segment_nodes[assigned_segment_ids==max_segment_id]
            mean_majority_dist = np.mean([dist_cost_mat[i,j] for i,j in list(zip(majority_nodes, [assignment[n] for n in majority_nodes]))])
            other_nodes = [n for n in segment_nodes if n not in majority_nodes]
            for node in other_nodes:
                if dist_cost_mat[node, assignment[node]] > 3 * mean_majority_dist: # cutoff is here
                    filtered_nodes.append(node)

            # remove outliers pointing to other segments without partners
            segment_nodes = [n for n in segment_nodes if n not in filtered_nodes]
            if len(segment_nodes) == 1:
                filtered_nodes.append(segment_nodes[0])
            else:
                for index, node in enumerate(segment_nodes):
                    seg_id = assigned_segment_ids[index]
                    if np.isnan(seg_id):
                        continue
                    elif seg_id != max_segment_id:
                        neighs = full_graph_m.neighbors(node)
                        count = np.sum([neigh in seg_to_nodes[seg_id] for neigh in neighs])
                        if count == 0: # if no neigh pointing to the same segment
                            filtered_nodes.append(node)
#     print('Number of filtered nodes is', len(filtered_nodes))

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

    # store the assignments
    linked_nodes.append(linked); terminated_nodes.append(terminated); initiated_nodes.append(initiated)

    # output stats
    avg_cost = cost_matrix[np.arange(number_m),assignment[:number_m]].sum() / number_m
    end = time.time()
    print('Number of nodes at frame {}, {} are {}, {}'.format(frame, frame+tracking_interval, number_m, number_n),
          '\nNumber of nodes linked, terminated at frame {}: {}, {}'.format(frame, len(linked), len(terminated)),
          '\nNumber of nodes initiated at frame {}: {}'.format(frame+tracking_interval, len(initiated)),
          '\nAverage cost: {:.0f}'.format(avg_cost),
          '\nTracking is done and takes {:.2f} s'.format(end - start))
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

    # record the assignment costs
    previous_costs += list(cost_matrix[np.arange(number_m),assignment[:number_m]])

### Save frame-to-frame tracking ###
suffix = ''
np.save(output_dir+'linked_nodes_each_frame'+suffix+'.npy', np.array(linked_nodes, dtype=object))
np.save(output_dir+'terminated_nodes_each_frame'+suffix+'.npy', np.array(terminated_nodes, dtype=object))
np.save(output_dir+'initiated_nodes_each_frame'+suffix+'.npy', np.array(initiated_nodes, dtype=object))

# each element in all_tracks is a track with 1) frame numbers; 2) node indices; 3) segment ids of the node, 4) frag ids of the node; 4) node coords; 5) node intensities; 6) node widths
terminated_tracks = np.array(terminated_tracks, dtype=object)
ongoing_tracks = np.array(ongoing_tracks, dtype=object)
all_tracks = np.concatenate([terminated_tracks, ongoing_tracks])
np.save(output_dir+'all_tracks.npy', np.array(all_tracks, dtype=object))
### tracking saved ###