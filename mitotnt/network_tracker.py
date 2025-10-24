import time
import warnings
import numpy as np
import pandas as pd
from tqdm.auto import tqdm, trange
from scipy.optimize import linear_sum_assignment as lap_solver
from fastdist import fastdist
from mitotnt.skeletonized_mito import SkeletonizedMito
from mitotnt.tracked_mito import TrackedMito

class NetworkTracker:
    """
    This class provides functionality to track skeletonized mitochondria
    in time-lapse imaging data. It uses distance and graph-based score to link
    mitochondrial nodes frame-to-frame and reconstruct trajectories
    under configurable tracking constraints.

    Parameters
    ----------
    segmented_mito : SkeletonizedMito
        Pre-processed mitochondrial data with skeletonized structures.
    frame_interval : float, optional
        Time between consecutive frames, in seconds (default is 1.0).
    start_frame : int, optional
        Index of the first frame to include in tracking (default is 0).
    end_frame : int, optional
        Index of the last frame to include in tracking. If None, uses
        the total number of frames from `segmented_mito` (default is None).
    tracking_interval : int, optional
        The frame number difference for the two frames to be tracked
        (default is 1 so track every frame).
    cutoff_num_neighbor : int, optional
        Maximum number of nearest neighbors to consider when making node
        assignments across frames (default is 10).
    cutoff_speed : float, optional
        Maximum allowed speed of mitochondria (microns/frame interval).
        If None, no speed cutoff is applied (default is None).
    graph_matching_depth : int, optional
        Depth of neighborhood considered for graph matching (default is 2).
    dist_exponent : int, optional
        Exponent applied to spatial distance when computing costs
        (default is 1).
    top_exponent : int, optional
        Exponent applied to topological distance when computing costs
        (default is 1).
    min_track_size : int, optional
        Minimum number of frames required for a valid track (default is 4).
    max_gap_size : int, optional
        Maximum number of consecutive missing frames allowed when
        performing gap closing tracks.
    block_size_factor : float, optional
        Values less than 1 allows using sliding blocks of cost matrix during gap closing
        to prevent memory overflow due to large number of tracks.
        (default is 1, close all tracks at the same time).

    Attributes
    ----------
    segmented_mito : SkeletonizedMito
        Reference to the input mitochondrial skeleton data.
    start_frame : int
        First frame index used for tracking.
    end_frame : int
        Last frame index used for tracking.
    frame_interval : float
        Time between frames.
    tracking_interval : int
        Step size in frames when tracking.
    cutoff_num_neighbor : int
        Maxium number of neighbors to consider.
    cutoff_speed : float
        Speed threshold for tracking.
    graph_matching_depth : int
        Graph depth used in matching.
    dist_exponent : int
        Distance exponent in matching cost.
    top_exponent : int
        Topology exponent in matching cost.
    min_track_size : int
        Minimum valid trajectory length.
    max_gap_size : int
        Maximum number of frames for an allowed gap.
    block_size_factor : float
        Factor for using sliding cost blocks during gap closing.
    """

    def __init__(self, segmented_mito: SkeletonizedMito, frame_interval: float = None,
                 start_frame: int = 0, end_frame: int = None, tracking_interval: int = 1,
                 cutoff_num_neighbor: int = 10, cutoff_speed: float = None,
                 graph_matching_depth: int = 2, dist_exponent: int = 1, top_exponent: int = 1,
                 min_track_size: int = 4, max_gap_size: int = 3, block_size_factor: float = 1.0):

        self.segmented_mito = segmented_mito
        self.start_frame = start_frame
        if end_frame is None:
            self.end_frame = self.segmented_mito.num_frames
        else:
            self.end_frame = end_frame
        if frame_interval is None:
            raise Exception()
        self.frame_interval = frame_interval
        self.tracking_interval = tracking_interval
        self.cutoff_num_neighbor = cutoff_num_neighbor
        self.cutoff_speed = cutoff_speed
        self.graph_matching_depth = graph_matching_depth
        self.dist_exponent = dist_exponent
        self.top_exponent = top_exponent
        self.min_track_size = min_track_size
        self.max_gap_size = max_gap_size
        self.block_size_factor = block_size_factor


    def reload_results(self):
        """
        Same as `run()` except reload previously saved results into `TrackedMito`.

        Returns
        ----------
        object
            `TrackedMito` object with network tracking results and metadata.
            `TrackedMito.node_tracks` contains tabular data with each tracked node at one frame
            as rows and the following columns:

            - ``frame_id`` (int): the frame number.
            - ``frame_node_id`` (int): node id with frame-wise indexing.
            - ``unique_node_id`` (int): node id shared by all the nodes in the same track at different frames. Each track is uniquely indexed throughout the whole trajectory. This contains the tracking information.
            - ``frame_seg_id`` (int): segment id for each mitochondrial segment (between non-degree-2 nodes) with frame-wise indexing.
            - ``frame_frag_id`` (int): fragment id for each mitochondrial fragment (connected component) with frame-wise indexing.
            - ``connected_unique_node_id`` (str): space-delimited `unique_node_id` for tracked neigboring nodes in the graph. Note that the topology may be different from static graphs due to absence of untracked nodes.
            - ``x``, ``y``, ``z``: coordinates for the node.
            - ``intensity``, ``width``: pixel intensity and tubular width for the node from MitoGraph.

        """

        try:
            node_tracks = pd.read_csv(self.segmented_mito.save_path + 'mito_node_tracks.csv')
            linked_nodes = np.load(self.segmented_mito.save_path + 'mito_linked_nodes.npy', allow_pickle=True)

        except:
            raise Exception('Could not locate saved results.')

        return TrackedMito(self.segmented_mito, self.frame_interval, self.start_frame, self.end_frame, self.tracking_interval,
                           node_tracks, linked_nodes)


    def run(self):
        """
        Perform frame-to-frame tracking of mitochondria using the parameters declared in `NetworkTracker`.

        Returns
        ----------
        object
            `TrackedMito` object with network tracking results and metadata.
            `TrackedMito.node_tracks` contains tabular data with each tracked node at one frame
            as rows and the following columns:

            - ``frame_id`` (int): the frame number.
            - ``frame_node_id`` (int): node id with frame-wise indexing.
            - ``unique_node_id`` (int): node id shared by all the nodes in the same track at different frames. Each track is uniquely indexed throughout the whole trajectory. This contains the tracking information.
            - ``frame_seg_id`` (int): segment id for each mitochondrial segment (between non-degree-2 nodes) with frame-wise indexing.
            - ``frame_frag_id`` (int): fragment id for each mitochondrial fragment (connected component) with frame-wise indexing.
            - ``connected_unique_node_id`` (str): space-delimited `unique_node_id` for tracked neigboring nodes in the graph. Note that the topology may be different from static graphs due to absence of untracked nodes.
            - ``x``, ``y``, ``z``: coordinates for the node.
            - ``intensity``, ``width``: pixel intensity and tubular width for the node from MitoGraph.

        """

        # store the data for all frames for easy access
        full_graphs_all_frames = self.segmented_mito.full_graphs
        segment_node_all_frames = self.segmented_mito.segment_nodes
        local_simple_graphs_all_frames= self.segmented_mito.local_simple_graphs

        # declare useful data holders
        linked_nodes, terminated_nodes, initiated_nodes = [], [], []
        terminated_tracks, ongoing_tracks = [], []

        num_frames = len(full_graphs_all_frames)
        if self.end_frame > num_frames:
            self.end_frame = num_frames - 1
            warnings.warn("The end frame specified is less than the number of frames in tracking inputs. "
                          "End frame has been changed to the maximum number of frames.")

        for frame in trange(self.start_frame, self.end_frame - self.tracking_interval, self.tracking_interval, desc="Frame-to-frame network tracking"):

            start = time.time()
            print(f"Start tracking frame {frame} and {frame + self.tracking_interval} ...")

            ### Load data ###

            # load full graph
            full_graph_m = full_graphs_all_frames[frame]
            full_graph_n = full_graphs_all_frames[frame + self.tracking_interval]
            cc_m, cc_n = full_graph_m.components(), full_graph_n.components()

            # get number of nodes and coordinates
            number_m, number_n = len(full_graph_m.vs), len(full_graph_n.vs)
            coords_m, coords_n = full_graph_m.vs['coordinate'], full_graph_n.vs['coordinate']

            # warnings for unusual mitograph outputs
            max_number = 5000
            if number_m > max_number or number_n > max_number:
                warnings.warn("The number of nodes is relatively large and may take longer time to process! "
                              "Recommend to crop a smaller region. Alternatively you can increase node_gap_size during generate_tracking_inputs.generate()")

            fluctuation_percent_max = 30
            fluctuation_percent = round(abs(number_n - number_m) / number_m, 3) * 100
            if fluctuation_percent > fluctuation_percent_max:
                warnings.warn(
                    f"The number of nodes changes by {fluctuation_percent}% between the two frames. "
                    "Please check for imaging or segmentation artifacts.")

            # get properties
            intensity_m, intensity_n = full_graph_m.vs['intensity'], full_graph_n.vs['intensity']
            width_m, width_n = full_graph_m.vs['width'], full_graph_n.vs['width']

            # load contracted graphs
            simple_graphs_m, simple_graphs_n = local_simple_graphs_all_frames[frame], \
            local_simple_graphs_all_frames[frame + self.tracking_interval]

            # load nodes for each segment
            all_segment_nodes_m, all_segment_nodes_n = segment_node_all_frames[frame], segment_node_all_frames[
                frame + self.tracking_interval]

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
            for segment_id, segment in enumerate(all_segment_nodes_m):  # segment consists of segment nodes
                for b in segment:
                    if b in branching_nodes_m:
                        node_to_segment_m[b] = np.nan
                    else:
                        node_to_segment_m[b] = segment_id
            node_to_segment_n = {}
            for segment_id, segment in enumerate(all_segment_nodes_n):  # segment consists of segment nodes
                for b in segment:
                    if b in branching_nodes_n:
                        node_to_segment_n[b] = np.nan
                    else:
                        node_to_segment_n[b] = segment_id
            ### Finish data loading ###

            ### Calculate distance cost matrix ###
            cost_start = time.time()

            coords_m_mat = np.array(coords_m)
            coords_n_mat = np.array(coords_n)
            dist_cost_mat = fastdist.matrix_to_matrix_distance(coords_m_mat, coords_n_mat, fastdist.euclidean, 'euclidean')

            min_dists = []

            # neighbor cutoff
            for i in range(number_m):
                row = dist_cost_mat[i, :]
                if len(row) > self.cutoff_num_neighbor:
                    neighbor_cutoff = sorted(row)[self.cutoff_num_neighbor]
                    row[row > neighbor_cutoff] = np.nan
                    dist_cost_mat[i, :] = row

                min_dists.append(np.nanmin(row))

            # displacement cutoff
            if self.cutoff_speed is None:
                valid_min_dists = [d for d in min_dists if not np.isnan(d)]
                if len(valid_min_dists) == 0:
                    disp_cutoff = np.inf  # no valid distances found
                else:
                    disp_cutoff = np.mean(min_dists) + 3 * np.std(min_dists)  # global estimate based on all nodes
            else:
                disp_cutoff = self.cutoff_speed * self.frame_interval

            dist_cost_mat[dist_cost_mat > disp_cutoff] = np.nan

            valid_node_pairs = np.argwhere(~np.isnan(dist_cost_mat))
            cost_end = time.time()
            # print('Distance cost matrix takes {:.2f} s'.format(cost_end - cost_start))
            ### Distance matrix complete ###

            ### Calculate topology cost matrix ###
            cost_start = time.time()

            topology_cost_mat = np.empty([number_m, number_n])
            topology_cost_mat[:] = np.nan

            for i, j in valid_node_pairs:
                topology_cost_mat[i, j] = _local_graph_comparison_score(self.graph_matching_depth, i, j, simple_graphs_m, simple_graphs_n)

            cost_end = time.time()
            # print('Topology cost matrix takes {:.2f} s'.format(cost_end - cost_start))
            ### Topology matrix complete ###

            ### Build final cost matrix ###
            # add pseudo-counts for zero scores
            dist_cost_mat += 0.01 * np.nanmax(dist_cost_mat.ravel())
            topology_cost_mat += 0.01 * np.nanmax(topology_cost_mat.ravel())

            # construct linking cost matrix based on three matrics and relative scaling
            cost_m_n = dist_cost_mat ** self.dist_exponent * topology_cost_mat ** self.top_exponent

            # construct termination cost matrix
            cost_m_m = np.empty([number_m, number_m])
            cost_m_m[:] = np.nan
            for i in range(number_m):
                row = cost_m_n[i, :]
                if np.isnan(row).all():
                    cost_m_m[i, i] = 0  # must be assigned to itself since all other nodes exceed max radius
                else:
                    min_cost = np.nanmin(row)
                    cost_m_m[i, i] = 3 * min_cost

            # construct initiation cost matrix
            cost_n_n = np.empty([number_n, number_n])
            cost_n_n[:] = np.nan
            for j in range(number_n):
                column = cost_m_n[:, j]
                if np.isnan(column).all():
                    cost_n_n[j, j] = 0  # must be assigned to itself since all other nodes exceed max radius
                else:
                    min_cost = np.nanmin(column)
                    cost_n_n[j, j] = 3 * min_cost

            # construct auxiliary block
            cost_n_m = cost_m_n.T.copy()
            cost_n_m[:] = np.nanmin(cost_m_n)  # this matrix is needed for LAP solver but not used for tracking

            # assemble into one matrix
            left_block = np.concatenate((cost_m_n, cost_n_n), axis=0)
            right_block = np.concatenate((cost_m_m, cost_n_m), axis=0)
            cost_matrix = np.concatenate((left_block, right_block), axis=1)
            cost_matrix[np.isnan(cost_matrix)] = np.inf  # blocking values
            ### Final matrix is complete ###

            ### Solve LAP ###
            assignment = lap_solver(cost_matrix)[1]

            ### Remove unrealistic tracking ###
            assigned_m = assignment[:number_m]  # find linked nodes at frame m, n

            linked_m, linked_n = [], []
            for i in range(len(assigned_m)):
                if assigned_m[i] < number_n:
                    linked_m.append(i)
                    linked_n.append(assigned_m[i])

            filtered_nodes = []
            for segment_id in range(len(all_segment_nodes_m)):  # segment consists of segment nodes

                segment_nodes_m = all_segment_nodes_m[segment_id]

                # find only linked nodes in the segment
                current_seg_nodes_m = np.array(
                    [node for node in segment_nodes_m if node in linked_m and node not in branching_nodes_m])

                if len(current_seg_nodes_m) == 0:
                    continue

                # useful mappings
                node_m_to_seg_n, seg_n_to_node_m = _get_mappings(assignment, current_seg_nodes_m, node_to_segment_n)

                # node count of each segment
                node_count_per_seg = {}
                for seg_id in seg_n_to_node_m.keys():
                    node_count_per_seg[seg_id] = len(seg_n_to_node_m[seg_id])

                max_segment_id = max(node_count_per_seg, key=node_count_per_seg.get)

                # remove extremely long arrows
                majority_nodes_m = seg_n_to_node_m[max_segment_id]
                linked_majority_nodes_m = [n for n in majority_nodes_m if assignment[n] < number_n]
                other_nodes_m = [n for n in current_seg_nodes_m if n not in linked_majority_nodes_m]
                linked_other_nodes_m = [n for n in other_nodes_m if assignment[n] < number_n]

                majority_dists = [dist_cost_mat[i, j] for i, j in zip(linked_majority_nodes_m, [assignment[n] for n in linked_majority_nodes_m]) if not np.isnan(dist_cost_mat[i, j])]
                if len(majority_dists) == 0:
                    mean_majority_dist = np.nan
                else:
                    mean_majority_dist = np.mean(majority_dists)

                for node in linked_other_nodes_m:
                    if dist_cost_mat[node, assignment[node]] > 3 * mean_majority_dist:  # cutoff is here
                        filtered_nodes.append(node)

                # correct crossing arrows to align with the overall direction of the segment motion
                current_seg_nodes_m = [n for n in current_seg_nodes_m if n not in filtered_nodes]
                node_m_to_seg_n, seg_n_to_node_m = _get_mappings(assignment, current_seg_nodes_m, node_to_segment_n)
                if len(current_seg_nodes_m) == 0:
                    continue

                for seg_id in seg_n_to_node_m.keys():

                    if np.isnan(seg_id):
                        continue

                    current_nodes = seg_n_to_node_m[seg_id]  # nodes point to the same segment
                    assigned_nodes = [assignment[node] for node in current_nodes]
                    current_coords = [full_graph_m.vs[node]['coordinate'] for node in current_nodes]
                    assigned_coords = [full_graph_n.vs[node]['coordinate'] for node in assigned_nodes]
                    current_coords = np.array(current_coords)
                    assigned_coords = np.array(assigned_coords)

                    if len(current_coords) < 3:
                        continue

                    linking_vectors = assigned_coords - current_coords

                    mean_vector = np.mean(linking_vectors, axis=0)

                    angles = []
                    concerted_node_idx, outlier_node_idx = [], []
                    for node_idx, linking_vector in enumerate(linking_vectors):
                        dot = np.dot(mean_vector, linking_vector)  # scalar
                        norms = np.linalg.norm(mean_vector) * np.linalg.norm(linking_vector)

                        if norms == 0:
                            continue

                        cos_theta = np.clip(dot / norms, -1.0, 1.0)  # keep in [-1, 1]
                        angle = np.arccos(cos_theta)

                        angles.append(angle)
                        angle_cutoff = min(3 * np.std(angles), np.pi / 4)

                        if angle <= angle_cutoff:
                            concerted_node_idx.append(node_idx)
                        else:
                            outlier_node_idx.append(node_idx)

                    if len(concerted_node_idx) == 0:
                        continue

                    reference_vector = np.mean(linking_vectors[concerted_node_idx], axis=0)

                    concerted_nodes_n = [assigned_nodes[n] for n in concerted_node_idx]
                    outlier_nodes_m = [current_nodes[n] for n in outlier_node_idx]
                    for outlier_node_m in outlier_nodes_m:
                        expected_position = full_graph_m.vs[outlier_node_m]['coordinate'] + reference_vector
                        distance = np.linalg.norm(assigned_coords - expected_position, axis=1)
                        closest_nodes = [assigned_nodes[n] for n in np.argsort(distance)[:2]]

                        for node_n in closest_nodes:
                            if node_n not in concerted_nodes_n and node_n != assignment[outlier_node_m]:
                                assignment[assignment == node_n] = number_n + 1  # terminate old assignment
                                assignment[outlier_node_m] = node_n  # correct assignment
                                concerted_nodes_n.append(node_n)  # avoid overwrite assignment
                                break

            # update assignment after filtering
            assignment_filtered = assignment.copy()

            for index_m in filtered_nodes:
                index_n = linked_n[linked_m.index(index_m)]
                # set the linked node to initiated
                assignment_filtered[number_m + index_n] = index_n
                # set the node to terminated
                assignment_filtered[index_m] = number_n + index_m

            ### Report assignments ###
            assigned_m, assigned_n = assignment_filtered[:number_m], assignment_filtered[number_m:]

            linked, terminated, initiated = [], [], []
            for i in range(len(assigned_m)):
                if assigned_m[i] < number_n:
                    linked.append([i, assigned_m[i]])  # first being index for frame t and second for frame t+tracking_interval
                else:
                    terminated.append(i)

            for i in range(len(assigned_n)):
                if assigned_n[i] < number_n:
                    initiated.append(i)

            linked = np.array(linked); terminated = np.array(terminated); initiated = np.array(initiated)
            linked_nodes.append(linked); terminated_nodes.append(terminated); initiated_nodes.append(initiated)

            dist_cost_assigned = [dist_cost_mat[a, b] for (a, b) in linked]

            end = time.time()

            # output stats
            print(f"Number of nodes at frame {frame} and {frame + self.tracking_interval} are {number_m}, {number_n}")
            print(f"Number of nodes linked, terminated at frame {frame}, initiated at frame {frame + self.tracking_interval}: {len(linked)}, {len(terminated)}, {len(initiated)}")

            max_linked = min(number_m, number_n)
            percent_linked_min = 70
            percent_linked = round(len(linked) / max_linked, 3) * 100

            if percent_linked < percent_linked_min:
                warnings.warn(
                    f"Only {percent_linked}% of the {max_linked} nodes are tracked. "
                    "This is likely due to large distance or inconsistent topology between the two frames.")

            # print(f"Mean speed for tracked nodes: {(np.nanmean(dist_cost_assigned) / self.frame_interval):2f} μm/s")

            if (np.mean(dist_cost_assigned) / self.frame_interval) >= 1.0:
                warnings.warn(
                    "The mean node speed is greater than 1 μm/s. "
                    "This is extremely fast and tracking may be unreliable!")

            print(f"Tracking for frame {frame} and {frame + self.tracking_interval} is complete and took {(end - start):.2f} s\n")
            ### Assignments reported ###

            ### Update tracks ###
            nodes_m, nodes_n = linked[:, 0].tolist(), linked[:, 1].tolist()
            tracks_to_remove = []

            if frame == self.start_frame:
                for i in range(len(nodes_m)):
                    # initiate with first two frames
                    ongoing_tracks.append([[frame, frame + self.tracking_interval],
                                           [nodes_m[i], nodes_n[i]],
                                           [node_to_segment_m[nodes_m[i]], node_to_segment_n[nodes_n[i]]],
                                           [cc_m.membership[nodes_m[i]], cc_n.membership[nodes_n[i]]],
                                           [coords_m[nodes_m[i]], coords_n[nodes_n[i]]],
                                           [intensity_m[nodes_m[i]], intensity_n[nodes_n[i]]],
                                           [width_m[nodes_m[i]], width_n[nodes_n[i]]], ])

            else:
                for idx, track in enumerate(ongoing_tracks):
                    # if terminated, remove track; else append linked node
                    if track[1][-1] in terminated:
                        terminated_tracks.append(track)
                        tracks_to_remove.append(idx)
                    else:
                        linked_index = nodes_m.index(track[1][-1])
                        linked_node = nodes_n[linked_index]
                        track[0].append(frame + self.tracking_interval)
                        track[1].append(linked_node)
                        track[2].append(node_to_segment_n[linked_node])
                        track[3].append(cc_n.membership[linked_node])
                        track[4].append(coords_n[linked_node])
                        track[5].append(intensity_n[linked_node])
                        track[6].append(width_n[linked_node])

            # delete terminated tracks from ongoing tracks
            ongoing_tracks = [t for i, t in enumerate(ongoing_tracks) if i not in tracks_to_remove]

            for init_node in initiated:
                ongoing_tracks.append([[frame + self.tracking_interval],
                                       [init_node],
                                       [node_to_segment_n[init_node]],
                                       [cc_n.membership[init_node]],
                                       [coords_n[init_node]],
                                       [intensity_n[init_node]],
                                       [width_n[init_node]]])

        linked_nodes = np.array(linked_nodes, dtype=object); terminated_nodes = np.array(terminated_nodes, dtype=object); initiated_nodes = np.array(initiated_nodes, dtype=object)

        # each element in all_tracks is a track with 1) frame numbers; 2) node indices; 3) segment ids of the node, 4) frag ids of the node; 4) node coords; 5) node intensities; 6) node widths
        terminated_tracks = np.array(terminated_tracks, dtype=object); ongoing_tracks = np.array(ongoing_tracks, dtype=object)
        all_tracks = np.concatenate([terminated_tracks, ongoing_tracks])

        # filter out too short tracks
        short_tracks = []
        for i in range(len(all_tracks)):
            if len(all_tracks[i][0]) < self.min_track_size:
                short_tracks.append(i)

        all_tracks = np.delete(all_tracks, short_tracks, axis=0)
        all_tracks = sorted(all_tracks, key=lambda t: t[0][0])  # sort by start frame
        num_tracks = len(all_tracks)

        ### Gap closing ###
        print('\nInitiate gap closing ...')

        # get track displacements
        all_track_disps = []
        for t in all_tracks:
            track_coords = t[4]  # use index for the coordinates
            all_track_disps.append(
                [np.linalg.norm(track_coords[t + 1] - track_coords[t]) for t in range(len(track_coords) - 1)])

        all_track_assignments = {}
        partition_start = 0
        iter_num = 1

        while partition_start < num_tracks:

            partition_size = int(self.block_size_factor * num_tracks)
            overlap_size = int(0.2 * partition_size)

            partition_end = partition_start + partition_size
            if partition_end > num_tracks:
                partition_end = num_tracks
            print('Block ' + str(iter_num) + ' index:', partition_start, 'to', partition_end)
            print('Computing cost terms for block ' + str(iter_num))

            track_cost_m_n = np.empty([partition_size, partition_size])
            track_cost_m_n[:] = np.nan

            for i in range(partition_start, partition_end):

                track_frames_m, track_nodes_m, track_coords_m = all_tracks[i][0], all_tracks[i][1], all_tracks[i][4]
                end_frame_m, end_node_m, end_coord_m = track_frames_m[-1], track_nodes_m[-1], track_coords_m[-1]

                simple_graphs_m = local_simple_graphs_all_frames[end_frame_m]
                disps_m = all_track_disps[i]

                for j in range(i+1, partition_end):  # no need to check index less than i because the start frame is sorted

                    track_frames_n, track_nodes_n, track_coords_n = all_tracks[j][0], all_tracks[j][1], all_tracks[j][4]
                    start_frame_n, start_node_n, start_coord_n = track_frames_n[0], track_nodes_n[0], track_coords_n[0]

                    gap_size = start_frame_n - end_frame_m - 1

                    # check only those within max gap size
                    if 1 <= gap_size <= self.max_gap_size:  # if gap == 1 it should have been linked before - so skip it

                        # load data
                        simple_graphs_n = local_simple_graphs_all_frames[start_frame_n]

                        # compute distance cutoff based on the two tracks
                        disps_n = all_track_disps[j]
                        comb_disps = disps_m + disps_n

                        dist_cutoff = (gap_size + 1) * (3 * np.std(comb_disps)) ** 2

                        # compute node-to-node distance
                        dist = end_coord_m - start_coord_n
                        dist_cost = np.sum(dist ** 2)

                        # filter by distance cutoff
                        if dist_cost > dist_cutoff:
                            continue

                        # compute topology cost
                        topology_cost = _local_graph_comparison_score(self.graph_matching_depth, end_node_m, start_node_n,
                                                                      simple_graphs_m, simple_graphs_n)

                        # assign g.c. cost
                        track_cost_m_n[i - partition_start, j - partition_start] = dist_cost ** self.dist_exponent * topology_cost ** self.top_exponent

            # construct termination cost matrix
            track_cost_m_m = np.empty([partition_size, partition_size])
            track_cost_m_m[:] = np.nan

            for i in range(partition_size):
                row = track_cost_m_n[i, :]

                if np.isnan(row).all():
                    track_cost_m_m[i, i] = 0  # must be assigned to itself since all other nodes exceed max radius
                else:
                    min_cost = np.nanmin(row)
                    track_cost_m_m[i, i] = 3 * min_cost

            # assemble into one matrix
            track_cost_matrix = np.concatenate((track_cost_m_n, track_cost_m_m), axis=1)
            track_cost_matrix[np.isnan(track_cost_matrix)] = np.inf

            # evaluate memory usage
            # print('Cost matrix memory usage: {:.1f} MB\n'.format(track_cost_matrix.nbytes / 1024 ** 2))

            # solve LAP and store linking results
            assignment = lap_solver(track_cost_matrix)[1]

            linked = []
            for i in range(len(assignment)):
                if assignment[i] < partition_size:
                    linked.append(
                        [i, assignment[i]])  # first being index for frame t and second for frame t+tracking_interval

            for pair in linked:
                all_track_assignments[partition_start + pair[0]] = partition_start + pair[
                    1]  # offset by start index of the partition

            # go to next partition
            if self.block_size_factor == 1:
                break
            else:
                partition_start = partition_start + partition_size - overlap_size
                iter_num += 1

        # convert dictionary to array and overwrite assignment by next partition's assignment for the overlapped region
        linked_tracks = np.zeros([len(all_track_assignments.keys()), 2], dtype=int)
        for i, a in enumerate(all_track_assignments.keys()):
            linked_tracks[i, 0] = a
            linked_tracks[i, 1] = all_track_assignments[a]

        # combine tracks for gap closing
        print('Start combining closed tracks ...')
        if linked_tracks.shape[0] > 0:
            linked_tracks_for_update = linked_tracks.copy()  # used to record appended tracks
            tracks_of_track = []  # list of linked tracks
            all_linked_tracks = []  # record tracks that are closed

            # recursive function for finding linked tracks
            def find_all_linked_tracks(tracks, track_id):
                for i in range(0, len(linked_tracks_for_update)):
                    if linked_tracks_for_update[i, 0] == track_id:  # find the first track
                        tracks.append(linked_tracks_for_update[i, 1])
                        linked_tracks_for_update[i, 0] = -1  # note that the track is already appended
                        find_all_linked_tracks(tracks, linked_tracks_for_update[
                            i, 1])  # go to the next track and find linked_tracks track

            # for each track find the series of linked tracks
            for index in range(len(linked_tracks)):
                track_id = linked_tracks[index, 0]
                if track_id in linked_tracks_for_update[:, 0]:
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
            sort_by_length = sorted(all_closed_tracks, key=lambda track: len(track[0]),
                                    reverse=True)  # first sort by size of track
            sort_by_start = sorted(sort_by_length, key=lambda track: track[0][0])  # then sort by start frame
            all_closed_tracks = sort_by_start

        else:
            all_closed_tracks = all_tracks

        print(f"Number of tracks and average track length before gap closing: {len(all_tracks)}, {np.mean([len(track[0]) for track in all_tracks]):.2f}")
        print(f"Number of tracks and average track length after gap closing: {len(all_closed_tracks)}, {np.mean([len(track[0]) for track in all_closed_tracks]):.2f}")

        print('\nSaving final node trajectory file. This might take a few minutes for large files.')

        records = []
        append = records.append 

        for track_id, track in enumerate(all_closed_tracks):
            track_frames, track_nodes, track_segs, track_frags, track_coords, track_ints, track_widths = track
            for f in range(len(track_frames)):
                x, y, z = track_coords[f]
                append((
                    track_frames[f], track_id, track_nodes[f], track_segs[f], track_frags[f],
                    x, y, z, track_ints[f], track_widths[f]
                ))

        tracks = pd.DataFrame.from_records(
            records,
            columns=[
                'frame_id', 'unique_node_id', 'frame_node_id', 'frame_seg_id', 'frame_frag_id',
                'x', 'y', 'z', 'intensity', 'width'
            ]
        )
        tracks.sort_values(['unique_node_id', 'frame_id'], inplace=True, ignore_index=True)


        # add connected_unique_nodes for tracked graph creation by skipping untracked nodes
        final_tracks = pd.DataFrame()
        for frame in range(self.start_frame, self.end_frame - self.tracking_interval, self.tracking_interval):

            full_graph = full_graphs_all_frames[frame]
            tracks_frame = tracks[tracks['frame_id'] == frame]

            tracked_frame_nodes = tracks_frame['frame_node_id'].astype('int').tolist()
            unique_nodes = tracks_frame['unique_node_id'].astype('int').tolist()
            frame_to_unique = {tracked_frame_nodes[i]: unique_nodes[i] for i in range(len(tracks_frame))}

            all_connected_unique_nodes = []
            for node in tracked_frame_nodes:
                connected_unique_nodes = _find_connected_unique_nodes(
                    node, full_graph, tracked_frame_nodes, frame_to_unique
                )
                all_connected_unique_nodes.append(connected_unique_nodes)

            tracks_frame.insert(2, 'connected_unique_node_id', [_list_to_str(a) for a in all_connected_unique_nodes])
            final_tracks = pd.concat([final_tracks, tracks_frame])  # accumulate tracks from each frame

        final_tracks = final_tracks[[
            'frame_id', 'unique_node_id', 'frame_node_id', 'frame_seg_id', 'frame_frag_id',
            'connected_unique_node_id', 'x', 'y', 'z', 'intensity', 'width'
        ]]

        # save data
        final_tracks.to_csv(self.segmented_mito.save_path + 'mito_node_tracks.csv', index=False)
        np.save(self.segmented_mito.save_path + 'mito_linked_nodes', linked_nodes)
        print('Tracking is complete and files are saved!')

        return TrackedMito(self.segmented_mito, self.frame_interval, self.start_frame, self.end_frame, self.tracking_interval,
                           final_tracks, linked_nodes)



def _list_to_str(list1):

    string = ""
    if len(list1) == 0:
        return string

    for i in range(len(list1) - 1):
        string += str(list1[i])
        string += " "
    string += str(list1[-1])

    return string


def _dissimilarity_score(sel_edge_len_m, sel_edge_len_n):
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
            cost_mat[i, j] = abs(sel_edge_len_m[i] - sel_edge_len_n[j]) / max(sel_edge_len_m[i], sel_edge_len_n[
                j])  # should never have two zeros

    # solve LAP to find minimum score
    assigned = lap_solver(cost_mat)[1]
    min_sum_cost = 0
    for i in range(num_edges):
        min_sum_cost += cost_mat[i, assigned[i]]

    return min_sum_cost


def _local_graph_comparison_score(depth, node_i, node_j, contracted_graphs_m, contracted_graphs_n):
    frag_m, frag_n = contracted_graphs_m[node_i].copy(), contracted_graphs_n[node_j].copy()
    root_m, root_n = frag_m.vs['index'].index(node_i), frag_n.vs['index'].index(node_j)

    node_mapping = {root_m: root_n}

    # iterate each level
    for n_level in range(depth):

        visited_nodes_m = node_mapping.keys()
        visited_nodes_n = node_mapping.values()

        if n_level > 0:
            last_level_m = frag_m.neighborhood(vertices=root_m, order=n_level - 1, mindist=n_level - 1)
            last_level_n = frag_n.neighborhood(vertices=root_n, order=n_level - 1, mindist=n_level - 1)
        else:
            last_level_m = []
            last_level_n = []

        # abort once the graph is fully mapped
        current_level_m = frag_m.neighborhood(vertices=root_m, order=n_level, mindist=n_level)
        if len(current_level_m) == 0:
            break

        for node_m in current_level_m:
            node_n = node_mapping[node_m]  # use mapping determined from last level

            neighbors_m = frag_m.neighbors(node_m)
            neighbors_n = frag_n.neighbors(node_n)

            # replace each cycle edge with two pseudo-edges of same lengths and add two pseudo-nodes
            for neigh in neighbors_m:
                if neigh in visited_nodes_m and neigh not in last_level_m:
                    dist = frag_m.es[frag_m.get_eid(node_m, neigh)]['distance']
                    frag_m.delete_edges(frag_m.get_eid(node_m, neigh))
                    frag_m.add_vertices(2)
                    frag_m.add_edges([[node_m, frag_m.vs[-2].index]], {'distance': dist})
                    frag_m.add_edges([[neigh, frag_m.vs[-1].index]], {'distance': dist})

            for neigh in neighbors_n:
                if neigh in visited_nodes_n and neigh not in last_level_n:
                    dist = frag_n.es[frag_n.get_eid(node_n, neigh)]['distance']
                    frag_n.delete_edges(frag_n.get_eid(node_n, neigh))
                    frag_n.add_vertices(2)
                    frag_n.add_edges([[node_n, frag_n.vs[-2].index]], {'distance': dist})
                    frag_n.add_edges([[neigh, frag_n.vs[-1].index]], {'distance': dist})

            # add pseudo-nodes and pseudo-edges of 0 at this level
            num_node_diff = len(neighbors_m) - len(neighbors_n)
            if num_node_diff >= 0:
                for it in range(num_node_diff):
                    frag_n.add_vertices(1)
                    frag_n.add_edges([[node_n, frag_n.vs[-1].index]], {'distance': 0})
            else:
                for it in range(-1 * num_node_diff):
                    frag_m.add_vertices(1)
                    frag_m.add_edges([[node_m, frag_m.vs[-1].index]], {'distance': 0})

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
            index_mapping_m = {i: neighbors_m[i] for i in range(len(neighbors_m))}
            index_mapping_n = {i: neighbors_n[i] for i in range(len(neighbors_n))}

            num_node = max(len(neighbors_m), len(neighbors_n))
            cost_mat = np.zeros((num_node, num_node))

            # fill cost matrix with dissimilarity scores of each node m and n
            for i in range(num_node):
                for j in range(num_node):
                    sel_edge_len_m = frag_m.es[frag_m.incident(index_mapping_m[i])]['distance']
                    sel_edge_len_n = frag_n.es[frag_n.incident(index_mapping_n[j])]['distance']
                    cost_mat[i, j] = _dissimilarity_score(sel_edge_len_m, sel_edge_len_n)

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

    index_mapping_m = {sorted(node_mapping.keys())[i]: i for i in range(len(node_mapping.keys()))}
    index_mapping_n = {sorted(node_mapping.values())[i]: i for i in range(len(node_mapping.values()))}

    for i in node_mapping.keys():
        for j in frag_m.neighbors(i):
            if j in node_mapping.keys():
                weighted_adj_mat_m[index_mapping_m[i], index_mapping_m[j]] = \
                    frag_m.es[frag_m.get_eid(i, j)]['distance']

    for i in node_mapping.keys():
        mapped_i = node_mapping[i]
        for mapped_j in frag_n.neighbors(mapped_i):
            if mapped_j in node_mapping.values():
                j = reverse_node_mapping[mapped_j]
                weighted_adj_mat_n[index_mapping_m[i], index_mapping_m[j]] = \
                    frag_n.es[frag_n.get_eid(mapped_i, mapped_j)]['distance']

    # calculate Euclidean distance between the two weighted adjacency matrices
    weight_diff = weighted_adj_mat_m - weighted_adj_mat_n

    score = np.sum((weight_diff.ravel()) ** 2)

    return score


def _get_mappings(assignment, segment_nodes, node_to_seg_mapping):
    node_m_to_seg_n = {}  # node at t to tracked segment at t+1
    seg_n_to_node_m = {}  # tracked segment at t+1 to all nodes at t

    for node in segment_nodes:
        linked_node = assignment[node]
        if linked_node in node_to_seg_mapping.keys():  # node to seg mapping for all seg
            linked_seg = node_to_seg_mapping[linked_node]
        else:
            linked_seg = np.nan

        node_m_to_seg_n[node] = linked_seg

        if linked_seg in seg_n_to_node_m.keys():
            seg_n_to_node_m[linked_seg] += [node]
        else:
            seg_n_to_node_m[linked_seg] = [node]

    return node_m_to_seg_n, seg_n_to_node_m


def _find_connected_unique_nodes(this_node, full_graph, tracked_frame_nodes, frame_to_unique, visited=None):
    if visited is None:
        visited = set()

    visited.add(this_node)
    connected_unique_nodes = []

    # Get all neighbors of this node
    neighs = full_graph.neighbors(this_node)

    for neigh in neighs:
        if neigh in visited:
            continue

        # If the neighbor is tracked in this frame, add its unique ID
        if neigh in tracked_frame_nodes:
            connected_unique_nodes.append(frame_to_unique[neigh])
        else:
            # Recursively search through untracked nodes
            deeper_nodes = _find_connected_unique_nodes(
                neigh, full_graph, tracked_frame_nodes, frame_to_unique, visited
            )
            connected_unique_nodes.extend(deeper_nodes)

    return connected_unique_nodes