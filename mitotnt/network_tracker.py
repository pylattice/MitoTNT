import time
import warnings
import numpy as np
import pandas as pd
from tqdm.auto import tqdm, trange
from scipy.optimize import linear_sum_assignment as lap_solver
from scipy.spatial import cKDTree
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import min_weight_full_bipartite_matching
import lap  # sparse Jonker-Volgenant (lapmod) for large sparse assignment problems
import multiprocessing as mp
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
                 min_track_size: int = 4, max_gap_size: int = 3, block_size_factor: float = 1.0,
                 num_workers: int = 20):

        self.segmented_mito = segmented_mito
        self.start_frame = start_frame
        if end_frame is None:
            self.end_frame = self.segmented_mito.num_frames
        else:
            self.end_frame = end_frame
        if frame_interval is None:
            raise ValueError("frame_interval (seconds per frame) is required and must be provided.")
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
        # Bound the worker count so the process pools stay sane on shared machines.
        self.num_workers = max(1, min(int(num_workers), 20))
        # Per-stage wall-clock accumulators, filled during run().
        self._timing = {}


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

        # per-stage timing accumulators (seconds)
        self._timing = {'candidates': 0.0, 'topology': 0.0, 'lap': 0.0,
                        'postfilter': 0.0, 'gap_closing': 0.0, 'save': 0.0}

        num_frames = len(full_graphs_all_frames)
        if self.end_frame > num_frames:
            warnings.warn(f"The end frame specified ({self.end_frame}) is greater than the number of frames "
                          f"in tracking inputs ({num_frames}). End frame has been clamped to {num_frames - 1}.")
            self.end_frame = num_frames - 1

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
            max_number = 50000
            if number_m > max_number or number_n > max_number:
                warnings.warn(
                    f"Large frame: up to {max(number_m, number_n)} nodes (> {max_number}); tracking will be slower. "
                    "To reduce the node count, increase node_gap_size in SkeletonizedMito.extract_graphs().")

            fluctuation_percent_max = 30
            fluctuation_percent = abs(number_n - number_m) / number_m * 100
            if fluctuation_percent > fluctuation_percent_max:
                warnings.warn(
                    f"The number of nodes changes by {fluctuation_percent:.1f}% between the two frames. "
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

            ### Calculate sparse candidate links (replaces dense distance matrix) ###
            # Avoid the dense (m x n) distance matrix and its per-row sort: a
            # KD-tree keeps only each m-node's nearest candidate n-nodes,
            # reproducing the neighbour + displacement cutoffs (modulo arbitrary
            # tie-breaking at exactly the cutoff distance) while staying sparse.
            cost_start = time.time()

            coords_m_mat = np.asarray(coords_m, dtype=float)
            coords_n_mat = np.asarray(coords_n, dtype=float)

            tree_n = cKDTree(coords_n_mat)

            # nearest-neighbour distance per m-node -> displacement cutoff
            nn_dist, _ = tree_n.query(coords_m_mat, k=1, workers=self.num_workers)
            nn_dist = np.atleast_1d(np.asarray(nn_dist, dtype=float))

            if self.cutoff_speed is None:
                if nn_dist.size == 0:
                    disp_cutoff = np.inf
                else:
                    disp_cutoff = np.mean(nn_dist) + 3 * np.std(nn_dist)  # global estimate based on all nodes
            else:
                disp_cutoff = self.cutoff_speed * self.frame_interval

            # Keep up to (cutoff_num_neighbor + 1) nearest n-nodes within disp_cutoff.
            # The +1 reproduces the original `sorted(row)[cutoff_num_neighbor]` threshold,
            # which retains the (cutoff_num_neighbor + 1) smallest distances.
            k_query = min(self.cutoff_num_neighbor + 1, number_n)
            cand_dist, cand_idx = tree_n.query(coords_m_mat, k=k_query, workers=self.num_workers)
            if k_query == 1:
                cand_dist = cand_dist[:, None]
                cand_idx = cand_idx[:, None]

            keep = cand_dist <= disp_cutoff
            cand_i = np.repeat(np.arange(number_m), k_query)[keep.ravel()]
            cand_j = cand_idx.ravel()[keep.ravel()].astype(np.int64)
            cand_d = cand_dist.ravel()[keep.ravel()]

            self._timing['candidates'] += time.time() - cost_start
            ### Candidate links complete ###

            ### Calculate topology cost for candidate pairs only ###
            cost_start = time.time()
            cand_topo = self._compute_topology_costs(cand_i, cand_j,
                                                     simple_graphs_m, simple_graphs_n)
            self._timing['topology'] += time.time() - cost_start
            ### Topology complete ###

            ### Build sparse cost terms and solve LAP ###
            cost_start = time.time()

            # add pseudo-counts for zero scores (reproduces dense `+= 0.01 * nanmax`)
            max_dist = cand_d.max() if cand_d.size else 0.0
            max_topo = cand_topo.max() if cand_topo.size else 0.0
            pdist = cand_d + 0.01 * max_dist
            ptopo = cand_topo + 0.01 * max_topo

            # linking cost for each candidate pair
            link_cost = pdist ** self.dist_exponent * ptopo ** self.top_exponent

            # sparse stand-in for the dense (pseudo-counted) distance matrix used
            # downstream during segment-based filtering / arrow correction
            dist_cost_mat = _SparseDist(cand_i, cand_j, pdist, number_m, number_n)

            # Solve the linear assignment problem on the sparse augmented bipartite
            # graph (links + termination + initiation + candidate-pattern auxiliary
            # block). This is mathematically equivalent to the dense formulation:
            # the auxiliary block uses constant cost C = min(link_cost) on exactly
            # the candidate pattern, and for every link (i, j) the symmetric edge
            # (birth_j -> death_i) exists, so any link set feasible in the dense
            # problem is feasible here at identical cost.
            assignment = _solve_lap_sparse(cand_i, cand_j, link_cost, number_m, number_n)

            self._timing['lap'] += time.time() - cost_start
            ### LAP solved ###

            _postfilter_start = time.time()

            ### Remove unrealistic tracking ###
            assigned_m = assignment[:number_m]  # find linked nodes at frame m, n

            linked_m, linked_n = [], []
            for i in range(len(assigned_m)):
                if assigned_m[i] < number_n:
                    linked_m.append(i)
                    linked_n.append(assigned_m[i])

            # Membership sets (linked_m is also used as an ordered list elsewhere
            # via .index(), so keep both). Avoids O(seg * |linked_m|) list scans.
            linked_m_set = set(linked_m)
            branching_nodes_m_set = set(branching_nodes_m)

            filtered_nodes = []
            for segment_id in range(len(all_segment_nodes_m)):  # segment consists of segment nodes

                segment_nodes_m = all_segment_nodes_m[segment_id]

                # find only linked nodes in the segment
                current_seg_nodes_m = np.array(
                    [node for node in segment_nodes_m if node in linked_m_set and node not in branching_nodes_m_set])

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

            # map each linked m-node to its n-node (linked_m entries are unique),
            # avoiding a linear linked_m.index() per filtered node
            linked_m_to_n = dict(zip(linked_m, linked_n))
            for index_m in filtered_nodes:
                index_n = linked_m_to_n[index_m]
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
            percent_linked = len(linked) / max_linked * 100

            if percent_linked < percent_linked_min:
                warnings.warn(
                    f"Only {percent_linked:.1f}% of the {max_linked} nodes are tracked. "
                    "This is likely due to large distance or inconsistent topology between the two frames.")

            if (np.mean(dist_cost_assigned) / self.frame_interval) >= 1.0:
                warnings.warn(
                    "The mean node speed is greater than 1 μm/s. "
                    "This is extremely fast and tracking may be unreliable!")

            print(f"Tracking for frame {frame} and {frame + self.tracking_interval} is complete and took {(end - start):.2f} s\n")
            ### Assignments reported ###

            ### Update tracks ###
            nodes_m, nodes_n = linked[:, 0].tolist(), linked[:, 1].tolist()
            tracks_to_remove = []
            # O(1) lookups for the per-track update below (nodes_m are unique
            # m-node indices; terminated is a node-index array): replaces a
            # linear nodes_m.index() and an array `in` scan per ongoing track.
            m_to_n = dict(zip(nodes_m, nodes_n))
            terminated_set = set(terminated.tolist())

            if frame == self.start_frame:
                for i in range(len(nodes_m)):
                    # initiate with first two frames
                    ongoing_tracks.append([[frame, frame + self.tracking_interval],
                                           [nodes_m[i], nodes_n[i]],
                                           [node_to_segment_m.get(nodes_m[i],np.nan),node_to_segment_n.get(nodes_n[i],np.nan)],
                                           [cc_m.membership[nodes_m[i]], cc_n.membership[nodes_n[i]]],
                                           [coords_m[nodes_m[i]], coords_n[nodes_n[i]]],
                                           [intensity_m[nodes_m[i]], intensity_n[nodes_n[i]]],
                                           [width_m[nodes_m[i]], width_n[nodes_n[i]]], ])

            else:
                for idx, track in enumerate(ongoing_tracks):
                    # if terminated, remove track; else append linked node
                    if track[1][-1] in terminated_set:
                        terminated_tracks.append(track)
                        tracks_to_remove.append(idx)
                    else:
                        linked_node = m_to_n[track[1][-1]]
                        track[0].append(frame + self.tracking_interval)
                        track[1].append(linked_node)
                        track[2].append(node_to_segment_n.get(linked_node,np.nan))
                        track[3].append(cc_n.membership[linked_node])
                        track[4].append(coords_n[linked_node])
                        track[5].append(intensity_n[linked_node])
                        track[6].append(width_n[linked_node])

            # delete terminated tracks from ongoing tracks (set membership)
            remove_set = set(tracks_to_remove)
            ongoing_tracks = [t for i, t in enumerate(ongoing_tracks) if i not in remove_set]

            for init_node in initiated:
                ongoing_tracks.append([[frame + self.tracking_interval],
                                       [init_node],
                                       [node_to_segment_n.get(init_node,np.nan)],
                                       [cc_n.membership[init_node]],
                                       [coords_n[init_node]],
                                       [intensity_n[init_node]],
                                       [width_n[init_node]]])

            # accumulate post-LAP filtering + track-update time for this frame pair
            self._timing['postfilter'] += time.time() - _postfilter_start

        linked_nodes = np.array(linked_nodes, dtype=object); terminated_nodes = np.array(terminated_nodes, dtype=object); initiated_nodes = np.array(initiated_nodes, dtype=object)

        # Each track is a list of 7 per-node fields: 1) frame numbers; 2) node
        # indices; 3) segment ids; 4) fragment ids; 5) coordinates; 6) intensities;
        # 7) widths. Build an explicit (num_tracks, 7) object array. Letting
        # np.array(..., dtype=object) infer the shape mis-collapses to 2-D/3-D when
        # tracks happen to have uniform-length sub-lists (e.g. a single frame pair),
        # which breaks downstream; an explicit (N, 7) shape is robust in all cases
        # and keeps each row indexable/`.tolist()`-able as the rest of run() expects.
        all_tracks_list = list(terminated_tracks) + list(ongoing_tracks)
        all_tracks = np.empty((len(all_tracks_list), 7), dtype=object)
        for _idx, _track in enumerate(all_tracks_list):
            all_tracks[_idx] = _track

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
        _gap_start = time.time()

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

            # Compute the finite (within-cutoff) gap-closing costs (parallelized),
            # then solve the block assignment on a sparse augmented matrix with
            # lap.lapmod instead of a dense [P, 2P] scipy solve (see
            # _solve_gap_lap_sparse). This avoids the dense P x P allocation and
            # scales to large track counts.
            entries = self._compute_gap_cost_entries(
                all_tracks, all_track_disps, local_simple_graphs_all_frames,
                partition_start, partition_end)

            for li, lj in _solve_gap_lap_sparse(entries, partition_start, partition_end):
                # offset local block indices back to global track indices
                all_track_assignments[partition_start + li] = partition_start + lj

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
            tracks_of_track = []  # list of linked-track chains
            all_linked_tracks = []  # record tracks that are closed

            # Each source track maps to exactly one target; follow source->target
            # links to build chains. An adjacency dict + a "consumed sources" set
            # replaces the recursive O(num_tracks^2) scan while reproducing it
            # exactly: a source row is consumed only once it is followed, so a
            # track can still be appended as the *target* of a later chain even
            # if it already headed its own chain.
            src_to_tgt = {int(s): int(t) for s, t in linked_tracks}
            consumed = set()  # sources whose link has been followed
            for index in range(len(linked_tracks)):
                track_id = int(linked_tracks[index, 0])
                if track_id in consumed:
                    continue
                chain = [track_id]
                current = track_id
                while current in src_to_tgt and current not in consumed:
                    target = src_to_tgt[current]
                    chain.append(target)
                    consumed.add(current)
                    current = target
                tracks_of_track.append(chain)
                all_linked_tracks += chain

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

            # add unclosed tracks back  (set membership: was O(num_tracks^2))
            all_linked_set = set(all_linked_tracks)
            for t in range(num_tracks):
                if t not in all_linked_set:
                    all_closed_tracks.append(all_tracks[t].tolist())

            # sort tracks
            sort_by_length = sorted(all_closed_tracks, key=lambda track: len(track[0]),
                                    reverse=True)  # first sort by size of track
            sort_by_start = sorted(sort_by_length, key=lambda track: track[0][0])  # then sort by start frame
            all_closed_tracks = sort_by_start

        else:
            all_closed_tracks = all_tracks

        self._timing['gap_closing'] += time.time() - _gap_start

        print(f"Number of tracks and average track length before gap closing: {len(all_tracks)}, {np.mean([len(track[0]) for track in all_tracks]):.2f}")
        print(f"Number of tracks and average track length after gap closing: {len(all_closed_tracks)}, {np.mean([len(track[0]) for track in all_closed_tracks]):.2f}")

        print('\nSaving final node trajectory file. This might take a few minutes for large files.')
        _save_start = time.time()

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
        self._timing['save'] += time.time() - _save_start
        print('Tracking is complete and files are saved!')

        # per-stage timing summary (helps locate remaining bottlenecks)
        print('\n=== Tracking stage timing (s) ===')
        for _stage, _dt in self._timing.items():
            print(f"  {_stage:<12}: {_dt:8.2f}")
        print(f"  {'TOTAL':<12}: {sum(self._timing.values()):8.2f}")

        return TrackedMito(self.segmented_mito, self.frame_interval, self.start_frame, self.end_frame, self.tracking_interval,
                           final_tracks, linked_nodes)


    def _compute_topology_costs(self, cand_i, cand_j, simple_graphs_m, simple_graphs_n):
        """
        Compute the local-graph topology comparison score for each candidate
        (m-node, n-node) pair. Every pair is independent, so for large workloads
        this fans out across up to `self.num_workers` processes.

        Parallelism uses the 'fork' start method so each worker inherits the two
        local-graph lists via copy-on-write — only the lightweight index chunks are sent to workers,
        never the graphs themselves. The score function copies each fragment
        before mutating it, so the shared graphs are read-only and fork-safe.
        """
        depth = self.graph_matching_depth
        n_pairs = len(cand_i)
        if n_pairs == 0:
            return np.empty(0, dtype=float)

        # Serial path for small workloads (avoids process-pool setup overhead).
        if self.num_workers <= 1 or n_pairs < 5000:
            out = np.empty(n_pairs, dtype=float)
            for k in range(n_pairs):
                out[k] = _local_graph_comparison_score(
                    depth, int(cand_i[k]), int(cand_j[k]), simple_graphs_m, simple_graphs_n)
            return out

        # Parallel path: publish graphs as globals BEFORE forking so workers
        # inherit them; pass only index chunks through the pool.
        global _WORKER_TOPO
        _WORKER_TOPO = {'m': simple_graphs_m, 'n': simple_graphs_n, 'depth': depth}
        try:
            n_chunks = min(n_pairs, self.num_workers * 4)
            split = np.array_split(np.arange(n_pairs), n_chunks)
            chunk_args = [(cand_i[c], cand_j[c]) for c in split]
            ctx = mp.get_context('fork')
            with ctx.Pool(self.num_workers) as pool:
                results = pool.map(_topo_worker, chunk_args)
        finally:
            _WORKER_TOPO = {}  # release references in the parent

        out = np.empty(n_pairs, dtype=float)
        for c, r in zip(split, results):
            out[c] = r
        return out


    def _compute_gap_cost_entries(self, all_tracks, all_track_disps,
                                  local_simple_graphs_all_frames,
                                  partition_start, partition_end):
        """
        Compute the gap-closing cost for candidate track pairs (i, j) with i in
        [partition_start, partition_end). Candidate targets are found with a
        per-start-frame KD-tree on track start coordinates instead of an
        O(num_tracks^2) all-pairs scan: only tracks whose start lies within a
        safe radius of source i's end point can satisfy the distance cutoff. The
        exact gap-size / distance / topology checks are still applied to every
        returned pair, so the result is identical to the brute-force version.
        Independent per source track, so it fans out across up to
        `self.num_workers` processes (mirrors `_compute_topology_costs`).
        Returns a flat list of (i, j, cost).
        """
        idxs = list(range(partition_start, partition_end))
        n = len(idxs)
        if n == 0:
            return []

        # Per-track endpoints (keyed by global track index), grouped by start frame.
        end_frame, end_coord, end_node = {}, {}, {}
        start_node, start_coord = {}, {}
        by_start_frame = {}
        max_disp = 0.0
        for idx in idxs:
            frames, nodes, coords = all_tracks[idx][0], all_tracks[idx][1], all_tracks[idx][4]
            end_frame[idx] = frames[-1]
            end_node[idx] = nodes[-1]
            end_coord[idx] = np.asarray(coords[-1], dtype=float)
            start_node[idx] = nodes[0]
            start_coord[idx] = np.asarray(coords[0], dtype=float)
            by_start_frame.setdefault(frames[0], []).append(idx)
            disps = all_track_disps[idx]
            if disps:
                dm = max(disps)
                if dm > max_disp:
                    max_disp = dm

        # One KD-tree of start coordinates per start frame; tree point k maps
        # back to global track index tree_idx[sf][k].
        trees, tree_idx = {}, {}
        for sf, members in by_start_frame.items():
            tree_idx[sf] = members
            trees[sf] = cKDTree(np.array([start_coord[i] for i in members], dtype=float))

        # Safe search radius: a pair's distance cutoff is
        # (gap+1)*(3*std(comb_disps))^2, and the std of non-negative displacements
        # is bounded by their maximum value, so the cutoff is <=
        # (max_gap+1)*(3*max_disp)^2 = radius^2. Querying within `radius` therefore
        # keeps every true candidate; the exact cutoff is re-checked per pair.
        radius = 3.0 * max_disp * np.sqrt(self.max_gap_size + 1)

        global _WORKER_GAP
        _WORKER_GAP = {'disps': all_track_disps, 'lsg': local_simple_graphs_all_frames,
                       'max_gap': self.max_gap_size, 'depth': self.graph_matching_depth,
                       'dist_exp': self.dist_exponent, 'top_exp': self.top_exponent,
                       'end_frame': end_frame, 'end_coord': end_coord, 'end_node': end_node,
                       'start_node': start_node, 'start_coord': start_coord,
                       'trees': trees, 'tree_idx': tree_idx, 'radius': radius,
                       'tracking_interval': self.tracking_interval}
        try:
            # Serial path for small blocks (avoids process-pool setup overhead).
            if self.num_workers <= 1 or n < 2000:
                return _gap_cost_worker(idxs)

            n_chunks = min(n, self.num_workers * 4)
            chunks = [idxs[k::n_chunks] for k in range(n_chunks)]
            chunks = [c for c in chunks if c]
            ctx = mp.get_context('fork')
            with ctx.Pool(self.num_workers) as pool:
                results = pool.map(_gap_cost_worker, chunks)
        finally:
            _WORKER_GAP = {}  # release references in the parent

        out = []
        for r in results:
            out.extend(r)
        return out


# Populated in the parent process before a fork-based pool is created; worker
# processes read the two local-graph lists from here via copy-on-write.
_WORKER_TOPO = {}


def _topo_worker(chunk):
    """Compute topology scores for one chunk of candidate pairs (forked worker)."""
    ci, cj = chunk
    gm, gn, depth = _WORKER_TOPO['m'], _WORKER_TOPO['n'], _WORKER_TOPO['depth']
    out = np.empty(len(ci), dtype=float)
    for t in range(len(ci)):
        out[t] = _local_graph_comparison_score(depth, int(ci[t]), int(cj[t]), gm, gn)
    return out


# Populated before a fork-based gap-closing pool is created; workers read the
# track data + local-graph lists from here via copy-on-write.
_WORKER_GAP = {}


def _gap_cost_worker(i_chunk):
    """Compute gap-closing costs for a chunk of source-track indices (forked
    worker), using the per-start-frame KD-trees to find candidate targets within
    the safe radius. Returns a list of (i, j, cost) for pairs passing the cutoffs."""
    G = _WORKER_GAP
    all_track_disps, lsg = G['disps'], G['lsg']
    max_gap, depth = G['max_gap'], G['depth']
    dist_exp, top_exp = G['dist_exp'], G['top_exp']
    end_frame, end_coord, end_node = G['end_frame'], G['end_coord'], G['end_node']
    start_node, start_coord = G['start_node'], G['start_coord']
    trees, tree_idx, radius = G['trees'], G['tree_idx'], G['radius']
    tracking_interval = G['tracking_interval']

    out = []
    for i in i_chunk:
        end_frame_m = end_frame[i]
        end_coord_m = end_coord[i]
        end_node_m = end_node[i]
        disps_m = all_track_disps[i]
        simple_graphs_m = lsg[end_frame_m]

        # Tracked frames are spaced by tracking_interval, so a target separated
        # by gap_size missing tracked-frames starts at
        # end_frame_m + (gap_size + 1) * tracking_interval. (gap_size == 0 would
        # already have been linked frame-to-frame, so gap_size >= 1.)
        for gap_size in range(1, max_gap + 1):
            start_frame_n = end_frame_m + (gap_size + 1) * tracking_interval
            tree = trees.get(start_frame_n)
            if tree is None:
                continue
            members = tree_idx[start_frame_n]
            for k in tree.query_ball_point(end_coord_m, radius):
                j = members[k]

                disps_n = all_track_disps[j]
                comb_disps = disps_m + disps_n
                dist_cutoff = (gap_size + 1) * (3 * np.std(comb_disps)) ** 2

                dist = end_coord_m - start_coord[j]
                dist_cost = np.sum(dist ** 2)
                if dist_cost > dist_cutoff:
                    continue

                topology_cost = _local_graph_comparison_score(depth, end_node_m, start_node[j],
                                                              simple_graphs_m, lsg[start_frame_n])
                out.append((i, j, dist_cost ** dist_exp * topology_cost ** top_exp))
    return out


def _solve_lap_sparse(cand_i, cand_j, link_cost, m, n):
    """
    Solve the births/deaths linear assignment problem on a SPARSE augmented
    bipartite graph, returning an `assignment` array of length (m + n) with the
    same semantics as the dense formulation it replaces:

        assignment[i]   for i in [0, m): linked n-node index if < n, else terminated
        assignment[m+j] for j in [0, n): < n  => n-node j initiated

    Layout of the (m+n) x (n+m) augmented matrix (identical to the dense code):
        rows [0, m)      = m-nodes        cols [0, n)      = n-nodes
        rows [m, m+n)    = birth_j         cols [n, n+m)    = death_i

    Edges:
        m_i  -> n_j        link cost            (candidate pattern)
        m_i  -> death_i    3 * min link in row  (always present)
        birth_j -> n_j     3 * min link in col  (always present)
        birth_j -> death_i constant C           (candidate pattern; feasibility filler)

    The death/birth diagonals alone guarantee a perfect matching exists. The
    constant-C auxiliary block on the candidate pattern makes this equivalent to
    the dense constant-`nanmin` block: #aux edges used == #links and each costs C,
    and for every link (i, j) the edge (birth_j -> death_i) exists, so any dense-
    feasible link set is feasible here at identical cost.
    """
    ncand = len(cand_i)

    # Degenerate empty problem (no nodes in either frame): nothing to assign.
    # lap.lapmod raises on a 0-dim matrix, so return an empty assignment.
    if m + n == 0:
        return np.empty(0, dtype=np.int64)

    row_min = np.full(m, np.inf)
    col_min = np.full(n, np.inf)
    if ncand > 0:
        np.minimum.at(row_min, cand_i, link_cost)
        np.minimum.at(col_min, cand_j, link_cost)
        C = float(link_cost.min())
    else:
        C = 0.0

    death = np.where(np.isfinite(row_min), 3.0 * row_min, 0.0)
    birth = np.where(np.isfinite(col_min), 3.0 * col_min, 0.0)

    rows, cols, vals = [], [], []
    if ncand > 0:
        # link edges
        rows.append(cand_i);            cols.append(cand_j);            vals.append(link_cost)
        # auxiliary edges (constant C on candidate pattern, transposed)
        rows.append(m + cand_j);        cols.append(n + cand_i);        vals.append(np.full(ncand, C))
    # death diagonal
    rows.append(np.arange(m));          cols.append(n + np.arange(m));  vals.append(death)
    # birth diagonal
    rows.append(m + np.arange(n));      cols.append(np.arange(n));      vals.append(birth)

    rows = np.concatenate(rows).astype(np.int64)
    cols = np.concatenate(cols).astype(np.int64)
    # Add a tiny constant so no stored weight is exactly zero (scipy sparse treats
    # explicit zeros as absent edges). A uniform shift adds (m+n)*eps to every
    # perfect matching, so it never changes which assignment is optimal.
    vals = np.concatenate(vals).astype(float) + 1e-12

    dim = m + n
    cost = coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()

    # Solve with lap.lapmod (sparse Jonker-Volgenant), which scales to large
    # sparse assignment problems where a dense solver does not. The death/birth
    # diagonals guarantee a complete (perfect-matching-feasible) sparse cost
    # matrix, which lapmod requires. lapmod returns x[i] = column assigned to
    # row i, which is exactly the `assignment` array semantics expected here.
    _, x, _ = lap.lapmod(dim,
                         cost.data.astype(np.float64),
                         cost.indptr.astype(np.int32),
                         cost.indices.astype(np.int32))
    return x.astype(np.int64)


def _solve_gap_lap_sparse(entries, partition_start, partition_end):
    """
    Sparse equivalent of the dense gap-closing assignment. The original built a
    dense [P, 2P] matrix (link block | per-track termination diagonal) and solved
    it with scipy.linear_sum_assignment; this reproduces the same optimum on a
    sparse augmented [2P, 2P] matrix solved by lap.lapmod, avoiding the dense
    O(P^2) allocation (P = number of tracks in the block).

    `entries` is the list of finite (i, j, cost) link costs with GLOBAL track
    indices; partition_start/partition_end define the block. Returns the linked
    pairs as LOCAL (source, target) block indices.

    Augmented layout (rows x cols, both 2P):
        rows [0,P)   sources    cols [0,P)   targets
        rows [P,2P)  births     cols [P,2P)  deaths
    Edges:
        source_i -> target_j   link cost             (candidate pattern)
        source_i -> death_i    3 * min link in row   (termination; always present)
        birth_j  -> target_j   0                     (free, unlinked target)
        birth_j  -> death_i    0                     (candidate pattern; filler)
    The death/birth diagonals guarantee a feasible perfect matching. With zero
    birth and auxiliary costs the augmented optimum equals the original
    rectangular optimum: every perfect matching uses exactly (#links) auxiliary
    edges, each costing 0, so the two objective values are identical.
    """
    P = partition_end - partition_start
    if P == 0:
        return []

    li = np.array([e[0] - partition_start for e in entries], dtype=np.int64)
    lj = np.array([e[1] - partition_start for e in entries], dtype=np.int64)
    lc = np.array([e[2] for e in entries], dtype=float)
    ncand = len(li)

    # termination cost per source = 3 * min link cost in its row (0 if no links)
    row_min = np.full(P, np.inf)
    if ncand:
        np.minimum.at(row_min, li, lc)
    death = np.where(np.isfinite(row_min), 3.0 * row_min, 0.0)

    rows, cols, vals = [], [], []
    if ncand:
        rows.append(li);               cols.append(lj);               vals.append(lc)
        rows.append(P + lj);           cols.append(P + li);           vals.append(np.zeros(ncand))
    rows.append(np.arange(P));         cols.append(P + np.arange(P)); vals.append(death)
    rows.append(P + np.arange(P));     cols.append(np.arange(P));     vals.append(np.zeros(P))

    rows = np.concatenate(rows).astype(np.int64)
    cols = np.concatenate(cols).astype(np.int64)
    # Uniform tiny shift so stored zeros survive COO->CSR (which drops explicit
    # zeros); a constant added to every weight leaves the optimal matching fixed.
    vals = np.concatenate(vals).astype(float) + 1e-12

    dim = 2 * P
    cost = coo_matrix((vals, (rows, cols)), shape=(dim, dim)).tocsr()
    _, x, _ = lap.lapmod(dim,
                         cost.data.astype(np.float64),
                         cost.indptr.astype(np.int32),
                         cost.indices.astype(np.int32))
    return [(i, int(x[i])) for i in range(P) if x[i] < P]


class _SparseDist:
    """
    Sparse stand-in for the dense (m x n) distance matrix. Stores the
    pseudo-counted distance for each candidate pair and returns NaN for
    non-candidate pairs, supporting the scalar `dist_cost_mat[i, j]` access used
    downstream during segment-based filtering and arrow correction.
    """
    __slots__ = ("_d", "shape")

    def __init__(self, cand_i, cand_j, values, m, n):
        self._d = {(int(i), int(j)): float(v)
                   for i, j, v in zip(cand_i.tolist(), cand_j.tolist(), values.tolist())}
        self.shape = (m, n)

    def __getitem__(self, key):
        i, j = key
        return self._d.get((int(i), int(j)), np.nan)


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
    
    # Early exit for two lone nodes with zero edges to avoid LAP solver on empty matrix
    if num_edges == 0:
        return 0.0
        
    cost_mat = np.zeros((num_edges, num_edges))

    # fill cost matrix
    for i in range(num_edges):
        for j in range(num_edges):
            # FIX: Handle the case where both edge lengths are exactly 0
            max_len = max(sel_edge_len_m[i], sel_edge_len_n[j])
            
            if max_len == 0:
                cost_mat[i, j] = 0.0
            else:
                cost_mat[i, j] = abs(sel_edge_len_m[i] - sel_edge_len_n[j]) / max_len 

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
                    # check if the edge still exists
                    eid = frag_m.get_eid(node_m, neigh, error=False)
                    if eid != -1:
                        dist = frag_m.es[eid]['distance']
                        frag_m.delete_edges(eid)
                        frag_m.add_vertices(2)
                        frag_m.add_edges([[node_m, frag_m.vs[-2].index]], {'distance': dist})
                        frag_m.add_edges([[neigh, frag_m.vs[-1].index]], {'distance': dist})

            for neigh in neighbors_n:
                if neigh in visited_nodes_n and neigh not in last_level_n:
                    # check if the edge still exists
                    eid = frag_n.get_eid(node_n, neigh, error=False)
                    if eid != -1:
                        dist = frag_n.es[eid]['distance']
                        frag_n.delete_edges(eid)
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

            # Update neighbor list to include pseudo-nodes
            neighbors_m = frag_m.neighbors(node_m)
            neighbors_n = frag_n.neighbors(node_n)

            # Use list comprehensions to exclude parents. 
            # Never use .remove() inside a loop as it skips elements.
            neighbors_m = [n for n in neighbors_m if n not in last_level_m]
            neighbors_n = [n for n in neighbors_n if n not in last_level_n]

            # Re-map indices to the filtered neighbors
            index_mapping_m = {idx: val for idx, val in enumerate(neighbors_m)}
            index_mapping_n = {idx: val for idx, val in enumerate(neighbors_n)}

            num_node = max(len(neighbors_m), len(neighbors_n))
            cost_mat = np.zeros((num_node, num_node))

            # fill cost matrix with dissimilarity scores
            for i in range(num_node):
                for j in range(num_node):
                    # Add bounds checking for mismatched neighbor counts.
                    # If an index is missing, treat it as a pseudo-edge with 0 distance.
                    if i in index_mapping_m:
                        sel_edge_len_m = frag_m.es[frag_m.incident(index_mapping_m[i])]['distance']
                    else:
                        sel_edge_len_m = [0.0]

                    if j in index_mapping_n:
                        sel_edge_len_n = frag_n.es[frag_n.incident(index_mapping_n[j])]['distance']
                    else:
                        sel_edge_len_n = [0.0]

                    cost_mat[i, j] = _dissimilarity_score(sel_edge_len_m, sel_edge_len_n) 

            # run LAP solver
            nc = lap_solver(cost_mat)[1]

            # get node correspondence between local graphs m and n
            for a in range(len(nc)):
                if a in index_mapping_m and nc[a] in index_mapping_n:
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
