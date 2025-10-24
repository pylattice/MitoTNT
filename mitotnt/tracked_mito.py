import numpy as np
import pandas as pd
import igraph as ig
from scipy.linalg import lstsq
from tqdm.auto import tqdm, trange
from mitotnt.skeletonized_mito import SkeletonizedMito

class TrackedMito:

    """
    This class stores the results of mitochondrial tracking based on
    skeletonized network representations.It also adds features such as diffusivity and remodeling events.

    Parameters
    ----------
    segmented_mito : SkeletonizedMito
        Pre-processed mitochondrial data with skeletonized structures.
    frame_interval : float
        Time between consecutive frames, in seconds.
    start_frame : int
        Index of the first frame used for tracking.
    end_frame : int
        Index of the last frame used for tracking.
    tracking_interval : int
        The frame number difference for the two frames to be tracked.
    node_tracks : pandas.DataFrame
        Table of tracked nodes across frames.
    linked_nodes : numpy.ndarray
        Array describing frame-to-frame node linkages produced by the
        tracking algorithm.

    Attributes
    ----------
    data_path : str
        Path to the directory containing the original segmented data.
    save_path : str
        Path where tracking results are stored.
    list_of_folders : list of str
        Frame-level segmentation folders.
    num_frames : int
        Total number of frames in the dataset.
    xy_pxl_size : float
        Pixel resolution in the x–y plane (microns).
    z_pxl_size : float
        Pixel resolution in the z axis (microns).
    full_graphs : list
        Frame-wise full skeleton graphs of mitochondria.
    segment_nodes : list
        Nodes of individual mitochondrial segments.
    simple_graphs : list
        Simplified skeleton graphs for each frame.
    local_simple_graphs : list
        Simplified skeleton graphs restricted to local neighborhoods for each node.
    start_frame : int
        First frame index used for tracking.
    end_frame : int
        Last frame index used for tracking.
    frame_interval : float
        Time between frames.
    tracking_interval : int
        Step size in frames for linking nodes.
    node_tracks : pandas.DataFrame
        Computed from `NetworkTracker.run()`.
        Tabular data with each tracked node at one frame as a row, and columns
        include 'frame_id', 'frame_node_id', 'unique_node_id', 'frame_seg_id',
        'frame_frag_id', 'connected_unique_node_id', 'x', 'y', 'z', 'intensity',
        'width'.
    linked_nodes : numpy.ndarray
        Computed from `NetworkTracker.run()`.
        Array of node linkages between frames.
    remodeling_events : dict or None
        Placeholder for detected mitochondrial remodeling events (e.g.,
        fusion, fission). None until computed.
    node_diffusivity : pandas.DataFrame or None
        Placeholder for diffusivity metrics at the node level. None until
        computed.
    segment_diffusivity : pandas.DataFrame or None
        Placeholder for diffusivity metrics at the segment level. None until
        computed.
    fragment_diffusivity : pandas.DataFrame or None
        Placeholder for diffusivity metrics at the fragment level. None until
        computed.
    """

    def __init__(self, segmented_mito:SkeletonizedMito, frame_interval: float,
                 start_frame: int, end_frame: int, tracking_interval: int,
                 node_tracks: pd.DataFrame, linked_nodes: np.array):

        # copy from SegmentedMito object
        self.data_path = segmented_mito.data_path
        self.save_path = segmented_mito.save_path
        self.list_of_folders = segmented_mito.list_of_folders
        self.num_frames = segmented_mito.num_frames
        self.xy_pxl_size = segmented_mito.xy_pxl_size
        self.z_pxl_size = segmented_mito.z_pxl_size
        self.full_graphs = segmented_mito.full_graphs
        self.segment_nodes = segmented_mito.segment_nodes
        self.simple_graphs = segmented_mito.simple_graphs
        self.local_simple_graphs = segmented_mito.local_simple_graphs

        self.start_frame = start_frame
        self.end_frame = end_frame
        self.frame_interval = frame_interval
        self.tracking_interval = tracking_interval
        self.node_tracks = node_tracks
        self.linked_nodes = linked_nodes

        self.remodeling_events = None
        self.node_diffusivity = None
        self.segment_diffusivity = None
        self.fragment_diffusivity = None


    def detect_fusion_fission(self, start_frame: int = None, end_frame: int = None, tracking_interval: int = None,
                              half_win_size: int = 4, min_tracked_frames: int = 2):
        """
        Detect mitochondrial fusion and fission events.

        Scans tracked mitochondrial nodes across frames to identify topological
        remodeling events (fusion and fission). A sliding temporal window is
        applied to ensure that events are stable over multiple frames and not
        due to transient noise.

        Parameters
        ----------
        start_frame : int, optional
            First frame index to analyze. If None, defaults to previous
            `start_frame` + `half_win_size`.
        end_frame : int, optional
            Last frame index to analyze. If None, defaults to previous
            `end_frame` - `half_win_size`.
        tracking_interval : int, optional
            Frame step size for event detection. If None, defaults to the object's
            `tracking_interval`.
        half_win_size : int, optional
            Half-width of the temporal window (in frames) used to confirm
            events (default is 4).
        min_tracked_frames : int, optional
            Minimum number of consecutive frames a node must be tracked
            to be considered in event detection (default is 2).

        Returns
        -------
        None
            Results updated to `remodeling_events` attribute as a
            pandas DataFrame and written to `remodeling_events.csv`.
            Each row is a detected fusion or fission event, with the following columns:
            'type', 'frame_id', 'frame_id_before', 'frame_id_after',
            'node_id_before', 'node_id_after', 'frag_id_before','frag_id_after', 'unique_node_id'.
        """

        if start_frame is None:
            start_frame = self.start_frame + half_win_size
            
        if end_frame is None:
            end_frame = self.end_frame - half_win_size

        if tracking_interval is None:
            tracking_interval = self.tracking_interval

        full_graphs_all_frames = self.full_graphs

        tracks = self.node_tracks
        tracks.loc[:,'frame_id'] = tracks.loc[:,'frame_id'].astype(int)
        tracks.loc[:,'frame_frag_id'] = tracks.loc[:,'frame_frag_id'].astype(int)

        if start_frame < half_win_size:
            raise Exception(f"start_frame must be >= half_win_size but start_frame is {start_frame} and half_win_size is {half_win_size}")
        if end_frame > self.end_frame - half_win_size:
            raise Exception(f"end_frame must be <= self.num_frames - half_win_size but end_frame is {end_frame} and self.num_frames - half_win_size is {self.num_frames-half_win_size}")

        all_fragments = []
        for frame in range(0, self.num_frames):
            full_graph = full_graphs_all_frames[frame]
            frags = full_graph.components()
            all_fragments.append(frags)

        event_list = []
        node_list = []

        print('Performing fusion and fission event detection ...')
        for current_frame in trange(start_frame, end_frame, tracking_interval, desc="Detecting events per sliding window"):

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
                        frag_list1, frag_list2 = _remove_untracked_entries(forward_track_frag_list[index_node], forward_track_frag_list[index_neigh])
                        this_edge['forward_jaccard'] = _overlap_score(frag_list1, frag_list2, min_frames=min_tracked_frames)

                        this_edge = frag_track_graph.es[-1]
                        this_edge['backward_node_list'] = (backward_track_node_list[index_node], backward_track_node_list[index_neigh])
                        this_edge['backward_frag_list'] = (backward_track_frag_list[index_node], backward_track_frag_list[index_neigh])
                        frag_list1, frag_list2 = _remove_untracked_entries(backward_track_frag_list[index_node], backward_track_frag_list[index_neigh])
                        this_edge['backward_jaccard'] = _overlap_score(frag_list1, frag_list2, min_frames=min_tracked_frames)

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
                fission_clusters = _find_site_nodes(frag_track_graph, fission_edges)
                fusion_clusters = _find_site_nodes(frag_track_graph, fusion_edges)

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
                                       'frame_id_before':_list_to_str([current_frame] * len(fission_next_frame_id[i])),
                                       'frame_id_after':_list_to_str(fission_next_frame_id[i]),
                                       'node_id_before':_list_to_str(fission_clusters[i]),
                                       'node_id_after':_list_to_str(fission_clusters_next_frame[i]),
                                       'frag_id_before':_list_to_str([frag_id] * len(fission_frags_next_frame[i])),
                                       'frag_id_after':_list_to_str(fission_frags_next_frame[i]),
                                       'unique_node_id':_list_to_str([frame_to_unique[n] for n in fission_clusters[i]])})
                # FUSION
                for i in range(len(fusion_clusters)):
                    event_list.append({'type':'fusion',
                                       'frame_id': current_frame-1,
                                       'frame_id_before':_list_to_str(fusion_last_frame_id[i]),
                                       'frame_id_after':_list_to_str([current_frame] * len(fusion_last_frame_id[i])),
                                       'node_id_before':_list_to_str(fusion_clusters_last_frame[i]),
                                       'node_id_after':_list_to_str(fusion_clusters[i]),
                                       'frag_id_before':_list_to_str(fusion_frags_last_frame[i]),
                                       'frag_id_after':_list_to_str([frag_id] * len(fusion_frags_last_frame[i])),
                                       'unique_node_id':_list_to_str([frame_to_unique[n] for n in fusion_clusters[i]])})

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
        event_list.to_csv(self.save_path+'remodeling_events.csv', index=False)
        self.remodeling_events = event_list

        print("Fusion/fission events recorded.")


    def compute_mitochondria_diffusivity(self, levels : list = None, max_tau: int = 5, tracked_ratio: float = 0.3, half_win_size: int = 10, selected_frames: list = None):
        """
        Wrapper to compute mitochondrial diffusivity at one or more structural levels:
        nodes, segments, or fragments.

        Parameters
        ----------
        levels : list of strings, optional
            Structural levels at which to compute diffusivity. Choose from any combination of "node", "segment", "fragment"
            (default is ["node", "segment", "fragment"]).
        max_tau : int, optional
            Maximum time lag (in frames) for MSD calculation (default=5).
        tracked_ratio : float, optional
            Minimum fraction of nodes tracked at a given lag (default=0.3).
            Used only for segment and fragment diffusivity.
        half_win_size : int, optional
            Half-width of the temporal window (in frames) around each center frame used for displacement calculations (default=10).
            Used only for segment and fragment diffusivity.
        selected_frames : list of int, optional
            Specific frame indices to use as centers for analysis. If None, frames are auto-selected at intervals of
            ``2 * half_win_size`` across the dataset. Used only for segment and fragment diffusivity.

        Returns
        -------
        None
            Results are stored in class attributes and written to CSV files.

            For node diffusivity, each row is a tracked node with the following columns:
                - ``unique_node_id`` (int): unique node track ID from TrackedMito.node_tracks
                - ``diffusivity`` (float): Estimated diffusivity coefficient D.
                - ``msd`` (float): Fitted mean squared displacement per frame.
                - ``r_squared`` (float): Goodness of fit (R²) for the linear regression.
                - ``num_points`` (int): Number of time lags included in the fit.

            For segment/fragment diffusivity, each row is a tracked segment/fragment with the following columns:
                - ``center_frame_id`` (int): Center frame used for the analysis.
                - ``seg_id/frag_id`` (int): Segment identifier within the frame.
                - ``diffusivity`` (float): Estimated diffusivity coefficient D.
                - ``msd`` (float): Fitted mean squared displacement per frame.
                - ``r_squared`` (float): Goodness of fit (R²) for the linear regression.
                - ``num_points`` (int): Number of time lags included in the fit.
        """

        if levels is None:
            levels = ["node", "segment", "fragment"]
        levels = [lvl.lower() for lvl in levels]

        for level in levels:
            if level == "node":
                self.compute_node_diffusivity(max_tau=max_tau)

            elif level == "segment":
                self.compute_segment_diffusivity(
                    max_tau=max_tau,
                    tracked_ratio=tracked_ratio,
                    half_win_size=half_win_size,
                    selected_frames=selected_frames,
                )

            elif level == "fragment":
                self.compute_fragment_diffusivity(
                    max_tau=max_tau,
                    tracked_ratio=tracked_ratio,
                    half_win_size=half_win_size,
                    selected_frames=selected_frames,
                )

            else:
                raise ValueError("levels must be a list containing 'node', 'segment', or 'fragment'")


    def compute_node_diffusivity(self, max_tau: int = 5):
        """
        Compute node-level diffusivity from tracked mitochondrial trajectories.

        For each tracked node, the time-averaged mean squared displacement (TA-MSD)
        is calculated over time lags up to `max_tau`. A linear least-squares fit of
        MSD vs. lag time is then used to estimate the diffusivity coefficient (D),
        assuming 3D diffusion (`MSD = 6Dt`).

        Parameters
        ----------
        max_tau : int, optional
            Maximum time lag (in frames) used to compute TA-MSD (default is 5).

        Returns
        -------
        None
            Results are stored in the attribute `self.node_diffusivity` as a
            pandas DataFrame and written to `node_diffusivity.csv`. Each row is
            a tracked node with the following columns:

            - ``unique_node_id`` (int): Unique node track ID from TrackedMito.node_tracks
            - ``diffusivity`` (float): Estimated diffusivity coefficient D.
            - ``msd`` (float): Fitted mean squared displacement per frame.
            - ``r_squared`` (float): Goodness of fit (R²) for the linear regression.
            - ``num_points`` (int): Number of time lags included in the fit.
        """

        tracks = self.node_tracks
        num_tracks = int(np.max(tracks['unique_node_id'])) + 1

        all_msd = []
        for track_id in trange(num_tracks, desc="Computing node diffusivity for all tracks"):
            track = tracks[tracks['unique_node_id']==track_id]
            frames = track['frame_id'].to_numpy()
            coords = track.loc[:,'x':'z'].to_numpy()
            coords = coords

            # Calculate TA-MSD
            track_msd = []
            for tau in range(1, max_tau):

                disp = []
                frame = 0
                next_frame = frame + tau

                while next_frame < frames[-1]:
                    if next_frame in frames and frame in frames:
                        node = frames.tolist().index(frame)
                        next_node = frames.tolist().index(next_frame)

                        disp.append(np.sum((coords[next_node]-coords[node])**2))

                    frame += 1 # next start frame
                    next_frame = frame + tau

                if len(disp) < 2:
                    break
                else:
                    track_msd.append(np.mean(disp))

            all_msd.append(track_msd)

        msd_matrix = np.zeros((num_tracks, max_tau))
        for track_id, track_msd in enumerate(all_msd):
            for i in range(max_tau):
                if i < len(track_msd):
                    msd_matrix[track_id,i] = track_msd[i]
                else:
                    msd_matrix[track_id,i] = np.nan

        node_diffusivity = []
        for n in range(num_tracks):
            eata_msd = [msd for msd in msd_matrix[n] if not np.isnan(msd)]
            eata_msd.insert(0,0)

            # choose the number of data points to fit
            if len(eata_msd) < 2:
                node_diffusivity.append({'unique_node_id':n, 'diffusivity':np.nan, 'msd':np.nan, 'r_squared':np.nan, 'num_points':1})
                continue
            elif len(eata_msd) > max_tau:
                n_points = max_tau
            else:
                n_points = len(eata_msd)

            all_tau = np.arange(0, n_points * self.frame_interval, self.frame_interval)[:,np.newaxis]
            slope, res, rnk, s = lstsq(all_tau[:n_points], eata_msd[:n_points])
            d = slope[0] / 6
            msd_per_frame = 6 * d * self.frame_interval

            # get r^2
            msd_mean = np.mean(eata_msd[:n_points])
            total_sum = np.sum([(i-msd_mean)**2 for i in eata_msd[:n_points]])

            if total_sum == 0:
                r_squared = np.nan
            else:
                r_squared = 1 - res/total_sum

            # store data
            node_diffusivity.append({'unique_node_id':n, 'diffusivity':d, 'msd':msd_per_frame, 'r_squared':r_squared, 'num_points':n_points})

        node_diffusivity = pd.DataFrame.from_dict(node_diffusivity)
        node_diffusivity.to_csv(self.save_path+'node_diffusivity.csv', index=False)
        self.node_diffusivity = node_diffusivity


    def compute_segment_diffusivity(self, max_tau: int = 5, tracked_ratio: float = 0.3, half_win_size: int = 10, selected_frames: list = None):
        """
        Compute segment-level diffusivity from tracked mitochondrial trajectories.

        For each segment in selected center frames, this method aggregates node-level
        displacements into segment-averaged displacements, computes time-averaged
        mean squared displacement (TA-MSD), and estimates diffusivity coefficients
        using a linear least-squares fit of MSD vs. lag time.

        Parameters
        ----------
        max_tau : int, optional
            Maximum time lag (in frames) for MSD calculation (default is 5).
        tracked_ratio : float, optional
            Minimum fraction of segment nodes that must be tracked at a given lag
            for the segment displacement to be considered (default is 0.3).
        half_win_size : int, optional
            Half-width of the temporal window (in frames) around each center frame
            used for displacement calculations (default is 10).
        selected_frames : list of int, optional
            Specific frame indices to use as centers for segment diffusivity
            analysis. If None, frames are auto-selected at intervals of
            ``2 * half_win_size`` across the dataset.

        Returns
        -------
        None
            Results are stored in the attribute `self.segment_diffusivity` as a
            pandas DataFrame and written to `segment_diffusivity.csv`.
            Each row is a tracked segment with the following columns:

                - ``center_frame_id`` (int): Center frame used for the analysis.
                - ``seg_id`` (int): Segment identifier within the frame.
                - ``diffusivity`` (float): Estimated diffusivity coefficient D.
                - ``msd`` (float): Fitted mean squared displacement per frame.
                - ``r_squared`` (float): Goodness of fit (R²) for the linear regression.
                - ``num_points`` (int): Number of time lags included in the fit.
        """

        if selected_frames is None:
            selected_frames = []
            num_selected_frames = int(self.end_frame / (2 * half_win_size))
            for i in range(num_selected_frames):
                selected_frames.append(half_win_size + 2 * i * half_win_size)
        else:
            for i in selected_frames:
                if i not in range(self.end_frame):
                    raise Exception('Selected frame index out of range?')

        # load inputs
        tracks = self.node_tracks
        segment_nodes_all_frames = self.segment_nodes

        seg_diffusivity = []

        # iterate the center frames
        for center_frame in tqdm(selected_frames, desc="Computing segment diffusivity around center frames"):
            segment_nodes = segment_nodes_all_frames[center_frame]
            num_segs = len(segment_nodes)

            frame_tracks = tracks[tracks.frame_id==center_frame] # find only tracks that intersects with the center frame
            unique_nodes = np.array(frame_tracks['unique_node_id'].tolist(), dtype=int)
            frame_nodes = np.array(frame_tracks['frame_node_id'].tolist(), dtype=int)
            frame_to_unique = {frame_nodes[i]:unique_nodes[i] for i in range(len(frame_nodes))}
            unique_to_frame = {unique_nodes[i]:frame_nodes[i] for i in range(len(unique_nodes))}

            all_msd, all_tau = [], []
            for seg_id in range(len(segment_nodes)):

                segment = segment_nodes[seg_id]
                seg_coords = np.zeros((len(segment), 2*half_win_size), dtype=object)
                seg_coords.fill(np.array([np.nan, np.nan, np.nan], dtype=object))

                for node_id, frame_node in enumerate(segment):
                    if frame_node in frame_to_unique.keys():
                        full_track = tracks[tracks.unique_node_id==frame_to_unique[frame_node]]

                        for i in range(len(full_track)):
                            frame = int(full_track.iloc[i]['frame_id'])
                            coord = full_track.iloc[i]['x':'z'].to_numpy()
                            arr_index = frame - (center_frame - half_win_size)

                            if arr_index < 2*half_win_size and arr_index >= 0:
                                seg_coords[node_id,arr_index] = coord
                            else:
                                break

                # Calculate TA-MSD
                seg_msd, seg_tau = [], []
                for tau in range(1, max_tau):

                    disp = []
                    frame = 0
                    next_frame = frame + tau

                    while next_frame < 2 * half_win_size:
                        coords_m = seg_coords[:,frame]
                        coords_n = seg_coords[:,next_frame]
                        vector_diff = coords_n - coords_m # get displacement vector for each node

                        disp_vectors = np.array([i for i in vector_diff if not pd.isnull(i).any()]) # collect the valid vectors that belong to tracked nodes
                        num_disp = disp_vectors.shape[0]

                        if num_disp >= tracked_ratio * len(segment): # if enough of the segment is tracked
                            # average all displacement vector to get segment vector
                            average_disp = (np.linalg.norm(np.mean(disp_vectors, axis=0))) ** 2
                            disp.append(average_disp)

                        frame += 1 # next start frame
                        next_frame = frame + tau

                    if len(disp) < 2:
                        break
                    else:
                        seg_msd.append(np.mean(disp))
                        seg_tau.append(tau)

                all_msd.append(seg_msd)
                all_tau.append(seg_tau)

            # print('Percent of segments tracked around center frame', center_frame, ':',
            #       (1 - np.sum([len(msd)==0 for msd in all_msd]) / len(all_msd)) * 100, '%')

            # compute diffusivity from MSD
            num_msd = len(all_msd)
            msd_matrix = np.zeros((num_msd, max_tau))
            for track_id, track_msd in enumerate(all_msd):
                for i in range(max_tau):
                    if i < len(track_msd):
                        msd_matrix[track_id,i] = track_msd[i]
                    else:
                        msd_matrix[track_id,i] = np.nan

            for seg_id in range(num_segs):
                eata_msd = [msd for msd in msd_matrix[seg_id] if not np.isnan(msd)]
                eata_msd.insert(0,0)

                # choose the number of data points to fit
                if len(eata_msd) <= 1:
                    seg_diffusivity.append({'center_frame_id': center_frame, 'seg_id':seg_id, 'diffusivity':np.nan, 'msd':np.nan, 'r_squared':np.nan, 'num_points':1})
                    continue
                elif len(eata_msd) > max_tau:
                    n_points = max_tau
                else:
                    n_points = len(eata_msd)

                all_tau = np.arange(0, n_points * self.frame_interval, self.frame_interval)[:,np.newaxis]
                slope, res, rnk, s = lstsq(all_tau[:n_points], eata_msd[:n_points])
                d = slope[0] / 6
                msd_per_frame = 6 * d * self.frame_interval

                # get r^2
                msd_mean = np.mean(eata_msd[:n_points])
                total_sum = np.sum([(i-msd_mean)**2 for i in eata_msd[:n_points]])

                if total_sum == 0:
                    r_squared = np.nan
                else:
                    r_squared = 1 - res/total_sum

                seg_diffusivity.append({'center_frame_id': center_frame, 'seg_id':seg_id, 'diffusivity':d, 'msd':msd_per_frame, 'r_squared':r_squared, 'num_points':n_points})

        seg_diffusivity = pd.DataFrame.from_dict(seg_diffusivity)
        seg_diffusivity.to_csv(self.save_path+'segment_diffusivity.csv', index=False)
        self.segment_diffusivity = seg_diffusivity


    def compute_fragment_diffusivity(self, max_tau: int = 5, tracked_ratio: float = 0.5, half_win_size: int = 10, selected_frames: list = None):
        """
        Compute fragment-level diffusivity from tracked mitochondrial trajectories.

        For each fragment in selected center frames, this method aggregates node-level
        displacements into fragment-averaged displacements, computes time-averaged
        mean squared displacement (TA-MSD), and estimates diffusivity coefficients
        using a linear least-squares fit of MSD vs. lag time.

        Parameters
        ----------
        max_tau : int, optional
            Maximum time lag (in frames) for MSD calculation (default is 5).
        tracked_ratio : float, optional
            Minimum fraction of fragment nodes that must be tracked at a given lag
            for the fragment displacement to be considered (default is 0.3).
        half_win_size : int, optional
            Half-width of the temporal window (in frames) around each center frame
            used for displacement calculations (default is 10).
        selected_frames : list of int, optional
            Specific frame indices to use as centers for fragment diffusivity
            analysis. If None, frames are auto-selected at intervals of
            ``2 * half_win_size`` across the dataset.

        Returns
        -------
        None
            Results are stored in the attribute `self.fragment_diffusivity` as a
            pandas DataFrame and written to `fragment_diffusivity.csv`.
            Each row is a tracked fragment with the following columns:

                - ``center_frame_id`` (int): Center frame used for the analysis.
                - ``seg_id`` (int): Fragment identifier within the frame.
                - ``diffusivity`` (float): Estimated diffusivity coefficient D.
                - ``msd`` (float): Fitted mean squared displacement per frame.
                - ``r_squared`` (float): Goodness of fit (R²) for the linear regression.
                - ``num_points`` (int): Number of time lags included in the fit.
        """

        if selected_frames is None:
            selected_frames = []
            num_selected_frames = int(self.end_frame / (2 * half_win_size))
            for i in range(num_selected_frames):
                selected_frames.append(half_win_size + 2 * i * half_win_size)
        else:
            for i in selected_frames:
                if i not in range(self.end_frame):
                    raise Exception('Selected frame index out of range?')

        # load inputs
        tracks = self.node_tracks
        full_graphs_all_frames = self.full_graphs

        frag_diffusivity = []

        # iterate the center frames
        for center_frame in tqdm(selected_frames, desc="Computing fragment diffusivity around center frames"):

            graph = full_graphs_all_frames[center_frame]
            fragment_nodes = graph.components()
            num_frags = len(fragment_nodes)

            frame_tracks = tracks[tracks.frame_id == center_frame]  # find only tracks that intersects with the center frame
            unique_nodes = np.array(frame_tracks['unique_node_id'].tolist(), dtype=int)
            frame_nodes = np.array(frame_tracks['frame_node_id'].tolist(), dtype=int)
            frame_to_unique = {frame_nodes[i]: unique_nodes[i] for i in range(len(frame_nodes))}
            unique_to_frame = {unique_nodes[i]: frame_nodes[i] for i in range(len(unique_nodes))}

            all_msd, all_tau = [], []
            for frag_id in range(len(fragment_nodes)):

                fragment = fragment_nodes[frag_id]
                frag_coords = np.zeros((len(fragment), 2 * half_win_size), dtype=object)
                frag_coords.fill(np.array([np.nan, np.nan, np.nan], dtype=object))

                for node_id, frame_node in enumerate(fragment):
                    if frame_node in frame_to_unique.keys():
                        full_track = tracks[tracks.unique_node_id == frame_to_unique[frame_node]]

                        for i in range(len(full_track)):
                            frame = int(full_track.iloc[i]['frame_id'])
                            coord = full_track.iloc[i]['x':'z'].to_numpy()
                            arr_index = frame - (center_frame - half_win_size)

                            if arr_index < 2 * half_win_size and arr_index >= 0:
                                frag_coords[node_id, arr_index] = coord
                            else:
                                break

                # Calculate TA-MSD
                frag_msd, frag_tau = [], []
                for tau in range(1, max_tau):

                    disp = []
                    frame = 0
                    next_frame = frame + tau

                    while next_frame < 2 * half_win_size:
                        coords_m = frag_coords[:, frame]
                        coords_n = frag_coords[:, next_frame]
                        vector_diff = coords_n - coords_m  # get displacement vector for each node

                        disp_vectors = np.array([i for i in vector_diff if not pd.isnull(
                            i).any()])  # collect the valid vectors that belong to tracked nodes
                        num_disp = disp_vectors.shape[0]

                        if num_disp >= tracked_ratio * len(fragment):  # if enough of the fragment is tracked
                            # average all displacement vector to get fragment vector
                            average_disp = (np.linalg.norm(np.mean(disp_vectors, axis=0))) ** 2
                            disp.append(average_disp)

                        frame += 1  # next start frame
                        next_frame = frame + tau

                    if len(disp) < 2:
                        break
                    else:
                        frag_msd.append(np.mean(disp))
                        frag_tau.append(tau)

                all_msd.append(frag_msd)
                all_tau.append(frag_tau)

            # print('Percent of fragments tracked around center frame', center_frame, ':',
            #       (1 - np.sum([len(msd) == 0 for msd in all_msd]) / len(all_msd)) * 100, '%')

            # compute diffusivity from MSD
            num_msd = len(all_msd)
            msd_matrix = np.zeros((num_msd, max_tau))
            for track_id, track_msd in enumerate(all_msd):
                for i in range(max_tau):
                    if i < len(track_msd):
                        msd_matrix[track_id, i] = track_msd[i]
                    else:
                        msd_matrix[track_id, i] = np.nan

            for frag_id in range(num_frags):
                eata_msd = [msd for msd in msd_matrix[frag_id] if not np.isnan(msd)]
                eata_msd.insert(0, 0)

                # choose the number of data points to fit
                if len(eata_msd) <= 1:
                    frag_diffusivity.append(
                        {'center_frame_id': center_frame, 'frag_id': frag_id, 'diffusivity': np.nan, 'msd': np.nan,
                         'r_squared': np.nan, 'num_points': 1})
                    continue
                elif len(eata_msd) > max_tau:
                    n_points = max_tau
                else:
                    n_points = len(eata_msd)

                all_tau = np.arange(0, n_points * self.frame_interval, self.frame_interval)[:, np.newaxis]
                slope, res, rnk, s = lstsq(all_tau[:n_points], eata_msd[:n_points])
                d = slope[0] / 6
                msd_per_frame = 6 * d * self.frame_interval

                # get r^2
                msd_mean = np.mean(eata_msd[:n_points])
                total_sum = np.sum([(i - msd_mean) ** 2 for i in eata_msd[:n_points]])

                if total_sum == 0:
                    r_squared = np.nan
                else:
                    r_squared = 1 - res / total_sum

                frag_diffusivity.append(
                    {'center_frame_id': center_frame, 'frag_id': frag_id, 'diffusivity': d, 'msd': msd_per_frame,
                     'r_squared': r_squared, 'num_points': n_points})

        frag_diffusivity = pd.DataFrame.from_dict(frag_diffusivity)
        frag_diffusivity.to_csv(self.save_path + 'fragment_diffusivity.csv', index=False)
        self.fragment_diffusivity = frag_diffusivity


def _list_to_str(list1):
    string = ""
    if len(list1) == 0:
        return string

    for i in range(len(list1)-1):
        string += str(list1[i])
        string += " "
    string += str(list1[-1])

    return string


def _remove_untracked_entries(frag_list1, frag_list2):
    new_list1, new_list2 = [], []
    for i in range(np.min([len(frag_list1), len(frag_list2)])):
        if np.isnan(frag_list1[i]) or np.isnan(frag_list2[i]):
            continue
        else:
            new_list1.append(frag_list1[i])
            new_list2.append(frag_list2[i])
    return new_list1, new_list2


def _overlap_score(list1, list2, min_frames):
    if len(list1) != len(list2):
        raise Exception('Two lists should be the same size!')

    num_frames = len(list1)
    if num_frames >= min_frames:
        score = np.sum(np.array(list1)==np.array(list2))
        return 1 - score / num_frames

    else: # less than min_frames
        return np.nan


def _find_site_nodes(graph, event_edges, max_edge_gap=5):
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
