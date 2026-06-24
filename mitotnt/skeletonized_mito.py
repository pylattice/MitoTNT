import glob
import io
import os
import numpy as np
import warnings
import pandas as pd
import igraph as ig
import multiprocessing as mp
from functools import partial
from scipy.spatial import cKDTree
from tqdm.auto import tqdm, trange
from mitotnt.mito_segmenter import MitoSegmenter
from pathlib import Path

class SkeletonizedMito:
    """
    This class stores metadata and skeletonized graph representations
    of mitochondria extracted from raw segmentation results produced by
    `MitoSegmenter`. It provides access to frame-level information and
    graph structures that can later be used for tracking.

    Parameters
    ----------
    mito_segmenter : MitoSegmenter
        Object that contains metadata.

    Attributes
    ----------
    data_path : str
        Path to the directory containing the segmented mitochondrial data.
    list_of_folders : list of str
        List of frame-level folders containing segmentation outputs.
    num_frames : int
        Total number of frames in the dataset.
    xy_pxl_size : float
        Pixel size in the x–y plane (microns).
    z_pxl_size : float
        Pixel size in the z axis (microns).
    save_path : str
        Path to where processed mitochondrial data will be saved.
    full_graphs : list or None
        Frame-wise full graphs with all nodes. None until computed.
    segment_nodes : list or None
        Nodes for each mitochondrial segment (branch). None until computed.
    simple_graphs : list or None
        Frame-wise simple graphs (full_graphs where degree-2 nodes are substituted by edges). None until computed.
    local_simple_graphs : list or None
        Simple graphs restricted to the local neighborhood of each node. They are used for computing topology cost when tracking. None until computed.
    """

    def __init__(self, mito_segmenter: MitoSegmenter):

        self.data_path = mito_segmenter.data_path
        self.list_of_folders = mito_segmenter.list_of_folders
        self.num_frames = mito_segmenter.num_frames
        self.xy_pxl_size = mito_segmenter.xy_pxl_size
        self.z_pxl_size = mito_segmenter.z_pxl_size

        self.save_path = os.path.join(str(Path(self.data_path).parent), "mitotnt", "") # add trailing slash is critical
        self.full_graphs = None
        self.segment_nodes = None
        self.simple_graphs = None
        self.local_simple_graphs = None

        if self.num_frames == 0:
            raise ValueError(f"No matching folders in {self.data_path!r}")


    def extract_graphs(self, overwrite: bool = False, node_gap_size: int = 0,
                       num_workers: int = 20):
        """
        Extract the graph representations of the previously segmented mitochondria used for tracking.

        Parameters
        ----------
        overwrite : bool
            Whether to overwrite existing graphs (default is False).
        node_gap_size : int, optional
            Number of nodes to skip when generating the full graphs (default is 0, keep all).
        num_workers : int, optional
            Max processes used to extract frames in parallel (default 20, capped
            at the CPU count and the number of frames). Each extraction stage is
            independent per frame, so this fans out across frames.

        Returns
        -------
        None
            The computed graphs and segment nodes are stored internally as `SkeletonizedMito` attributes.
        """

        if os.path.exists(self.save_path+'extracted_graphs.npz') and not overwrite:
            self._load_graphs()
            print("Graphs have already been extracted. Reload previous data.\n")
        else:
            self.extract_full_graphs_and_segment_nodes(node_gap_size=node_gap_size,
                                                       num_workers=num_workers)
            self.extract_simple_graphs(num_workers=num_workers)
            self.extract_local_simple_graphs(num_workers=num_workers)
            self._save_graphs()


    def extract_full_graphs_and_segment_nodes(self, node_gap_size: int = 0,
                                              num_workers: int = 20):
        # Per-frame work is independent -> fan out across frames. Each worker
        # reads its frame's MitoGraph files and builds the full graph + segment
        # nodes; results come back in frame order.
        worker = partial(_extract_full_graph_one, node_gap_size=node_gap_size)
        results = _map_frames(worker, self.list_of_folders, num_workers,
                              "Extracting full graphs and segment nodes")

        all_full_graphs = [r[0] for r in results]
        all_segment_nodes = [r[1] for r in results]
        self.full_graphs = np.array(all_full_graphs, dtype=object)
        self.segment_nodes = np.array(all_segment_nodes, dtype=object)


    def extract_simple_graphs(self, num_workers: int = 20):
        results = _map_frames(_extract_simple_graph_one, self.list_of_folders,
                              num_workers, "Extracting simple graphs")
        self.simple_graphs = np.array(results, dtype=object)


    def extract_local_simple_graphs(self, num_workers: int = 20):

        if self.full_graphs is None:
            raise Exception("No full graphs to extract. Run extract_full_graphs_and_segment_nodes() first.")

        # Publish the (read-only) full graphs as a module global BEFORE forking
        # so workers inherit them via copy-on-write; only the frame index is sent
        # through the pool. induced_subgraph()/_contract_edges() operate on copies,
        # so the shared graphs are never mutated (fork-safe).
        global _WORKER_FULL_GRAPHS
        _WORKER_FULL_GRAPHS = list(self.full_graphs)
        try:
            results = _map_frames(_extract_local_simple_one,
                                  list(range(len(_WORKER_FULL_GRAPHS))),
                                  num_workers, "Extracting local simple graphs")
        finally:
            _WORKER_FULL_GRAPHS = []  # release references in the parent

        self.local_simple_graphs = np.array(results, dtype=object)
    
    
    def _save_graphs(self):
        # save all inputs as a single compressed .npz file
        os.makedirs(self.save_path, exist_ok=True)
        path = self.save_path + 'extracted_graphs.npz'
        data = {'full_graphs': self.full_graphs,
                'segment_nodes': self.segment_nodes,
                'simple_graphs': self.simple_graphs,
                'local_simple_graphs': self.local_simple_graphs}
        # Compressed: the per-node local_simple_graphs dominate file size and
        # compress well. np.load reads compressed npz transparently.
        np.savez_compressed(path, **data)
    
    
    def _load_graphs(self):
        try:
            data = np.load(self.save_path+'extracted_graphs.npz', allow_pickle=True)
        except:
            raise Exception(f"No extracted_graphs.npz found under {self.save_path!r}.")
    
        self.full_graphs = data['full_graphs']
        self.segment_nodes = data['segment_nodes']
        self.simple_graphs = data['simple_graphs']
        self.local_simple_graphs = data['local_simple_graphs']


# --------------------------------------------------------------------------- #
# Parallel per-frame extraction helpers                                         #
# --------------------------------------------------------------------------- #

# Populated in the parent process before a fork-based pool is created; worker
# processes read the full-graph list from here via copy-on-write.
_WORKER_FULL_GRAPHS = []


def _map_frames(worker, items, num_workers, desc):
    """
    Apply `worker` to each item (one per frame), preserving order. Runs in a
    fork-based process pool for real parallelism on CPU-bound igraph work
    (threads would be serialized by the GIL). Falls back to a serial loop for a
    single worker or a single frame.
    """
    n = len(items)
    nw = max(1, min(int(num_workers), os.cpu_count() or 1, n))
    if nw <= 1 or n <= 1:
        return [worker(it) for it in tqdm(items, desc=desc)]

    ctx = mp.get_context('fork')
    with ctx.Pool(nw) as pool:
        return list(tqdm(pool.imap(worker, items), total=n, desc=desc))


def _extract_full_graph_one(folder, node_gap_size=0):
    """Build the full graph + segment nodes for ONE frame folder.

    Construction is fully batched: vertices, vertex attributes and edges are
    accumulated in Python lists and committed to the igraph object in a single
    add_vertices / attribute-assignment / add_edges call each. This avoids the
    per-node add_vertices(1)/add_edge() pattern, which is O(n^2) in python-igraph
    (every incremental call rebuilds internal structures) and dominated runtime
    on large frames. The resulting graph is identical to the incremental build:
    same vertex ids/attributes, same edges (added in the same order), and the
    same per-segment distances (computed from the .txt coordinates, exactly the
    values the incremental version read back from the vertices it had just set)."""
    if len(glob.glob(folder + '/*.coo')) == 1 and len(glob.glob(folder + '/*.gnet')) == 1 and len(glob.glob(folder + '//*.txt')) == 1:
        simple_graph_coords = _read_array(glob.glob(folder + '/*.coo')[0])
        segment_node_data = _read_table(glob.glob(folder + '/*.txt')[0], delimiter='\t')
    else:
        raise Exception(f"{folder!r} has none/duplicate MitoGraph outputs.")

    simple_graph_coords = _round_coord(simple_graph_coords)

    # KD-tree for fast nearest-node lookup (replaces per-endpoint O(N)
    # linear search; identical nearest-neighbour / argmin semantics)
    coord_tree = cKDTree(simple_graph_coords)

    num_base = len(simple_graph_coords)

    # base (degree != 2) node attributes; endpoints get overwritten below
    coords = list(simple_graph_coords)
    intensity = [np.nan] * num_base
    width = [np.nan] * num_base

    line_ids = np.unique(segment_node_data['line_id'])
    frame_segment_nodes = []

    # group rows by line_id once (replaces per-line O(N) boolean filter)
    line_groups = {k: v.reset_index() for k, v in segment_node_data.groupby('line_id')}

    edges = []
    edge_dists = []
    next_index = num_base  # next bulk-node id to assign

    for line in line_ids:
        line_nodes = line_groups[line]
        end_index = len(line_nodes) - 1

        # per-line numpy views (positional indexing avoids slow pandas .loc)
        xyz = line_nodes.loc[:, 'x':'z'].to_numpy()
        pint = line_nodes['pixel_intensity'].to_numpy()
        pwid = line_nodes['width_(um)'].to_numpy()

        # find index of network nodes in .coo based on (rounded) end coords
        index_end_a = _coord_to_node(coord_tree, _round_coord(xyz[0]))
        index_end_b = _coord_to_node(coord_tree, _round_coord(xyz[end_index]))

        # overwrite endpoint attributes with the (unrounded) .txt values
        coords[index_end_a] = xyz[0]
        intensity[index_end_a] = pint[0]
        width[index_end_a] = pwid[0]
        coords[index_end_b] = xyz[end_index]
        intensity[index_end_b] = pint[end_index]
        width[index_end_b] = pwid[end_index]

        # walk the segment, emitting bulk nodes/edges; distances use the .txt
        # coords (== the coordinates the incremental version had just assigned)
        last_node = index_end_a
        last_coord = xyz[0]
        segment_nodes = [index_end_a]
        sel_node_index = range(0, end_index, 1 + node_gap_size)
        if len(sel_node_index) == 0:
            continue
        last_sel = sel_node_index[-1]

        for index in sel_node_index:
            # add bulk nodes
            if index > 0 and index < last_sel:
                current_node = next_index
                next_index += 1
                coords.append(xyz[index])
                intensity.append(pint[index])
                width.append(pwid[index])

                edges.append((last_node, current_node))
                edge_dists.append(np.linalg.norm(last_coord - xyz[index]))
                last_node = current_node
                last_coord = xyz[index]

                segment_nodes.append(current_node)

            # link last bulk node to the other network node
            if index == last_sel:
                edges.append((last_node, index_end_b))
                edge_dists.append(np.linalg.norm(last_coord - xyz[end_index]))

                segment_nodes.append(index_end_b)  # get the node index of another end
                frame_segment_nodes.append(segment_nodes)  # finish this segment

    # commit the graph in bulk
    full_graph = ig.Graph()
    full_graph.add_vertices(next_index)
    full_graph.vs['index'] = list(range(next_index))  # index == vertex id for all nodes
    full_graph.vs['coordinate'] = coords
    full_graph.vs['intensity'] = intensity
    full_graph.vs['width'] = width
    if edges:
        full_graph.add_edges(edges, {'distance': edge_dists})

    return full_graph, frame_segment_nodes


def _extract_simple_graph_one(folder):
    """Build the simple graph for ONE frame folder.

    Batched construction (see _extract_full_graph_one): node attributes and the
    full edge list are committed in one shot instead of per-edge add_edge()
    calls, which are O(E^2) in python-igraph. Output is identical."""
    if len(glob.glob(folder+'/*.coo')) == 1 and len(glob.glob(folder+'/*.gnet')) == 1 and len(glob.glob(folder+'//*.txt')) == 1:
        simple_graph_coords = _read_array(glob.glob(folder+'/*.coo')[0])
        edge_list = _read_array(glob.glob(folder+'/*.gnet')[0], skiprows=1)
        segment_node_data = _read_table(glob.glob(folder+'/*.txt')[0], delimiter='\t')
    else:
        raise Exception(f"{folder!r} has none/duplicate MitoGraph outputs.")

    simple_graph_coords = _round_coord(simple_graph_coords)

    # KD-tree for fast nearest-node lookup (see extract_full_graphs)
    coord_tree = cKDTree(simple_graph_coords)

    num_base = len(simple_graph_coords)
    coords = list(simple_graph_coords)
    intensity = [np.nan] * num_base
    width = [np.nan] * num_base

    # create all the network nodes
    line_ids = np.unique(segment_node_data['line_id'])

    # group rows by line_id once (replaces per-line O(N) boolean filter)
    line_groups = {k: v.reset_index() for k, v in segment_node_data.groupby('line_id')}

    for line in line_ids:
        line_nodes = line_groups[line]
        end_index = len(line_nodes) - 1

        xyz = line_nodes.loc[:, 'x':'z'].to_numpy()
        pint = line_nodes['pixel_intensity'].to_numpy()
        pwid = line_nodes['width_(um)'].to_numpy()

        # find index of network nodes in .coo based on (rounded) end coords
        index_end_a = _coord_to_node(coord_tree, _round_coord(xyz[0]))
        index_end_b = _coord_to_node(coord_tree, _round_coord(xyz[end_index]))

        coords[index_end_a] = xyz[0]
        intensity[index_end_a] = pint[0]
        width[index_end_a] = pwid[0]
        coords[index_end_b] = xyz[end_index]
        intensity[index_end_b] = pint[end_index]
        width[index_end_b] = pwid[end_index]

    simple_graph = ig.Graph()
    simple_graph.add_vertices(num_base)
    simple_graph.vs['index'] = list(range(num_base))
    simple_graph.vs['coordinate'] = coords
    simple_graph.vs['intensity'] = intensity
    simple_graph.vs['width'] = width

    if len(edge_list):
        edges = edge_list[:, :2].astype(int)
        simple_graph.add_edges([tuple(e) for e in edges], {'distance': edge_list[:, 2]})

    return simple_graph


def _extract_local_simple_one(frame_idx):
    """Build the per-node local simple graphs for ONE frame, reading the frame's
    full graph from the inherited _WORKER_FULL_GRAPHS global. Verbatim per-frame
    body of extract_local_simple_graphs()."""
    frame_full_graphs = _WORKER_FULL_GRAPHS[frame_idx]

    # load graphs of nodes and edges
    total_num_nodes = len(frame_full_graphs.vs)

    # get fragments
    all_frags = frame_full_graphs.components()

    # contract edges and update frag and return new root index
    frame_simple_graphs_per_node = []
    for node_index in range(total_num_nodes):
        # use full graph node id
        frag = frame_full_graphs.induced_subgraph(all_frags[all_frags.membership[node_index]])

        # find fragment graph node id
        root = frag.vs['index'].index(node_index)

        # contract edges to extract simple graph around the node
        simple_graph = _contract_edges(frag, root)
        frame_simple_graphs_per_node.append(simple_graph)

    return frame_simple_graphs_per_node


def _contract_edges(frag, root):

    segment_node_data = []
    edge_weights = []
    all_segments = []
    last_node = -1

    for i, node in enumerate(frag.dfsiter(root)):

        n = node.index
        degree = node.degree()

        if n != root:
            if degree == 2:

                # first node on a new segment after concluding a segment
                if last_node == -1:
                    segment_node_data.append(n)
                    last_node = n

                # sometimes in large graph a new segment is visited without reaching a degree!=2 node
                else:
                    # this may fail when reach end of one segment and jump to the start of another segment
                    try:
                        weight = frag.es[frag.get_eid(n, last_node)]['distance']

                    # when it fails just start another bulk
                    except Exception:
                        if len(segment_node_data) != 0:
                            all_segments.append([segment_node_data, sum(edge_weights)])
                            segment_node_data = []
                            edge_weights = []

                            segment_node_data.append(n)
                            last_node = n

                    # when the two nodes are on bulk we can just append distance and node
                    else:
                        edge_weights.append(weight)
                        segment_node_data.append(n)
                        last_node = n

                # conclude the segment if this is the last node traversed
                if i == len(frag.vs) - 1:
                    if len(segment_node_data) != 0:
                        all_segments.append([segment_node_data, sum(edge_weights)])
                        segment_node_data = []
                        edge_weights = []
                        last_node = -1

            # conclude the segment when reached a terminal or branching point
            else:
                if len(segment_node_data) != 0:
                    all_segments.append([segment_node_data, sum(edge_weights)])
                    segment_node_data = []
                    edge_weights = []
                    last_node = -1

    # add edges and delete nodes
    edge_nodes = []
    for f in all_segments:
        nodes = f[0]
        weight = f[1]

        if len(nodes) == 1:
            ends = frag.neighbors(nodes[0])
        else:
            neighs_a = frag.neighbors(nodes[0])
            neighs_b = frag.neighbors(nodes[-1])
            end_a = [n for n in neighs_a if n not in nodes]
            end_b = [n for n in neighs_b if n not in nodes]
            ends = end_a + end_b

        if len(ends) != 2:
            raise Exception('Invalid pairs to connect')
        else:
            # add edges that connect network nodes to bulk nodes
            weight += frag.es[frag.get_eid(nodes[0], ends[0])]['distance']
            weight += frag.es[frag.get_eid(nodes[-1], ends[1])]['distance']
            frag.add_edge(ends[0], ends[1], distance=weight)
            edge_nodes = edge_nodes + nodes

    frag.delete_vertices(edge_nodes)
    return frag

def _round_coord(coord, decimals=3):
    result = np.round(np.array(coord), decimals)
    return result


def _coord_to_node(coord_tree, coord):
    # nearest network node to `coord`. Uses a prebuilt KD-tree query instead of
    # an O(N) linear scan; this returns the same argmin-of-Euclidean index.
    _, min_node = coord_tree.query(np.asarray(coord, dtype=float))
    return int(min_node)


def _read_array(path, skiprows=0):
    """
    Read a whitespace-delimited numeric file into a 2-D float array.

    The whole file is pulled in with ONE sequential read, then parsed in memory.
    This is a drop-in, value-identical replacement for `np.loadtxt(path,
    skiprows=...)` that avoids issuing many tiny read() syscalls, which can be
    far slower than a single bulk read on uncached network filesystems.
    """
    with open(path, 'rb') as f:
        raw = f.read()
    return pd.read_csv(io.BytesIO(raw), sep=r'\s+', skiprows=skiprows,
                       header=None, dtype=np.float64).to_numpy()


def _read_table(path, delimiter='\t'):
    """Read a delimited table with header into a DataFrame, via one bulk read
    (see `_read_array` for why this matters on uncached network filesystems)."""
    with open(path, 'rb') as f:
        raw = f.read()
    return pd.read_csv(io.BytesIO(raw), delimiter=delimiter)

