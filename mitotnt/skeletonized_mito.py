import glob
import os
import numpy as np
import pandas as pd
import igraph as ig
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


    def extract_graphs(self, overwrite: bool = False, node_gap_size: int = 0):
        """
        Extract the graph representations of the previously segmented mitochondria used for tracking.

        Parameters
        ----------
        overwrite : bool
            Whether to overwrite existing graphs (default is False).
        node_gap_size : int, optional
            Number of nodes to skip when generating the full graphs (default is 0, keep all).

        Returns
        -------
        None
            The computed graphs and segment nodes are stored internally as `SkeletonizedMito` attributes.
        """

        if os.path.exists(self.save_path+'extracted_graphs.npz') and not overwrite:
            self._load_graphs()
            print("Graphs have already been extracted. Reloaded previous data.")
        else:
            self.extract_full_graphs_and_segment_nodes(node_gap_size=node_gap_size)
            self.extract_simple_graphs()
            self.extract_local_simple_graphs()
            self._save_graphs()


    def extract_full_graphs_and_segment_nodes(self, node_gap_size: int = 0):

        all_full_graphs = []
        all_segment_nodes = []
        for folder in tqdm(self.list_of_folders, desc="Extracting full graphs and segment nodes"):

            full_graph = ig.Graph()

            if len(glob.glob(folder + '/*.coo')) == 1 and len(glob.glob(folder + '/*.gnet')) == 1 and len(glob.glob(folder + '//*.txt')) == 1:
                raw_coords = np.loadtxt(glob.glob(folder + '/*.coo')[0])
                bulk_nodes = pd.read_csv(glob.glob(folder + '/*.txt')[0], delimiter='\t')
            else:
                raise Exception(f"{folder!r} has none/duplicate MitoGraph outputs.")

            coords = round_coord(raw_coords)

            # create all the network nodes
            full_graph.add_vertices(len(coords))

            line_ids = np.unique(bulk_nodes['line_id'])
            frame_segment_nodes = []

            for line in line_ids:
                line_nodes = bulk_nodes[bulk_nodes['line_id'] == line]
                line_nodes = line_nodes.reset_index()
                end_index = len(line_nodes) - 1

                # get branching and terminal nodes to contruct graph
                coord_end_a = round_coord(line_nodes.loc[0, 'x':'z'])
                coord_end_b = round_coord(line_nodes.loc[end_index, 'x':'z'])

                # find index of network nodes in .coo based on coords
                index_end_a = coord_to_node(coords, coord_end_a)
                index_end_b = coord_to_node(coords, coord_end_b)

                node_end_a = full_graph.vs[index_end_a]
                node_end_a['index'] = node_end_a.index
                node_end_a['coordinate'] = line_nodes.loc[0, 'x':'z'].to_numpy()
                node_end_a['intensity'] = line_nodes.loc[0, 'pixel_intensity']
                node_end_a['width'] = line_nodes.loc[0, 'width_(um)']

                node_end_b = full_graph.vs[index_end_b]
                node_end_b['index'] = node_end_b.index
                node_end_b['coordinate'] = line_nodes.loc[end_index, 'x':'z'].to_numpy()
                node_end_b['intensity'] = line_nodes.loc[end_index, 'pixel_intensity']
                node_end_b['width'] = line_nodes.loc[end_index, 'width_(um)']

                # get node id of the nodes on same segment, start with one end
                last_node = index_end_a
                segment_nodes = [index_end_a]
                sel_node_index = range(0, end_index, 1 + node_gap_size)

                for index in sel_node_index:
                    # add bulk nodes
                    if index > 0 and index < sel_node_index[-1]:
                        # add vertex and vertex attributes
                        full_graph.add_vertices(1)
                        bulk_node = full_graph.vs[-1]
                        bulk_node['index'] = bulk_node.index
                        bulk_node['coordinate'] = line_nodes.loc[index, 'x':'z'].to_numpy()
                        bulk_node['intensity'] = line_nodes.loc[index, 'pixel_intensity']
                        bulk_node['width'] = line_nodes.loc[index, 'width_(um)']

                        # add edge and edge attributes
                        current_node = len(full_graph.vs) - 1
                        dist = np.linalg.norm(
                            full_graph.vs[last_node]['coordinate'] - full_graph.vs[current_node]['coordinate'])
                        full_graph.add_edge(last_node, current_node, distance=dist)
                        last_node = current_node

                        # add segment node index
                        segment_nodes.append(current_node)

                    # link last bulk node to the other network node
                    if index == sel_node_index[-1]:
                        dist = np.linalg.norm(
                            full_graph.vs[last_node]['coordinate'] - full_graph.vs[index_end_b]['coordinate'])
                        full_graph.add_edge(last_node, index_end_b, distance=dist)

                        # add segment node index
                        segment_nodes.append(index_end_b)  # get the node index of another end
                        frame_segment_nodes.append(segment_nodes)  # finish this segment

            full_graph.simplify(combine_edges='sum')  # remove self-loops, and combine multiple edges if needed

            all_full_graphs.append(full_graph)
            all_segment_nodes.append(frame_segment_nodes)

        self.full_graphs = np.array(all_full_graphs, dtype=object)
        self.segment_nodes = np.array(all_segment_nodes, dtype=object)


    def extract_simple_graphs(self):

        all_simple_graphs = []
        for folder in tqdm(self.list_of_folders, desc="Extracting simple graphs"):

            simple_graph = ig.Graph()

            if len(glob.glob(folder+'/*.coo')) == 1 and len(glob.glob(folder+'/*.gnet')) == 1 and len(glob.glob(folder+'//*.txt')) == 1:
                raw_coords = np.loadtxt(glob.glob(folder+'/*.coo')[0])
                edge_list = np.loadtxt(glob.glob(folder+'/*.gnet')[0], skiprows=1)
                bulk_nodes = pd.read_csv(glob.glob(folder+'/*.txt')[0], delimiter='\t')
            else:
                raise Exception(f"{folder!r} has none/duplicate MitoGraph outputs.")

            coords = round_coord(raw_coords)

            simple_graph.add_vertices(len(coords))

            # create all the network nodes
            line_ids = np.unique(bulk_nodes['line_id'])

            for line in line_ids:
                line_nodes = bulk_nodes[bulk_nodes['line_id']==line]
                line_nodes = line_nodes.reset_index()
                end_index = len(line_nodes) - 1

                # get branching and terminal nodes to contruct graph
                coord_end_a = round_coord(line_nodes.loc[0, 'x':'z'])
                coord_end_b = round_coord(line_nodes.loc[end_index, 'x':'z'])

                # find index of network nodes in .coo based on coords
                index_end_a = coord_to_node(coords, coord_end_a)
                index_end_b = coord_to_node(coords, coord_end_b)

                node_end_a = simple_graph.vs[index_end_a]
                node_end_a['index'] = node_end_a.index
                node_end_a['coordinate'] = line_nodes.loc[0, 'x':'z'].to_numpy()
                node_end_a['intensity'] = line_nodes.loc[0, 'pixel_intensity']
                node_end_a['width'] = line_nodes.loc[0, 'width_(um)']

                node_end_b = simple_graph.vs[index_end_b]
                node_end_b['index'] = node_end_b.index
                node_end_b['coordinate'] = line_nodes.loc[end_index, 'x':'z'].to_numpy()
                node_end_b['intensity'] = line_nodes.loc[end_index, 'pixel_intensity']
                node_end_b['width'] = line_nodes.loc[end_index, 'width_(um)']

            for edge in edge_list:
                node_end_a, node_end_b, distance = int(edge[0]), int(edge[1]), edge[2]
                simple_graph.add_edge(node_end_a, node_end_b, distance=distance)

            all_simple_graphs.append(simple_graph)

        self.simple_graphs = np.array(all_simple_graphs, dtype=object)


    def extract_local_simple_graphs(self):
    
        all_local_simple_graphs = []
        try:
            all_full_graphs = self.full_graphs
        except:
            raise Exception("No full graphs to extract. Run extract_full_graphs_and_segment_nodes() first.")
    
        for frame_full_graphs in tqdm(all_full_graphs, desc="Extracting local simple graphs"):
    
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
                simple_graph = contract_edges(frag, root)
                frame_simple_graphs_per_node.append(simple_graph)
    
            all_local_simple_graphs.append(frame_simple_graphs_per_node)
    
        self.local_simple_graphs = np.array(all_local_simple_graphs, dtype=object)
    
    
    def _save_graphs(self):
        # save all inputs as a single compressed .npz file
        os.makedirs(self.save_path, exist_ok=True)
        path = self.save_path + 'extracted_graphs.npz'
        data = {'full_graphs': self.full_graphs,
                'segment_nodes': self.segment_nodes,
                'simple_graphs': self.simple_graphs,
                'local_simple_graphs': self.local_simple_graphs}
        np.savez(path, **data)
    
    
    def _load_graphs(self):
        try:
            data = np.load(self.save_path+'extracted_graphs.npz', allow_pickle=True)
        except:
            raise Exception(f"No extracted_graphs.npz found under {self.save_path!r}.")
    
        self.full_graphs = data['full_graphs']
        self.segment_nodes = data['segment_nodes']
        self.simple_graphs = data['simple_graphs']
        self.local_simple_graphs = data['local_simple_graphs']


def round_coord(coord, decimals=3):
    result = np.round(np.array(coord), decimals)
    return result


def coord_to_node(all_coords, coord):
    dist = [np.linalg.norm(disp) for disp in all_coords - coord]
    return np.argmin(dist)


def contract_edges(frag, root):

    bulk_nodes = []
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
                    bulk_nodes.append(n)
                    last_node = n

                # sometimes in large graph a new segment is visited without reaching a degree!=2 node
                else:
                    # this may fail when reach end of one segment and jump to the start of another segment
                    try:
                        weight = frag.es[frag.get_eid(n, last_node)]['distance']

                    # when it fails just start another bulk
                    except Exception:
                        if len(bulk_nodes) != 0:
                            all_segments.append([bulk_nodes, sum(edge_weights)])
                            bulk_nodes = []
                            edge_weights = []

                            bulk_nodes.append(n)
                            last_node = n

                    # when the two nodes are on bulk we can just append distance and node
                    else:
                        edge_weights.append(weight)
                        bulk_nodes.append(n)
                        last_node = n

                # conclude the segment if this is the last node traversed
                if i == len(frag.vs) - 1:
                    if len(bulk_nodes) != 0:
                        all_segments.append([bulk_nodes, sum(edge_weights)])
                        bulk_nodes = []
                        edge_weights = []
                        last_node = -1

            # conclude the segment when reached a terminal or branching point
            else:
                if len(bulk_nodes) != 0:
                    all_segments.append([bulk_nodes, sum(edge_weights)])
                    bulk_nodes = []
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
    frag.simplify(combine_edges='sum')

    return frag
