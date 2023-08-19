import os
import numpy as np
import pandas as pd
import igraph as ig
from tqdm.notebook import trange

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

def generate(data_dir, input_dir,
             start_frame, end_frame,
             node_gap_size=0):

    ### read mitograph outputs to create 1) classic graphs, 2) full-resolution graphs with node and edge attributes, 3) segment nodes, 4) classic graphs per node ###
    
    # 1) classic graphs
    all_classic_graphs = []
    for frame in trange(start_frame, end_frame+1, desc='Classic graphs in progress'):
        
        classic_graph = ig.Graph()

        network_nodes = np.loadtxt(data_dir+'frame_'+str(frame)+'/frame_'+str(frame)+'.coo')
        
        edge_list = np.loadtxt(data_dir+'frame_'+str(frame)+'/frame_'+str(frame)+'.gnet', skiprows=1)
        
        bulk_nodes = pd.read_csv(data_dir+'frame_'+str(frame)+'/frame_'+str(frame)+'.txt', delimiter='\t')
        
        coords = round_coord(network_nodes)

        classic_graph.add_vertices(len(coords))
        
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

            node_end_a = classic_graph.vs[index_end_a]
            node_end_a['index'] = node_end_a.index
            node_end_a['coordinate'] = line_nodes.loc[0, 'x':'z'].to_numpy()
            node_end_a['intensity'] = line_nodes.loc[0, 'pixel_intensity']
            node_end_a['width'] = line_nodes.loc[0, 'width_(um)']

            node_end_b = classic_graph.vs[index_end_b]
            node_end_b['index'] = node_end_b.index
            node_end_b['coordinate'] = line_nodes.loc[end_index, 'x':'z'].to_numpy()
            node_end_b['intensity'] = line_nodes.loc[end_index, 'pixel_intensity']
            node_end_b['width'] = line_nodes.loc[end_index, 'width_(um)']
        
        for edge in edge_list:
            node_end_a, node_end_b, distance = int(edge[0]), int(edge[1]), edge[2]
            classic_graph.add_edge(node_end_a, node_end_b, distance=distance)

        all_classic_graphs.append(classic_graph)
    
    all_classic_graphs = np.array(all_classic_graphs, dtype=object)

    
    # 2) full-resolution graphs with node and edge attributes + 3) segment nodes
    all_full_graphs = []
    all_segment_nodes = []
    for frame in trange(start_frame, end_frame+1, desc='Full-resolution graphs in progress'):
        
        full_graph = ig.Graph()

        raw_coords = np.loadtxt(data_dir+'frame_'+str(frame)+'/frame_'+str(frame)+'.coo')
        
        bulk_nodes = pd.read_csv(data_dir+'frame_'+str(frame)+'/frame_'+str(frame)+'.txt', delimiter='\t')
        
        coords = round_coord(raw_coords)

        # create all the network nodes
        full_graph.add_vertices(len(coords))

        line_ids = np.unique(bulk_nodes['line_id'])
        frame_segment_nodes = []

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
                    dist = np.linalg.norm(full_graph.vs[last_node]['coordinate'] - full_graph.vs[current_node]['coordinate'])
                    full_graph.add_edge(last_node, current_node, distance=dist)
                    last_node = current_node

                    # add segment node index
                    segment_nodes.append(current_node)

                # link last bulk node to the other network node
                if index == sel_node_index[-1]:
                    dist = np.linalg.norm(full_graph.vs[last_node]['coordinate'] - full_graph.vs[index_end_b]['coordinate'])
                    full_graph.add_edge(last_node, index_end_b, distance=dist)

                    # add segment node index
                    segment_nodes.append(index_end_b) # get the node index of another end
                    frame_segment_nodes.append(segment_nodes) # finish this segment

        full_graph.simplify(combine_edges='sum') # remove self-loops, and combine multiple edges if needed

        all_full_graphs.append(full_graph)
        all_segment_nodes.append(frame_segment_nodes)

    all_full_graphs = np.array(all_full_graphs, dtype=object)
    all_segment_nodes = np.array(all_segment_nodes, dtype=object)

    # 4) classic graphs around each node
    all_classic_graphs_per_node = []
    for frame in trange(start_frame, end_frame+1, desc='Classic graphs per node in progress'):

        # load graphs of nodes and edges
        full_graph = all_full_graphs[frame]
        coords = full_graph.vs['coordinate']
        total_num_nodes = len(full_graph.vs)

        # get fragments
        all_frags = full_graph.components()

        # contract edges and update frag and return new root index
        frame_classic_graphs_per_node = []
        for node_index in range(total_num_nodes):

            # use full graph node id
            frag = full_graph.induced_subgraph(all_frags[all_frags.membership[node_index]])

            # find fragment graph node id
            root = frag.vs['index'].index(node_index)

            # call function for edge contraction
            classic_graph = contract_edges(frag, root)
            frame_classic_graphs_per_node.append(classic_graph)

        all_classic_graphs_per_node.append(frame_classic_graphs_per_node)

    all_classic_graphs_per_node = np.array(all_classic_graphs_per_node, dtype=object)

    ### save all inputs as a single compressed .npz file ###
    path = input_dir+'tracking_inputs.npz'
    data = {}
    data['classic_graphs'] = all_classic_graphs
    data['full_graphs'] = all_full_graphs
    data['segment_nodes'] = all_segment_nodes
    data['classic_graphs_per_node'] = all_classic_graphs_per_node

    np.savez(path, **data)