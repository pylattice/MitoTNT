#     ### Save tracks and connected nodes using frame indexing ###
#     new_tracks = pd.DataFrame()
#     for frame in trange(start_frame, end_frame, tracking_interval):

#         full_graph = full_graph_all_frames[frame]
#         cc = full_graph.components()

#         segment_nodes = segment_node_all_frames[frame]
#         branching_nodes = []
#         for i in range(len(full_graph.vs)):
#             if full_graph.vs[i].degree() > 2:
#                 branching_nodes.append(i)
#         node_to_segment = {}
#         for segment_id, segment in enumerate(segment_nodes): # segment consists of of segment nodes
#             for b in segment:
#                 if b in branching_nodes:
#                     node_to_segment[b] = np.nan
#                 else:
#                     node_to_segment[b] = segment_id

#         coords = full_graph.vs['coordinate']
#         intensities, widths = full_graph.vs['intensity'], full_graph.vs['width']
#         tracks_frame = tracks[tracks['frame_id'] == frame]

#         # add new rows for untracked nodes
#         untracked = []
#         tracked_frame_node_id = tracks_frame['frame_node_id'].astype('int').tolist()
#         for node in range(len(coords)):
#             if node not in tracked_frame_node_id:
#                 coord = coords[node]
#                 x, y, z = coord[0], coord[1], coord[2]
#                 untracked.append({'frame_id': frame, 'unique_node_id': 'untracked',
#                                   'frame_node_id': node, 'frame_seg_id': node_to_segment[node], 'frame_frag_id': cc.membership[node],
#                                   'x': x, 'y': y, 'z': z,
#                                   'intensity':intensities[node], 'width':widths[node]})

#         untracked_tracks = pd.DataFrame.from_dict(untracked)
#         tracks_frame = pd.concat([tracks_frame, untracked_tracks])

#         # add new column with neighbor nodes for each node
#         all_connected_frame_nodes = []
#         tracked_frame_node_id = tracks_frame['frame_node_id'].astype('int').tolist()
#         for node in tracked_frame_node_id:
#             neighs = full_graph.neighbors(int(node))
#             all_connected_frame_nodes.append(neighs)

#         tracks_frame.insert(2, 'connected_frame_node_id', [list_to_str(a) for a in all_connected_frame_nodes])

#         # accumulate tracks from each frame
#         new_tracks = pd.concat([new_tracks, tracks_frame])

#     # reorder the columns
#     new_tracks = new_tracks[['frame_id', 'unique_node_id', 'frame_node_id', 'frame_seg_id', 'frame_frag_id', 'connected_frame_node_id', 'x', 'y', 'z', 'intensity', 'width']]
#     new_tracks.to_csv(output_dir+'tracks_with_connectivity_frame_index.csv', index=False)
#     ### Done saving ###