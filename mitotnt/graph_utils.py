import pandas as pd 
import numpy as np
import networkx as nx
from collections import defaultdict
from mayavi import mlab

class GraphObject:

	def __init__(self, nx_graph, positionDict, intensityDict, widthDict, 
			     fragmentDict=None, segmentDict=None, isConnected=False):
		
		'''
		Abstract data structure for mitochondrial temporal graph modelling
	
		Attributes
		----------
		nx_graph : networkx.Graph()
			The NetworkX graph representing the mitochondrial network
		positionDict : dict
			A dictionary mapping each node to a position tuple (x, y, z)
		intensityDict: dict
			A dictionary mapping each node to an intensity value
		widthDict: dict
			A dictionary mapping each node to a width value
		isConnected : bool, optional
			Bool indicating if the graph is connected
		'''

		self.graph = nx_graph
		self.positionDict = positionDict
		self.intensityDict = intensityDict
		self.widthDict = widthDict
		self.length = len(self.graph.nodes())
		self.isConnected = isConnected
		if fragmentDict is not None:
			self.fragmentDict = fragmentDict
		if segmentDict is not None:
			self.segmentDict = segmentDict

	@classmethod
	def createFromSubgraph(cls, subgraph, graphObject, isConnected=True):
		
		'''
		Constructor overload for creating an instance from a subgraph of another GraphObject
		Parameters
		----------
		subgraph: networkx.Graph()
			The subgraph from which the GraphObject is to be created
		graphObject: GraphObject
			The original GraphObject 
		Returns
		----------
		A GraphObject created from the subgraph
		'''
		subgraphNodes = subgraph.nodes()
		subgraphPositionDict = {k:graphObject.positionDict[k] for k in subgraphNodes}
		subgraphIntensityDict = {k:graphObject.intensityDict[k] for k in subgraphNodes}
		subgraphWidthDict = {k:graphObject.widthDict[k] for k in subgraphNodes}
		return cls(subgraph, subgraphPositionDict, subgraphIntensityDict, subgraphWidthDict, isConnected=isConnected)

class GraphLoader:

	def __init__(self, path, mode):

		'''
		A utility class for loading a GraphObject and processing it
		Attributes:
		----------
		path: str
			A path to a csv file with the data
		mode: str
			'frame' or 'unique'
		'''

		self.csv = pd.read_csv(path)
		self.mode = mode
		self.csv.loc[self.csv.unique_node_id == 'untracked', 'unique_node_id'] = '-1'
		for column in ['frame_node_id', 'unique_node_id', 'frame_id']:
			self.csv[column] = self.csv[column].apply(lambda x: int(float(x)))
		for column in ['x', 'y', 'z']:
			self.csv[column] = self.csv[column].apply(lambda x: float(x))
		self.csv['connected_{}_node_id'.format(self.mode)] = self.csv['connected_{}_node_id'.format(self.mode)].apply(self.convertToInt)
		self.csv['edgeList'] = self.csv[['{}_node_id'.format(self.mode), 'connected_{}_node_id'.format(self.mode)]].apply(self.getEdgeList, axis=1)
		self.csv['nodePosition'] = list(zip(self.csv.x, self.csv.y, self.csv.z))
		self.n_neighbours = defaultdict(list)
		self.n_unique_neighbours = defaultdict(set)
		self.consistent_graph = None
		self.max_fragment_index = 0

	def createFrameGraph(self, frame=0):

		'''
		Create the mitochondrial network for one frame
		Parameters
		----------
		frame: int
			The frame number for which the graph is to be created
		Returns
		----------
		A GraphObject for the given frame
		'''

		frame_df = self.csv[self.csv['frame_id']==frame]
		positionDict = dict(zip(frame_df['{}_node_id'.format(self.mode)], frame_df['nodePosition']))
		intensityDict = dict(zip(frame_df['{}_node_id'.format(self.mode)], frame_df['intensity']))
		widthDict = dict(zip(frame_df['{}_node_id'.format(self.mode)], frame_df['width']))
		fragmentDict = dict(zip(frame_df['{}_node_id'.format(self.mode)], frame_df['frame_frag_id']))
		try:
			segmentDict = dict(zip(frame_df['{}_node_id'.format(self.mode)], frame_df['frame_seg_id']))
		except:
			segmentDict = None
		edge_list = []
		edge_list.extend([tuple(item) for nodeList in frame_df['edgeList'] for item in nodeList])
		G = nx.Graph()
		G.add_nodes_from(list(frame_df['{}_node_id'.format(self.mode)]))
		G.add_edges_from(edge_list)
		return GraphObject(G, positionDict, intensityDict, widthDict, fragmentDict=fragmentDict, segmentDict=segmentDict)

	def convertToInt(self, x):
		try:
			return np.array(x.split(' '), dtype = int)
		except:
			return []

	def getEdgeList(self, x):
		a, b = [x[0]], x[1]
		return np.transpose([np.tile(a, len(b)), np.repeat(b, len(a))])

	def drawGraph(self, graphObject, withLabels=True, node_color='blue', scale=0.001, label=None, edge_color='blue', fill_between=True):
		if fill_between:
			max_node = max(graphObject.graph.nodes())
			for pair in list(graphObject.graph.edges()):
				point1 = graphObject.positionDict[pair[0]]
				point2 = graphObject.positionDict[pair[1]]
				dist = np.sqrt(np.sum(np.array(point1 - np.array(point2))**2))
				if dist > 1:
					npoints = dist//1 + 1
					x_list = list(np.linspace(point1[0], point2[0], npoints.astype(int)))
					y_list = list(np.linspace(point1[1], point2[1], npoints.astype(int)))
					for i in range(npoints.astype(int)):
						graphObject.graph.add_node(max_node+1)
						graphObject.positionDict[max_node+1] = (x_list[i], y_list[i], 0)
						graphObject.intensityDict[max_node+1] = np.mean(list(graphObject.intensityDict.values()))
						max_node += 1

		nx.draw_networkx(graphObject.graph, pos = {k:v[:2] for (k,v) in graphObject.positionDict.items()}, with_labels=withLabels,
			node_size = [s*scale for (k,s) in graphObject.intensityDict.items()], 
			font_size=18, node_color=node_color, edge_color=edge_color, label=label)

	def drawGraph3D(self, graphObject, voxel_size, scale, color, draw_axes=False, draw_lines=False, fill_between=False):
		
		'''
		Draw a 3D fragment plot using Mayavi
		Parameters
		----------
		graphObject: GraphObject
			The GraphObject for the given frame
		color: tuple
			A tuple (r, g, b) of floats
		draw_axes: bool, optional
			Whether to draw a box around the cell
		draw_lines: bool, optional
			Ignore
		fill_between: bool, optional
			Whether to fill in gaps in fragments
		'''

		node_xyz = np.array(list(graphObject.positionDict.values()))
		pts = mlab.points3d(node_xyz[:,0], node_xyz[:,1], node_xyz[:,2],
                    color=color,
                    scale_factor=scale,
                    scale_mode='none',
                    colormap='Blues',
                    resolution=20)

		if draw_lines:
			pts.mlab_source.dataset.lines = np.array(graphObject.graph.edges())
			tube = mlab.pipeline.tube(pts, tube_radius=0.1)
			mlab.pipeline.surface(tube, color=(0.2, 0.2, 0.2))

		if fill_between:
			fill_pts = [[], [], []]
			for pair in list(graphObject.graph.edges()):
				point1 = graphObject.positionDict[pair[0]]
				point2 = graphObject.positionDict[pair[1]]
				dist = np.sqrt(np.sum(np.array(point1 - np.array(point2))**2))
				if dist > 1.5 * voxel_size:
					npoints = dist//1 + 1
					fill_pts[0].extend(list(np.linspace(point1[0], point2[0], npoints.astype(int))))
					fill_pts[1].extend(list(np.linspace(point1[1], point2[1], npoints.astype(int))))
					fill_pts[2].extend(list(np.linspace(point1[2], point2[2], npoints.astype(int))))

			fill_pts = mlab.points3d(np.array(fill_pts[0]), np.array(fill_pts[1]), np.array(fill_pts[2]),
                    color=color,
                    scale_factor=scale,
                    scale_mode='none',
                    colormap='Blues',
                    resolution=20)

		if draw_axes:
			mlab.outline(color = (0.2, 0.2, 0.2))

	def computeFragmentStatistics(self, frame_graph_obj_1, frame_graph_obj_2, delta=1, min_track_size=3):
		if not frame_graph_obj_1.isConnected:
			tracked_components = self.trackFragmentsBipartite(frame_graph_obj_1, frame_graph_obj_2)
		else:
			tracked_components = {0:(frame_graph_obj_1, frame_graph_obj_2)}
		MeanRGArray = []
		MeanDisplacementArray = []
		TotalDistanceArray= []
		for fragment in tracked_components.values():
			if fragment[0].length >= min_track_size and fragment[1].length >= min_track_size:
				fragment1PositionArray = np.array(list(fragment[0].positionDict.values()))
				fragment2PositionArray = np.array(list(fragment[1].positionDict.values()))
				fragment1Centroid = np.nanmean(fragment1PositionArray, axis=0)
				fragment2Centroid = np.nanmean(fragment2PositionArray, axis=0)
				MeanDisplacementArray.append(np.sqrt(np.sum(fragment1Centroid - fragment2Centroid)**2))
				TotalDistanceArray.append(np.mean([np.sqrt(np.sum((fragment[0].positionDict[node] + \
					fragment2Centroid - fragment1Centroid) - fragment[1].positionDict[node])**2) \
					for node in set(fragment[0].graph.nodes()).intersection(set(fragment[1].graph.nodes()))]))
				fragment1RG = np.sqrt(np.sum(((fragment1PositionArray - fragment1Centroid)**2).reshape(3, -1))/fragment[0].length)
				MeanRGArray.append(fragment1RG)
		frame_labels = ['Radius of Gyration', 'Mean Displacement', 
						'Mean Fragment Intensity', 'Mean Fragment Width',
						'Density', 'Average Degree', 'Average Shortest Path Length',
						'Network Efficiency', 'Fragment Size', 'Temporal Intersection']
		MeanIntensityNode = sum([(fragment[0].intensityDict[node]) \
			for node in list(set(fragment[0].graph.nodes()))])/(len(list(set(fragment[0].graph.nodes()))) + np.finfo(float).eps)
		MeanWidthNode = sum([(fragment[0].widthDict[node]) \
			for node in list(set(fragment[0].graph.nodes()))])/(len(list(set(fragment[0].graph.nodes()))) + np.finfo(float).eps)
		density = nx.density(fragment[0].graph)
		averageDegree = sum(dict(nx.degree(fragment[0].graph)).values())/(len(fragment[0].graph) + np.finfo(float).eps)
		connectedComponents = sorted(nx.connected_components(fragment[0].graph), key=len, reverse=True)
		averagePath = sum([nx.average_shortest_path_length(fragment[0].graph.subgraph(c).copy()) for c in connectedComponents])/(len(connectedComponents) + np.finfo(float).eps)
		eff = nx.global_efficiency(fragment[0].graph)
		temporalIntersection = (len(list(set(fragment[0].graph) & set(fragment[1].graph))) + 1)/(len(list(set(fragment[0].graph))) + 1)

		return (frame_labels, np.array([np.nanmean(MeanRGArray), np.nanmean(MeanDisplacementArray), 
			    MeanIntensityNode, MeanWidthNode, density, averageDegree, averagePath,
			    eff, fragment[0].length, temporalIntersection]))

	def trackFragmentsBipartite(self, frame_graph_obj_1, frame_graph_obj_2, threshold=0.2, maintain_fused=False):
		connectedComponents1 = {'A' + str(i): c for i, c in enumerate(sorted(nx.connected_components(frame_graph_obj_1.graph), key=len, reverse=True))}
		connectedComponents2 = {'B' + str(i): c for i, c in enumerate(sorted(nx.connected_components(frame_graph_obj_2.graph), key=len, reverse=True))}
		connectedComponentsCombined = {**connectedComponents1, **connectedComponents2}
		if self.consistent_graph:
			assert self.consistent_graph == connectedComponents1
		bipartite = nx.Graph()
		bipartite.add_nodes_from([n for n in connectedComponents1.keys()], bipartite=0)
		bipartite.add_nodes_from([n for n in connectedComponents2.keys()], bipartite=1)
		bipartite.add_edges_from([(i, j) for i in connectedComponents1.keys() for j in \
			connectedComponents2.keys() if ((len(connectedComponents1[i].intersection(connectedComponents2[j]))/max(len(connectedComponents1[i]), len(connectedComponents2[j]))) > threshold or len(connectedComponents1[i].intersection(connectedComponents2[j])) > 10)])
		CC = {i: c for i, c in enumerate(sorted(nx.connected_components(bipartite), key=len, reverse=True))}
		
		def __flatten__(t):
			return [item for sublist in t for item in sublist]

		try:
			CC_index = {k: np.array(__flatten__([[self.id_dict.get(f, np.nan)]*len(connectedComponentsCombined[f]) if self.id_dict.get(f, np.nan) < len(CC.values()) else [np.nan] for f in v])) for k, v in CC.items()}
		except AttributeError as e:
			CC_index = {k: np.array([k for f in v]) for k, v in CC.items()}
		CC_reindex = {}
		for key, index_list in CC_index.items():
			try:
				CC_reindex[key] = np.argmax(np.bincount(index_list[~np.isnan(index_list)].astype(int)))
			except ValueError as e:
				try:
					if not self.max_fragment_index:
						self.max_fragment_index = max(CC_reindex.values())+1
					else:
						self.max_fragment_index += 1
					CC_reindex[key] = self.max_fragment_index
				except: pass
		new_CC = {}
		for k, v in CC_reindex.items():
			if maintain_fused:
				try:
					new_CC[v].update(CC[k])
				except:
					new_CC[v] = CC[k]
			else:
				if new_CC.get(v) is not None:
					if sum([len(connectedComponentsCombined[f]) for f in new_CC[v]]) > sum([len(connectedComponentsCombined[f]) for f in CC[k]]):
						max_index = max(new_CC.keys())
						new_CC[max_index + 1] = CC[k]
					else:
						max_index = max(new_CC.keys())
						new_CC[max_index + 1] = new_CC[v]
						new_CC[v] = CC[k]
				else:
					new_CC[v] = CC[k]
		CC = new_CC
		self.id_dict = {'A' + f[1:]:k for k, v in CC.items() for f in v if f[0]=='B'}
		A = {k:set().union(*[connectedComponents1[n] for n in CC[k] if n[0] == 'A']) for k in sorted(CC.keys())}
		B = {k:set().union(*[connectedComponents2[n] for n in CC[k] if n[0] == 'B']) for k in sorted(CC.keys())}
		self.consistent_graph = {'A' + k[1:]:v for k,v in connectedComponents2.items()}
		return {k:(GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(A.get(k))).copy(), frame_graph_obj_1), 
			    GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(B.get(k))).copy(), frame_graph_obj_2)) for k, v in sorted(CC.items(), key=lambda x: len(x[1]), reverse=True)}

	def trackMaximumVote(self, frame_graph_obj_1, frame_graph_obj_2, inverted=False, threshold=5, level='segment'):
		
		if level == 'segment':
			connectedComponents1 = {i: set([x for x in frame_graph_obj_1.segmentDict if frame_graph_obj_1.segmentDict[x] == f]) \
									for i, f in enumerate(frame_graph_obj_1.segmentDict)}
			connectedComponents2 = {i: set([x for x in frame_graph_obj_2.segmentDict if frame_graph_obj_2.segmentDict[x] == f]) \
								for i, f in enumerate(frame_graph_obj_2.segmentDict)}
		elif level == 'fragment':
			connectedComponents1 = {i: set([x for x in frame_graph_obj_1.fragmentDict if frame_graph_obj_1.fragmentDict[x] == f]) \
									for i, f in enumerate(frame_graph_obj_1.fragmentDict)}
			connectedComponents2 = {i: set([x for x in frame_graph_obj_2.fragmentDict if frame_graph_obj_2.fragmentDict[x] == f]) \
								for i, f in enumerate(frame_graph_obj_2.fragmentDict)}
		else:
			raise ValueError("Only implemented for 'segment' or 'fragment'")
		
		connectedComponents1 = {x[0]:x[1] for x in sorted(connectedComponents1.items(), key=lambda k: len(k[1]), reverse=True)}
		connectedComponents1 = {i: connectedComponents1[k] for i, k in enumerate(connectedComponents1) if len(connectedComponents1[k])}
		
		connectedComponents2 = {x[0]:x[1] for x in sorted(connectedComponents2.items(), key=lambda k: len(k[1]), reverse=True)}
		connectedComponents2 = {i: connectedComponents2[k] for i, k in enumerate(connectedComponents2) if len(connectedComponents2[k])}

		superFragments1 = []
		if self.consistent_graph is not None:
			superFragments1 = self.consistent_graph
		else:
			superFragments1 = connectedComponents1
		cost_matrix = np.array([[(len(f1 & c2))/max(len(f1), len(c2)) for c2 in connectedComponents2.values()] for f1 in superFragments1.values()])
		superFragmentsDict = defaultdict(set)
		for i, f in enumerate(connectedComponents2.values()):
			if len(f) < threshold:
				continue
			maxAlignment = np.max(cost_matrix[:, i])
			if maxAlignment > 0.2:
				maxAlignmentIndex = np.argmax(cost_matrix[:, i])
				if superFragmentsDict.get(list(superFragments1.keys())[maxAlignmentIndex]):
					self.max_fragment_index = max(max(superFragments1.keys()), self.max_fragment_index) + 1
					if len(superFragmentsDict.get(list(superFragments1.keys())[maxAlignmentIndex])) > len(f):
						superFragmentsDict[self.max_fragment_index].update(f)
						superFragments1[self.max_fragment_index] = set()
					else:
						superFragmentsDict[self.max_fragment_index] = superFragmentsDict[list(superFragments1.keys())[maxAlignmentIndex]]
						superFragments1[self.max_fragment_index] = set()
						superFragmentsDict[list(superFragments1.keys())[maxAlignmentIndex]] = f
				else:
					superFragmentsDict[list(superFragments1.keys())[maxAlignmentIndex]].update(f)
			else:
				self.max_fragment_index = max(max(superFragments1.keys()), self.max_fragment_index) + 1
				superFragmentsDict[self.max_fragment_index].update(f)
				superFragments1[self.max_fragment_index] = set()
		self.consistent_graph = superFragmentsDict
		superFragments2 = superFragmentsDict
		if inverted:
			return list(zip(np.array([GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(c)).copy(), frame_graph_obj_2)
			 for c in superFragments2], dtype = object), np.array([GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(c)).copy(), frame_graph_obj_1) 
			for c in superFragments1], dtype = object)))
		return {k:(GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(superFragments1.get(k))).copy(), frame_graph_obj_1), 
			    GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(superFragments2.get(k))).copy(), frame_graph_obj_2)) for k, v in superFragments2.items()}