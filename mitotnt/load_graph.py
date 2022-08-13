import pandas as pd 
import numpy as np
import math
from itertools import product, combinations
import networkx as nx
from networkx.algorithms import community, global_efficiency
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
import glob
from tqdm import tqdm 
from collections import defaultdict
from scipy.optimize import linear_sum_assignment
from operator import itemgetter
from collections import Counter
from mayavi import mlab
import vtk
from scipy.stats import multivariate_normal, entropy
import time
# from mpl_toolkits.mplot3d import Axes3D

class GraphObject:

	def __init__(self, nx_graph, positionDict, intensityDict, widthDict, isConnected=False):
		
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
		edge_list = []
		edge_list.extend([tuple(item) for nodeList in frame_df['edgeList'] for item in nodeList])
		G = nx.Graph()
		G.add_nodes_from(list(frame_df['{}_node_id'.format(self.mode)]))
		G.add_edges_from(edge_list)
		return GraphObject(G, positionDict, intensityDict, widthDict)

	def createTemporalDynamicGraph(self, frame_graph_obj_1, frame_graph_obj_2, timestep, directed=False):
		start_time = time.time()
		if timestep == 0:
			self.temporal_graph_nodes = [str(timestep) + '_' + str(node) \
				for node in frame_graph_obj_1.graph.nodes()]
			self.intensityDict = {str(timestep) + '_' + str(k):v for (k, v) in frame_graph_obj_1.intensityDict.items()} 
			self.positionDict = {str(timestep) + '_' + str(k):v for (k, v) in frame_graph_obj_1.positionDict.items()} 
			self.widthDict = {str(timestep) + '_' + str(k):v for (k, v) in frame_graph_obj_1.widthDict.items()} 
		self.temporal_graph_nodes.extend([str(timestep + 1) + '_' + str(node) \
			for node in frame_graph_obj_2.graph.nodes()])
		# print('Node List Time', time.time() - start_time)
		start_time = time.time()
		try:
			self.temporal_graph_edges.extend([(str(timestep) + '_' + str(node1), 
				str(timestep+1) + '_' + str(node2)) for node1 in frame_graph_obj_1.graph.nodes() \
			for node2 in frame_graph_obj_2.graph.nodes() if (node2 in frame_graph_obj_1.graph.adj[node1] or node1==node2)])
		except Exception as e:
			self.temporal_graph_edges = [(str(timestep) + '_' + str(node1), 
				str(timestep+1) + '_' + str(node2)) for node1 in frame_graph_obj_1.graph.nodes() \
			for node2 in frame_graph_obj_2.graph.nodes() if (node2 in frame_graph_obj_1.graph.adj[node1] or node1==node2)]
		# print('Edge List Time', time.time() - start_time)
		start_time = time.time()
		if directed:
			G = nx.DiGraph()
		else:
			G = nx.Graph()
		G.add_nodes_from(self.temporal_graph_nodes)
		G.add_edges_from(self.temporal_graph_edges)
		intensityDict2 = {str(timestep+1) + '_' + str(k):v for (k, v) in frame_graph_obj_2.intensityDict.items()}
		self.intensityDict = {**self.intensityDict, **intensityDict2}
		positionDict2 = {str(timestep+1) + '_' + str(k):v for (k, v) in frame_graph_obj_2.positionDict.items()}
		self.positionDict = {**self.positionDict, **positionDict2}
		widthDict2 = {str(timestep+1) + '_' + str(k):v for (k, v) in frame_graph_obj_2.widthDict.items()}
		self.widthDict = {**self.widthDict, **widthDict2}
		temporalGraph = GraphObject(G, self.positionDict, self.intensityDict, self.widthDict)
		# print(positionDict)
		# print('GraphObject Time', time.time() - start_time)
		start_time = time.time()
		# self.drawGraph(temporalGraph, withLabels=False)
		# plt.show()
		return temporalGraph

	def computeTemporalMetricsDP(self, frame_graph_obj_1, frame_graph_obj_2, timestep, damaged=False, centrality_dict=None, local=False, local_nodes = [], random=False):
		
		if random:
			for node in centrality_dict.keys():
				centrality_dict[node] = np.random.uniform()

		# Add Temporal Centrality
		if damaged:
			centrality_threshold = np.percentile(list(centrality_dict.values()), 95)

		if timestep == 0:
			if local:
				self.shortest_path_dict = {str(timestep) + '_' + str(node):set() for node in frame_graph_obj_1.graph.nodes()}
				for node in local_nodes:
					self.shortest_path_dict[str(timestep) + '_' + str(node)] = {node}
			else:
				self.shortest_path_dict = {str(timestep) + '_' + str(node):{node} for node in frame_graph_obj_1.graph.nodes()}

		frame_reachable_dict = {}

		for node in frame_graph_obj_2.graph.nodes():
			if damaged:
				if centrality_dict[str(timestep+1) + '_' + str(node)] > centrality_threshold:
					self.shortest_path_dict[str(timestep+1) + '_' + str(node)] = {node}
					frame_reachable_dict[node] = {node}
					continue
			try:
				for neighbour in frame_graph_obj_1.graph.adj[node]:
					try:
						self.shortest_path_dict[str(timestep+1) + '_' + str(node)].update(self.shortest_path_dict[str(timestep) + '_' + str(neighbour)])
						frame_reachable_dict[node].update(self.shortest_path_dict[str(timestep) + '_' + str(neighbour)])
					except:
						self.shortest_path_dict[str(timestep+1) + '_' + str(node)] = \
						self.shortest_path_dict[str(timestep) + '_' + str(node)].union(self.shortest_path_dict[str(timestep) + '_' + str(neighbour)])
						frame_reachable_dict[node] = \
						self.shortest_path_dict[str(timestep) + '_' + str(node)].union(self.shortest_path_dict[str(timestep) + '_' + str(neighbour)])
				if len(frame_graph_obj_1.graph.adj[node]) == 0:
					if not local:
						self.shortest_path_dict[str(timestep+1) + '_' + str(node)] = {node}
						frame_reachable_dict[node] = {node}
					else:
						self.shortest_path_dict[str(timestep+1) + '_' + str(node)] = set()
						frame_reachable_dict[node] = set()
			except KeyError as e:
				# print(frame_graph_obj_1.graph.adj)
				# print(timestep, str(timestep+1) + '_' + str(node))
				if not local:
					self.shortest_path_dict[str(timestep+1) + '_' + str(node)] = {node}
					frame_reachable_dict[node] = {node}
				else:
					self.shortest_path_dict[str(timestep+1) + '_' + str(node)] = set()
					frame_reachable_dict[node] = set()
				

		return self.shortest_path_dict, frame_reachable_dict

	def computeTemporalRandomWalk(self, frame_graph_obj_1, frame_graph_obj_2, timestep):
		
		self.num_token_dict = defaultdict(lambda: 0)

		if timestep == 0:
			self.token_dict = {node:1 for node in frame_graph_obj_1.graph.nodes()}
			self.num_token_dict.update({node:1 for node in frame_graph_obj_1.graph.nodes()})

		for node in frame_graph_obj_1.graph.nodes():
			try:
				if self.token_dict[node]:
					neighbour_nodes = set(frame_graph_obj_1.graph.adj[node].keys()).intersection(set(frame_graph_obj_2.graph.nodes()))
					if len(neighbour_nodes):
						neighbour = np.random.choice(list(neighbour_nodes))
						try:
							self.token_dict[neighbour] += 1
						except:
							self.token_dict[neighbour] = 1
						try:
							self.num_token_dict[neighbour] += 1
						except:
							self.num_token_dict[neighbour] = 1
						self.token_dict[node] -= 1
			except KeyError as e:
				self.token_dict[node] = 0
				self.num_token_dict[node] = 0

		frame_reachable_dict = {k:self.num_token_dict[k] for k in frame_graph_obj_2.graph.nodes()}

		return self.num_token_dict, frame_reachable_dict

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

	def drawGraph3D(self, graphObject, color, scale, draw_axes=False, draw_lines=False, fill_between=False):
		
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
		rand = np.random.rand(*node_xyz[:, 0].shape)*np.finfo(float).eps
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

		if draw_axes:
			# mlab.axes(color = (0.2, 0.2, 0.2))
			mlab.outline(color = (0.2, 0.2, 0.2))

		if fill_between:
			fill_pts = [[], [], []]
			for pair in list(graphObject.graph.edges()):
				point1 = graphObject.positionDict[pair[0]]
				point2 = graphObject.positionDict[pair[1]]
				dist = np.sqrt(np.sum(np.array(point1 - np.array(point2))**2))
				if dist > 0.2:
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


	def computeStatistics(self, frame_graph):
		connectedComponents = sorted(nx.connected_components(frame_graph), key=len, reverse=True)
		# nNodes = len(frame_graph)
		density = nx.density(frame_graph)
		# averageDegree = sum(dict(nx.degree(frame_graph)).values())/(len(frame_graph) + np.finfo(float).eps)
		averagePath = sum([nx.average_shortest_path_length(frame_graph.subgraph(c).copy()) for c in connectedComponents])/(len(connectedComponents) + np.finfo(float).eps)
		eff = global_efficiency(frame_graph)
		assortativity = nx.degree_assortativity_coefficient(frame_graph)
		# largestCC = len(max(nx.connected_components(frame_graph), key=len))
		# nCC = len(connectedComponents)
		# betweenCentrality = sum(nx.betweenness_centrality(frame_graph).values())/(len(frame_graph) + np.finfo(float).eps)
		frame_labels = ['Density', 'Efficiency', 'Assortativity']
		return(frame_labels, np.array([density, eff, assortativity]))

	def computeNodeStatistics(self, frame_graph):

		for node in frame_graph.graph.nodes():
			self.n_neighbours[node].extend(list(frame_graph.graph.edges(node)))
			self.n_unique_neighbours[node].update(frame_graph.graph.edges(node))

		return(self.n_neighbours, self.n_unique_neighbours)

	def computeSimilarity(self, frame_graph_obj_1, frame_graph_obj_2, delta):
		frame_graph_1, positionDict_1 = frame_graph_obj_1.graph, frame_graph_obj_1.positionDict
		frame_graph_2, positionDict_2 = frame_graph_obj_2.graph, frame_graph_obj_2.positionDict
		frame_1_degree = dict(nx.degree(frame_graph_1))
		frame_2_degree = dict(nx.degree(frame_graph_2))
		frame_1_neighbour_degree = dict(nx.average_neighbor_degree(frame_graph_1))
		frame_2_neighbour_degree = dict(nx.average_neighbor_degree(frame_graph_2))
		connectedComponents1 = sorted(nx.connected_components(frame_graph_1), key=len, reverse=True)
		connectedComponents2 = sorted(nx.connected_components(frame_graph_2), key=len, reverse=True)
		frame_graph_1_edges = frame_graph_1.edges()
		frame_graph_2_edges = frame_graph_2.edges()
		# print(len(list(set(frame_graph_1_edges))), len(list(set(frame_graph_2_edges))), len(list(set(frame_graph_1_edges) & set(frame_graph_2_edges))))
		temporalIntersection = (len(list(set(frame_graph_1_edges) & set(frame_graph_2_edges))) + 1)/(len(list(set(frame_graph_1_edges))) + 1)
		MeanEuclideanDistanceNode = sum([np.sqrt(np.sum((np.array(positionDict_1[node]) - np.array(positionDict_2[node]))**2)) \
			for node in list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))])/(len(list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))) + np.finfo(float).eps)
		MeanDegreeDiffNode = sum([np.abs(frame_1_degree[node] - frame_2_degree[node]) \
			for node in list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))])/(len(list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))) + np.finfo(float).eps)
		MeanNeighbourDegreeDiffNode = sum([np.abs(frame_1_neighbour_degree[node] - frame_2_neighbour_degree[node]) \
			for node in list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))])/(len(list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))) + np.finfo(float).eps)
		frame_labels = ['Temporal Intersection', 
		'Mean Euclidean Distance Per Subnode',
		'Mean Degree Difference Per Subnode',
		'Mean Neighbour Degree Difference Per Subnode']
		return(frame_labels, np.array([temporalIntersection, MeanEuclideanDistanceNode, MeanDegreeDiffNode, MeanNeighbourDegreeDiffNode]))

	def computeFlux(self, frame_graph_obj_1, frame_graph_obj_2, delta):
		frame_graph_1, intensityDict_1, widthDict_1 = frame_graph_obj_1.graph, frame_graph_obj_1.intensityDict, frame_graph_obj_1.widthDict
		frame_graph_2, intensityDict_2, widthDict_2 = frame_graph_obj_2.graph, frame_graph_obj_2.intensityDict, frame_graph_obj_2.widthDict

		MeanIntensityDiffNode = sum([(intensityDict_1[node] - intensityDict_2[node]) \
			for node in list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))])/(len(list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))) + np.finfo(float).eps)
		MeanWidthDiffNode = sum([widthDict_1[node] - widthDict_2[node] \
			for node in list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))])/(len(list(set(frame_graph_1.nodes()) & set(frame_graph_2.nodes()))) + np.finfo(float).eps)
		MeanIntensityNode = sum([(intensityDict_1[node]) \
			for node in list(set(frame_graph_1.nodes()))])/(len(list(set(frame_graph_1.nodes()))) + np.finfo(float).eps)
		MeanWidthNode = sum([(widthDict_1[node]) \
			for node in list(set(frame_graph_1.nodes()))])/(len(list(set(frame_graph_1.nodes()))) + np.finfo(float).eps)
		frame_labels = ['Mean Width Difference Per Subnode, Delta = {}'.format(delta), 'Log of Mean Intensity Difference Per Subnode, Delta = {}'.format(delta)]
		return (frame_labels, np.array([MeanWidthDiffNode, np.log(np.abs(MeanIntensityDiffNode))]))

	def computeFragmentStatistics(self, frame_graph_obj_1, frame_graph_obj_2, delta=1, fragmentThresh=5):
		if not frame_graph_obj_1.isConnected:
			tracked_components = self.trackFragmentsBipartite(frame_graph_obj_1, frame_graph_obj_2)
		else:
			tracked_components = {0:(frame_graph_obj_1, frame_graph_obj_2)}
		MeanRGArray = []
		MeanDisplacementArray = []
		TotalDistanceArray= []
		for fragment in tracked_components.values():
			if fragment[0].length >= fragmentThresh and fragment[1].length >= fragmentThresh:
				fragment1PositionArray = np.array(list(fragment[0].positionDict.values()))
				fragment2PositionArray = np.array(list(fragment[1].positionDict.values()))
				fragment1MassArray = np.array(list(fragment[0].intensityDict.values()))
				fragment1Centroid = np.nanmean(fragment1PositionArray, axis=0)
				fragment2Centroid = np.nanmean(fragment2PositionArray, axis=0)
				MeanDisplacementArray.append(np.sqrt(np.sum(fragment1Centroid - fragment2Centroid)**2))
				TotalDistanceArray.append(np.mean([np.sqrt(np.sum((fragment[0].positionDict[node] + \
					fragment2Centroid - fragment1Centroid) - fragment[1].positionDict[node])**2) \
					for node in set(fragment[0].graph.nodes()).intersection(set(fragment[1].graph.nodes()))]))
				fragment1RG = np.sqrt(np.sum(((fragment1PositionArray - fragment1Centroid)**2).reshape(3, -1))/fragment[0].length)
				MeanRGArray.append(fragment1RG)
		frame_labels = ['Radius of Gyration', 'Mean Displacement per Fragment', 
		'Mean Fragment Intensity', 'Mean Fragment Width',
		'Density', 'Average Degree', 'Average Shortest Path Length',
		'Network Efficiency', 'Fragment Length', 'Temporal Intersection']
		MeanIntensityNode = sum([(fragment[0].intensityDict[node]) \
			for node in list(set(fragment[0].graph.nodes()))])/(len(list(set(fragment[0].graph.nodes()))) + np.finfo(float).eps)
		MeanWidthNode = sum([(fragment[0].widthDict[node]) \
			for node in list(set(fragment[0].graph.nodes()))])/(len(list(set(fragment[0].graph.nodes()))) + np.finfo(float).eps)
		density = nx.density(fragment[0].graph)
		averageDegree = sum(dict(nx.degree(fragment[0].graph)).values())/(len(fragment[0].graph) + np.finfo(float).eps)
		connectedComponents = sorted(nx.connected_components(fragment[0].graph), key=len, reverse=True)
		averagePath = sum([nx.average_shortest_path_length(fragment[0].graph.subgraph(c).copy()) for c in connectedComponents])/(len(connectedComponents) + np.finfo(float).eps)
		eff = global_efficiency(fragment[0].graph)
		temporalIntersection = (len(list(set(fragment[0].graph) & set(fragment[1].graph))) + 1)/(len(list(set(fragment[0].graph))) + 1)
		# assortativity = nx.degree_assortativity_coefficient(fragment[0].graph)
		return (frame_labels, np.array([np.nanmean(MeanRGArray), np.nanmean(MeanDisplacementArray), 
			MeanIntensityNode, MeanWidthNode, density, averageDegree, averagePath,
			eff, fragment[0].length, temporalIntersection]))

	def computeCommunities(self, frame_graph_obj_1):
		# not useful on disconnected graphs
		communities = nx.algorithms.community.greedy_modularity_communities(frame_graph_obj_1.graph)
		communities = [(GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(c)).copy(), frame_graph_obj_1)) for c in communities]
		color = iter(cm.gist_rainbow(np.linspace(0, 1, len(communities))))
		for item in communities:
			c = next(color).reshape(1,-1)
			self.drawGraph(item, color = c, withLabels=False, scale=0.01)
		print(communities)

	def computeDiffusionMetrics(self, frame_graph_obj_1, scale=5e-3):
		# Maybe compare with a baseline distribution, not uniform - these distributions are really far from uniform anyway
		means = np.array([(x, y, z) for (x, y, z) in frame_graph_obj_1.positionDict.values()])
		cs = np.array(list(frame_graph_obj_1.intensityDict.values()))
		covariances = np.array([np.eye(3)*intensity*scale/cs.max() if intensity > 0 else np.eye(3)*cs.mean()*scale/cs.max() for intensity in frame_graph_obj_1.intensityDict.values()])
		# x, y = np.mgrid[means[:, 0].min():means[:, 0].max():1, means[:, 1].min():means[:, 1].max():1]
		# print(frame_graph_obj_1.graph.nodes())
		x, y, z = np.mgrid[means[:, 0].min():means[:, 0].max():5, means[:, 1].min():means[:, 1].max():5, means[:, 2].min():means[:, 2].max():5]
		pos = np.stack((x, y, z), axis=-1)
		# print(pos.shape)
		mixture = np.zeros(pos.shape[:3])
		for i in range(means.shape[0]):
			normal = multivariate_normal(mean=means[i], cov=covariances[i])
			mixture = mixture + normal.pdf(pos)
		mixture = mixture/np.sum(mixture)
		# plt.contourf(x, y, mixture)
		# plt.colorbar()
		# # plt.hist(mixture.ravel(), bins=20)
		# plt.show()
		base = (np.ones(mixture.shape)/mixture.ravel().shape[0]).ravel().shape[0]
		# print(entropy(mixture.ravel(), base=base), entropy((np.ones(mixture.shape)/mixture.ravel().shape[0]).ravel(), base=base))
		frame_labels = ['Entropy']
		return (frame_labels, np.array([entropy(mixture.ravel(), base=base)]))

	def trackFragmentsLinearAssignment(self, frame_graph_obj_1, frame_graph_obj_2):
		connectedComponents1 = [c for c in sorted(nx.connected_components(frame_graph_obj_1.graph), key=len, reverse=True)]
		connectedComponents2 = [c for c in sorted(nx.connected_components(frame_graph_obj_2.graph), key=len, reverse=True)]
		try:
			connectedComponents1 = list(np.array(connectedComponents1[:len(self.col_ind_prev)])[self.col_ind_prev])
			minLength = min(len(connectedComponents1), len(connectedComponents2), len(self.col_ind_prev))
		except AttributeError:
			minLength = min(len(connectedComponents1), len(connectedComponents2))
		cost_matrix = np.array([[-len(c1 & c2)/max(len(c1), len(c2)) for c2 in connectedComponents2[:minLength]] for c1 in connectedComponents1[:minLength]])
		row_ind, col_ind = linear_sum_assignment(cost_matrix)
		self.col_ind_prev = col_ind
		return list(zip(np.array([GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(c)).copy(), frame_graph_obj_1) 
			for c in connectedComponents1], dtype = object)[row_ind], 
			np.array([GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(c)).copy(), frame_graph_obj_2)
			 for c in connectedComponents2], dtype = object)[col_ind]))

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
				# print(index_list)
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
		# print(CC)
		self.id_dict = {'A' + f[1:]:k for k, v in CC.items() for f in v if f[0]=='B'}
		A = {k:set().union(*[connectedComponents1[n] for n in CC[k] if n[0] == 'A']) for k in sorted(CC.keys())}
		B = {k:set().union(*[connectedComponents2[n] for n in CC[k] if n[0] == 'B']) for k in sorted(CC.keys())}
		self.consistent_graph = {'A' + k[1:]:v for k,v in connectedComponents2.items()}
		# print(self.consistent_graph)
		return {k:(GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(A.get(k))).copy(), frame_graph_obj_1), 
			GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(B.get(k))).copy(), frame_graph_obj_2)) for k, v in sorted(CC.items(), key=lambda x: len(x[1]), reverse=True)}

	def trackFragmentsMaximumVote(self, frame_graph_obj_1, frame_graph_obj_2, inverted=False, threshold=5):
		connectedComponents1 = {i: c for i, c in enumerate(sorted(nx.connected_components(frame_graph_obj_1.graph), key=len, reverse=True))}
		connectedComponents2 = {i: c for i, c in enumerate(sorted(nx.connected_components(frame_graph_obj_2.graph), key=len, reverse=True))}
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
			# print(self.max_fragment_index)
		self.consistent_graph = superFragmentsDict
		superFragments2 = superFragmentsDict
		if inverted:
			return list(zip(np.array([GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(c)).copy(), frame_graph_obj_2)
			 for c in superFragments2], dtype = object), np.array([GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(c)).copy(), frame_graph_obj_1) 
			for c in superFragments1], dtype = object)))
		return {k:(GraphObject.createFromSubgraph(frame_graph_obj_1.graph.subgraph(list(superFragments1.get(k))).copy(), frame_graph_obj_1), 
			GraphObject.createFromSubgraph(frame_graph_obj_2.graph.subgraph(list(superFragments2.get(k))).copy(), frame_graph_obj_2)) for k, v in superFragments2.items()}

	def plotDynamicStats(self, frame_range, seed = 0, return_ = False):
		np.random.seed(seed)
		for frame in tqdm(frame_range):
			FrameGraphObject = self.createFrameGraph(frame=frame)
			frame_graph, positionDict = FrameGraphObject.graph, FrameGraphObject.positionDict
			labels, frame_statistics = self.computeStatistics(frame_graph)
			try:
				self.dynamicStatistics = np.concatenate((self.dynamicStatistics, np.expand_dims(frame_statistics, -1)), axis=1)
			except Exception as e:
				self.dynamicStatistics = np.expand_dims(frame_statistics, -1)

		if return_:
			return self.dynamicStatistics, labels

		colors = cm.rainbow(np.random.rand(len(labels)))

		for i, label in enumerate(labels):
			plt.subplot(len(labels)//2, 2, i+1)
			plt.title(labels[i])
			plt.plot(frame_range, self.dynamicStatistics[i], color=colors[i])

	def plotTemporalSimilarity(self, frame_range, delta = 1, seed = 0, return_ = False, color = 'r', cumulative=False):
		np.random.seed(seed)
		for frame in tqdm(frame_range):
			FrameGraphObject = self.createFrameGraph(frame=frame)
			frame_graph, positionDict = FrameGraphObject.graph, FrameGraphObject.positionDict
			if cumulative:
				FrameGraphObjectPrev = self.createFrameGraph(frame=frame_range[0])
				frame_graph_prev, positionDict_prev = FrameGraphObjectPrev.graph, FrameGraphObjectPrev.positionDict
			else:
				FrameGraphObjectPrev = self.createFrameGraph(frame=frame-delta)
				frame_graph_prev, positionDict_prev = FrameGraphObjectPrev.graph, FrameGraphObjectPrev.positionDict
			labels, frame_statistics = self.computeSimilarity(FrameGraphObject, 
				FrameGraphObjectPrev, delta=delta)
			try:
				self.dynamicSimilarity = np.concatenate((self.dynamicSimilarity, np.expand_dims(frame_statistics, -1)), axis=1)
			except Exception as e:
				self.dynamicSimilarity = np.expand_dims(frame_statistics, -1)

		if return_:
			return self.dynamicSimilarity, labels

		for i, label in enumerate(labels):
			plt.title(labels[i])
			plt.plot(frame_range, self.dynamicSimilarity[i].flatten(), color=color)

	def plotFragmentStatistics(self, frame_range, delta = 1, seed = 0, return_ = False, color = 'r', cumulative=False):
		np.random.seed(seed)

		# You can only plot this at a Delta of 1, since you need continuous timestep
		# from one period to the next, i.e., (0, 10) and (10, 20)
		assert delta==1
		for frame in tqdm(frame_range):
			FrameGraphObject = self.createFrameGraph(frame=frame)
			frame_graph, positionDict = FrameGraphObject.graph, FrameGraphObject.positionDict
			if cumulative:
				FrameGraphObjectPrev = self.createFrameGraph(frame=frame_range[0])
			else:
				FrameGraphObjectPrev = self.createFrameGraph(frame=frame-delta)
			labels, frame_statistics = self.computeFragmentStatistics(FrameGraphObjectPrev, 
				FrameGraphObject, delta=delta, fragmentThresh=5)
			try:
				self.fragmentStatistics = np.concatenate((self.fragmentStatistics, np.expand_dims(frame_statistics, -1)), axis=1)
			except Exception as e:
				self.fragmentStatistics = np.expand_dims(frame_statistics, -1)

		if return_:
			return self.fragmentStatistics, labels

		for i, label in enumerate(labels):
			plt.title(labels[i])
			plt.plot(frame_range, self.fragmentStatistics[i].flatten(), color=color)