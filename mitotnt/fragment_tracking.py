import pandas as pd 
import numpy as np
import math
import os
from collections import defaultdict
from itertools import product
import networkx as nx
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm
from matplotlib.collections import PolyCollection
from matplotlib.legend_handler import HandlerTuple	
from matplotlib.colors import to_rgb
from matplotlib import rcParams
import matplotlib.font_manager as font_manager
import glob
import imageio
import cv2
import time
from mayavi import mlab
import matplotlib as mpl
from tqdm import tqdm 
import json
import re
from scipy.optimize import linear_sum_assignment
from scipy.stats import pearsonr
from .load_graph import *
from .utils import *
import vtk
import seaborn as sns
from networkx.algorithms import community, global_efficiency
# vtk.vtkObject.GlobalWarningDisplayOff() 

def plotCascade(rootpath, start_frame=0, timesteps=1, delta=1, frag_range=None, 
	scale=1, plot_metrics=False, savepath='', mode='', show=False, view_params=None, get_frag_ids=True):

	'''
	Create the mitochondrial network for each frame and then plot it in 3D

	Parameters
	----------
	rootpath: str
		The path of the root data directory
	cell_id: str
		The ID of the cell to process
	timesteps: int, optional
		The number of timesteps to process
	n_fragments: int, optional
		The number of fragments to plot
	start_frame: int, optional
		The frame at which to start plotting
	delta: int, optional
		The frame skip
	scale: float, optional
		The scale of the nodes in the visualization
	plot_metrics: bool, optional
		Whether to compute and save metrics
	savepath: str, optional
		The directory to save images
	mode: str, optional
		Ignore
	'''
	start_time = time.time()
	GLG = GraphLoader(rootpath, 'unique')
	metricListFrame = []
	metricListFragment = []
	fragmentMetrics = defaultdict(lambda: defaultdict(list)) if mode =='classify' else []

	graph_object_frame_0 = GLG.createFrameGraph(frame=0)
	node_reach_dict_time = []

	if not os.path.exists(savepath):
		os.makedirs(savepath)

	for t in tqdm(range(timesteps)):
		graph_object_frame_1 = GLG.createFrameGraph(frame=start_frame+(t)*delta)
		graph_object_frame_2 = GLG.createFrameGraph(frame=start_frame+(t+1)*delta)
		tracked_components = GLG.trackFragmentsBipartite(graph_object_frame_1, graph_object_frame_2)
		if t == 0:
			tracked_reference = tracked_components

		if get_frag_ids:
			fragmentMetrics = saveIndividualFragmentMetrics(GLG, tracked_components, t,
															fragmentMetrics = fragmentMetrics, mode=mode)

		SaveForGIF3D(GLG, graph_object_frame_2, tracked_components, t, 
					 frag_range=frag_range, scale=scale, savepath=savepath, tracked_reference=tracked_reference, show=show, view_params=view_params)
	
	metricDf = pd.DataFrame(fragmentMetrics, columns = ['frame_id', 'unique_frag_id', 'frame_frag_id', 'fragment_unique_node_ids', 'weighted_centroid_x', 'weighted_centroid_x', 'weighted_centroid_x'])

	metricDf.to_csv(savepath+'fragment_tracks.csv', ignore_index=True)
		
def plotSnapshot(rootpath, start_frame=0, timesteps=1, delta=1, frag_range=None, 
	scale=1, plot_metrics=False, savepath='', mode='classify', show=False, view_params=None):

	'''
	Create the mitochondrial network for each frame and then plot it in 3D

	Parameters
	----------
	rootpath: str
		The path of the root data directory
	cell_id: str
		The ID of the cell to process
	timesteps: int, optional
		The number of timesteps to process
	n_fragments: int, optional
		The number of fragments to plot
	start_frame: int, optional
		The frame at which to start plotting
	delta: int, optional
		The frame skip
	scale: float, optional
		The scale of the nodes in the visualization
	plot_metrics: bool, optional
		Whether to compute and save metrics
	savepath: str, optional
		The directory to save images
	mode: str, optional
		Ignore
	'''
	
	GLG = GraphLoader(rootpath, 'unique')
	metricListFrame = []
	metricListFragment = []
	fragmentMetrics = defaultdict(lambda: defaultdict(list)) if mode =='classify' else []

	graph_object_frame_0 = GLG.createFrameGraph(frame=0)
	node_reach_dict_time = []

	if not os.path.exists(savepath):
		os.makedirs(savepath)

	graph_object_frame_1 = GLG.createFrameGraph(frame=start_frame+(timesteps)*delta)
	graph_object_frame_2 = GLG.createFrameGraph(frame=start_frame+(timesteps+1)*delta)
	tracked_components = GLG.trackFragmentsBipartite(graph_object_frame_1, graph_object_frame_2)

	tracked_reference = tracked_components

	SaveForGIF3D(GLG, graph_object_frame_2, tracked_components, timesteps, 
				 frag_range=frag_range, scale=scale, savepath=savepath, tracked_reference=tracked_reference, show=show, view_params=view_params)

def SaveForGIF3D(G, frame_graph, tracked_components, timestep, frag_range, scale, 
				 savepath, tracked_reference=None, fill_between=True, show=False, view_params=None):

	'''
	Save the 3D fragment plot using Mayavi

	Parameters
	----------
	G: GraphLoader
		The GraphLoader parent object
	frame_graph: GraphObject
		The GraphObject data structure for the given frame
	tracked_components: int, optional
		The number of timesteps to process
	n_fragments: int, optional
		The number of fragments to plot
	frag_range: python range object
		Utility parameter
	savepath: str
		The directory to save images
	tracked_reference: GraphObject, optional
		reference GraphObject
	'''

	np.random.seed(7)
	mlab.figure(1, bgcolor=(1, 1, 1), size=(600, 400))

	unicolor = (0.9, 0.5, 0.5)
	color = [(np.random.rand(), np.random.rand(), np.random.rand()) for n in np.linspace(0, 1, max(tracked_components.keys()))]
	for frag_index in frag_range:
		try:
			G.drawGraph3D(tracked_components[frag_index][1], color=color[frag_index], scale=scale, fill_between=fill_between)
		except Exception:
			pass

	if not os.path.exists(savepath):
		os.makedirs(savepath)

	if view_params is None:
		mlab.view()
	else:
		mlab.view(view_params[0], view_params[1], view_params[2])
	mlab.savefig(savepath + 'frame_{}.png'.format(timestep), magnification=2)
	mlab.close(all=True)

	if show:
		vis = plt.imread(savepath + 'frame_{}.png'.format(timestep))
		plt.axis('off')
		plt.imshow(vis)

def saveIndividualFragmentMetrics(G, tracked_components, timestep, fragmentMetrics, mode='classify'):
	labels = ['Radius of Gyration', 'Mean Displacement per Fragment', 
			  'Mean Fragment Intensity', 'Mean Fragment Width',
			  'Density', 'Average Degree', 'Average Shortest Path Length',
			  'Network Efficiency', 'Fragment Length', 'Temporal Intersection']

	# Add Fragment Temporal Intersection etc.
	if mode == 'classify':
		for frag_index in range(max(tracked_components.keys())):
			try:
				_, metrics = G.computeFragmentStatistics(tracked_components[frag_index][0], tracked_components[frag_index][1])
				for i, label in enumerate(labels):
					fragmentMetrics[frag_index][label].append(metrics[i])
			except KeyError:
				for i, label in enumerate(labels):
					fragmentMetrics[frag_index][label].append(np.nan)
			try:
				fragmentMetrics[frag_index]['Track Length'] = [np.argwhere(np.isnan(fragmentMetrics[frag_index]['Density']))[0]]
			except IndexError:
				fragmentMetrics[frag_index]['Track Length'] = [len(fragmentMetrics[frag_index]['Density'])]
	
	else:
		if timestep == 0:
			tracked_components = {k: v for k, v in sorted(tracked_components.items(), key = lambda x: x[1][0].length, reverse=True)}
			for i, frag_index in enumerate(list(tracked_components.keys())):
				if tracked_components[frag_index][0].length > 5:
					weighted_centroid = np.sum(np.array([np.array(tracked_components[frag_index][0].positionDict[k])*tracked_components[frag_index][0].intensityDict[k] for k in 
					tracked_components[frag_index][0].positionDict]), axis=0)/np.sum(np.array([tracked_components[frag_index][0].intensityDict[k] for k in 
					tracked_components[frag_index][0].positionDict]))
					fragmentMetrics.append([timestep, frag_index, i, list(tracked_components[frag_index][0].graph.nodes()), weighted_centroid[0], weighted_centroid[1], weighted_centroid[2]])
		tracked_components = {k: v for k, v in sorted(tracked_components.items(), key = lambda x: x[1][1].length, reverse=True)}
		for i, frag_index in enumerate(list(tracked_components.keys())):
			weighted_centroid = np.nansum(np.array([np.array(tracked_components[frag_index][1].positionDict[k])*tracked_components[frag_index][1].intensityDict[k] for k in 
					tracked_components[frag_index][1].positionDict]), axis=0)/np.nansum(np.array([tracked_components[frag_index][1].intensityDict[k] for k in 
					tracked_components[frag_index][1].positionDict]))
			if np.isnan(weighted_centroid).any():
				weighted_centroid = [0., 0., 0.]
			fragmentMetrics.append([timestep+1, frag_index, i, list(tracked_components[frag_index][1].graph.nodes()), weighted_centroid[0], weighted_centroid[1], weighted_centroid[2]])

	return fragmentMetrics