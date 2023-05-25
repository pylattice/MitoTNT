import pandas as pd 
import numpy as np
import os
from collections import defaultdict
from mayavi import mlab
from tqdm import trange 
from .graph_utils import *


def track_from_nodes(track_dir, output_dir, 
	num_frames=5, start_frame=0, frame_step=1, 
	level='fragment', frag_range=range(100),
	voxel_size=0.2, scale=0.5, view_params=None,
	compute_fragment_metrics=False, save_image=False):

	'''
	Create the mitochondrial network for each frame and then plot it in 3D
	Parameters
	----------
	track_dir: str
		the file for node tracks
	output_dir: str
		the directory for saving output files and images
	level: str ('fragment' or 'segment')
		whether to track on fragment or segment levels based on node tracks
    voxel_size: float
        The length of voxel diagonal in the unit of coordinates
	scale: float, optional
		The scale of the nodes in the visualization
	view_params: tuple of size 3, optional
		mlab.view() parameters
	compute_fragment_metrics: bool
		whether to compute the list of graph metrics
	save_image: bool
		whether to save snapshots of tracking
	'''

	GLG=GraphLoader(track_dir, 'unique')
	fragmentMetrics=[]

	if not os.path.exists(output_dir):
		os.makedirs(output_dir)

	for t in trange(num_frames):
		graph_object_frame_1=GLG.createFrameGraph(frame=start_frame+(t)*frame_step)
		graph_object_frame_2=GLG.createFrameGraph(frame=start_frame+(t+1)*frame_step)
		
		# if level == 'fragment':
		# 	tracked_components=GLG.trackFragmentsBipartite(graph_object_frame_1, graph_object_frame_2)
		# 	if t == 0:
		# 		tracked_reference=tracked_components

		# 	fragmentMetrics=saveIndividualFragmentMetrics(GLG, tracked_components, t,
		# 												  fragmentMetrics=fragmentMetrics,
		# 												  frag_range=frag_range, output_dir=output_dir, compute_fragment_metrics=compute_fragment_metrics)

		tracked_components=GLG.trackMaximumVote(graph_object_frame_1, graph_object_frame_2, level=level)
		if t == 0:
			tracked_reference=tracked_components

		fragmentMetrics=saveIndividualFragmentMetrics(GLG, tracked_components, t,
													  fragmentMetrics=fragmentMetrics,
													  frag_range=frag_range, output_dir=output_dir, compute_fragment_metrics=compute_fragment_metrics)
		
		if save_image:
			Save3D(GLG, graph_object_frame_2, tracked_components, t, frag_range=frag_range, output_dir=output_dir, 
				   tracked_reference=tracked_reference, voxel_size=voxel_size, scale=scale, view_params=view_params)
	
	labels=['Radius of Gyration', 'Mean Displacement', 
		    'Mean Fragment Intensity', 'Mean Fragment Width',
		    'Density', 'Average Degree', 'Average Shortest Path Length',
		    'Network Efficiency', 'Fragment Size', 'Temporal Intersection']

	if not compute_fragment_metrics:
		metricDf=pd.DataFrame(fragmentMetrics, columns=['frame', 'unique_frag_id', 'frame_frag_id', 'included_unique_node_id', 'weighted_centroid.x', 'weighted_centroid.y', 'weighted_centroid.z'])
	else:
		metricDf=pd.DataFrame(fragmentMetrics, columns=['frame', 'unique_frag_id', 'frame_frag_id', 'included_unique_node_id', 'weighted_centroid.x', 'weighted_centroid.y', 'weighted_centroid.z'] + labels)
	
	metricDf.to_csv(output_dir+level+'_tracks.csv', index=False)
		
def Save3D(G, frame_graph, tracked_components, timestep, frag_range, output_dir, 
		plot_reference=True, tracked_reference=None, voxel_size=0.2, scale=0.5, draw_axes=False, fill_between=True, view_params=('auto', 'auto', 'auto')):

	'''
	Save the 3D fragment plot using Mayavi
	Parameters
	----------
	G: GraphLoader
		The GraphLoader parent object
	frame graph: GraphObject
		The GraphObject data structure for the given frame
	tracked_components: int, optional
		The number of num_frames to process
	frag_range: python range object
		the 
	output_dir: str
		The directory to save images
	tracked_reference: GraphObject, optional
		Reference GraphObject
	'''

	mlab.figure(1, bgcolor=(1, 1, 1), size=(600, 400))

	color=[(np.random.rand(), np.random.rand(), np.random.rand()) for n in np.linspace(0, 1, max(tracked_components.keys()))]
	if plot_reference:
		G.drawGraph3D(frame_graph, voxel_size=voxel_size, scale=scale, draw_axes=draw_axes, color=(0.8, 0.8, 0.8), fill_between=fill_between)
	for frag_index in frag_range:
		try:
			G.drawGraph3D(tracked_components[frag_index][1], voxel_size=voxel_size, scale=scale, color=color[frag_index], fill_between=fill_between)
		except Exception:
			pass

	mlab.view(view_params[0], view_params[1], view_params[2])
	mlab.savefig(output_dir + '/frame_{}.png'.format(timestep), magnification=2)
	mlab.close(all=True)

def saveIndividualFragmentMetrics(G, tracked_components, timestep, fragmentMetrics, frag_range, output_dir, min_fragment_size=5, compute_fragment_metrics=False):
	
	labels=['Radius of Gyration', 'Mean Displacement', 
			'Mean Fragment Intensity', 'Mean Fragment Width',
			'Density', 'Average Degree', 'Average Shortest Path Length',
			'Network Efficiency', 'Fragment Size', 'Temporal Intersection']

	if timestep==0:
		tracked_components={k: v for k, v in sorted(tracked_components.items(), key=lambda x: x[1][0].length, reverse=True)}
		for i, frag_index in enumerate(list(tracked_components.keys())):
			if tracked_components[frag_index][0].length > min_fragment_size:
				weighted_centroid=np.sum(np.array([np.array(tracked_components[frag_index][0].positionDict[k])*tracked_components[frag_index][0].intensityDict[k] for k in 
				tracked_components[frag_index][0].positionDict]), axis=0)/np.sum(np.array([tracked_components[frag_index][0].intensityDict[k] for k in 
				tracked_components[frag_index][0].positionDict]))
				fragmentMetrics.append([timestep, frag_index, i, list(tracked_components[frag_index][0].graph.nodes()), weighted_centroid[0], weighted_centroid[1], weighted_centroid[2]])

	else:
		tracked_components={k: v for k, v in sorted(tracked_components.items(), key=lambda x: x[1][1].length, reverse=True)}
		for i, frag_index in enumerate(list(tracked_components.keys())):
			weighted_centroid=np.nansum(np.array([np.array(tracked_components[frag_index][1].positionDict[k])*tracked_components[frag_index][1].intensityDict[k] for k in tracked_components[frag_index][1].positionDict]), axis=0)\
			/ np.nansum(np.array([tracked_components[frag_index][1].intensityDict[k] for k in tracked_components[frag_index][1].positionDict]))
			if np.isnan(weighted_centroid).any():
				weighted_centroid=[0., 0., 0.]
			fragmentMetrics.append([timestep+1, frag_index, i, list(tracked_components[frag_index][1].graph.nodes()), weighted_centroid[0], weighted_centroid[1], weighted_centroid[2]])

			if compute_fragment_metrics:
				try:
					_, metrics=G.computeFragmentStatistics(tracked_components[frag_index][0], tracked_components[frag_index][1], min_fragment_size)
					for i, label in enumerate(labels):
						fragmentMetrics[-1].append(metrics[i])
				except KeyError:
					for i, label in enumerate(labels):
						fragmentMetrics[-1].append(np.nan)

	return fragmentMetrics