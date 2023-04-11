# MitoTNT:boom:: Mitochondrial Temporal Network Tracking
**MitoTNT** is a Python-based pipeline for the tracking, visualization, and dynamic analysis of 4D mitochondrial network data.

It is built upon mitochondria segmentation provided by MitoGraph, and visualization engine provided by ChimeraX.  

MitoTNT is written by Zichen (Zachary) Wang (ziw056@ucsd.edu), with the help from people in the [Johannes Schöneberg lab](https://www.schoeneberglab.org/) at UCSD.

# Installation
On github, download this [repository](https://github.com/pylattice/mitoTNT). Alternative clone the repository on terminal: 
`git clone https://github.com/pylattice/mitoTNT.git`  

### Software requirements
- **[Jupyter Notebook](https://jupyter.org/)** or get it from **[Anaconda](https://www.anaconda.com/products/distribution)**

- **[MitoGraph](https://github.com/vianamp/MitoGraph/)** for mitochondria segmentation. Available for MacOS. Windows and Linux binary can be built.

- **[ChimeraX](https://www.cgl.ucsf.edu/chimerax/)** for tracking visualization

### Python dependencies
We will create a conda environment that automatically installs all the required dependencies.
1. Open anaconda prompt on Windows or just terminal on MacOS/Linux
2. Go to the root directory of MitoTNT repository
3. Create the enviroment using the provided .yml file: `conda env create --name mitotnt --file=mitotnt_env.yml`

To use MitoTNT, first activate the environmnet we created with `conda activate mitotnt`, and then open notebook with `jupyter notebook`.

# Data Preparation
**To perform tracking, mitochondria 3D image stacks needs to be segmented using MitoGraph for all time points.**
Example dataset has been provided under `test_data/mitograph`, which are MitoGraph-segmented. You can skip this section if using the example dataset.

## 1. Install MitoGraph
MitoGraph can be installed [here](https://github.com/vianamp/MitoGraph/#how-to-install)

## 2. Save 3D image stacks in individual directories
If your data is 4D image stacks, you will need to save it into 3D image stacks for individual timepoints.
Each 3D image stack also needs to be placed in its own folder.
**We have provided a utility script `prepare_tif_for_mitograph.ipynb ` in the helper_scripts directory for this purpose.**

Example directory structure:

- frame_0/frame_0.tif

- frame_1/frame_1.tif  
  ......
 
## 3. Run mitochondria segmentation
To run MitoGraph on **command line** for one snapshot:  
`MitoGraph -xy lateral_pixel_size -z axial_pixel_size -path tif_dir`  
Note: you need to replace `lateral_pixel_size` and `axial_pixel_size` with microscope parameters, otherwise the segmentation is empty

To segment multiple z-stacks for a range of timepoints simultaneously:  
`for frame in frame*; do MitoGraph -xy lateral_pixel_size -z axial_pixel_size -path "$frame"; done`

**However, processing frames one by one can take long time. We have also provided a utility script `run_MitoGraph_parallel.ipynb` in the helper_scripts directory that can speed up the process 10-20X based on your local machine.**

Segmentation outputs will be saved under each timepoint folder.

Please find the description for the output files [here](https://github.com/vianamp/MitoGraph/#mitograph-outputs).

## 4. Check segmentation results
We recommend using [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) to open `.tif` and `.vtk` files for visual inspection of the segmentation results.

More information for using MitoGraph and selecting parameters can be found in the publication [here](https://doi.org/10.1016/j.ab.2018.02.022).
# Network Tracking
Open the terminal and type `jupyter notebook` to start Jupyter Notebook.  
Find the root directory of MitoTNT and open `mitotnt_tracking_pipeline.ipynb` for mitochondrial network tracking.

## 1. Generate inputs for tracking
**In this section we will process the raw data into a format that is used for the subsequent tracking.**

First specify the directories we will use:

- `work_dir`: the directory where data will be processed and stored. For test data, you can use the directory of `test_data` on your machine

- `data_dir`: the directory where MitoGraph segmented mitochondria is stored. For test data, this is `test_data/mitograph`.

- `input_dir`: the directory where the processed inputs used for tracking will be stored. For test data, you can use `test_data/tracking_inputs`. This is an empty folder that will be created.

After specifying the folders, we need to set a few parameters:

- `start_frame`, `end_frame`: the range of frames to process.

- `node_gap_size`: the number of nodes to skip when creating full-resolution graphs from mitograph `.gnet` files. Default to 0 (use all nodes).

All processed inputs will be saved as a single compressed `.npz` file in `input_dir`.

To start making inputs:
```
generate_tracking_inputs.generate()
```

## 2. Frame-to-frame tracking
**In this section we will perform node assignments for each consecutive frames.**

In addition to the directories declared above, we will create `output_dir` to store the tracking outputs. For test data, you can use `test_data/tracking_outputs`.

Additional parameters needed for frame-to-frame tracking:
- `tracking_interval`: the frame interval between the two frames to be tracked. Default to 1 (every consecutive frame).  

Distance threshold is computed as the minimum of two values:  
1) the distance to the N-th closest neighbor where N is given by `cutoff_num_neighbor`, default to 10;  
2) the frame interval times the maximum allowed speed given by `cutoff_speed`, default to `None` and estimated from the distance matrix.

- `graph_matching_depth`: the maximum level used for graph comparison. Default to 2 (usually sufficient).
- `dist_exponent`, `top_exponent`: the final cost term is given by D<sup>dist_exponent</sup> x T<sup>top_exponent</sup>, where D, T are the distance and topology costs respectively. Default both to 1 (equal weighting).


To run tracking:
```
network_tracking.frametoframe_tracking()
```

## 3. Gap closing
**In this section we attempt to merge tracks that are mistakenly terminated during frame-to-frame tracking.**

Additional parameters need to be set:

- `min_track_size`: the minimum number of frames for the tracks to be kept. Default to 5.

- `max_gap_size`: the maximum number of frames for which gap closing is allowed. Default to 3. Value of 1 indicates no gap closing.

- `block_size_factor`: values less than 1 allows the sliding block implementation of gap closing to prevent memory overflow due to large cost matrix. The size of the block is given by `block_size_factor` * total number of tracks to be closed. Default to 1 (close all tracks at once).

The final node trajectories are saved in `final_node_tracks.csv` file.
Each row is one node at one time point. 
Each column is an attribute of the given node, described below.

To use gap closing:
```
network_tracking.gap_closing()
```

## 4. Evaluate output

The final node trajectories are saved in final_node_tracks.csv file.  
Each row is one node at one time point.  
Each column is an attribute of the given node, described below.  

### Columns

- `frame_id`: frame number of the node.

- `frame_node_id`: node id at the given frame. Each frame has its own indexing.

- `unique_node_id`: node id shared by all the nodes in the same track at different frames. Each track is uniquely indexed throughout the whole trajectory. This is essetially the tracking information.

- `frame_seg_id`: segment id for all the nodes in the same segment. Each frame has its own indexing. The branching points are not assigned.

- `frame_frag_id`: fragment id for all the nodes in the same connected component. Each frame has its own indexing.

- `connected_unique_node_id`: unique_node_id for neigboring nodes in the graph. This has all the topology information.

- `x`, `y`, `z`: coordinates for the node.

- `intensity`, `width`: pixel intensity and tubular width for the node given by MitoGraph.

# Visualization
**In this section we will visualize the tracked mitochondrial networks in ChimeraX**

**Please first download [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)**

We need to specify the directory to save visualization files

- `vis_dir`: store `.cxc` commands to load in ChimeraX. You can use `work_dir/chimerax_visualization/` for example.

- `vis_data_dir`: store `.cmap ` and `.bild` files created for each frame and used for visualization. You can use `vis_dir/data/` for example.

## 1. Transform .tif to match MitoGraph coordinates (optional)
Because MitoGraph does coordinate transformation, original `.tif` files need to be transformed.
This is only needed if you want to show fluorescence cloud when visualizing tracking.

- `voxel_size`: provide the voxel_size same as inputs for MitoGraph segmentation, in the format of 'x_size,y_size,z_size'.  
For example, `voxel_size='0.2,0.2,0.4'` refers to lateral pixel size 0.2 μm and axial pixel size 0.4 μm.

```
tracking_visualization.generate_transformed_tif()
```
**After the output generate_transformed_tif.cxc file is created, load it in ChimeraX to save the transformed images.**

## 2. Create ChimeraX rendering of the skeleton (optional but recommended)
We can use MitoGraph-generated `*skeleton.vtk` files for visualizing skeleton, but this is not ideal because it has fixed width and color.\
Alternatively here, we can render the skeleton using BILD format in ChimeraX. This allows us to set the skeleton sizes, node sizes and color. However, it also takes much longer to load in ChimeraX.
- `skeleton_colors`: a list of two colors to render for current and next frames. We use blue for current frame and red for next frame. [See more colors.](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html)
- `skeleton_size`: diameter of the skeleton that connects nodes.
- `node_size`: diameter of the spheres that make up the nodes. If `node_size`= `skeleton_size`, the nodes are not visible (but still required to fill the gaps along the skeletons).
```
tracking_visualization.generate_chimerax_skeleton()
```

## 3. Create ChimeraX rendering of tracking vectors
We will use the frame-to-frame node assignments to draw the tracking vectors for two frames.

- `arrow_color`: color of the tracking arrows. Default to black.

- `arrow_size`: diameter of the arrow stem.
```
tracking_visualization.generate_tracking_arrows()
```
## 4. Visualize network tracking in ChimeraX
Now we can combine the visualization files created above to visualize the tracking of timelapse data.
- `show_tif`: if true include fluorescence cloud in background. 

- `tif_colors`: color of fluorescence cloud. See [colors](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html).

- `threshold_level`: if you want to change the contrast/thickness of the fluorescence cloud, use the [level argument](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/volume.html) in ChimeraX.

- `use_chimerax_skeleton`:  
if true use BILD format skeleton which is more flexible but slower to load;  
if false use mitograph-generated .vtk files of fixed color and size (not recommended but quicker).

- `skeleton_colors`: same as in step 2.

```
tracking_visualization.visualize_tracking()
```
**Open chimerax_visualization/visualize_tracking.cxc in ChimeraX. This may take some time. Click Home -> Backgound -> White to see it better.**
# Detect Remodeling Events
**In this section we will detect nodes that undergo fusion or fission events based on the tracking results.**

This is done using a sliding-window approach to identify nodes that undergo persistent structural changes as opposed to transient segmentation differences.  
First, the fragment indices for each node are recorded for the `half_win_size frames` before and after the current frame, to form the fragment list.  
Second, for each network edge, the fragment lists for the connected nodes are compared.  
Finally, Fission will be declared if the fragment lists before the current frame are strictly identical, as well as the fragment lists after the current frame are strictly non-overlapping. 
Since fusion events can be considered as fission events reversed in time, the opposite criterion is used for fusion detection. 

Note because of the sliding window approach:  
`start_frame` must be >= `half_win_size` and `end_frame` must be <= total number of frames - `half_win_size`

Please specify:

- `stride_size`: step size for sliding the window in number of frames. Default to 1 (to detect events happening in every frame).

- `half_win_size`: size of the half sliding window in number of frames. The higher the value the stricter the requirement for calling fusion/fission. Default to 4.

- `min_tracked_frames`: minimum number of frames that are tracked in both half window in order to declare an event. Default to 2.
To detect remodeling events:
```
detect_fusion_fission.detect()
```

The remodeling events are saved under `tracking_output/remodeling_events.csv`  
Multiple remodeling nodes located in proximity (less than 5 edges away) are grouped into a single fission/fusion site.
Columns in the output:

- `type`: fusion or fission

- `frame_id`: a single frame number for describing when the event happens

- `frame_id_before`: the frame numbers before the event for each detected node (may be different due to gap closing)

- `frame_id_after`: the frame numbers after the event for each detected node (may be different due to gap closing)

- `node_id_before`: the `frame_node_id` at corresponding frame before the event for each detected node

- `node_id_after`: the `frame_node_id` at corresponding frame after the event for each detected node

- `frag_id_before`: the `frame_frag_id` at corresponding frame before the event for each detected node

- `frag_id_after`: the `frame_frag_id` at corresponding frame after the event for each detected node

- `unique_node_id`: the `unique_node_id` for each detected node

`frame_node_id`, `frame_frag_id`, `unique_node_id` are as defined in the node tracking outputs
# Post-tracking Analysis
## Motility Measurements
**In this section we will use the tracking results to compute diffusivity at three levels of description and visualize motility in space.**

The diffusision coeffients can be computed for 1) nodes, 2) segments, 3) fragments, in the order of higher level of coarse graining.  

We need to specify the directory for save post analysis results:

- `analysis_dir`: umbrella directory for post analysis. A good choice is `work_dir+'post_analysis/'`. 

- `analy_motility_dir`: directory for saving motility measurements. This is usually a subfolder in `analysis_dir`, for example `analysis_dir+'motility/'`.  

Information needed for computing MSD vs. time delay curve:

- `max_tau`: maximum number of frames/datapoints used for linear regression. Default to 5.

Additional information for computing segment and fragment level diffusivity:

- `selected_frames`: because segment and fragments undergo constant remodeling, we need to select the frames at which the segment and fragments are evaluated. This is recommended to be separated by 2x half window size (see below). Each frame will be visualized separately.

- `half_win_size `: the time window size (frames) before and after the selected center frames for collecting track coordinates. Default to 10.

- `tracked_ratio`: the minimum ratio of tracked nodes in each segment/fragment to be qualified for calculating diffusivity. Default to 0.5.  

To compute diffusivity, use:
```
post_analysis.compute_node_diffusivity()

post_analysis.compute_segment_diffusivity()

post_analysis.compute_fragment_diffusivity()
```

The data is saved in `analy_motility_dir/*diffusivity.csv` files. Each row is a node/segment/fragment. The columns are explained below.
#### Columns

 - `center_frame_id`: selected frame for determining the segment/fragment diffusivity. N/A for node diffusivity.
 
 - `unique_node_id` as in `final_node_tracks.csv`. `seg_id`, `frag_id` are specific to `center_frame_id`.

 - `diffusivity`: slope of MSD vs. time delay curve divided by 6 to account for 3D random walk.

 - `msd`: MSD per frame, euqal to 6 x diffusivity x frame interval

 - `r_squared`: coefficient of determination for the linear regression.

 - `num_points`: number of points in MSD vs. time delay curve used for linear regression.

## Motility Visualization
**Now we can visualize the computed motility measurements by mapping diffusivity values onto corresponding mitochondrial skeleton nodes with different coloring.**

- `node_size`: the size scale of the nodes
- `selected_frames`: the frames for which the skeleton is to be used for visualization. Each frame will produce its own visualization files. All tracks that include this frame will be shown.

To color diffusivity for the tracks that include the selected frames, use:
```
post_analysis.map_node_motility_onto_surface()

post_analysis.map_segment_motility_onto_surface()

post_analysis.map_fragment_motility_onto_surface()
```
