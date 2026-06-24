# User Guide

## Data Preparation

**Before tracking, all 3D image stacks of mitochondria need to be segmented using MitoGraph.**

### 1. Convert 4D movie to 3D stacks
If your data is a 4D hyperstack, you will need to convert it into individual 3D image stacks.
We have provided a utility script `convert_to_tiff.ipynb` under `helper_scripts` directory for this purpose.

### 2. Extract single cell 3D z-stacks for all timepoints
MitoGraph and MitoTNT usually work with mitochondria in a single cell. This helps to avoid tracking mitochondria across different cells and reduces the computational time significantly.
In order to better organize the MitoGraph outputs, we also put the z-stack for each timepoint into a separate folder.

Directory structure:

- test_data/mitograph/frame_000/frame_000.tif
- test_data/mitograph/frame_001/frame_001.tif  
  . . .
- test_data/mitograph/frame_089/frame_089.tif  
  
We have included a utility script `crop_image_using_ImageJ_ROI.ipynb` under the `helper_scripts` directory to read ImageJ ROI files to extract single cell z-stacks for each timepoint. 
The user needs to draw rectangular boxes in ImageJ and save as `.roi` selection files under a ROI folder.
 
### 3. Run mitochondria segmentation
After installing MitoGraph, confirm it by running it on command line. To run MitoGraph on all the timepoints in parallel, 
we provide a python wrapper function:
```python
mito_segmenter = mitotnt.MitoSegmenter('D:/GitHub/MitoTNT/test_data/mitograph', 
                                       xy_pxl_size=0.145, z_pxl_size=0.145)

mito_segmenter.run_mitograph_cpu(max_workers=10)
```

MitoGraph outputs will be saved under each folder. Please find the documentation for the output files [here](https://github.com/vianamp/MitoGraph/#mitograph-outputs).

### 4. Check segmentation results
We recommend using [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) or [ParaView](https://www.paraview.org/) to visually inspect the quality of segmentation against your microscopy data.

To create ChimeraX-ready scripts, run:
```python
mito_segmenter.visualize_segmented_data()
```
Load file `mitograph/check_mitograph.cxc` in ChimeraX to visualize raw data and mitograph-segmented surfaces.
MitoGraph parameters may be fine-tuned if needed. More information can be found in this publication [here](https://doi.org/10.1016/j.ab.2018.02.022).

Example 4D dataset is provided under `test_data/mitograph`, along with MitoGraph outputs that can be used for the next section.

---


## Network Tracking

### 1. Generate inputs for tracking
First we need to compile the MitoGraph outputs into a format that is convenient for the subsequent tracking. We will pass the `mito_segmenter`
object in the previous section to the `SkeletonizedMito()` and extract various graph representations. 
See *API Reference* for a list of extracted data.

There are two types of graphs, full graphs that include all voxel-level mitochondrial nodes, and 
simple graphs that substitutes degree-2 nodes on linear segments with edges (typically used for network analyses). 
To extract graphs from MitoGraph results:
```python
skeletonized_mito = mitotnt.SkeletonizedMito(mito_segmenter)
skeletonized_mito.extract_graphs()
```
You can also opt to skip nodes while contructing full graphs for faster processing:
```python
skeletonized_mito.extract_graphs(node_gap_size=2)
```

Note that the `skeletonized_mito.simple_graphs` are igraph representation of segmented mitochondrial network 
that can be used for downstream network analyses.


### 2. Mitochondrial network tracking
Now we can perform the actual network tracking!

First the algorithm wants to eliminate nodes that are too far to be even considered.
Distance threshold is computed as the minimum of two values:  
1) the distance to the N-th closest neighbor where N is given by `cutoff_num_neighbor`, default to 10;  
2) the frame interval times the maximum allowed speed given by `cutoff_speed`, default to None and estimated from the distance matrix.

Then we want to construct a cost matrix composed of both distance and topology costs:

- `graph_matching_depth`: the maximum graph depth used for graph comparison (default is 2, usually sufficient).
- `dist_exponent`, `top_exponent`: the final cost term is given by D<sup>dist_exponent</sup> × T<sup>top_exponent</sup>, 
where D, T are the distance and topology costs respectively (default is 1, equal weighting).

After frame-to-frame tracking across the entire movie, we also try to merge tracks with untracked gaps.

- `min_track_size`: the minimum number of frames required for a track to be kept (default is 4).
- `max_gap_size`: the maximum number of frames for which gap closing is allowed (default is 3).

To run tracking, specify your movie's frame interval in seconds:
```python
network_tracker = mitotnt.NetworkTracker(skeletonized_mito, frame_interval=3.3)
```

By default, all frames as in the original MitoGraph folder will be tracked, you can also choose to track
a fraction of it, e.g. skip time when scope drifts:
```python
network_tracker = mitotnt.NetworkTracker(skeletonized_mito, frame_interval=3.3, start_frame=20, end_frame=40)
```

You can also choose to skip frames while tracking to speed up processing or look at slow dynamics:
```python
network_tracker = mitotnt.NetworkTracker(skeletonized_mito, frame_interval=3.3, tracking_interval=4)
```

Now run tracking, this can take some time:
```python
tracked_mito = network_tracker.run()
```

If the data has already been written to the disk, and you want to continue the run:
```python
tracked_mito = network_tracker.reload_results()
```

### Evaluate output
The tracking outputs are returned as a `TrackedMito` object and locally saved as `mitotnt/mito_node_tracks.csv`.  
Each row is tracked node at one frame with the following columns:

- `frame_id`: the frame number.
- `frame_node_id` node id with frame-wise indexing.
- `unique_node_id`: node id shared by all the nodes in the same track at different frames. Each track is uniquely indexed throughout the whole trajectory. This contains the tracking information.
- `frame_seg_id`: segment id for each mitochondrial segment (between non-degree-2 nodes) with frame-wise indexing.
- `frame_frag_id`: fragment id for each mitochondrial fragment (connected component) with frame-wise indexing.
- `connected_unique_node_id`: space-delimited `unique_node_id` for tracked neighboring nodes in the graph. Note that the topology may be different from static graphs due to absence of untracked nodes.
- `x`, `y`, `z`: coordinates for the node
- `intensity`, `width`: pixel intensity and tubular width for the node from MitoGraph.
---

## Visualization

**In this section we will visualize the tracked mitochondrial networks in ChimeraX.**  
**Please first download [ChimeraX](https://www.cgl.ucsf.edu/chimerax/).**

Create a `Visualizer` object by feeding a `TrackedMito` object.

```python
visualization = mitotnt.Visualizer(tracked_mito)
```

### 1. Create ChimeraX rendering of the skeleton (optional but recommended)
We can use MitoGraph-generated `*skeleton.vtk` files for visualizing skeleton, but this is not ideal because it has fixed width and color.  
Alternatively here, we can render the skeleton using BILD format in ChimeraX. This allows us to set the skeleton sizes, node sizes and color. 
However, it also takes longer to load.

- `skeleton_colors`: a list of two colors to render for current and next frames (default is blue for current frame and red for next frame). [See more colors.](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html)
- `skeleton_size`: radius of the cylinder that connects nodes (default is 0.03).
- `node_size`: radius of the nodes (default is 0.03, same as `skeleton_size` and the nodes are not highlighted).
```python
visualization.generate_chimerax_skeleton()
```

### 2. Create ChimeraX rendering of tracking vectors
We will use the frame-to-frame node assignments to draw the tracking vectors for two frames.

- `arrow_color`: color of the tracking arrows ( default is 'black').
- `arrow_size`: Radius of arrow shafts (default is 0.03). Arrowheads are scaled proportionally.
```python
visualization.generate_tracking_arrows()
```

### 3. Visualize network tracking in ChimeraX
Now we can combine the visualization files created above to visualize the tracking of timelapse data.

- `show_tif`: whether to include fluorescence cloud in background (default is True).
- `tif_colors`: List of two ChimeraX color names for rendering consecutive TIFF frames (defaults to `['deep sky blue', 'tomato']`). See [colors](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html) for more options.
- `threshold_level`: Threshold setting for volume rendering (default is '', use ChimeraX defaults). See the [level argument](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/volume.html) in ChimeraX.
- `use_chimerax_skeleton`: If True, load skeletons exported in ChimeraX BILD format. If False, load original Mitograph VTK skeletons (default is True).
- `skeleton_colors`: same as previous step.

```python
visualization.visualize_tracking()
```

**Open `mitotnt/visualization/visualize_tracking.cxc` in ChimeraX. This may take some time. Click `Home -> Background -> White` enhance visibility.**

---


## Post Tracking Analyses

### 1. Detect network remodeling events
**In this section we will detect nodes that undergo fusion or fission events based on the tracking data.**

This is done using a sliding-window approach to identify nodes that undergo persistent structural changes as opposed to transient segmentation differences.  
First, the fragment indices for each node are recorded for `half_win_size` frames before and after the current frame, to form the fragment list.  
Second, for each network edge, the fragment lists for the connected nodes are compared.  
Finally, fission will be declared if the fragment lists before the current frame are strictly identical, as well as the fragment lists after the current frame are strictly non-overlapping. 
Since fusion events can be considered as fission events reversed in time, the opposite criterion is used for fusion detection.  
Multiple remodeling nodes located in proximity (less than 5 edges away) are grouped into a single fission/fusion site.  


- `tracking_interval`: step size for the sliding window in frames (default is None, inherits the value used during tracking).
- `half_win_size`: size of the half sliding window in number of frames (default is 4). The higher the value the stricter the requirement for calling fusion/fission.
- `min_tracked_frames`: minimum number of frames that are tracked in both half window in order to declare an event (default is 2).

To detect remodeling events:
```python
tracked_mito.detect_fusion_fission()
```

Results are updated to `tracked_mito.remodeling_events` as a pandas DataFrame and written to `mitotnt/remodeling_events.csv`.

Each row is a detected fusion or fission event, with the following columns (numbers are space-delimited):

- `type`: fusion or fission.
- `frame_id`: a single frame number for describing when the event happens.
- `frame_id_before`: the frame numbers before the event for each detected node (maybe different due to gap closing).
- `frame_id_after`: the frame numbers after the event for each detected node (maybe different due to gap closing).
- `node_id_before`: the `frame_node_id` at corresponding frame before the event for each detected node.
- `node_id_after`: the `frame_node_id` at corresponding frame after the event for each detected node.
- `frag_id_before`: the `frame_frag_id` at corresponding frame before the event for each detected node.
- `frag_id_after`: the `frame_frag_id` at corresponding frame after the event for each detected node.
- `unique_node_id`: the `unique_node_id` for each detected node.
`frame_node_id`, `frame_frag_id`, `unique_node_id` are as defined in the node tracking outputs.


### 2. Motility Measurements

**In this section we will use the tracking outputs to compute diffusivity at three levels in mitochondria.**

The motility can be computed for 1) nodes, 2) segments, 3) fragments, in the order of more coarse-graining.  
A node is the most fundamental tracked unit — basically a single voxel along the mitochondrial skeleton at a given frame.
A segment, or branch is a continuous edge in the skeleton graph, composed of all the nodes between branching or terminal nodes .
A fragment is a connected component of the mitochondrial network (i.e. a whole mitochondrion, whether a small blob or an entire tubular network).

Motility as measured by diffusivity coefficient is computed by tracking how far mitochondria (nodes, segments, or fragments) move over increasing time lags and calculating
the mean squared displacement (MSD). The slope of MSD versus time is then divided by 6 (for 3D motion) to estimate the diffusivity coefficient.

Parameters for computing MSD vs. time delay curve:

- `max_tau`: maximum number of frames/datapoints used for linear fitting (default is 5).

Additional parameters for computing segment and fragment level diffusivity:

- `tracked_ratio`: the minimum ratio of tracked nodes in each segment/fragment to be qualified for calculating diffusivity (default is 0.3).
- `selected_frames`: the frames that segment and fragment topology references are taken. This is needed due to the constant network remodeling through
fusion and fission, thus changing segments/fragments. 
If not given, frames are auto-selected at intervals of 2×`half_win_size` across the movie.
- `half_win_size `: the time window size (frames) before and after the selected center frames for checking topology (default is 10).

To compute diffusivity, use:
```python
tracked_mito.compute_mitochondria_diffusivity()
```

Or choose to only compute node diffusivity:
```python
tracked_mito.compute_mitochondria_diffusivity(levels=['node'])
```

Results are stored in the attributes `tracked_mito.node_diffusivity`, `tracked_mito.segment_diffusivity`, `tracked_mito.fragment_diffusivity`.
The corresponding .csv files are also saved onto disk.

Each row is a tracked node/segment/fragment with the following columns:

- `unique_node_id`: Unique node track ID from mito_node_tracks.csv (only for node_diffusivity).
- `center_frame_id`: Center frame used for the analysis (only for segment_diffusivity and fragment_diffusivity).
- `seg_id/frag_id`: Segment identifier within the frame (only for segment_diffusivity and fragment_diffusivity).
- `diffusivity`: Estimated diffusivity coefficient D.
- `msd`: Fitted mean squared displacement per frame.
- `r_squared`: Goodness of fit (R²) for the linear regression.
- `num_points`: Number of time lags included in the fit.


**The remodeling events and motility data are now yours to explore. Feel free to use them in your own research and plot figures in your desired style.**

---


## Feature Extraction

**In this section we will aggregate the per-frame results into tidy tables for statistics and plotting.**

`FeatureExtractor` collects the outputs produced in the previous steps (the extracted graphs, the diffusivity CSVs, and the remodeling events) into three feature tables. Create it from a `TrackedMito` object, or directly from a results directory:

```python
feature_extractor = mitotnt.FeatureExtractor(tracked_mito)   # or FeatureExtractor(save_path='.../mitotnt/')
```

Run feature extraction after the post-tracking analyses, since `compute_diffusivity` reads the `*_diffusivity.csv` files and `compute_remodeling_rates` reads `remodeling_events.csv`.

### 1. Graph (topology) features
Per-frame network topology computed from the extracted graphs: node/segment/fragment counts, segment lengths and widths, fragment length/diameter/tortuosity, graph density, global efficiency, clustering coefficient, betweenness, and branching ratios.
```python
feature_extractor.compute_graph_features(save_csv=True)
```
Results are stored in `feature_extractor.graph_features` and written to `mitotnt/graph_features.csv`.

### 2. Motility features
Diffusivity values gathered from the node/segment/fragment diffusivity CSVs, keeping only well-fit tracks (R² > 0.8 and at least 3 points).
```python
feature_extractor.compute_diffusivity(save_csv=True)
```
Results are stored in `feature_extractor.motility_features` and written to `mitotnt/motility_features.csv`.

### 3. Remodeling rates
Per-frame fusion and fission rates, normalized by the number of nodes in each frame.
```python
feature_extractor.compute_remodeling_rates(save_csv=True)
```
Results are stored in `feature_extractor.remodeling_features` and written to `mitotnt/remodeling_features.csv`.

### 4. Plot feature distributions
Each computed feature category can be plotted as histograms, violin plots, or box plots:
```python
feature_extractor.plot_features_as_histogram()
feature_extractor.plot_features_as_violinplot()
feature_extractor.plot_features_as_boxplot()
```
