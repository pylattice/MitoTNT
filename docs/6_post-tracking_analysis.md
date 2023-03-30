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
