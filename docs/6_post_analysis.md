# Post Analysis
## Motility Measurements
**In this section we will use the tracking results to compute diffusivity at three levels of description and visualize motility in space.**

We need to specify the directory for save post analysis results:

- `analysis_dir`: umbrella directory for post analysis. A good choice is `work_dir+'post_analysis/'`. 

- `analy_motility_dir`: directory for saving motility measurements. This is usually a subfolder in `analysis_dir`, for example `analysis_dir+'motility/'`.

Information needed for computing MSD vs. time delay curve:

- `frame_interval`: frame interval for the data in seconds.

- `max_tau`: maximum number of frames/datapoints used for linear regression. Default to 5.

The diffusision coeffients can be computed for 1) nodes, 2) segments, 3) fragments, in the order of higher level of coarse graining.  
To compute diffusivity, use:
```
post_analysis.compute_node_diffusivity()

post_analysis.compute_segment_diffusivity()

post_analysis.compute_fragment_diffusivity()
```
