# Detect Remodeling Events
**In this section we will detect nodes that undergo fusion or fission events based on the tracking results.**

This is done using a sliding-window approach to identify nodes that undergo persistent structural changes as opposed to transient segmentation differences.  
First, the fragment indices for each node are recorded for the half_win_size frames before and after the current frame, to form the fragment list.  
Second, for each network edge, the fragment lists for the connected nodes are compared.  
Finally, Fission will be declared if the fragment lists before the current frame are strictly identical, as well as the fragment lists after the current frame are strictly non-overlapping. 
Since fusion events can be considered as fission events reversed in time, the opposite criterion is used for fusion detection. 

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
