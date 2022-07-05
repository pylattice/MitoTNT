# Visualization
**In this section we will visualize the tracked mitochondrial networks in ChimeraX**

**Please first download [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)**

We need to specify the directory to save visualization files

- `vis_dir`: store `.cxc` commands to load in ChimeraX. You can use `work_dir/chimerax_visualization/` for example.

- `vis_data_dir`: store `.cmap ` and `.bild` files created for each frame and used for visualization. You can use `vis_dir/data/` for example.

## 1. Transform .tif to match MitoGraph coordinates (optional)
Because MitoGraph does coordinate transformation, original `.tif` files need to be transformed.  
This is only needed if you want to show fluorescence cloud when visualizing tracking.
```
tracking_visualization.generate_transformed_tif(data_dir, vis_dir, vis_data_dir,
                                                start_frame, end_frame)
```

## 2. Create ChimeraX rendering of the skeleton (optional)
We can use MitoGraph-generated `*skeleton.vtk` files for visualizing skeleton, but this is not ideal because it has fixed width and color.  
Alternatively here, we can render the skeleton using BILD format in ChimeraX. This allows us to set the skeleton sizes, node sizes and color.

- `skeleton_colors`: a list of colors to render. Typically two colors are needed in order to differentiate current and next frames.  

We use blue for current frame and red for next frame.

- `skeleton_size`: diameter of the cynlinder that connects nodes. Default to 0.2.

- `node_size`: diameter of the spheres that make up the nodes.  
If `node_size`= `skeleton_size`, the nodes are not obvious (but needed to fill the holes between skeletons). Default to 0.2.
```
tracking_visualization.generate_chimerax_skeleton()
```

## 3. Create ChimeraX rendering of tracking vectors
We will use the frame-to-frame node assignments to draw the tracking vectors for two frames.

- `arrow_color`: color of the tracking arrows. Default to black.

- `arrow_size`: diameter of the arrow head. Default to 0.3.
```
tracking_visualization.generate_tracking_arrows()
```

## 4. Visualize network tracking in ChimeraX
Now we can combine the visualization files created above to visualize the tracking of timeseries data

- `show_tif`: if true include fluorescence cloud in background

- `use_chimerax_skeleton`: # if true use BILD format skeleton which is more flexible but slower to load, if false use mitograph-generated .vtk files of fixed color and size
```
tracking_visualization.visualize_tracking()
```
