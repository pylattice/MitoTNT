# Visualization
**In this section we will visualize the tracked mitochondrial networks in ChimeraX**

**Please first download [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)**

We need to specify the directory to save visualization files

- `vis_dir`: store `.cxc` commands to load in ChimeraX. You can use `work_dir/chimerax_visualization/` for example.

- `vis_data_dir`: store `.cmap ` and `.bild` files created for each frame and used for visualization. You can use `vis_dir/data/` for example.

## 1. Create ChimeraX rendering of the skeleton (optional but recommended)
We can use MitoGraph-generated `*skeleton.vtk` files for visualizing skeleton, but this is not ideal because it has fixed width and color.\
Alternatively here, we can render the skeleton using BILD format in ChimeraX. This allows us to set the skeleton sizes, node sizes and color. However, it also takes much longer to load in ChimeraX.
- `skeleton_colors`: a list of two colors to render for current and next frames. We use blue for current frame and red for next frame. [See more colors.](https://www.cgl.ucsf.edu/chimerax/docs/user/commands/colornames.html)
- `skeleton_size`: diameter of the skeleton that connects nodes.
- `node_size`: diameter of the spheres that make up the nodes. If `node_size`= `skeleton_size`, the nodes are not visible (but still required to fill the gaps along the skeletons).
```
tracking_visualization.generate_chimerax_skeleton()
```

## 2. Create ChimeraX rendering of tracking vectors
We will use the frame-to-frame node assignments to draw the tracking vectors for two frames.

- `arrow_color`: color of the tracking arrows. Default to black.

- `arrow_size`: diameter of the arrow stem.
```
tracking_visualization.generate_tracking_arrows()
```
## 3. Visualize network tracking in ChimeraX
Now we can combine the visualization files created above to visualize the tracking of timelapse data.
- `show_tif`: if true include fluorescence cloud in background. 

- `voxel_size`: voxel_size same as that used for MitoGraph segmentation, in the order of x,y,z.  
For example, `voxel_size='0.2,0.2,0.4'` refers to lateral pixel size 0.2 μm and axial pixel size 0.4 μm.

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
