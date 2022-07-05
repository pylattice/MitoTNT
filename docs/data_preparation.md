# Data Preparation
**To perform tracking, mitochondria 3D image stacks needs to be segmented using MitoGraph for all time points.**

## 1. Install MitoGraph
MitoGraph can be installed here:
https://github.com/vianamp/MitoGraph

## 2. Save 3D image stacks in individual directories
If your data is 4D image stacks, you will need to save it as 3D image stacks for individual timepoints. This is usually done in ImageJ or Python.  
Each 3D image stack needs to be placed in its own folder. This can be easily done using Bash commands.  

Example directory structure:
- frame_0/frame_0.tif

- frame_1/frame_1.tif  
  ......
 
## 3. Run mitochondria segmentation
To run MitoGraph on **command line**:  
`MitoGraph -xy lateral_pixel_size -z axial_pixel_size -adaptive 10 -path tif_dir`  
if adaptive thresholding is desired

To segment multiple z-stacks for a range of timepoints simultaneously:  
`for frame in frame*; do -xy lateral_pixel_size -z axial_pixel_size -adaptive 10 -path "$frame"; done`


Segmentation outputs will be saved under each timepoint folder.

Please find the description for the output files in [MitoGraph doc](https://github.com/vianamp/MitoGraph/).

## 4. Check segmentation results
We recommend using [ChimeraX](https://www.cgl.ucsf.edu/chimerax/)  to open `.tif` and `.vtk` files for visual inspection of the segmentation results.

More information for using MitoGraph and selecting parameters can be found [here](https://doi.org/10.1016/j.ab.2018.02.022).
