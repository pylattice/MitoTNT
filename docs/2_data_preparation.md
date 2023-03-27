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

Segmentation outputs will be saved under each timepoint folder.

Please find the description for the output files [here](https://github.com/vianamp/MitoGraph/#mitograph-outputs).

## 4. Check segmentation results
We recommend using [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) to open `.tif` and `.vtk` files for visual inspection of the segmentation results.

More information for using MitoGraph and selecting parameters can be found in the publication [here](https://doi.org/10.1016/j.ab.2018.02.022).
