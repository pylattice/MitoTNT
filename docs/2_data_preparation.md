# Data Preparation
**To perform tracking, mitochondria 3D image stacks needs to be segmented using MitoGraph for all time points.**
Example dataset has been provided under `test_data/mitograph`, which are MitoGraph-segmented. You can skip this section if using the example dataset.

## 1. Cell segmentation
MitoGraph and MitoTNT usually work with mitochondria in a single cell. This helps to avoid tracking mitochondria across different cells and reduces the computational time significantly. For this reason, single cell z-stacks need to be prepared. We have included a sample cell segmentation scripts `crop_image_using_ImageJ_ROI` under `helper_scripts` directory, where it can read ImageJ ROI files to extract single cell patches for each timepoint. The users need to draw rectangular boxes in ImageJ and save as `.roi` selection files under a ROI folder for each movie.

## 2. Save 3D image stacks in individual directories
If your data is 4D image stacks, you will need to save it into 3D image stacks for individual timepoints.
Each 3D image stack also needs to be placed in its own folder.
We have provided a utility script `convert_to_tiff.ipynb ` under `helper_scripts` directory for this purpose.

Example directory structure:

- frame_0/frame_0.tif

- frame_1/frame_1.tif  
  ......
 
## 3. Run mitochondria segmentation
To run MitoGraph on **command line** for one snapshot:  
`MitoGraph -xy lateral_pixel_size -z axial_pixel_size -path tif_dir`  
Note: you need to replace `lateral_pixel_size` and `axial_pixel_size` with microscope parameters, otherwise the segmentation is empty

To segment multiple z-stacks for a number of timepoints simultaneously:  
`for frame in frame*; do MitoGraph -xy lateral_pixel_size -z axial_pixel_size -path "$frame"; done`

Segmentation outputs will be saved under each timepoint folder.

Please find the description for the output files [here](https://github.com/vianamp/MitoGraph/#mitograph-outputs).

## 4. Check segmentation results
We recommend using [ChimeraX](https://www.cgl.ucsf.edu/chimerax/) or [ParaView](https://www.paraview.org/) to open `.tif` and `.vtk` files for visual inspection of the segmentation results.

More information for using MitoGraph and parameters can be found in the publication [here](https://doi.org/10.1016/j.ab.2018.02.022).
