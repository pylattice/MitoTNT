import glob
import os
import concurrent.futures
from tqdm.auto import tqdm, trange

class MitoSegmenter:
    """
    Scans a directory for subfolders matching a given pattern
    (e.g., `"frame_000"`, `"frame_001"`), stores metadata about
    the dataset, and records pixel resolution information.

    Parameters
    ----------
    data_path : str
        Path to the root directory containing frame subfolders of raw tif.
    pattern : str, optional
        Substring pattern used to identify frame folders
        (default is "frame").
    xy_pxl_size : float, optional
        Pixel size in the x–y plane, in microns (default is 1.0).
    z_pxl_size : float, optional
        Pixel size in the z axis, in microns (default is 1.0).

    Attributes
    ----------
    data_path : str
        Path to the root data directory.
    list_of_folders : list of str
        Sorted list of frame-level subfolders matching the given pattern.
    num_frames : int
        Number of frames detected in the dataset.
    xy_pxl_size : float
        Pixel resolution in the x–y plane.
    z_pxl_size : float
        Pixel resolution in the z axis.
    """

    def __init__(self, data_path: str, pattern: str = "frame", xy_pxl_size: float = 1.0, z_pxl_size: float = 1.0):

        self.data_path = data_path
        self.list_of_folders = [os.path.normpath(f) for f in sorted(glob.glob(f"{self.data_path}/*{pattern}*"))]
        self.num_frames = len(self.list_of_folders)
        self.xy_pxl_size = xy_pxl_size
        self.z_pxl_size = z_pxl_size

        if self.num_frames == 0:
            raise ValueError(f"No folders found matching '{pattern}' in {data_path!r}")


    def run_mitograph_cpu(self, overwrite: bool = False, extra_params: str = "", max_workers: int = 4):
        """
        Run MitoGraph on each frame of the dataset using CPU in parallel with a progress bar.

        Parameters
        ----------
        overwrite : bool
            Whether to overwrite existing MitoGraph outputs (default is False).
        extra_params : str
            Additional MitoGraph parameters (default is "").
        max_workers : int
            Number of threads to use in parallel (default is 4).

        Returns
        -------
        None
            The MitoGraph outputs are written to each frame subfolder containing `.tif`.
        """

        def process_folder(folder):
            if len(glob.glob(folder + '/*.txt')) > 0 and not overwrite:
                return f"{folder!r} has already been processed"
            try:
                os.system(
                    f"MitoGraph -xy {str(self.xy_pxl_size)} -z {str(self.z_pxl_size)} -path {folder} {extra_params}")
                return f"MitoGraph completed for folder {folder!r}"
            except Exception as e:
                return f"MitoGraph failed for folder {folder!r}: {e}"

        results = []
        with concurrent.futures.ThreadPoolExecutor(max_workers=max_workers) as executor:
            # Use as_completed so we can update tqdm as each finishes
            futures = {executor.submit(process_folder, folder): folder for folder in self.list_of_folders}
            for future in tqdm(concurrent.futures.as_completed(futures), total=len(futures), desc="Running MitoGraph"):
                results.append(future.result())

        # Print results in order they were submitted (not finished)
        for folder, res in zip(self.list_of_folders, results):
            print(res)


    def visualize_segmented_data(self, frames: list = None):
        """
        Create ChimeraX script to visualize the segmented mitochondrial surface and the raw microscopy data.
        This is critical to evaluate the quality of mitochondria segmentation.

        Parameters
        ----------
        frames : list of int
            The frames to visualize (default is [0, -1], the first and last frames).

        Returns
        -------
        None
            check_mitograph.cxc is saved and can be opened in ChimeraX for visualization.
        """

        if frames is None:
            frames = [0, -1]

        voxel_size = str(self.xy_pxl_size) + ',' + str(self.xy_pxl_size) + ',' + str(self.z_pxl_size)

        file_dir = self.data_path+'/check_mitograph.cxc'
        if os.path.exists(file_dir):
            os.remove(file_dir)
        script = open(file_dir, 'x')

        try:
            commands = []

            idx = 1
            commands.append('close\n')
            for frame in frames:

                # load tif
                commands.append('open \"' + sorted(glob.glob(self.list_of_folders[frame] + '/*.tif'))[0] + '\"\n')
                commands.append('volume #' + str(idx) + ' voxelSize ' + voxel_size + '\n')
                commands.append('volume #' + str(idx) + ' color white style image' + '\n')
                commands.append('volume flip #' + str(idx) + ' axis y\n')
                commands.append('close #' + str(idx) + '\n')
                commands.append('rename #' + str(idx + 1) + ' id #' + str(idx) + '\n')

                # load .vtk surface
                commands.append('open \"' + sorted(glob.glob(self.list_of_folders[frame] + '/*mitosurface.vtk'))[0] + '\"\n')

                # combine the models
                commands.append('rename #' + str(idx) + '-' + str(idx + 1) + ' id ' + str(idx) + '\n')

                idx += 1

            commands.append('view # 1 clip false' + '\n')

            if len(frames) > 1:
                commands.append('mseries slider #1-' + str(idx - 1) + '\n')

            script.writelines(commands)

        except Exception as e:
            print(f"An error occurred: {e}")
            script.close()

        print(f'Load file {self.data_path}/check_mitograph.cxc in ChimeraX to visualize raw data and mitograph-segmented surfaces.')
