import numpy as np
import os
import matplotlib.pyplot as plt
import glob
from mitotnt.tracked_mito import TrackedMito

class Visualizer:
    """
    Prepare ChimeraX files for visualizing tracking.

    Parameters
    ----------
    tracked_mito : TrackedMito
        Contains segmented mitochondria and network tracking results.

    Attributes
    ----------
    tracked_mito : TrackedMito
        Reference to the input `TrackedMito` object.
    visual_path : str
        Path to the directory where visualization scripts are saved (default is `TrackedMito.save_path` + 'visualization/').
    visual_data_path : str
        Path to the directory where intermediate visualization data is saved (default is `TrackedMito.save_path` + 'visualization/data/').
    """

    def __init__(self, tracked_mito: TrackedMito):

        self.tracked_mito = tracked_mito
        self.visual_path = os.path.join(self.tracked_mito.save_path, "visualization", "")
        self.visual_data_path = os.path.join(self.tracked_mito.save_path, "visualization", "data", "")
        
        os.makedirs(self.visual_path, exist_ok=True)
        os.makedirs(self.visual_data_path, exist_ok=True)


    def generate_chimerax_skeleton(self, start_frame: int = 0, end_frame: int = 5,
                                   skeleton_colors: list = None, skeleton_size: float=0.03, node_size: float=0.03):
        """
            Export mitochondrial network skeletons to ChimeraX BILD format.

            Generates `.bild` files representing mitochondrial skeleton graphs
            for a specified frame range. Each skeleton is rendered as cylinders
            (edges) and spheres (nodes), which can be directly loaded into
            ChimeraX for 3D visualization.

            Parameters
            ----------
            start_frame : int, optional
                Index of the first frame to export (default is 0).
            end_frame : int, optional
                Index of the last frame to export (default is 5).
            skeleton_colors : list of str, optional
                List of colors (ChimeraX color names or abbreviations) used for rendering
                skeletons. Each color is applied sequentially across frames. If None,
                defaults to blue for current frame and red for next frame.
            skeleton_size : float, optional
                Cylinder radius for skeleton edges (default is 0.03).
            node_size : float, optional
                Sphere radius for nodes (default is 0.03).

            Returns
            -------
            None
                BILD files are written to the directory `self.visual_data_path`.
        """

        tracking_interval = self.tracked_mito.tracking_interval

        if skeleton_colors is None:
            skeleton_colors = ['b', 'r']
            
        full_graphs_all_frames = self.tracked_mito.full_graphs

        for color in skeleton_colors:
            for frame in range(start_frame, end_frame + tracking_interval, tracking_interval):

                frame_name = str(frame).zfill(3)
                file_dir = self.visual_data_path + 'frame_' + frame_name + '_chimerax_skeleton_' + color + '.bild'
                if os.path.exists(file_dir):
                    os.remove(file_dir)

                if frame < len(full_graphs_all_frames):

                    bild = open(file_dir, 'x')
                    commands = ['.color ' + color + '\n']

                    full_graph = full_graphs_all_frames[frame]

                    for edge in full_graph.es:
                        coord_m = _coord_to_str(full_graph.vs[edge.source]['coordinate'])
                        coord_n = _coord_to_str(full_graph.vs[edge.target]['coordinate'])
                        commands.append('.cylinder ' + coord_m + coord_n + str(skeleton_size) + '\n')

                    for node in full_graph.vs:
                        coord = _coord_to_str(node['coordinate'])
                        commands.append('.sphere ' + coord + str(node_size) + '\n')

                    bild.writelines(commands)
                    bild.close()

                else:
                    print('End frame is invalid. Using the last frame.')
                    break

        print('Generated network skeleton in BILD format')


    def generate_tracking_arrows(self, start_frame: int = 0, end_frame: int = 5,
                                 arrow_color: str='black', arrow_size: float = 0.03):
        """
        Export tracking arrows between frames to ChimeraX BILD format.

        Generates `.bild` files that visualize temporal linkages of tracked
        mitochondrial nodes as arrows between consecutive frames. Each arrow
        connects the coordinate of a node in frame `t` to its linked node in
        frame `t + tracking_interval`.

        Parameters
        ----------
        start_frame : int, optional
            Index of the first frame to export arrows for (default is 0).
        end_frame : int, optional
            Index of the last frame to export arrows for (default is 5).
        arrow_color : str, optional
            Color of the arrows (ChimeraX color name or abbreviation, default is 'black').
        arrow_size : float, optional
            Radius of arrow shafts (default is 0.03). Arrowheads are scaled proportionally.

        Returns
        -------
        None
            BILD files are written to the directory `self.visual_data_path`.
        """
        tracking_interval = self.tracked_mito.tracking_interval

        full_graphs_all_frames = self.tracked_mito.full_graphs
        linked_nodes = self.tracked_mito.linked_nodes

        for fid, frame in enumerate(range(start_frame, end_frame-tracking_interval, tracking_interval)):

            arrows = ['.color ' + arrow_color + '\n']

            frame_name = str(frame).zfill(3)
            file_dir = self.visual_data_path + 'frame_' + frame_name + '_arrows.bild'
            if os.path.exists(file_dir):
                os.remove(file_dir)

            if frame + tracking_interval < len(full_graphs_all_frames):

                bild = open(file_dir, "x")

                coords_m = full_graphs_all_frames[frame].vs['coordinate']
                coords_n = full_graphs_all_frames[frame + tracking_interval].vs['coordinate']

                # create linking vectors for frame m,n
                linked = linked_nodes[fid]
                for i in range(len(linked)):
                    start, end = linked[i, 0], linked[i, 1]
                    start_coord = coords_m[start]
                    end_coord = coords_n[end]

                    comparison = (start_coord == end_coord)
                    if comparison.all():
                        end_coord = coords_n[end] + np.random.normal(0, arrow_size, 3)

                    arrows.append('.arrow ' + _coord_to_str(start_coord) + _coord_to_str(end_coord) + str(arrow_size) + ' ' + str(arrow_size * 2) + ' 0.6\n')

                bild.writelines(arrows)
                bild.close()

            else:
                print('End frame is invalid. Using the last frame.')
                break

        print('Generated tracking arrows in BILD format')


    def visualize_tracking(self, start_frame: int = 0, end_frame: int = 5,
                           show_tif: bool = True, tif_colors = None, threshold_level = '',
                           use_chimerax_skeleton: bool = True, skeleton_colors = None):
        """
            Generate a ChimeraX command script for visualizing mitochondrial tracking.

            Creates a `.cxc` script that loads skeletons, arrows, and optional TIFF volumes
            into ChimeraX for side-by-side visualization of tracked mitochondrial
            networks. Skeletons can be displayed either as constructed ChimeraX BILD skeletons
            or as raw Mitograph VTK files.

            Parameters
            ----------
            start_frame : int, optional
                Index of the first frame to visualize (default is 0).
            end_frame : int, optional
                Index of the last frame to visualize (default is 5).
            show_tif : bool, optional
                If True, include raw TIFF volumes in the visualization (default is True).
            tif_colors : list of str, optional
                List of two ChimeraX color names for rendering consecutive TIFF frames.
                If None, defaults to ['deep sky blue', 'tomato'].
            threshold_level : str, optional
                Threshold setting for volume rendering (default is '', use ChimeraX defaults).
            use_chimerax_skeleton : bool, optional
                If True, load skeletons exported in ChimeraX BILD format.
                If False, load original Mitograph VTK skeletons (default is True).
            skeleton_colors : list of str, optional
                List of two ChimeraX color names for consecutive skeletons. If None,
                defaults to blue for current frame and red for next frame.

            Returns
            -------
            None
                A ChimeraX command script `visualize_tracking.cxc` is written to
                `self.visual_path`.
        """

        voxel_size = str(self.tracked_mito.xy_pxl_size)+','+str(self.tracked_mito.xy_pxl_size)+','+str(self.tracked_mito.z_pxl_size)

        tracking_interval = self.tracked_mito.tracking_interval

        if tif_colors is None:
            tif_colors = ['deep sky blue', 'tomato']

        if skeleton_colors is None:
            skeleton_colors = ['b', 'r']

        file_dir = self.visual_path + 'visualize_tracking.cxc'
        if os.path.exists(file_dir):
            os.remove(file_dir)
        bild = open(file_dir, 'x')

        commands = ['close\n']
        idx = 1

        for frame in range(start_frame, end_frame-tracking_interval, tracking_interval):

            if frame + tracking_interval < len(self.tracked_mito.full_graphs):

                frame_m = str(frame).zfill(3)
                frame_n = str(frame + tracking_interval).zfill(3)

                # load arrow
                commands.append('open \"' + self.visual_data_path + 'frame_' + frame_m + '_arrows.bild\"\n')

                if use_chimerax_skeleton:
                    # load chimerax skeleton
                    commands.append(
                        'open \"' + self.visual_data_path + 'frame_' + frame_m + '_chimerax_skeleton_' + skeleton_colors[0] + '.bild\"' + '\n')
                    commands.append(
                        'open \"' + self.visual_data_path + 'frame_' + frame_n + '_chimerax_skeleton_' + skeleton_colors[1] + '.bild\"' + '\n')

                else:
                    # load mitograph skeleton
                    commands.append(
                        'open \"' + sorted(glob.glob(self.tracked_mito.list_of_folders[int(frame_m)]+'/*skeleton.vtk'))[0] + '\"\n')
                    commands.append(
                        'open \"' + sorted(glob.glob(self.tracked_mito.list_of_folders[int(frame_n)]+'/*skeleton.vtk'))[0] + '\"\n')

                # load tif
                if show_tif:
                    # frame_m
                    commands.append('open \"' + sorted(glob.glob(self.tracked_mito.list_of_folders[int(frame_m)]+'/*.tif'))[0] + '\"\n')
                    commands.append('volume #' + str(idx + 3) + ' voxelSize ' + voxel_size + '\n')
                    commands.append('volume flip #' + str(idx + 3) + ' axis y\n')
                    commands.append('close #' + str(idx + 3) + '\n')
                    commands.append('rename #' + str(idx + 4) + ' id #' + str(idx + 3) + '\n')
                    commands.append('volume #' + str(idx + 3) + ' color ' + tif_colors[0] + ' style image ' + threshold_level + '\n')

                    # frame_n
                    commands.append('open \"' + sorted(glob.glob(self.tracked_mito.list_of_folders[int(frame_n)]+'/*.tif'))[0] + '\"\n')
                    commands.append('volume #' + str(idx + 4) + ' voxelSize ' + voxel_size + '\n')
                    commands.append('volume flip #' + str(idx + 4) + ' axis y\n')
                    commands.append('close #' + str(idx + 4) + '\n')
                    commands.append('rename #' + str(idx + 5) + ' id #' + str(idx + 4) + '\n')
                    commands.append('volume #' + str(idx + 4) + ' color ' + tif_colors[1] + ' style image ' + threshold_level + '\n')

                # combine the models
                if show_tif:
                    commands.append('rename #' + str(idx) + '-' + str(idx + 4) + ' id ' + str(idx) + '\n')
                else:
                    commands.append('rename #' + str(idx) + '-' + str(idx + 2) + ' id ' + str(idx) + '\n')

                idx += 1

            else:
                print('End frame is invalid. Using the last frame.')
                break

        if end_frame - start_frame > 1:
            commands.append('mseries slider #1-' + str(idx - 1))

        bild.writelines(commands)
        bild.close()
        print(f'Load file {self.visual_path}visualize_tracking.cxc in ChimeraX to visualize tracking.')


    def map_mitochondria_motility(self, levels: list = None, node_size: float = 0.3, selected_frames: list = None):
        """
        Wrapper to map mitochondrial motility (diffusivity) to 3D visualization
        at one or more structural levels: nodes, segments, or fragments.

        Parameters
        ----------
        levels : list of strings, optional
            Structural levels at which to map motility. Choose from any combination of "node", "segment", "fragment"
            (default is ["node", "segment", "fragment"]).
        node_size : float, optional
            Radius of spheres used to represent nodes (default=0.3).
        selected_frames : list of int
            Frame indices to visualize. Must be explicitly provided.

        Returns
        -------
        None
            `.bild` visualization files are written to `self.visual_path`.
        """

        if levels is None:
            levels = ["node", "segment", "fragment"]
        if selected_frames is None:
            raise Exception("Please specify the frames to visualize.")

        levels = [lvl.lower() for lvl in levels]

        for level in levels:
            if level == "node":
                self.map_mitochondria_node_motility(node_size=node_size, selected_frames=selected_frames)

            elif level == "segment":
                self.map_mitochondria_segment_motility(node_size=node_size, selected_frames=selected_frames)

            elif level == "fragment":
                self.map_mitochondria_fragment_motility(node_size=node_size, selected_frames=selected_frames)

            else:
                raise ValueError("levels must be a list containing 'node', 'segment', or 'fragment'")


    def map_mitochondria_node_motility(self, node_size: float = 0.3, selected_frames = None):
        """
            Map node-level mitochondrial motility (diffusivity) to 3D visualization.

            For each selected frame, node-level diffusivity values are extracted from
            `self.tracked_mito.node_diffusivity` and mapped onto the skeleton graph.
            Nodes are rendered as spheres with colors representing normalized diffusivity.
            The output is exported in ChimeraX BILD format for visualization.

            Parameters
            ----------
            node_size : float, optional
                Radius of spheres used to represent nodes (default is 0.3).
            selected_frames : list of int
                Frame indices to visualize. Must be explicitly provided.

            Returns
            -------
            None
                `.bild` files are written to the directory `self.visual_path`.
        """
        if selected_frames is None:
            raise Exception('Please specify the frames to visualize.')

        full_graphs_all_frames = self.tracked_mito.full_graphs
        tracks = self.tracked_mito.node_tracks
        node_diffusivity_df = self.tracked_mito.node_diffusivity

        for center_frame in selected_frames:
            graph = full_graphs_all_frames[center_frame]
            num_nodes = len(graph.vs)

            node_diffusivity = np.empty(num_nodes)
            node_diffusivity[:] = np.nan

            frame_tracks = tracks[tracks.frame_id==center_frame]
            frame_track_ids = [int(n) for n in frame_tracks.unique_node_id.tolist()]
            unique_nodes = np.array(frame_tracks['unique_node_id'].tolist(), dtype=int)
            frame_nodes = np.array(frame_tracks['frame_node_id'].tolist(), dtype=int)
            frame_to_unique = {frame_nodes[i]:unique_nodes[i] for i in range(len(frame_nodes))}
            unique_to_frame = {unique_nodes[i]:frame_nodes[i] for i in range(len(unique_nodes))}

            for track_id in frame_track_ids:

                d = node_diffusivity_df.loc[track_id].diffusivity
                r_squared = node_diffusivity_df.loc[track_id].r_squared

                frame_node_index = unique_to_frame[track_id]

                if r_squared >= 0.8:
                    node_diffusivity[frame_node_index] = d

            print(f"{np.sum(~np.isnan(node_diffusivity))} nodes are mapped out of total {num_nodes} nodes\n")

            # get normalized diffusivity
            d_max = np.nanpercentile(node_diffusivity, 90)
            d_normalized = node_diffusivity / d_max

            # make .bild file
            file_dir = self.visual_path+'map_node_motility_frame_'+str(center_frame).zfill(3)+'.bild'
            if os.path.exists(file_dir):
                os.remove(file_dir)
            bild = open(file_dir, "x")
            commands = []

            try:
                coords = graph.vs['coordinate']
                for idx in range(num_nodes):
                    commands.append('.color ' + _color_motility(d_normalized[idx]) + '\n')
                    commands.append('.sphere ' + _coord_to_str(coords[idx]) + str(node_size) + '\n')
                bild.writelines(commands)
                bild.close()
            except:
                bild.close()


    def map_mitochondria_segment_motility(self, node_size: float = 0.3, selected_frames = None):
        """
            Map segment-level mitochondrial motility (diffusivity) to 3D visualization.

            For each selected frame, segment-level diffusivity values are extracted from
            `self.tracked_mito.segment_diffusivity` and mapped onto the skeleton graph.
            Each segment is rendered by coloring its constituent nodes, represented as
            spheres, with colors reflecting normalized diffusivity values. The output
            is exported in ChimeraX BILD format for visualization.

            Parameters
            ----------
            node_size : float, optional
                Radius of spheres used to represent segment nodes (default is 0.3).
            selected_frames : list of int
                Frame indices to visualize. Must be explicitly provided.

            Returns
            -------
            None
                `.bild` files are written to the directory `self.visual_path`.
        """

        if selected_frames is None:
            raise Exception('Please specify the frames to visualize.')

        full_graphs_all_frames = self.tracked_mito.full_graphs
        all_segment_nodes = self.tracked_mito.segment_nodes
        seg_diffusivity_df = self.tracked_mito.segment_diffusivity

        for center_frame in selected_frames:
            graph = full_graphs_all_frames[center_frame]
            segment_nodes = all_segment_nodes[center_frame]
            num_segs = len(segment_nodes)

            seg_diffusivity = seg_diffusivity_df[seg_diffusivity_df['center_frame_id']==center_frame].diffusivity

            print(f"{np.sum(~np.isnan(seg_diffusivity))} segments are mapped out of total {num_segs} segments\n")

            # get normalized diffusivity
            d_max = np.nanpercentile(seg_diffusivity, 90)
            d_normalized = seg_diffusivity / d_max

            # make .bild file
            file_dir = self.visual_path+'map_segment_motility_frame_'+str(center_frame).zfill(3)+'.bild'
            if os.path.exists(file_dir):
                os.remove(file_dir)
            bild = open(file_dir, "x")
            commands = []

            try:
                coords = graph.vs['coordinate']
                for seg_id, seg in enumerate(segment_nodes):
                    commands.append('.color ' + _color_motility(d_normalized[seg_id]) + '\n')
                    for node in seg:
                        commands.append('.sphere ' + _coord_to_str(coords[node]) + str(node_size) + '\n')
                bild.writelines(commands)
                bild.close()
            except:
                bild.close()


    def map_mitochondria_fragment_motility(self, node_size: float = 0.3, selected_frames = None):
        """
            Map fragment-level mitochondrial motility (diffusivity) to 3D visualization.

            For each selected frame, fragment-level diffusivity values are extracted from
            `self.tracked_mito.fragment_diffusivity` and mapped onto the skeleton graph.
            Each fragment is rendered by coloring its constituent nodes, represented as
            spheres, with colors reflecting normalized diffusivity values. The output
            is exported in ChimeraX BILD format for visualization.

            Parameters
            ----------
            node_size : float, optional
                Radius of spheres used to represent fragment nodes (default is 0.3).
            selected_frames : list of int
                Frame indices to visualize. Must be explicitly provided.

            Returns
            -------
            None
                `.bild` files are written to the directory `self.visual_path`.
        """

        if selected_frames is None:
            raise Exception('Please specify the frames to visualize.')

        full_graphs_all_frames = self.tracked_mito.full_graphs
        frag_diffusivity_df = self.tracked_mito.fragment_diffusivity

        for center_frame in selected_frames:
            graph = full_graphs_all_frames[center_frame]
            fragment_nodes = graph.components()
            num_frags = len(fragment_nodes)

            frag_diffusivity = frag_diffusivity_df[frag_diffusivity_df['center_frame_id']==center_frame].diffusivity

            print(f"{np.sum(~np.isnan(frag_diffusivity))} fragments are mapped out of total {num_frags} fragments\n")

            # get normalized diffusivity
            d_max = np.nanpercentile(frag_diffusivity, 90)
            d_normalized = frag_diffusivity / d_max

            # make .bild file
            file_dir = self.visual_path+'map_fragment_motility_frame_'+str(center_frame).zfill(3)+'.bild'
            if os.path.exists(file_dir):
                os.remove(file_dir)
            bild = open(file_dir, "x")
            commands = []

            try:
                coords = graph.vs['coordinate']
                for frag_id, frag in enumerate(fragment_nodes):
                    commands.append('.color ' + _color_motility(d_normalized[frag_id]) + '\n')
                    for node in frag:
                        commands.append('.sphere ' + _coord_to_str(coords[node]) + str(node_size) + '\n')
                bild.writelines(commands)
                bild.close()
            except:
                bild.close()
                
                
def _time_to_rgb(time, cmap):
    cmap = plt.get_cmap(cmap)
    color = cmap(time)
    return str(color[0]) + ' ' + str(color[1]) + ' ' + str(color[2])


def _coord_to_str(coord):
    string = ''
    for s in coord:
        string = string + str(np.round(s, 3)) + ' '
    return string


def _color_motility(diffusivity):
    if np.isnan(diffusivity):
        return '1.0 1.0 1.0' # white
    else:
        cmap = plt.get_cmap('coolwarm')
        color = cmap(diffusivity)
    return str(color[0])+' '+str(color[1])+' '+str(color[2])