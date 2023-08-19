import numpy as np
import os
import pandas as pd
import matplotlib.pyplot as plt
from tqdm.notebook import trange

def time_to_rgb(time, cmap):
    cmap = plt.get_cmap(cmap)
    color = cmap(time)
    return str(color[0])+' '+str(color[1])+' '+str(color[2])

def coord_to_str(coord):
    string = ''
    for s in coord:
        string = string + str(np.round(s,3)) + ' '
    return string


skeleton_colors = ['b','r']
def generate_chimerax_skeleton(input_dir, vis_dir, vis_data_dir,
                               start_frame, end_frame, tracking_interval,
                               skeleton_colors, skeleton_size=0.02, node_size=0.02):

    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']

    print('Generate network skeleton in BILD format ...')
    for color in skeleton_colors:
        print('Color', color)
        for frame in trange(start_frame, end_frame+tracking_interval, tracking_interval):

            fid = str(frame)
            file_dir = vis_data_dir+'frame_'+fid+'_chimerax_skeleton_'+color+'.bild'
            if os.path.exists(file_dir):
                os.remove(file_dir)
            bild = open(file_dir, 'x')
            commands = []

            commands.append('.color '+color+'\n')

            full_graph = full_graph_all_frames[frame]
            for edge in full_graph.es:
                coord_m = coord_to_str(full_graph.vs[edge.source]['coordinate'])
                coord_n = coord_to_str(full_graph.vs[edge.target]['coordinate'])
                commands.append('.cylinder '+coord_m+coord_n+str(skeleton_size)+'\n')

            for node in full_graph.vs:
                coord = coord_to_str(node['coordinate'])
                commands.append('.sphere '+coord+str(node_size)+'\n')

            bild.writelines(commands)
            bild.close()
    print('Done')

        
def generate_tracking_arrows(input_dir, output_dir, vis_data_dir,
                              start_frame, end_frame, tracking_interval,
                              arrow_color='black', arrow_size=0.02):

    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']

    outputs = np.load(output_dir+'frametoframe_tracking_outputs.npz', allow_pickle=True)
    linked_nodes = outputs['linked_nodes']

    print('Generate tracking arrows in BILD format ...')
    for frame in trange(start_frame, end_frame, tracking_interval):
        arrows = []
        arrows.append('.color '+arrow_color+'\n')

        file_dir = vis_data_dir+'frame_'+str(frame)+'_arrows.bild'
        if os.path.exists(file_dir):
            os.remove(file_dir)
        bild = open(file_dir, "x")

        coords_m = full_graph_all_frames[frame].vs['coordinate']
        coords_n = full_graph_all_frames[frame+tracking_interval].vs['coordinate']

        # create linking vectors for frame m,n
        linked = linked_nodes[frame]
        for i in range(len(linked)):
            start, end = linked[i,0], linked[i,1]
            start_coord = coords_m[start]
            end_coord = coords_n[end]

            comparison = (start_coord == end_coord)
            if comparison.all() == True:
                end_coord = coords_n[end] + np.random.normal(0, arrow_size, 3)

            arrows.append('.arrow '+coord_to_str(start_coord)+coord_to_str(end_coord)+str(arrow_size)+' '+str(arrow_size*2)+' 0.6\n')

        bild.writelines(arrows)
        bild.close()

    print('Done')


def visualize_tracking(data_dir, input_dir, vis_dir, vis_data_dir,
                       start_frame, end_frame, tracking_interval, 
                       show_tif, voxel_size, tif_colors, threshold_level, 
                       use_chimerax_skeleton, skeleton_colors):

    file_dir = vis_dir+'Visualize tracking.cxc'
    if os.path.exists(file_dir):
        os.remove(file_dir)
    script = open(file_dir, 'x')
    commands = []

    idx = 1
    commands.append('close\n')
    for frame in range(start_frame, end_frame, tracking_interval):

        frame_m = str(frame)
        frame_n = str(frame+tracking_interval)

        # load arrow
        commands.append('open \"'+vis_data_dir+'frame_'+frame_m+'_arrows.bild\"\n')

        if use_chimerax_skeleton:
            # load chimerax skeleton
            commands.append('open \"'+vis_data_dir+'frame_'+frame_m+'_chimerax_skeleton_'+skeleton_colors[0]+'.bild\"'+'\n')
            commands.append('open \"'+vis_data_dir+'frame_'+frame_n+'_chimerax_skeleton_'+skeleton_colors[1]+'.bild\"'+'\n')
            
 
        else:
            # load mitograph skeleton
            commands.append('open \"'+data_dir+'frame_'+frame_m+'/frame_'+frame_m+'_skeleton.vtk\"\n')
            commands.append('open \"'+data_dir+'frame_'+frame_n+'/frame_'+frame_n+'_skeleton.vtk\"\n')

        # load tif
        if show_tif:
            # frame_m
            commands.append('open \"'+data_dir+'frame_'+frame_m+'/frame_'+frame_m+'.tif\"\n')
            commands.append('volume #'+str(idx+3)+' voxelSize '+voxel_size+'\n')
            commands.append('volume flip #'+str(idx+3)+' axis y\n')
            commands.append('close #'+str(idx+3)+'\n')
            commands.append('rename #'+str(idx+4)+' id #'+str(idx+3)+'\n')
            commands.append('volume #'+str(idx+3)+' color '+tif_colors[0]+' style image '+threshold_level+'\n')
            
            # frame_n
            commands.append('open \"'+data_dir+'frame_'+frame_n+'/frame_'+frame_n+'.tif\"\n')
            commands.append('volume #'+str(idx+4)+' voxelSize '+voxel_size+'\n')
            commands.append('volume flip #'+str(idx+4)+' axis y\n')
            commands.append('close #'+str(idx+4)+'\n')
            commands.append('rename #'+str(idx+5)+' id #'+str(idx+4)+'\n')
            commands.append('volume #'+str(idx+4)+' color '+tif_colors[1]+' style image '+threshold_level+'\n')
            
            
        # combine the models
        if show_tif:
            commands.append('rename #'+str(idx)+'-'+str(idx+4)+' id '+str(idx)+'\n')
        else:
            commands.append('rename #'+str(idx)+'-'+str(idx+2)+' id '+str(idx)+'\n')
            
        idx += 1
    
    if end_frame - start_frame > 1:
        commands.append('mseries slider #1-'+str(idx-1))

    script.writelines(commands)
    script.close()

    print('Load file', file_dir, '\nin ChimeraX to visualize tracking')
