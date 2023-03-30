import time, os
import numpy as np
import pandas as pd
import igraph as ig
import pickle
import matplotlib.pyplot as plt
from scipy.linalg import lstsq
from tqdm.notebook import trange

def compute_node_diffusivity(input_dir, output_dir, analy_motility_dir,
                             frame_interval, max_tau):

    tracks = pd.read_csv(output_dir+'final_node_tracks.csv')
    num_tracks = int(np.max(tracks['unique_node_id'])) + 1

    print('Computing MSD for tracks ...')
    all_msd = []
    for track_id in trange(num_tracks):
        track = tracks[tracks['unique_node_id']==track_id]
        frames = track['frame_id'].to_numpy()
        coords = track.loc[:,'x':'z'].to_numpy()
        coords = coords

        # Calculate TA-MSD
        track_msd = []
        for tau in range(1, max_tau):

            disp = []
            frame = 0
            next_frame = frame + tau

            while next_frame < frames[-1]:
                if next_frame in frames and frame in frames:
                    node = frames.tolist().index(frame)
                    next_node = frames.tolist().index(next_frame)

                    disp.append(np.sum((coords[next_node]-coords[node])**2))

                frame += 1 # next start frame
                next_frame = frame + tau

            if len(disp) < 2:
                break
            else:
                track_msd.append(np.mean(disp))

        all_msd.append(track_msd)

    msd_matrix = np.zeros((num_tracks, max_tau))
    for track_id, track_msd in enumerate(all_msd):
        for i in range(max_tau):
            if i < len(track_msd):
                msd_matrix[track_id,i] = track_msd[i]
            else:
                msd_matrix[track_id,i] = np.nan

    print('Fitting linear curve and computing diffusion coefficient ...')
    node_diffusivity = []
    for n in trange(num_tracks):
        eata_msd = [msd for msd in msd_matrix[n] if not np.isnan(msd)]
        eata_msd.insert(0,0)

        # choose the number of data points to fit
        if len(eata_msd) < 2:
            node_diffusivity.append({'unique_node_id':n, 'diffusivity':np.nan, 'msd':np.nan, 'r_squared':np.nan, 'num_points':1})
            continue
        elif len(eata_msd) > max_tau:
            n_points = max_tau
        else:
            n_points = len(eata_msd)

        all_tau = np.arange(0, n_points*frame_interval, frame_interval)[:,np.newaxis]
        slope, res, rnk, s = lstsq(all_tau[:n_points], eata_msd[:n_points])
        d = slope[0] / 6
        msd_per_frame = 6*d*frame_interval

        # get r^2
        msd_mean = np.mean(eata_msd[:n_points])
        total_sum = np.sum([(i-msd_mean)**2 for i in eata_msd[:n_points]])

        if total_sum == 0:
            r_squared = np.nan
        else:
            r_squared = 1 - res/total_sum

        # store data
        node_diffusivity.append({'unique_node_id':n, 'diffusivity':d, 'msd':msd_per_frame, 'r_squared':r_squared, 'num_points':n_points})

    node_diffusivity = pd.DataFrame.from_dict(node_diffusivity)
    node_diffusivity.to_csv(analy_motility_dir+'node_diffusivity.csv', index=False)

    print('Data saved at '+analy_motility_dir+'node_diffusivity.csv')


def compute_segment_diffusivity(input_dir, output_dir, analy_motility_dir,
                                frame_interval, max_tau, tracked_ratio, half_win_size, selected_frames):
    # load inputs
    tracks = pd.read_csv(output_dir+'final_node_tracks.csv')
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    segment_node_all_frames = inputs['segment_nodes']

    seg_diffusivity = []

    # iterate the center frames
    for center_frame in selected_frames:
        segment_nodes = segment_node_all_frames[center_frame]
        num_segs = len(segment_nodes)

        frame_tracks = tracks[tracks.frame_id==center_frame] # find only tracks that intersects with the center frame
        unique_nodes = np.array(frame_tracks['unique_node_id'].tolist(), dtype=int)
        frame_nodes = np.array(frame_tracks['frame_node_id'].tolist(), dtype=int)
        frame_to_unique = {frame_nodes[i]:unique_nodes[i] for i in range(len(frame_nodes))}
        unique_to_frame = {unique_nodes[i]:frame_nodes[i] for i in range(len(unique_nodes))}

        all_msd, all_tau = [], []
        print('Computing MSD for tracks around frame', center_frame, '...')
        for seg_id in trange(len(segment_nodes)):

            segment = segment_nodes[seg_id]
            seg_coords = np.zeros((len(segment), 2*half_win_size), dtype=object)
            seg_coords.fill(np.array([np.nan, np.nan, np.nan], dtype=object))

            for node_id, frame_node in enumerate(segment):
                if frame_node in frame_to_unique.keys():
                    full_track = tracks[tracks.unique_node_id==frame_to_unique[frame_node]]

                    for i in range(len(full_track)):
                        frame = int(full_track.iloc[i]['frame_id'])
                        coord = full_track.iloc[i]['x':'z'].to_numpy()
                        arr_index = frame - (center_frame - half_win_size)

                        if arr_index < 2*half_win_size and arr_index >= 0:
                            seg_coords[node_id,arr_index] = coord
                        else:
                            break

            # Calculate TA-MSD
            seg_msd, seg_tau = [], []
            for tau in range(1, max_tau):

                disp = []
                frame = 0
                next_frame = frame + tau

                while next_frame < 2 * half_win_size:
                    coords_m = seg_coords[:,frame]
                    coords_n = seg_coords[:,next_frame]
                    vector_diff = coords_n - coords_m # get displacement vector for each node

                    disp_vectors = np.array([i for i in vector_diff if not pd.isnull(i).any()]) # collect the valid vectors that belong to tracked nodes
                    num_disp = disp_vectors.shape[0]

                    if num_disp >= tracked_ratio * len(segment): # if enough of the segment is tracked
                        # average all displacement vector to get segment vector
                        average_disp = (np.linalg.norm(np.mean(disp_vectors, axis=0))) ** 2
                        disp.append(average_disp)

                    frame += 1 # next start frame
                    next_frame = frame + tau

                if len(disp) < 2:
                    break
                else:
                    seg_msd.append(np.mean(disp))
                    seg_tau.append(tau)

            all_msd.append(seg_msd)
            all_tau.append(seg_tau)

        print('Percent of segments tracked at frame', center_frame, ':', (1 - np.sum([len(msd)==0 for msd in all_msd]) / len(all_msd)) * 100, '%\n')

        # compute diffusivity from MSD
        num_msd = len(all_msd)
        msd_matrix = np.zeros((num_msd, max_tau))
        for track_id, track_msd in enumerate(all_msd):
            for i in range(max_tau):
                if i < len(track_msd):
                    msd_matrix[track_id,i] = track_msd[i]
                else:
                    msd_matrix[track_id,i] = np.nan

        print('Fitting linear curve and computing diffusion coefficient at frame', center_frame, '...')
        for seg_id in trange(num_segs):
            eata_msd = [msd for msd in msd_matrix[seg_id] if not np.isnan(msd)]
            eata_msd.insert(0,0)

            # choose the number of data points to fit
            if len(eata_msd) <= 1:
                seg_diffusivity.append({'center_frame_id': center_frame, 'seg_id':seg_id, 'diffusivity':np.nan, 'msd':np.nan, 'r_squared':np.nan, 'num_points':1})
                continue
            elif len(eata_msd) > max_tau:
                n_points = max_tau
            else:
                n_points = len(eata_msd)

            all_tau = np.arange(0, n_points*frame_interval, frame_interval)[:,np.newaxis]
            slope, res, rnk, s = lstsq(all_tau[:n_points], eata_msd[:n_points])
            d = slope[0] / 6
            msd_per_frame = 6*d*frame_interval

            # get r^2
            msd_mean = np.mean(eata_msd[:n_points])
            total_sum = np.sum([(i-msd_mean)**2 for i in eata_msd[:n_points]])

            if total_sum == 0:
                r_squared = np.nan
            else:
                r_squared = 1 - res/total_sum

            seg_diffusivity.append({'center_frame_id': center_frame, 'seg_id':seg_id, 'diffusivity':d, 'msd':msd_per_frame, 'r_squared':r_squared, 'num_points':n_points})

    seg_diffusivity = pd.DataFrame.from_dict(seg_diffusivity)
    seg_diffusivity.to_csv(analy_motility_dir+'segment_diffusivity.csv', index=False)

    print('Data saved at '+analy_motility_dir+'segment_diffusivity.csv')



def compute_fragment_diffusivity(input_dir, output_dir, analy_motility_dir,
                                 frame_interval, max_tau, tracked_ratio, half_win_size, selected_frames):
    # load inputs
    tracks = pd.read_csv(output_dir+'final_node_tracks.csv')
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']

    frag_diffusivity = []

    # iterate the center frames
    for center_frame in selected_frames:
        graph = full_graph_all_frames[center_frame]
        fragment_nodes = graph.components()
        num_frags = len(fragment_nodes)

        frame_tracks = tracks[tracks.frame_id==center_frame] # find only tracks that intersects with the center frame
        unique_nodes = np.array(frame_tracks['unique_node_id'].tolist(), dtype=int)
        frame_nodes = np.array(frame_tracks['frame_node_id'].tolist(), dtype=int)
        frame_to_unique = {frame_nodes[i]:unique_nodes[i] for i in range(len(frame_nodes))}
        unique_to_frame = {unique_nodes[i]:frame_nodes[i] for i in range(len(unique_nodes))}

        all_msd, all_tau = [], []
        print('Computing MSD for tracks around frame', center_frame, '...')
        for frag_id in trange(len(fragment_nodes)):

            fragment = fragment_nodes[frag_id]
            frag_coords = np.zeros((len(fragment), 2*half_win_size), dtype=object)
            frag_coords.fill(np.array([np.nan, np.nan, np.nan], dtype=object))

            for node_id, frame_node in enumerate(fragment):
                if frame_node in frame_to_unique.keys():
                    full_track = tracks[tracks.unique_node_id==frame_to_unique[frame_node]]

                    for i in range(len(full_track)):
                        frame = int(full_track.iloc[i]['frame_id'])
                        coord = full_track.iloc[i]['x':'z'].to_numpy()
                        arr_index = frame - (center_frame - half_win_size)

                        if arr_index < 2*half_win_size and arr_index >= 0:
                            frag_coords[node_id,arr_index] = coord
                        else:
                            break

            # Calculate TA-MSD
            frag_msd, frag_tau = [], []
            for tau in range(1, max_tau):

                disp = []
                frame = 0
                next_frame = frame + tau

                while next_frame < 2 * half_win_size:
                    coords_m = frag_coords[:,frame]
                    coords_n = frag_coords[:,next_frame]
                    vector_diff = coords_n - coords_m # get displacement vector for each node

                    disp_vectors = np.array([i for i in vector_diff if not pd.isnull(i).any()]) # collect the valid vectors that belong to tracked nodes
                    num_disp = disp_vectors.shape[0]

                    if num_disp >= tracked_ratio * len(fragment): # if enough of the fragment is tracked
                        # average all displacement vector to get fragment vector
                        average_disp = (np.linalg.norm(np.mean(disp_vectors, axis=0))) ** 2
                        disp.append(average_disp)

                    frame += 1 # next start frame
                    next_frame = frame + tau

                if len(disp) < 2:
                    break
                else:
                    frag_msd.append(np.mean(disp))
                    frag_tau.append(tau)

            all_msd.append(frag_msd)
            all_tau.append(frag_tau)

        print('Percent of fragments tracked at frame', center_frame, ':', (1 - np.sum([len(msd)==0 for msd in all_msd]) / len(all_msd)) * 100, '%\n')

        # compute diffusivity from MSD
        num_msd = len(all_msd)
        msd_matrix = np.zeros((num_msd, max_tau))
        for track_id, track_msd in enumerate(all_msd):
            for i in range(max_tau):
                if i < len(track_msd):
                    msd_matrix[track_id,i] = track_msd[i]
                else:
                    msd_matrix[track_id,i] = np.nan

        print('Fitting linear curve and computing diffusion coefficient at frame', center_frame, '...')
        for frag_id in trange(num_frags):
            eata_msd = [msd for msd in msd_matrix[frag_id] if not np.isnan(msd)]
            eata_msd.insert(0,0)

            # choose the number of data points to fit
            if len(eata_msd) <= 1:
                frag_diffusivity.append({'center_frame_id': center_frame, 'frag_id':frag_id, 'diffusivity':np.nan, 'msd':np.nan, 'r_squared':np.nan, 'num_points':1})
                continue
            elif len(eata_msd) > max_tau:
                n_points = max_tau
            else:
                n_points = len(eata_msd)

            all_tau = np.arange(0, n_points*frame_interval, frame_interval)[:,np.newaxis]
            slope, res, rnk, s = lstsq(all_tau[:n_points], eata_msd[:n_points])
            d = slope[0] / 6
            msd_per_frame = 6*d*frame_interval

            # get r^2
            msd_mean = np.mean(eata_msd[:n_points])
            total_sum = np.sum([(i-msd_mean)**2 for i in eata_msd[:n_points]])

            if total_sum == 0:
                r_squared = np.nan
            else:
                r_squared = 1 - res/total_sum

            frag_diffusivity.append({'center_frame_id': center_frame, 'frag_id':frag_id, 'diffusivity':d, 'msd':msd_per_frame, 'r_squared':r_squared, 'num_points':n_points})

    frag_diffusivity = pd.DataFrame.from_dict(frag_diffusivity)
    frag_diffusivity.to_csv(analy_motility_dir+'fragment_diffusivity.csv', index=False)

    print('Data saved at '+analy_motility_dir+'fragment_diffusivity.csv')

def coord_to_str(coord):
    string = ''
    for s in coord:
        string = string + str(np.round(s,3)) + ' '
    return string

def color_motility(diffusivity):
    if np.isnan(diffusivity):
        return '1.0 1.0 1.0' # white
    else:
        cmap = plt.get_cmap('coolwarm')
        color = cmap(diffusivity)
    return str(color[0])+' '+str(color[1])+' '+str(color[2])

def map_node_motility_onto_surface(input_dir, output_dir, analy_motility_dir,
                                   node_size, selected_frames):

    # load data
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']

    tracks = pd.read_csv(output_dir+'final_node_tracks.csv')

    node_diffusivity_df = pd.read_csv(analy_motility_dir+'node_diffusivity.csv')

    for center_frame in selected_frames:
        graph = full_graph_all_frames[center_frame]
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

        print('{} nodes are mapped out of total {} nodes\n'.format(np.sum(~np.isnan(node_diffusivity)), num_nodes))

        # get normalized diffusivity
        d_max = np.nanpercentile(node_diffusivity, 90)
        d_normalized = node_diffusivity / d_max

        # make .bild file
        file_dir = analy_motility_dir+'map_node_motility_frame_'+str(center_frame)+'.bild'
        if os.path.exists(file_dir):
            os.remove(file_dir)
        bild = open(file_dir, "x")
        commands = []

        try:
            coords = graph.vs['coordinate']

            for idx in range(num_nodes):
                commands.append('.color '+color_motility(d_normalized[idx])+'\n')
                commands.append('.sphere '+coord_to_str(coords[idx])+str(node_size)+'\n')

            bild.writelines(commands)
            bild.close()
        except:
            bild.close()

        print('Load file', file_dir, '\nin ChimeraX to visualize motility map')


def map_segment_motility_onto_surface(input_dir, output_dir, analy_motility_dir,
                                      node_size, selected_frames):

    # load data
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']
    all_segment_nodes = inputs['segment_nodes']

    seg_diffusivity_df = pd.read_csv(analy_motility_dir+'segment_diffusivity.csv')

    for center_frame in selected_frames:
        graph = full_graph_all_frames[center_frame]

        segment_nodes = all_segment_nodes[center_frame]
        num_segs = len(segment_nodes)

        seg_diffusivity = seg_diffusivity_df[seg_diffusivity_df['center_frame_id']==center_frame].diffusivity


        print('{} segments are mapped out of total {} segments\n'.format(np.sum(~np.isnan(seg_diffusivity)), num_segs))

        # get normalized diffusivity
        d_max = np.nanpercentile(seg_diffusivity, 90)
        d_normalized = seg_diffusivity / d_max

        # make .bild file
        file_dir = analy_motility_dir+'map_segment_motility_frame_'+str(center_frame)+'.bild'
        if os.path.exists(file_dir):
            os.remove(file_dir)
        bild = open(file_dir, "x")
        commands = []

        try:
            coords = graph.vs['coordinate']

            for seg_id, seg in enumerate(segment_nodes):
                commands.append('.color '+color_motility(d_normalized[seg_id])+'\n')
                for node in seg:
                    commands.append('.sphere '+coord_to_str(coords[node])+str(node_size)+'\n')

            bild.writelines(commands)
            bild.close()
        except:
            bild.close()

        print('Load file', file_dir, '\nin ChimeraX to visualize motility map')


def map_fragment_motility_onto_surface(input_dir, output_dir, analy_motility_dir,
                                       node_size, selected_frames):

    # load data
    inputs = np.load(input_dir+'tracking_inputs.npz', allow_pickle=True)
    full_graph_all_frames = inputs['full_graphs']

    frag_diffusivity_df = pd.read_csv(analy_motility_dir+'fragment_diffusivity.csv')

    for center_frame in selected_frames:
        graph = full_graph_all_frames[center_frame]

        fragment_nodes = graph.components()
        num_frags = len(fragment_nodes)

        frag_diffusivity = frag_diffusivity_df[frag_diffusivity_df['center_frame_id']==center_frame].diffusivity

        print('{} fragments are mapped out of total {} fragments\n'.format(np.sum(~np.isnan(frag_diffusivity)), num_frags))

        # get normalized diffusivity
        d_max = np.nanpercentile(frag_diffusivity, 90)
        d_normalized = frag_diffusivity / d_max

        # make .bild file
        file_dir = analy_motility_dir+'map_fragment_motility_frame_'+str(center_frame)+'.bild'
        if os.path.exists(file_dir):
            os.remove(file_dir)
        bild = open(file_dir, "x")
        commands = []

        try:
            coords = graph.vs['coordinate']

            for frag_id, frag in enumerate(fragment_nodes):
                commands.append('.color '+color_motility(d_normalized[frag_id])+'\n')
                for node in frag:
                    commands.append('.sphere '+coord_to_str(coords[node])+str(node_size)+'\n')

            bild.writelines(commands)
            bild.close()
        except:
            bild.close()

        print('Load file', file_dir, '\nin ChimeraX to visualize motility map')