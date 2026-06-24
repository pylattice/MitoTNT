import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

class FeatureExtractor:
    """
    Aggregate per-frame mitochondrial features from saved MitoTNT results into
    tidy (Feature, Value) tables and plot their distributions.

    Reads the artifacts written under the tracking ``save_path``
    (``extracted_graphs.npz``, the ``*_diffusivity.csv`` files, and
    ``remodeling_events.csv``) and produces three feature tables:
    ``graph_features``, ``motility_features``, and ``remodeling_features``.

    Parameters
    ----------
    tracked_mito : TrackedMito, optional
        A `TrackedMito` object; its ``save_path`` is used to locate results.
        Preferred, consistent with the rest of the API.
    save_path : str, optional
        Path to the directory containing MitoTNT results (must end with a
        path separator). Used only if `tracked_mito` is not given.

    Attributes
    ----------
    save_path : str
        Directory the result files are read from / written to.
    graph_features : pandas.DataFrame or None
        Per-frame graph/topology features. None until computed.
    motility_features : pandas.DataFrame or None
        Diffusivity features by level. None until computed.
    remodeling_features : pandas.DataFrame or None
        Per-frame fusion/fission rates. None until computed.
    """

    def __init__(self, tracked_mito=None, save_path: str = None):

        if tracked_mito is not None:
            self.save_path = tracked_mito.save_path
        elif save_path is not None:
            self.save_path = save_path
        else:
            raise ValueError("Provide either a TrackedMito object or a save_path.")

        self.graph_features = None
        self.motility_features = None
        self.remodeling_features = None

    def compute_graph_features(self, selected_frames: list = None, save_csv: bool = False):
        mitograph_data = np.load(self.save_path + 'extracted_graphs.npz', allow_pickle=True)
        full_graphs_all_frames = mitograph_data['full_graphs']
        classic_graphs_all_frames = mitograph_data['simple_graphs']
        segment_nodes_all_frames = mitograph_data['segment_nodes']

        num_frames = len(full_graphs_all_frames)
        if selected_frames is None:
            selected_frames = range(num_frames)

        # Use a list to collect DataFrames for efficiency
        all_frames_list = []

        for frame_id in selected_frames:
            feature_list, value_list = [], []

            full_graphs = full_graphs_all_frames[frame_id]
            classic_graphs = classic_graphs_all_frames[frame_id]
            fragments = classic_graphs.components()
            segment_nodes = segment_nodes_all_frames[frame_id]

            # Check if graph has edges/distances to avoid crashes
            try:
                if 'distance' not in classic_graphs.es.attributes():
                    continue
            except:
                continue

            # --- Node and Segment Counts ---
            feature_list.append('total_node_count')
            value_list.append(full_graphs.vcount())

            feature_list.append('total_segment_count')
            value_list.append(classic_graphs.ecount())

            # --- Segment Lengths (Vectorized addition) ---
            lengths = classic_graphs.es['distance']
            feature_list.extend(['segment_length'] * len(lengths))
            value_list.extend(lengths)

            # --- Segment Widths ---
            for seg in segment_nodes:
                seg_widths = np.array(full_graphs.vs[seg]['width'])
                feature_list.append('segment_sum_width')
                value_list.append(seg_widths.sum())
                feature_list.append('segment_average_width')
                value_list.append(seg_widths.mean())

            # --- Fragment Features ---
            feature_list.append('total_fragment_count')
            value_list.append(len(fragments))

            for frag_nodes in fragments:
                g = classic_graphs.induced_subgraph(frag_nodes)
                if g.vcount() >= 2:
                    frag_length = np.sum(g.es['distance'])
                    frag_diameter = g.diameter(weights='distance')
                    
                    feature_list.append('fragment_length')
                    value_list.append(frag_length)
                    feature_list.append('fragment_diameter')
                    value_list.append(frag_diameter)
                    
                    if frag_diameter > 0:
                        feature_list.append('fragment_tortuosity')
                        value_list.append(frag_length / frag_diameter)

                    # Topological Features
                    if g.vcount() >= 5 and g.ecount() >= 4:
                        feature_list.append('graph_density')
                        value_list.append(g.density())
                        
                        feature_list.append('graph_efficiency')
                        value_list.append(_global_efficiency(g))
                        
                        feature_list.append('graph_clustering_coefficient')
                        value_list.append(g.transitivity_undirected())

                        node_betweenness = np.array(g.betweenness())
                        feature_list.append('graph_max_betweenness')
                        value_list.append(node_betweenness.max())
                        feature_list.append('graph_mean_betweenness')
                        value_list.append(node_betweenness.mean())

                        feature_list.append('fragment_branching_index')
                        value_list.append(g.ecount() / g.vcount())

                        degrees = np.array(g.degree())
                        num_branchpoint = np.sum(degrees >= 3)
                        num_endpoint = np.sum(degrees == 1)
                        if num_endpoint > 0:
                            feature_list.append('fragment_branchpoint_to_endpoint_ratio')
                            value_list.append(num_branchpoint / num_endpoint)

            # Create a DataFrame for this frame and store it
            if feature_list:
                frame_df = pd.DataFrame({'Feature': feature_list, 'Value': value_list})
                frame_df['Frame'] = frame_id
                all_frames_list.append(frame_df)

        if all_frames_list:
            self.graph_features = pd.concat(all_frames_list, ignore_index=True)
        else:
            self.graph_features = pd.DataFrame(columns=['Feature', 'Value', 'Frame'])

        if save_csv and not self.graph_features.empty:
            self.graph_features.to_csv(self.save_path + 'graph_features.csv', index=False)

        return


    def compute_diffusivity(self, levels: list = ['node', 'segment', 'fragment'], save_csv: bool = False):

        created = False
        data = None
        feature_list, value_list = [], []

        for level in levels:
            track_diffusivity_df = pd.read_csv(self.save_path+level.lower()+'_diffusivity.csv')
            track_diffusivity_df['level'] = level

            if not created:
                data = track_diffusivity_df
                created = True
            else:
                data = pd.concat([data, track_diffusivity_df])
    
        # Derive the feature name from the level so this stays correct for any
        # subset or ordering of `levels` (previously paired a hardcoded feature
        # list against levels[idx], which mislabeled / IndexError'd otherwise).
        for level in levels:
            select_values = data[(data.r_squared > 0.8) & (data.num_points >= 3) & (data.level == level)]['diffusivity'].tolist()
            feature_list += len(select_values) * [f'{level}_diffusivity']
            value_list += select_values

        self.motility_features = pd.DataFrame({'Feature': feature_list, 'Value': value_list})
    
        if save_csv:
            self.motility_features.to_csv(self.save_path+'motility_features.csv', encoding='utf-8', index=False)

        return


    def compute_remodeling_rates(self, save_csv: bool = False):

        created = False
        data = None
        feature_list, value_list = [], []

        mitograph_data = np.load(self.save_path+'extracted_graphs.npz', allow_pickle=True)
        full_graphs_all_frames = mitograph_data['full_graphs']
        remodeling_events = pd.read_csv(self.save_path+'remodeling_events.csv')

        # No events detected (e.g. too few frames for the detection window) ->
        # nothing to compute; np.min/max of an empty column would be NaN and
        # range(NaN, NaN) raises. Emit an empty feature table instead.
        if remodeling_events.empty:
            self.remodeling_features = pd.DataFrame({'Feature': [], 'Value': []})
            if save_csv:
                self.remodeling_features.to_csv(self.save_path+'remodeling_features.csv', index=False)
            return

        for frame in range(int(np.min(remodeling_events['frame_id'])), int(np.max(remodeling_events['frame_id']))):

            full_graph = full_graphs_all_frames[frame]
            num_nodes = len(full_graph.vs)
            
            frame_events = remodeling_events[remodeling_events['frame_id']==frame]
            
            frame_fusions = frame_events[frame_events['type']=='fusion']
            num_fusions = 0
            for frags in frame_fusions['frag_id_before'].values:
                # frags is a space-joined string; an empty event field would make
                # ''.split(' ') -> [''] and int('') crash. Count only non-empty.
                if isinstance(frags, str) and frags.strip():
                    num_fusions += 1

            frame_fissions = frame_events[frame_events['type']=='fission']
            num_fissions = 0
            for frags in frame_fissions['frag_id_after'].values:
                if isinstance(frags, str) and frags.strip():
                    num_fissions += 1

            feature_list.append('fusion_rate')
            value_list.append(num_fusions/num_nodes)

            feature_list.append('fission_rate')
            value_list.append(num_fissions/num_nodes)

        self.remodeling_features = pd.DataFrame({'Feature': feature_list, 'Value': value_list})
    
        if save_csv:
            self.remodeling_features.to_csv(self.save_path+'remodeling_features.csv', encoding='utf-8', index=False)

        return


    def plot_features_as_violinplot(self, feature_category: list = None,
                                font_scale: float = 1.0,
                                specify_ylim: dict = None):

        if feature_category is None:

            feature_category = ['graph_features', 'motility_features', 'remodeling_features']

        for category in feature_category:

            if category == 'graph_features':
                df = self.graph_features
            elif category == 'motility_features':
                df = self.motility_features
            elif category == 'remodeling_features':
                df = self.remodeling_features
            else:
                raise ValueError(
                    f"feature_category must be one of: graph_features, motility_features, remodeling_features. Got: {category}"
                )

            for feature in df['Feature'].unique():

                plt.figure(figsize=(8,8))
                sns.set_theme(style='whitegrid', font_scale=font_scale)
                ax = plt.subplot(1, 1, 1)

                current_feature = df[df['Feature'] == feature]
                sns.violinplot(ax=ax, data=current_feature, y='Value', inner='box', saturation=0.8, width=0.8, cut=0)

                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title(feature)
                
                try:
                    ax.set_ylim(specify_ylim[feature])
                except:
                    pass

                plt.tight_layout()
                plt.show()

        return


    def plot_features_as_boxplot(self, feature_category: list = None,
                                 font_scale: float = 1.0,
                                 specify_ylim: dict = None):

        if feature_category is None:
            feature_category = ['graph_features', 'motility_features', 'remodeling_features']

        for category in feature_category:

            if category == 'graph_features':
                df = self.graph_features
            elif category == 'motility_features':
                df = self.motility_features
            elif category == 'remodeling_features':
                df = self.remodeling_features
            else:
                raise ValueError(
                    f"feature_category must be one of: graph_features, motility_features, remodeling_features. Got: {category}"
                )

            for feature in df['Feature'].unique():

                plt.figure(figsize=(8, 8))
                sns.set_theme(style='whitegrid', font_scale=font_scale)
                ax = plt.subplot(1, 1, 1)

                current_feature = df[df['Feature'] == feature]

                sns.boxplot(ax=ax, data=current_feature, y='Value', width=0.5, whis=1.5)

                ax.set_xlabel('')
                ax.set_ylabel('')
                ax.set_title(feature)

                try:
                    ax.set_ylim(specify_ylim[feature])
                except:
                    pass

                plt.tight_layout()
                plt.show()

        return


    def plot_features_as_histogram(self, feature_category: list = None,
                                   font_scale: float = 1.0,
                                   specify_ylim: dict = None,
                                   bins: int = 30,
                                   kde: bool = True):

        if feature_category is None:
            feature_category = ['graph_features', 'motility_features', 'remodeling_features']

        for category in feature_category:

            if category == 'graph_features':
                df = self.graph_features
            elif category == 'motility_features':
                df = self.motility_features
            elif category == 'remodeling_features':
                df = self.remodeling_features
            else:
                raise ValueError(
                    f"feature_category must be one of: graph_features, motility_features, remodeling_features. Got: {category}"
                )

            for feature in df['Feature'].unique():

                plt.figure(figsize=(8, 6))
                sns.set_theme(style='whitegrid', font_scale=font_scale)
                ax = plt.subplot(1, 1, 1)


                current_feature = df[df['Feature'] == feature]

                sns.histplot(data=current_feature, x='Value', bins=bins, kde=kde, stat='probability', alpha=0.8)

                ax.set_xlabel('')
                ax.set_ylabel('Probability')
                ax.set_title(feature)

                try:
                    ax.set_ylim(specify_ylim[feature])
                except:
                    pass

                plt.tight_layout()
                plt.show()

        return


def _global_efficiency(graph):
    n = graph.vcount()
    if n < 2:
        return 0
    
    # distances() returns a matrix; convert to numpy array immediately
    dist_matrix = np.array(graph.distances(weights='distance'))
    
    # Mask the diagonal (self-distances) and infinite distances (unreachable nodes)
    # We only care about 1/d for i != j where d > 0 and d != inf
    with np.errstate(divide='ignore'):
        inv_distances = 1.0 / dist_matrix
    
    # Replace infinity (from 1/0) and NaN with 0
    inv_distances[np.isinf(inv_distances)] = 0
    inv_distances[np.isnan(inv_distances)] = 0
    
    # Efficiency is the sum of inverse distances normalized by n(n-1)
    return np.sum(inv_distances) / (n * (n - 1))
