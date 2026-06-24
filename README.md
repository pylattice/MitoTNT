## Installation
### Install using PyPI
```bash
pip install mitotnt
```

Download the [repository](https://github.com/pylattice/MitoTNT) or `git clone https://github.com/pylattice/mitoTNT.git` to use the test data and notebook.


### Supporting Software

- **[MitoGraph](https://github.com/vianamp/MitoGraph/)** for mitochondria segmentation and skeletonization.  
Note that the original MitoGraph is currently only available on MacOS and Windows, but you can build it from source on Linux. Please contact ziw056@ucsd.edu if you need the Linux build.
We have also developed a Python-based GPU-accelerated MitoGraph that is cross-platform. This version has been validated against the original CPU implementation.

- **[ChimeraX](https://www.cgl.ucsf.edu/chimerax/)** for the visualization of MitoGraph segmentation and MitoTNT tracking


## Quick Start

MitoTNT is a Python toolkit for analyzing the temporal dynamics of mitochondrial networks.  
It integrates with **MitoGraph** for per-frame segmentation, then performs **tracking, remodeling event detection, motility analysis, feature extraction, and ChimeraX visualization**.

This page walks through a typical workflow step by step.

---

## Fast and scalable

The core pipeline has been re-engineered to process large 4D datasets efficiently while producing results numerically equivalent to the original implementation:

- **Multiprocessing** — graph extraction, topology-cost computation, and gap closing fan out across CPU cores (fork-based and copy-on-write, so the large graph objects are shared rather than pickled to each worker).
- **Batched graph construction** — skeleton graphs are assembled with single bulk `add_vertices`/`add_edges` calls instead of per-node incremental edits, turning the previously O(N²) extraction into near-linear time (≈20× faster at ~140k nodes/frame).
- **KD-tree spatial search** — `scipy.spatial.cKDTree` replaces brute-force O(N²) nearest-neighbor and coordinate-matching scans for candidate links and gap closing.
- **Sparse linear assignment** — frame-to-frame and gap-closing linking are solved on sparse augmented cost matrices with the Jonker–Volgenant solver (`lap.lapmod`) instead of allocating dense matrices, drastically reducing memory on large frames.
- **Vectorized lookups** — dict/set-based mappings throughout remove repeated linear scans and per-frame DataFrame re-filtering.

In practice the full pipeline (tracking → remodeling → motility) processes a **90-frame movie with ~1,400 nodes per frame in ~85 s on a multi-core workstation**, and scales to substantially larger datasets. Worker counts are configurable via the `num_workers` argument.

---

## 1. Prepare your data

MitoTNT expects your raw time-lapse movie as a sequence of 3D image stacks (TIFF).  
Each frame in time should be placed in its own folder, with a single `.tif` inside:

```
frame_000/frame_000.tif
frame_001/frame_001.tif
frame_002/frame_002.tif
...

```

This format makes it easy for MitoGraph to write outputs back into each frame’s folder.  
Make sure you have MitoGraph installed and accessible in your PATH.

---

## 2. Segment mitochondria with MitoGraph

MitoGraph takes each 3D stack and generates:
- Surface meshes `.vtk` for visualization.
- Supporting `.gnet`, `.coo`, and `.txt` files describing network structure.

Run it through MitoTNT:

```python
import mitotnt

seg = mitotnt.MitoSegmenter('D:/GitHub/MitoTNT/test_data/mitograph', 
                            xy_pxl_size=0.145, z_pxl_size=0.145) # add image pixel sizes in microns/pixel

seg.run_mitograph_cpu(max_workers=10) # specify the number of threads
```

* `xy_pxl_size` and `z_pxl_size` are **microns per pixel**, important for scaling.
* By default, results are saved inside each frame folder.

You can preview segmentation results in **ChimeraX**:

```python
seg.visualize_segmented_data()
```
This generates a `.cxc` script for visual quality check of raw TIFFs alongside segmented surfaces.


## 3. Extract skeleton graphs

To prepare for tracking, you need to build graphs from MitoGraph outputs:

```python
skel = mitotnt.SkeletonizedMito(seg)
skel.extract_graphs()
```
These graphs capture mitochondria location and connectivity and are the inputs to tracking.

---

## 4. Track mitochondrial networks

Tracking links network fragments across timepoints using both spatial and network information:

```python
tracker = mitotnt.NetworkTracker(skel, frame_interval=3.3)  # frame interval in seconds per frame
tracked = tracker.run()
```

* The algorithm maps nodes across frames, building trajectories for each mitochondrial fragment.
* Results are written to: `mitotnt/mito_node_tracks.csv` (all node trajectories).

If you’ve already run once, reload results without recomputing to continue with the next sections:

```python
tracked = tracker.reload_results()
```

You can also restrict tracking to a subset of frames, or skip frames (useful for slower processes or faster processing).

```python
tracker = mitotnt.NetworkTracker(skel, frame_interval=3.3,
                                 start_frame=10, end_frame=30,
                                 tracking_interval=2)
```
---

## 5. Visualize in ChimeraX

MitoTNT can produce ready-to-load visualization scripts for **UCSF ChimeraX**:

```python
viz = mitotnt.Visualizer(tracked)
viz.generate_chimerax_skeleton()  # skeleton structure
viz.generate_tracking_arrows()  # tracking arrows
viz.visualize_tracking()  # build full .cxc script
```

This creates a script:

```
mitotnt/visualization/visualize_tracking.cxc
```

Open it in ChimeraX to interactively explore:

* **Raw fluorescence cloud** rendered for context.
* **Network skeletons** across time.
* **Tracking arrows** showing mitochondrial network movement.
---

## 6. Analyze remodeling and motility

Once tracking is complete, MitoTNT can quantify key mitochondrial behaviors:

```python
# Detect network remodeling (fusion & fission events)
tracked.detect_fusion_fission()

# Compute motility at three levels:
# node (voxel), segment (linear branch), fragment (connected component)
tracked.compute_mitochondria_diffusivity() 
```

Each analysis writes a CSV file:

* `remodeling_events.csv` → list of detected fusions/fissions.
* `node_diffusivity.csv`, `segment_diffusivity.csv`, `fragment_diffusivity.csv` → quantitative mobility measures.

These files can be imported for downstream analyses and plotting.

---

## 7. Extract and plot features

`FeatureExtractor` aggregates the per-frame results from the previous steps into tidy tables for statistics and plotting. Run it after step 6, since it reads the diffusivity and remodeling outputs.

```python
fe = mitotnt.FeatureExtractor(tracked)  # or FeatureExtractor(save_path='.../mitotnt/')

fe.compute_graph_features(save_csv=True)     # topology: lengths, widths, tortuosity, branching, efficiency, betweenness, ...
fe.compute_diffusivity(save_csv=True)        # node/segment/fragment diffusivity (well-fit tracks only)
fe.compute_remodeling_rates(save_csv=True)   # per-frame fusion and fission rates

# Plot distributions for any computed category
fe.plot_features_as_histogram()   # also: plot_features_as_violinplot(), plot_features_as_boxplot()
```

Each step writes a CSV:

* `graph_features.csv` → network topology features per frame.
* `motility_features.csv` → diffusivity values by level.
* `remodeling_features.csv` → fusion/fission rates per frame.

---


## Summary

1. Organize your movie into per-frame 3D stacks.
2. Run MitoGraph via `MitoSegmenter` to obtain skeletons.
3. Extract graphs with `SkeletonizedMito`.
4. Track fragments across frames using `NetworkTracker`.
5. Visualize skeletons and dynamics in ChimeraX with `Visualizer`.
6. Analyze remodeling events and motility with `TrackedMito`.
7. Aggregate and plot features with `FeatureExtractor`.
