# Installation
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
