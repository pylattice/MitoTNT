# Installation
On github, download this [repository](https://github.com/pylattice/mitoTNT). Alternative clone the repository on terminal: 
`git clone https://github.com/pylattice/mitoTNT.git`  

### Software requirements
- **[Jupyter Notebook](https://jupyter.org/)** or get it from **[Anaconda](https://www.anaconda.com/products/distribution)**

- **[MitoGraph](https://github.com/vianamp/MitoGraph/)** for mitochondria segmentation. Available for MacOS. Windows and Linux binary can be built.

- **[ChimeraX](https://www.cgl.ucsf.edu/chimerax/)** for tracking visualization

### Python dependencies

We will create a conda environment that automatically install all the required dependencies:
1. Open anaconda prompt on Windows or just terminal on MacOS/Linux
2. Go to the root directory of MitoTNT repository
3. Create the enviroment using the .yml file  
``
conda env create --name mitotnt --file=mitotnt_env.yml
``

Before using MitoTNT, do `conda activate mitotnt` and then open notebook with `jupyter notebook`.
