# Installation
On github, download this [repository](https://github.com/pylattice/mitoTNT).  

Or clone the repository on terminal:  
`git clone https://github.com/pylattice/mitoTNT.git`  

### Software requirements
- **[Jupyter Notebook](https://jupyter.org/)** or get it from **[Anaconda](https://www.anaconda.com/products/distribution)**

- **[MitoGraph](https://github.com/vianamp/MitoGraph/)** for mitochondria segmentation. Available for MacOS and Linux (needs to compile from source on Linux).

- **[ChimeraX](https://www.cgl.ucsf.edu/chimerax/)** for tracking visualization

### Python dependencies

We will create a conda environment that automatically install all the required dependencies:
Open anaconda prompt on Windows or just terminal on MacOS/Linux:
1. Go to the root directory of mitotnt repository
2. Create the enviroment using the .yml file
``
conda env create --name mitotnt --file=mitotnt_env.yml
``

Before using MitoTNT, do 
