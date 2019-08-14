# mcdiff

Monte Carlo method to obtain diffusion constants

mcdiff is a Python program for obtaining diffusion and free energy profiles from molecular dynamics trajectories.
It uses Bayesian Analysis and Monte-Carlo importance sampling to determine those profiles and in turn calculate membrane
permeabilities via the Inhomogeneous Solubility Diffusion model.

mcdiff uses the following third-party libraries

- numpy >= 1.16.0
- scipy >= 1.2.0
- matplotlib >= 2.0.2

## Getting Started
### System Dependencies

- Linux or MacOS 
- python >= 3.5

### Install the Requirements 

Installation via anaconda

```
conda install "numpy>=1.16.0" "scipy>=1.2.0" "matplotlib>=2.0.2"
```

or pip

```
pip install "numpy>=1.16.0" "scipy>=1.2.0" "matplotlib>=2.0.2"
```

### Download and Install mcdiff

```
git clone https://github.com/annekegh/mcdiff.git
cd mcdiff
python setup.py install
```

### Executables and Examples

The setup process installs the following terminal commands 

```run-mcdiff```

The main command for running the Bayesian Analysis. For seeing the available options, 
type `run-mcdiff --help`. For examples on how to create and format the transition matrix,
take a look into the examples directory.

```mcdiff-infty```

A command to extrapolate profiles to infinite lag time.

```plotresults```

A command to plot the final profiles.

### License

MIT License. 
See LICENSE file.

### Authors
- An Ghysels, University of Ghent
- Gerhard Hummer, MPI for Biophysics

