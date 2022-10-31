![Website](https://img.shields.io/website?down_color=red&down_message=offline&label=scperturb.org&up_message=online&url=http%3A%2F%2Fprojects.sanderlab.org%2Fscperturb%2F)
![PyPI - Downloads](https://img.shields.io/pypi/dm/scperturb?label=PyPI%20downloads)
![GitHub issues](https://img.shields.io/github/issues-raw/sanderlab/scperturb)

# scPerturb: A resource and a python tool for single-cell perturbation data

## Where to find the data
The datasets are available to download on [scperturb.org](https://scperturb.org/) and alternatively on Zenodo:
- [RNA data](https://zenodo.org/record/7041849)
- [ATAC data](https://zenodo.org/record/7058382)

## scperturb for python
A python package to compute E-distances in single-cell perturbation data and perform E-tests.

### Install
Just install via pip:

```
pip install scperturb
```

### Usage example

Check out [this notebook](https://github.com/sanderlab/scPerturb/blob/master/package/notebooks/e-distance.ipynb) for a tutorial.
Basic usage is:
```
# E-distances
estats = edist(adata, obs_key='perturbation')
# E-distances to a specific group (e.g. 'control')
estats_control = estats.loc['control']
# E-test for difference to control
df = etest(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean', control='control', alpha=0.05, runs=100)
```

## Instructions to run the code to recreate the manuscript figures
- install conda if necessary (also check out mamba, it's way faster)
- run "conda env create -f sc_env.yaml" to create a new conda environment with all the necessary packages to run the code (**you will need this**)
- activate the environment with "conda activate sc_env"
- download all the datasets from [scperturb.org](https://scperturb.org/)
- in "config.yaml", change paths, especially the one to the directory where the data was downloaded to
- run the notebooks in "notebooks" to produce the figures and tables found in the paper / supplement and the website (each saved to a different folder)


