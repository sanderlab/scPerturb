![Website](https://img.shields.io/website?down_color=red&down_message=offline&label=scperturb.org&up_message=online&url=http%3A%2F%2Fprojects.sanderlab.org%2Fscperturb%2F)
![GitHub issues](https://img.shields.io/github/issues-raw/sanderlab/scperturb)
![GitHub last commit](https://img.shields.io/github/last-commit/sanderlab/scperturb)

# scPerturb: A resource and a python tool for single-cell perturbation data
For the publication see: [Peidli, S., Green, T.D., et al. Nature Methods (2024)](https://www.nature.com/articles/s41592-023-02144-y).

## Where to find the data
The datasets are available to download on [scperturb.org](https://scperturb.org/) (where you can also find an [interactive table](http://projects.sanderlab.org/scperturb/datavzrd/scPerturb_vzrd_v1/dataset_info/index_1.html) of all included datasets). The latest versions are available in full on Zenodo, depending on the modality you are interested in:
- [RNA data](https://zenodo.org/records/13350497)
- [ATAC data](https://zenodo.org/record/7058382)

## scperturb for python (integrates with scanpy)
![PyPI - Downloads](https://img.shields.io/pypi/dm/scperturb?label=PyPI%20downloads)

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

## scperturbR for R (integrates with Seurat)
[![R-CMD-check](https://github.com/sanderlab/scPerturb/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/sanderlab/scPerturb/actions/workflows/R-CMD-check.yaml)

We wrote an R version of scperturb that works with Seurat objects. You can find it as [scperturbR on CRAN](https://cran.r-project.org/package=scperturbR). A basic usage vignette is WIP.
Install using:
```
install.packages('scperturbR')
```

## Reproducibility
Instructions to run the code to reproduce the figures and tables in the paper and supplement:
- install conda if necessary (also check out mamba, it's way faster)
- run "conda env create -f sc_env.yaml" to create a new conda environment with all the necessary packages to run the code (**you will need this**)
- activate the environment with "conda activate sc_env"
- download all the datasets from [scperturb.org](https://scperturb.org/)
- in "config.yaml", change paths, especially the one to the directory where the data was downloaded to
- run the notebooks in "notebooks" to produce the figures and tables found in the paper / supplement and the website (each saved to a different folder)


