# Code repository for scPerturb

## Instructions to run the code
- install conda if necessary (also check out mamba, it's way faster)
- run "conda env create -f sc_env.yaml" to create a new conda environment with all the necessary packages to run the code (**you will need this**)
- activate the environment with "conda activate sc_env"
- download all the datasets from [scperturb.org](https://scperturb.org/)
- in "config.yaml", change paths, especially the one to the directory where the data was downloaded to
- run the notebooks in "notebooks" to produce the figures and tables found in the paper / supplement and the website (each saved to a different folder)

## Where to find the data
The datasets are available to download on [scperturb.org](https://scperturb.org/) and alternatively on zenodo (in case the website does not work):
- [RNA data](https://zenodo.org/record/7041849)
- [ATAC data](https://zenodo.org/record/7058382)
