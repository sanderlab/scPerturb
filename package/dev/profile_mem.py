# run as
# mprof run profile_mem.py
# mprof plot
# or
# python -m memory_profiler profile_mem.py
import scanpy as sc
import matplotlib.pyplot as pl
import anndata as ad
import pandas as pd
import numpy as np
import seaborn as sns
print('Scanpy version:', sc.__version__)

import sys
sys.path.append("..")
from src.scperturb import *


from memory_profiler import profile
@profile
def exec_fun():
    estats = edist(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean')

if __name__ == '__main__':
    # !wget https://zenodo.org/record/7041849/files/DatlingerBock2021.h5ad
    adata = sc.read('../notebooks/DatlingerBock2021.h5ad')

    if 'processed' in adata.uns.keys():
        print('The dataset is already processed. Skipping processing...')
    else:
        adata.layers['counts'] = adata.X.copy()

        # basic qc and pp
        sc.pp.filter_cells(adata, min_counts=1000)
        sc.pp.normalize_per_cell(adata)
        sc.pp.filter_genes(adata, min_cells=50)
        sc.pp.log1p(adata)

        # high class imbalance
        #adata = equal_subsampling(adata, 'perturbation', N_min=50)
        sc.pp.filter_genes(adata, min_cells=3)  # sanity cleaning

        # select HVGs
        n_var_max = 2000  # max total features to select
        sc.pp.highly_variable_genes(adata, n_top_genes=n_var_max, subset=False, flavor='seurat_v3', layer='counts')
        sc.pp.pca(adata, use_highly_variable=True)
        sc.pp.neighbors(adata)

        adata.uns['processed'] = True
    print('starting profiling...')
    import cProfile
    cProfile.run("estats = edist(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean')", 'cpstats.prof')
    exec_fun()