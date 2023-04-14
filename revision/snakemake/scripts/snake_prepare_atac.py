from pathlib import Path
import scanpy as sc
import numpy as np
import sys
# project utils
sys.path.insert(1, '../../')
from utils import read_from_singles

path = Path(snakemake.input[0]).parents[0]
adata = read_from_singles(path)
adata.var_names = adata.var_names.astype(str)
sc.pp.filter_genes(adata, min_cells=0)
if np.max(adata.X) > 100:
    print('logging data')
    sc.pp.log1p(adata)
if adata.n_vars > 2000:
    print('using HVG')
    sc.pp.highly_variable_genes(adata, n_top_genes=2000)
if adata.n_vars > 50:
    print('using pca')
    sc.pp.pca(adata)
    adata.obsm[f'X_pca_{snakemake.wildcards.feature}'] = adata.obsm['X_pca'].copy()
else:
    print('using direct')
    sc.pp.scale(adata)
    adata.obsm[f'X_raw_{snakemake.wildcards.feature}'] = adata.X.copy()

adata.write(snakemake.output[0])