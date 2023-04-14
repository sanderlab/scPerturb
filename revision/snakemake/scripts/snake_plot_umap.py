import scanpy as sc
import scvelo as scv
import matplotlib.pyplot as pl
import sys
# project utils
sys.path.insert(1, '../../')
from utils import pseudo_bulk

adata = sc.read(snakemake.input[0])
title = snakemake.wildcards.dataset+' single-cell umap'
if len(adata.obs.perturbation.unique()) > 20:
    adata = pseudo_bulk(adata, ['perturbation'])
    sc.pp.normalize_per_cell(adata)
    sc.pp.log1p(adata)
    sc.pp.pca(adata)
    sc.pp.neighbors(adata)
    adata.obs['perturbation'] = adata.obs['perturbation'].astype('category')
    title = snakemake.wildcards.dataset+' pseudobulk umap'
sc.tl.umap(adata)
scv.pl.scatter(adata, color='perturbation', show=False, dpi=120, legend_loc='right margin')
pl.savefig(snakemake.output[0], bbox_inches='tight', dpi=120)
pl.close()