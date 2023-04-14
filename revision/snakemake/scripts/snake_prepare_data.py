import scanpy as sc
import sys
# project utils
sys.path.insert(1, '../../')
from utils import equal_subsampling, random_split

adata = sc.read(snakemake.input[0])
adata.layers['counts'] = adata.X.copy()

# basic qc and pp
sc.pp.filter_cells(adata, min_counts=snakemake.params.min_counts)
sc.pp.normalize_per_cell(adata)
sc.pp.filter_genes(adata, min_cells=snakemake.params.min_cells)
sc.pp.log1p(adata)
print('Shape after filtering: ', adata.shape)

# subsample against high class imbalance
adata = equal_subsampling(adata, 'perturbation', N_min=snakemake.params.min_perturb_size)
sc.pp.filter_genes(adata, min_cells=3)  # sanity cleaning
print('Shape after equal-subsampling: ', adata.shape)

# select HVGs
n_var_max = 2000  # max total features to select
sc.pp.highly_variable_genes(adata, n_top_genes=snakemake.params.n_hvgs,
                            subset=False, flavor='seurat_v3', layer='counts')
sc.pp.pca(adata, use_highly_variable=True, n_comps=snakemake.params.n_pcs)
sc.pp.neighbors(adata)

# annotate split for evaluation
random_split(adata)

adata.write(snakemake.output[0])
print('Control in perturbation column:', 'control' in adata.obs.perturbation)