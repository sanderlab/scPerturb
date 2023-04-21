import scanpy as sc
import sys
# project utils
sys.path.insert(1, '../../')
from utils import equal_subsampling, random_split

print('Reading data...')
adata = sc.read(snakemake.input[0])

# basic qc and pp
print('Shape before filtering: ', adata.shape)
sc.pp.filter_cells(adata, min_counts=snakemake.params.min_counts)
adata = equal_subsampling(adata, 'perturbation', N_min=snakemake.params.min_perturb_size)
adata.layers['counts'] = adata.X.copy()
print('Normalizing counts...')
sc.pp.normalize_per_cell(adata)
sc.pp.filter_genes(adata, min_cells=snakemake.params.min_cells)
print('Applying log1p...')
sc.pp.log1p(adata)
print('Shape after filtering: ', adata.shape)

# select HVGs
print('Selecting highly variable genes...')
n_var_max = 2000  # max total features to select
sc.pp.highly_variable_genes(adata, n_top_genes=snakemake.params.n_hvgs,
                            subset=False, flavor='seurat_v3', layer='counts')
print('Running PCA...')
sc.pp.pca(adata, use_highly_variable=True, n_comps=snakemake.params.n_pcs)

# annotate split for evaluation
random_split(adata)

print('Writing to file...')
adata.write(snakemake.output[0])
print('Control in perturbation column:', 'control' in adata.obs.perturbation)