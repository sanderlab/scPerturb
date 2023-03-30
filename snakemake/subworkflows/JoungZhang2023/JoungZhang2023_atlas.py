import pandas as pd
import scanpy as sc
import numpy as np
import sys
import re
import h5py

from scipy.sparse import csr_matrix, vstack

# Custom functions
sys.path.insert(1, '../')
from utils import *

### TF atlas ###
print('TF atlas')
# adata = sc.read_h5ad(snakemake.input['atlas'])
f = h5py.File(snakemake.input['atlas'], "r")
# We load the original dense data per chunk and convert each to sparse (csr) right on the spot
# This avoids loading all data as dense into memory (too big)
chunk_size = 10000
iterations = int(np.ceil(f['X'].shape[0] / chunk_size) + 1)
chunks = [csr_matrix(f['X'][i*chunk_size:(i+1)*chunk_size]) for i in tqdm(range(iterations))]
X = vstack(chunks)
print('X.shape', X.shape)
# sparsity is high (~95% of entries are empty)
print(f'Sparsity: {np.round(100* X.count_nonzero() / np.product(X.shape), 2)}%')
# obs
keys = [k for k in f['obs'].keys() if not k.startswith('__')]
obs = pd.DataFrame(np.vstack([np.array(f['obs'][k]) for k in keys]).T, columns=keys)
obs.set_index('_index', inplace=True)
obs.index = obs.index.astype(str)
obs.index.name='cell_barcode'
# var
keys = [k for k in f['var'].keys() if not k.startswith('__')]
var = pd.DataFrame(np.vstack([np.array(f['var'][k]) for k in keys]).T, columns=keys)
var.set_index('_index', inplace=True)
var.index = var.index.astype(str)
var.index.name='gene_symbol'

# assemble
adata = sc.AnnData(X=X, obs=obs, var=var)

# set metadata
# TODO: Perturbation here is just a number. Find out what it means!!!
adata.obs = adata.obs.rename({'n_genes': 'ngenes', 'n_counts': 'ncounts', 'TF': 'perturbation'}, axis=1)
adata.obs['disease'] = "healthy"
adata.obs['cancer'] = False
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "hESCs"
adata.obs["celltype"] = 'hESCs'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'ORF overexpression'
adata.obs['nperts'] = 1
annotate_qc(adata, species='human')
adata.obs.index.name = 'cell_barcode'

assert_annotations(adata)
# anndata can not handle these
adata.obs.perturbation = adata.obs.perturbation.astype(str)
adata.obs.batch = adata.obs.batch.astype(str)
adata.obs.louvain = adata.obs.louvain.astype(str)
adata.obs.ncounts = adata.obs.ncounts.astype(int)
adata.obs.ngenes = adata.obs.ngenes.astype(int)
adata.obs.percent_mito = adata.obs.percent_mito.astype(float)
adata.write(snakemake.output['atlas'], compression='gzip')
