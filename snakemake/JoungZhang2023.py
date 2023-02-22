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

def adapt(adata):
    adata.obs.rename({'n_genes': 'ngenes', 'n_counts': 'ncounts', 'TF': 'perturbation'}, axis=1)
    adata.obs['disease'] = "healthy"
    adata.obs['cancer'] = False
    adata.obs['tissue_type']="cell_line"
    adata.obs["cell_line"] = "hESCs"
    adata.obs["celltype"] = 'hESCs'
    adata.obs['organism'] = 'human'
    adata.obs['nperts'] = [p.count('+') + 1 - p.count('GFP') for p in adata.obs.perturbation]
    family_genes = np.logical_or(*[adata.var_names.str.startswith(k) for k in ['RPS', 'RPL']])
    adata.obs['percent_ribo'] = np.ravel(np.sum(adata[:, family_genes].X, axis=1) * 100 / np.sum(adata.X, axis=1))
    adata.obs.index.name = 'cell_barcode'
    return adata

### Combinatorial screen ###
adata = sc.read_h5ad(snakemake.input['combinatorial'])
# their adata.X was scaled and their adata.raw.X was saved as dense...
adata.X = csr_matrix(adata.raw.X) 
adata.obs['perturbation'] = ['+'.join(sorted(re.sub(r'TFORF[0-9]{4}-', ' ', s).replace(' ','').split(','))) for s in adata.obs.perturbation]  # cleaner names
adata = adapt(adata)
assert_annotations(adata)
adata.write(snakemake.output['combinatorial'], compression='gzip')

### TF atlas ###
# adata = sc.read_h5ad(snakemake.input['atlas'])
f = h5py.File(snakemake.input['atlas'], "r")
# We load the original dense data per chunk and convert each to sparse (csr) right on the spot
# This avoids loading all data as dense into memory (too big)
chunk_size = 10000
iterations = int(np.ceil(f['X'].shape[0] / 10000) + 1)
chunks = [csr_matrix(f['X'][i*chunk_size:(i+1)*chunk_size]) for i in tqdm(range(iterations))]
X = vstack(chunks)
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
adata = sc.AnnData(X=X, obs=obs, var=var)
adata.obs['perturbation'] = [re.sub(r'TFORF[0-9]{4}-', ' ', s).replace(' ','') for s in adata.obs.perturbation]  # cleaner names
adata = adapt(adata)
assert_annotations(adata)
adata.write(snakemake.output['atlas'], compression='gzip')
