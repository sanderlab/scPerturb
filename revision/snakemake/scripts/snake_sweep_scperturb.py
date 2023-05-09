import scanpy as sc
import pandas as pd
import numpy as np
# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import etest

# Load data
adata = sc.read(snakemake.input[0])
control = 'control'
groupby = 'perturbation'

# basic qc
sc.pp.filter_cells(adata, min_counts=1000)
adata.layers['counts'] = adata.X.copy()

def strict_equal_subsampling(adata, obs_key, N):
    '''
    Subsample to same class sizes. Classes given by obs_key pointing to categorical in adata.obs.
    '''
    counts = adata.obs[obs_key].value_counts()
    counts = counts[counts >= N]  # remove classes with less than N samples
    # subsample indices per group defined by obs_key
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=False) for group in counts.index]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()

ncounts = int(snakemake.wildcards['ncounts'])
ncells = int(snakemake.wildcards['ncells'])
npcs = int(snakemake.wildcards['npcs'])
nhvgs = int(snakemake.wildcards['nhvgs'])

tdata = adata.copy()
# subsample cells
tdata = strict_equal_subsampling(tdata, 'perturbation', N=ncells)
# subsample counts
tdata = sc.pp.downsample_counts(tdata, total_counts=ncounts*ncells, copy=True)
# preprocess
sc.pp.normalize_per_cell(tdata)
sc.pp.filter_genes(tdata, min_cells=10)
sc.pp.log1p(tdata)
# select nHVGs
sc.pp.highly_variable_genes(tdata, n_top_genes=nhvgs)
# select nPCs
sc.pp.pca(tdata, n_comps=npcs)

# compute E-statistics
# ed = edist_to_control(tdata, 'perturbation', obsm_key='X_pca', control='control')
et = etest(tdata, runs=int(snakemake.params.iterations), n_jobs=int(snakemake.threads))

# export results
et.to_csv(snakemake.output[0])