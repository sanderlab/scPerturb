import scanpy as sc
import matplotlib.pyplot as pl
import numpy as np
import pandas as pd
import seaborn as sns
import sys
# scperturb functions
sys.path.insert(1, '../../package/src/')
from scperturb import edist_to_control

def leftsided_chebyshev_nodes(N):
    # Takes negative chebyshev nodes of first kind, forces 0 to be included, then adds 1 to all.
    return np.polynomial.chebyshev.chebpts1(N*2+1)[:(N+1)] + 1

adata = sc.read(snakemake.input[0])
control = 'control'
groupby = 'perturbation'

sigs = {}
deltas = {}
tests = {}

# Define sample points ("nodes")
ns_ = leftsided_chebyshev_nodes(snakemake.params.iterations)
ns_ = np.array(ns_ * adata.n_vars, dtype=int)
ns_ = np.unique(ns_)
ns = ns_[ns_>0]
n_max = np.max(ns)

eds = {}
for n in ns:
    tdata = adata.copy()
    try:
        sc.pp.highly_variable_genes(tdata, n_top_genes=int(n), flavor='seurat_v3', layer='counts')
    except:
        print('Failed to compute HVGs for n=', n)
        continue
    sc.pp.pca(tdata, use_highly_variable=True)
    eds[n] = edist_to_control(tdata, groupby, control='control', verbose=False)

sdf = pd.concat([
    pd.concat(eds)
], axis=1).reset_index()
sdf.columns = ['n', groupby, 'edist']

# write to file
sdf.to_csv(snakemake.output[0])

# plot
with sns.axes_style('whitegrid'):
    fig, ax = pl.subplots(figsize=(6,4), dpi=120)
lp = sns.lineplot(data=sdf, x='n', y='edist', hue='perturbation', ax=ax, marker='o')
lp.legend_.remove()
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('Number of top genes ordered by variance used for PCA')
ax.set_ylabel('E-distance to unperturbed')
ax.axvline(2000, linestyle='--', zorder=-1, color='grey', linewidth=3)
ax.set_title(f'E-distance and Feature Selection\ndataset: {snakemake.wildcards.dataset}')
pl.savefig(snakemake.output[1], bbox_inches='tight')
pl.close()