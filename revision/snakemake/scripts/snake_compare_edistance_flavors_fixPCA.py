import scanpy as sc
import matplotlib.pyplot as pl
import pandas as pd
import numpy as np
import seaborn as sns
import sys
from tqdm import tqdm

# scperturb import
sys.path.insert(1, '../../package/src/')
from scperturb import edist_to_control

def leftsided_chebyshev_nodes(N):
    # Takes negative chebyshev nodes of first kind, forces 0 to be included, then adds 1 to all.
    return np.polynomial.chebyshev.chebpts1(N*2+1)[:(N+1)] + 1

def subsampling(adata, obs_key, N):
    '''
    Subsample to same class sizes. Classes given by obs_key pointing to categorical in adata.obs.
    Discards any class with less than N entries.
    '''
    counts_per_class = adata.obs[obs_key].value_counts()
    selected_classes = counts_per_class.index[counts_per_class>=N]
    if len(selected_classes) == 0:
        return None
    # subsample indices per group defined by obs_key
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=False) for group in selected_classes]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()

adata = sc.read_h5ad(snakemake.input[0])
adata.layers['counts'] = adata.X.copy()
# basic qc and pp
sc.pp.filter_cells(adata, min_counts=1000)
sc.pp.normalize_per_cell(adata)
sc.pp.filter_genes(adata, min_cells=50)
sc.pp.log1p(adata)
# select HVGs
sc.pp.highly_variable_genes(adata, n_top_genes=2000, subset=False, flavor='seurat_v3', layer='counts')
sc.pp.pca(adata, use_highly_variable=True)

# at each N, subsample multiple times to compute the expected value of e-distance
min_cells = 10
res = []
nodes = leftsided_chebyshev_nodes(snakemake.params.N_eval_points)
eval_points = np.array(snakemake.params.N_max * nodes, dtype=int) + min_cells  # Stützstellen
for N in tqdm(eval_points):
    Z = []
    for i in range(snakemake.params.N_resample):
        sdata = subsampling(adata, 'perturbation', N=N)
        if sdata is None:
            continue

        edtc_0 = edist_to_control(sdata, flavor=0, verbose=False)  # uncorrected
        edtc_1 = edist_to_control(sdata, flavor=1, verbose=False)  # naive_correction
        edtc_2 = edist_to_control(sdata, flavor=2, verbose=False)  # rizzo_correction
        edtc_3 = edist_to_control(sdata, flavor=3, verbose=False)  # delta_correction

        tab = pd.concat([edtc_0, edtc_1, edtc_2, edtc_3],axis=1)
        tab.columns = ['uncorrected', 'naive_correction', 'rizzo_correction', 'delta_correction']
        tab['N_cells'] = N
        tab['resample_id'] = i
        Z.append(tab)
    if len(Z)>0:
        res.append(Z)
    
# merge
df = pd.concat([item for sublist in res for item in sublist], axis=0).reset_index()
df['log10_cells'] = np.log10(df.N_cells)  # for better color range
df.to_csv(snakemake.output[0])

# expected value under resampling
means = df.groupby(['perturbation', 'N_cells']).mean().reset_index()
sdevs = df.groupby(['perturbation', 'N_cells']).std().reset_index()

# plot with error bars (std)
flavors = ['uncorrected', 'naive_correction', 'rizzo_correction', 'delta_correction']
M = len(flavors)
with sns.axes_style('whitegrid'):
    fig, axs = pl.subplots(1, M, figsize=[4*M, 3], dpi=100)
for ax, key in zip(axs, flavors):
    for p in df.perturbation.unique():
        mask = means.perturbation==p
        ax.errorbar(means[mask].N_cells, means[mask][key], yerr=sdevs[mask][key], linewidth=2, capsize=3)
    ax.set_xscale('log')
    ax.set_title(f'{key} (fixed PCA) (N={snakemake.params.N_resample})')
pl.tight_layout()
pl.savefig(snakemake.output[1], bbox_inches='tight')
pl.close()
