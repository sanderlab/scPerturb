import scanpy as sc
import matplotlib.pyplot as pl
import pandas as pd
import numpy as np
import seaborn as sns
import sys
from tqdm import tqdm
from sklearn.metrics import pairwise_distances

# scperturb import
sys.path.insert(1, '../../package/src/')
from scperturb import edist_to_control, equal_subsampling

def simul(mu=(2,-2), S=(5,3), D=30, N=1000, iterations=50, correction_factor=True):
    from itertools import product
    D = [D] if type(D) is int else D
    N = [N] if type(N) is int else N
    evals = [x for x in product(D, N)]
    res = []
    for d, n in tqdm(evals):
        factor = (n/(n-1)) if correction_factor else 1
        for i in range(iterations):
            X = np.random.normal(mu[0], S[0], size=(n, d))
            Y = np.random.normal(mu[1], S[1], size=(n, d))
            dxy = np.mean(pairwise_distances(X,Y, metric='sqeuclidean'))
            sx = np.mean(pairwise_distances(X,X, metric='sqeuclidean'))*factor
            sy = np.mean(pairwise_distances(Y,Y, metric='sqeuclidean'))*factor
            ed = 2*dxy-sx-sy
            res.append([ed, dxy, sx, sy, i, d, n])
    return pd.DataFrame(res, columns=['edistance', 'delta', 'sigma_x', 'sigma_y', 
                                      'iteration', 'dimension', 'sample_size'])

# Show that dimension has no impact
df = simul(mu=(2,-2), S=(5,3), D=np.arange(1,21), N=200, iterations=100)
# Note that e.g. PCA will try to represent to overall variance in the data, so
# this will look different for PCA. Each dimension added in this simulated data
# has the same variance, so the e-distance scales linearly with dimension.
pl.errorbar(np.arange(1,21), df.groupby('dimension').mean().edistance, yerr=df.groupby('dimension').std().edistance)
pl.xlabel('dimension')
pl.ylabel('edistance')
pl.title('E-distance linearly scales with dimension\nassuming each dimension has same variance')
pl.savefig(snakemake.output[0], dpi=300, bbox_inches='tight')
pl.close()

# Show that sample size has an impact when not correcting
N = [2,5,20,50,200,500,2000]
df_low = simul(mu=(2,-2), S=(5,3), D=30, N=N, iterations=50, correction_factor=False)
df_high = simul(mu=(2,-2), S=(10,50), D=30, N=N, iterations=50, correction_factor=False)
# since PCA will try to represent to overall variance in the data, this is not bad at all!
keys = ['edistance', 'delta', 'sigma_x']
fig, axs = pl.subplots(1, len(keys), figsize=(6*len(keys), 4), dpi=100)
for ax, key in zip(axs, keys):
    ax.errorbar(N, df_low.groupby('sample_size').mean()[key], 
                yerr=df_low.groupby('sample_size').std()[key], 
                label=f'Low Variances {(5,3)}')
    ax.errorbar(N, df_high.groupby('sample_size').mean()[key], 
                yerr=df_high.groupby('sample_size').std()[key], 
                label=f'High Variances {(10,50)}')
    ax.set_xlabel('Number of samples')
    ax.set_ylabel(key)
    ax.set_title(f'{key}')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
pl.tight_layout()
pl.savefig(snakemake.output[1], dpi=300, bbox_inches='tight')
pl.close()

# Show that sample size has NO impact when not correcting
N = [2,5,20,50,200,500,2000]
df_low = simul(mu=(2,-2), S=(5,3), D=30, N=N, iterations=50, correction_factor=True)
df_high = simul(mu=(2,-2), S=(10,50), D=30, N=N, iterations=50, correction_factor=True)
# since PCA will try to represent to overall variance in the data, this is not bad at all!
keys = ['edistance', 'delta', 'sigma_x']
fig, axs = pl.subplots(1, len(keys), figsize=(6*len(keys), 4), dpi=100)
for ax, key in zip(axs, keys):
    ax.errorbar(N, df_low.groupby('sample_size').mean()[key], 
                yerr=df_low.groupby('sample_size').std()[key], 
                label=f'Low Variances {(5,3)}')
    ax.errorbar(N, df_high.groupby('sample_size').mean()[key], 
                yerr=df_high.groupby('sample_size').std()[key], 
                label=f'High Variances {(10,50)}')
    ax.set_xlabel('Number of samples')
    ax.set_ylabel(key)
    ax.set_title(f'{key}')
    ax.set_yscale('log')
    ax.set_xscale('log')
    ax.legend()
pl.tight_layout()
pl.savefig(snakemake.output[2], dpi=300, bbox_inches='tight')
pl.close()
