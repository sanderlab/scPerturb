import scanpy as sc
import matplotlib.pyplot as pl
import sys
import numpy as np

# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import edist
# project utils
sys.path.insert(1, '../../')
from utils import cluster_matrix, plot_heatmap

adata = sc.read(snakemake.input[0])
tab = edist(adata, 'perturbation', obsm_key='X_pca', flavor=1)
tab.to_csv(snakemake.output[1])

# plot
tab = (1/tab).replace([np.inf, -np.inf], 0, inplace=False)
tab = cluster_matrix(tab, 'both')
plot_heatmap(tab, snakemake.wildcards.dataset+' edistances')
pl.savefig(snakemake.output[0], bbox_inches='tight', dpi=120)
pl.close()