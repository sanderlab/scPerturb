import scanpy as sc
import sys
# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import edist_to_control

adata = sc.read(snakemake.input[0])
tab = edist_to_control(adata, 'perturbation', obsm_key='X_pca', control='control', flavor=1)
tab.to_csv(snakemake.output[0])
