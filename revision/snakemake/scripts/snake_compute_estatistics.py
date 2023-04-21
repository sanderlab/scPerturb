import scanpy as sc
import numpy as np
import sys

# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import etest

adata = sc.read(snakemake.input[0])
df = etest(adata, runs=snakemake.params.iterations, n_jobs=snakemake.threads)

# export table
df.to_csv(snakemake.output[0])