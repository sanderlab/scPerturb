import scanpy as sc
import matplotlib.pyplot as pl
import seaborn as sns
import numpy as np
import sys

# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import etest

adata = sc.read(snakemake.input[0])
df = etest(adata, runs=snakemake.params.iterations)
df['log10_edist'] = np.log10(np.clip(df.edist, 0, np.infty))

# export table
df.to_csv(snakemake.output[0])