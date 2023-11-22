import pandas as pd
import scanpy as sc
import numpy as np
import sys
from tqdm import tqdm
from pathlib import Path

from scipy.io import mmread
from scipy.sparse import csr_matrix

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

suffix = 'whole_tx'
X = csr_matrix(mmread(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}_count_matrix.mtx'))
obs = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}.Barcodes.tsv', index_col=0, names=['cell_barcode'])
var = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}.Genes.tsv', index_col=0, names=['gene_symbol'])
adata = sc.AnnData(X.T, obs, var)
suffix = 'nascent_tx'
X = csr_matrix(mmread(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}_count_matrix.mtx'))
obs = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}.Barcodes.tsv', index_col=0, names=['cell_barcode'])
var = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}.Genes.tsv', index_col=0, names=['gene_symbol'])
ndata = sc.AnnData(X.T, obs, var)
suffix = 'sgRNA'
X = csr_matrix(mmread(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}_count_matrix.mtx'))
obs = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}.Barcodes.tsv', index_col=0, names=['cell_barcode'])
var = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_{suffix}.Genes.tsv', index_col=0, names=['gene_symbol'])
sdata = sc.AnnData(X.T, obs, var)

# check alignment
assert all(adata.var_names==ndata.var_names)
assert all(adata.obs_names==ndata.obs_names)
assert all(adata.obs_names==sdata.obs_names)

adata.layers['nascent_counts'] = ndata.X.copy()
adata.obsm['sgRNA_counts'] = pd.DataFrame(sdata.X.A, index=sdata.obs_names, columns=sdata.var_names, dtype=int)
full_obs = pd.read_csv(TEMPDIR / 'XuCao2023/' / f'GSM6752591_on_target_cell_metadata.csv', index_col=0)
assert all(adata.obs_names==full_obs.index)
adata.obs = full_obs.copy()

# harmonize metadata
adata.obs.rename({'UMI_counts': 'ncounts', 'target_genes': 'perturbation', 'target': 'guide_id'}, axis=1, inplace=True)
adata.obs.perturbation.replace('NO-TARGET', 'control', inplace=True)
cols = ['perturbation', 'ncounts', 'nascent_UMI_counts', 'nascent_ratio', 'guide_id', 'gRNA_UMI_counts', 'nascent_MT_ratio', 'Cell_cycle_phase', 'whole_exon_ratio', 'new_exon_ratio']
adata.obs = adata.obs[cols]
adata.obs['target'] = adata.obs['perturbation'].astype(str)

adata.obs['perturbation_type'] = 'CRISPRi'
adata.obs['disease'] = "healthy"
adata.obs['cancer'] = True
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "HEK293"
adata.obs["celltype"] = 'embryonic kidney cells'
adata.obs['organism'] = 'human'
adata.obs['nperts'] = [p.count('_')+1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
annotate_qc(adata, species='human')
adata.obs.index.name = 'cell_barcode'

assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')