import pandas as pd
import scanpy as sc
import numpy as np
import sys

from scipy.io import mmread
from scipy.sparse import csr_matrix, vstack
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

X = mmread(TEMPDIR / 'SunshineHein2023/matrix.mtx')
obs = pd.read_csv(TEMPDIR / 'SunshineHein2023/barcodes.tsv.gz', index_col=0, sep='\t', names=['cell_barcode'])
var = pd.read_csv(TEMPDIR / 'SunshineHein2023/features.tsv.gz', index_col=1, sep='\t', names=['ensembl_id', 'gene_symbol', 'feature_type'])
ids = pd.read_csv(TEMPDIR / 'SunshineHein2023/cell_identities.csv', index_col=0)

adata = sc.AnnData(csr_matrix(X).T, pd.concat([obs, ids], axis=1), var)
adata.var.drop('feature_type', axis=1, inplace=True)  # trivial
adata.var_names_make_unique()

# move non-gene features to obsm
group1 = adata.var.index.str.startswith('SCV_')
non_genes = list(adata.var_names[group1])
adata.obsm['SCOV_expression'] = pd.DataFrame(adata[:, non_genes].X.A, index=adata.obs_names, columns=non_genes)
adata = adata[:, ~group1].copy()

group2 = adata.var.index.str.startswith('lenti_')
non_genes = list(adata.var_names[group2])
adata.obsm['lentivirus_capture'] = pd.DataFrame(adata[:, non_genes].X.A, index=adata.obs_names, columns=non_genes)
adata = adata[:, ~adata.var.index.str.startswith('lenti_')].copy()

# harmonize metadata
adata.obs['perturbation'] = ['_'.join(np.unique([y.split('_')[0] for y in x.split(';')])).replace('non-targeting', 'control').replace('control_','').replace('_control','') if type(x)==str else None for x in adata.obs.guide_identity]
adata.obs = adata.obs.rename({'guide_identity': 'guide_id', 
                              'number_of_guides': 'nperts', 
                             }, axis=1)
adata.obs['perturbation_type'] = 'CRISPR-cas9'
adata.obs['disease'] = "lung adenocarcinoma and SARS-CoV-2"
adata.obs['cancer'] = True
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "Calu-3"
adata.obs["celltype"] = 'lung epithelial cells'
adata.obs['organism'] = 'human'

annotate_qc(adata)
assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')


