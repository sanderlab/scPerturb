import pandas as pd
import scanpy as sc
import numpy as np
import sys

from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

tab_DM = pd.read_csv(TEMPDIR / f'LiangWang2023/OSCAR_DM_expression_matrix.csv', index_col=0)
tab_EM = pd.read_csv(TEMPDIR / f'LiangWang2023/OSCAR_EM_expression_matrix.csv', index_col=0)
meta_df = pd.read_csv(TEMPDIR / f'LiangWang2023/OSCAR_metadata.csv', index_col=0)

tab = pd.concat([tab_DM, tab_EM], axis=1)

adata = sc.AnnData(csr_matrix(tab.T.values))
adata.obs_names = tab.columns
adata.var_names = tab.index
adata.obs = meta_df.loc[adata.obs_names]

# Harmonize Metadata
adata.obs_names = [x.split('-')[0] for x in adata.obs_names]
adata.var.index.name='gene_symbol'
adata.obs.drop(['cell'], axis=1, inplace=True)
adata.obs.index.name = 'cell_barcode'
adata.obs.rename({
    'barcode': 'guide_id',
    'sgrna': 'guide_sequence',
    'gene': 'perturbation',
    'type': 'medium'
}, axis=1, inplace=True)
adata.obs.perturbation[adata.obs.perturbation=='Non-Targeting'] = 'control'
adata.obs.medium.replace({'DM': 'Differentiation Medium', 'EM': 'Expansion Medium'}, inplace=True)
adata.obs['perturbation_type'] = 'CRISPR-cas9'
adata.obs['nperts']= 1 - adata.obs.perturbation.str.count('control')
adata.obs['organism'] = 'mouse'
adata.obs['tissue_type'] = 'organoid'
adata.obs['disease'] = 'healthy'
adata.obs['celltype'] = 'hepatocyte'
adata.obs['cancer'] = False
adata.obs['cell_line'] = 'intrahepatic cholangiocyte organoids'
annotate_qc(adata)
assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')
