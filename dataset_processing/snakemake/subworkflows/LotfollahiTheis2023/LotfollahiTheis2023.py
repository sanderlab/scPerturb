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

obs = pd.read_csv(TEMPDIR / f'LotfollahiTheis2023/GSE206741_cell_metadata.tsv', index_col=0, sep='\t')
var = pd.read_csv(TEMPDIR / f'LotfollahiTheis2023/GSE206741_gene_metadata.tsv', index_col=0, sep='\t')
X = csr_matrix(mmread(TEMPDIR / f'LotfollahiTheis2023/GSE206741_count_matrix.mtx'))

adata = sc.AnnData(X.T, obs, var)

adata.var.set_index('gene_short_name', inplace=True, drop=False)
adata.var.columns = ['ensembl_id']
adata.var.index.name = 'gene_symbol'

adata.obs['perturbation'] = ['_'.join(np.sort([d1,d2])).replace('\xa0','') for d1, d2 in zip(adata.obs.Drug1, adata.obs.Drug2)]
adata.obs['perturbation'] = [x.replace('DMSO_', '').replace('_DMSO', '').replace('DMSO', 'control') for x in adata.obs.perturbation]
adata.obs.index.name = 'cell_barcode'
adata.obs.rename({
    'n.umi': 'ncounts', 
}, axis=1, inplace=True)
adata.obs.drop(['sample', 'Drug1', 'Drug2'], axis=1, inplace=True)
adata.obs = adata.obs[['perturbation', 'Size_Factor', 'ncounts', 'RT_well', 'Well']]
adata.obs['nperts'] = [p.count('_')+1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
adata.obs['perturbation_type'] = 'drug'
adata.obs['disease'] = "lung adenocarcinoma"
adata.obs['cancer'] = True
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "A549"
adata.obs["celltype"] = 'lung epthelial cells'
adata.obs['organism'] = 'human'
annotate_qc(adata, species='human')
assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')


