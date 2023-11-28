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

def read_count_matrix(key):
    X = csr_matrix(mmread(TEMPDIR / f'WesselsSatija2023/{key}.matrix.mtx').T)
    obs = pd.read_csv(TEMPDIR / f'WesselsSatija2023/{key}.barcodes.tsv', sep='\t', index_col=0, names=['cell_barcode'])
    var = pd.read_csv(TEMPDIR / f'WesselsSatija2023/{key}.features.tsv', sep='\t', index_col=1, names=['ENSEMBL_ID', 'gene_symbol', 'feature_type'])
    return X, obs, var

adatas = []
for j in ['1', '2', '3', '4']:
    X, obs, var = read_count_matrix(f'THP1-CaRPool-seq_and_HEK293FTstabRNA.GEXGDO{j}')
    adata = sc.AnnData(X, obs=obs, var=var)
    adata.obs.index = [x[:-2] for x in obs.index]
    adata.obs.index = [f'L{j}_{x}' for x in adata.obs.index]
    
    # Add ADT
    X_ADT, obs_ADT, var_ADT = read_count_matrix(f'THP1-CaRPool-seq_and_HEK293FTstabRNA.ADT{j}')
    var_ADT = var_ADT.set_index('ENSEMBL_ID')
    var_ADT.index.name='ADT_id'
    var_ADT.feature_type = 'Surface Protein'
    obs_ADT.index = [f'L{j}_{x}' for x in obs_ADT.index]
    ADT = pd.DataFrame(X_ADT.A, index=obs_ADT.index, columns=var_ADT.index)
    diff = adata.obs.index[~adata.obs.index.isin(ADT.index)]  #  2 idxs...
    for d in diff:
        ADT.loc[d]=0
    adata.obsm['Surface_protein']=ADT.loc[adata.obs_names]

    # Add HTO
    X_HTO, obs_HTO, var_HTO = read_count_matrix(f'THP1-CaRPool-seq_and_HEK293FTstabRNA.HTO{j}')
    var_HTO = var_HTO.set_index('ENSEMBL_ID')
    var_HTO.index.name='HTO_id'
    var_HTO.feature_type = 'Hashtag Oligo'
    obs_HTO.index = [f'L{j}_{x}' for x in obs_HTO.index]
    HTO = pd.DataFrame(X_HTO.A, index=obs_HTO.index, columns=var_HTO.index)
    diff = adata.obs.index[~adata.obs.index.isin(HTO.index)]  #  2 idxs...
    for d in diff:
        HTO.loc[d]=0
    adata.obsm['Sample_Tags']=HTO.loc[adata.obs_names]
    
    adata.var_names_make_unique()  # Otherwise error
    adatas.append(adata)

# merge and add metadata
adata = sc.concat(adatas)
tab = pd.read_csv(TEMPDIR / 'WesselsSatija2023/GSE213957_THP1-CaRPool-seq.metadata.tsv', sep='\t')
tab.index = [x[:-2] for x in tab.index]
adata = adata[tab.index].copy()
adata.obs = tab

# Harmonize
adata.obs.drop(['S.Score', 'G2M.Score'], axis=1, inplace=True)
adata.obs.rename({
    'nCount_RNA': 'ncounts', 
    'nFeature_RNA': 'ngenes', 
    'percent.mt': 'percent_mito',
    'nCount_HTO': 'ncounts_HTO',
    'nCount_ADT': 'ncounts_ADT',
    'TenX.Lane': '10X_lane',
    'CRISPR.Array': 'CRISPRcas13_Array',
    'GenePair': 'perturbation',
    'Guides': 'guides',
    'Phase': 'Cell_cycle_phase'
}, axis=1, inplace=True)
adata.obs.perturbation[adata.obs.perturbation=='NT_NT'] = 'control'
adata.obs.perturbation = [x.replace('NT_','').replace('_NT', '') for x in adata.obs.perturbation]
adata.obs['nperts'] = [p.count('_')+1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
adata.obs['perturbation_type'] = 'CRISPR-cas13'
adata.obs['disease'] = "leukemia"
adata.obs['cancer'] = True
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "THP-1"
adata.obs["celltype"] = 'monocytes'
adata.obs['organism'] = 'human'
annotate_qc(adata, species='human')
adata.obs.index.name = 'cell_barcode'
assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')


