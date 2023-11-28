import pandas as pd
import scanpy as sc
import numpy as np
import sys
from tqdm import tqdm
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

def process_adata(key):
    # build adata
    adata = sc.read(TEMPDIR / f'YaoCleary2023/{key}_temp.h5ad')
    adata.obs = adata.obs.rename({'Total_RNA_count': 'ncounts', 'Total_unique_genes': 'ngenes', 'Biological_replicate': 'replicate',
                      'Percent_mitochondrial_reads': 'percent_mito', 'Guides': 'full_guides', 
                      'Guides_collapsed_by_gene': 'guides', 'Total_number_of_guides': 'nguides'}, axis=1)
    adata.var.index.name = 'gene_symbol'
    # adata.var = adata.var.set_index('features').drop('_index', axis=1)
    # h5py does not support sparse matrices yet...
    # tab = pd.read_csv(TEMPDIR / f'YaoCleary2023/{key}_perturbations.txt', sep='\t').T
    # stab = tab.astype(pd.SparseDtype("int", 0))  # make sparse
    # adata.obsm['barcodes'] = stab

    # harmonize metadata
    adata.obs['perturbation'] = [x.replace('safe-targeting', 'control').replace('non-targeting', 'control') for x in adata.obs.guides]
    adata.obs['perturbation'] = [x.replace('--', '_') for x in adata.obs['perturbation']]
    adata.obs.perturbation = [k.replace('control_', '').replace('_control', '') for k in adata.obs.perturbation]  # collapse controls
    adata.obs['perturbation_type'] = 'CRISPR-cas9' if '_KO_' in key else 'CRISPRi'
    if 'replicate' in adata.obs.columns:
        cols = ['perturbation', 'replicate', 'ncounts', 'ngenes', 'nguides', '10X_channel', 'percent_mito', 'guides', 'full_guides', 'S_score', 'G2M_score', 'Cell_cycle_phase', 'perturbation_type']
    else:
        cols = ['perturbation', 'ncounts', 'ngenes', 'nguides', '10X_channel', 'percent_mito', 'guides', 'full_guides', 'S_score', 'G2M_score', 'Cell_cycle_phase', 'perturbation_type']
    adata.obs = adata.obs[cols]
    adata.obs['disease'] = "leukemia"
    adata.obs['cancer'] = True
    adata.obs['tissue_type']="cell_line"
    adata.obs["cell_line"] = "THP-1"
    adata.obs["celltype"] = 'monocytes'
    adata.obs['organism'] = 'human'
    adata.obs['nperts'] = [p.count('_')+1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
    annotate_qc(adata, species='human')
    adata.obs.index.name = 'cell_barcode'
    return adata

adatas = {}
for key in ['GSM6858447_KO_conventional', 'GSM6858449_KD_conventional', 'GSM6858448_KO_cell_pooled', 'GSM6858450_KD_guide_pooled']:
    print(key)
    adata = process_adata(key=key)
    assert_annotations(adata)
    adatas[key] = adata
adata = sc.concat(adatas, label='dataset')
adata.write(snakemake.output[0], compression='gzip')