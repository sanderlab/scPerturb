import pandas as pd
import scanpy as sc
import sys
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

# load main dataset
adata = sc.read(snakemake.input[0])

# obs
adata.obs.set_index('Barcode', drop=True, inplace=True)
adata.obs.index.name = 'cell_barcode'
adata.obs['organism'] = 'Mus musculus'
adata.obs['disease'] = 'healthy'
adata.obs['cancer'] = False
adata.obs['perturbation_type'] = 'CRISPR-cas9'
adata.obs['tissue_type'] = 'primary'
adata.obs['tissue'] = 'brain'
adata.obs.rename({
    'per_gene': 'perturbation',
    'nCount_RNA':'ncounts',
    'nFeature_RNA': 'ngenes',
    'gRNAs': 'guide_id'
}, axis=1, inplace=True)
adata.obs.perturbation = adata.obs.perturbation.astype(str)
adata.obs.perturbation[adata.obs.perturbation=='Safe_H'] = 'control'
adata.obs['nperts'] = [1-p.count('control') if type(p)==str else None for p in adata.obs.perturbation]

# var
adata.var.index.name = 'gene_symbol'
annotate_qc(adata, species='mouse')
assert_annotations(adata)

adata.write(snakemake.output[0])