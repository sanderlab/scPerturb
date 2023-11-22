import pandas as pd
import scanpy as sc
import sys
from pathlib import Path

# Custom functions
sys.path.insert(1, '../')
from utils import *

# WARNING:
# The raw count file was obtained through personal communication with the authors.
# It is not publicly available.
WDIR = Path(snakemake.config['WDIR'])

adata = sc.read(WDIR / 'CRISPR_transcriptome_with_sgRNA_raw_all_data.h5ad')
adata.var = adata.var.drop('feature_types', axis=1)  # trivial

# move non-gene features to obsm
is_non_gene = ~adata.var.gene_ids.str.startswith('ENSG')
non_genes = list(adata.var_names[is_non_gene])
adata.obsm['reporters'] = pd.DataFrame(adata[:, non_genes].X.A, index=adata.obs_names, columns=non_genes)
adata = adata[:, adata.var.gene_ids.str.startswith('ENSG')].copy()

# harmonize metadata
adata.obs['perturbation'] = adata.obs.sgRNA_gene_identity
adata.obs = adata.obs.rename({'sgRNA_gene_identity': 'guide_id', 'n_count': 'ncounts', 'ChromHMM': 'ChromHMM_chromatin_state',
                              'CRISPR_num_features': 'nguides_detected', 'CRISPR_feature_call': 'guides_call', 'CRISPR_num_umis': 'guidewise_counts', 
                              'CRISPR_umis_sum': 'guide_ncounts'
                             }, axis=1)
adata.obs = adata.obs.drop(['Basal_count_bulk'], axis=1)  # irrelevant for general audience
adata.obs.perturbation = adata.obs.perturbation.replace({'Non-Targeting': 'control'})  # Non-targeting = control
adata = adata[adata.obs.perturbation!='nan'].copy()  # barcode undetermined
adata.obs = adata.obs[['perturbation', 'nguides_detected', 'guides_call', 'guidewise_counts', 'guide_ncounts', 'guide_id', 'ChromHMM_chromatin_state', 'ncounts']] # reorder
adata.obs = adata.obs[['perturbation', 'nguides_detected', 'guides_call', 'guidewise_counts', 'guide_ncounts', 'guide_id', 'ChromHMM_chromatin_state', 'ncounts']] # reorder
adata.obs['perturbation_type'] = 'CRISPRa'
adata.obs['disease'] = "healthy"
adata.obs['cancer'] = False
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "hPSCs"
adata.obs["celltype"] = 'stem cells'
adata.obs['organism'] = 'human'
adata.obs.perturbation = [x.replace('-', '_') for x in adata.obs.perturbation]  # convention for double perturbations
adata.obs['nperts'] = [p.count('_')+1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
annotate_qc(adata, species='human')
adata.obs.index.name = 'cell_barcode'

assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')