import pandas as pd
import scanpy as sc
import sys

from tqdm import tqdm
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

# merge
adatas = []
for f in tqdm(snakemake.input):
    adata = sc.read(f)
    adatas.append(adata)
adata = sc.concat(adatas, axis=0)

# Obs
adata.obs.rename({
    'nCount_RNA': 'ncounts', 
    'nFeature_RNA': 'ngenes',
    'nCount_HTO': 'ncounts_tags',
    'biological_replicate_number': 'bio_replicate',
    'sample': 'perturbation',
    'rep': 'bio_replicate'
}, axis=1, inplace=True)
adata.obs.drop(['orig.ident'], axis=1, inplace=True)
adata.obs.perturbation = adata.obs.perturbation.astype(str)
adata.obs['perturbation'][adata.obs['perturbation']=='PBS'] = 'control'
adata.obs['nperts'] = [1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
adata.obs['perturbation_type'] = 'cytokines'
adata.obs['disease'] = "healthy"
adata.obs['cancer'] = False
adata.obs['tissue_type']="primary"
adata.obs['organism'] = 'mouse'
annotate_qc(adata, species='mouse')
assert_annotations(adata)

adata.write(snakemake.output[0], compression='gzip')
print('Done.')