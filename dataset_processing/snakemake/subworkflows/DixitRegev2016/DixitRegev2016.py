import pandas as pd
import scanpy as sc
import numpy as np
import sys
import os

from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

adatas = [sc.read_h5ad(file) for file in snakemake.input]
adata = sc.concat(adatas, index_unique='-')
adata.obs['tissue_type']='cell_line'
adata.obs['organism'] = 'human'

obs = adata.obs.copy()
# obs['grna_lenient']=[x.replace(' + ',';') for x in obs['grna_lenient']]
# obs['guide_id']= obs['perturbation'].copy()

# Intergenic is proper control
obs.perturbation = obs.perturbation.astype(str)
obs.perturbation = [x.replace('_', '-').replace('p-sg', '').replace('p-', '').replace(' + ', '_') for x in obs.perturbation]  # formatting
obs.perturbation = ['control' if ('INTERGENIC' in x) and ('_' not in x) else x for x in obs.perturbation]  # annotate control
obs.perturbation = ['_'.join(np.unique([y.split('-')[0] if 'INTERGENIC' not in y else 'INTERGENIC' for y in x.split('_')])) for x in obs.perturbation]  # collapse guides
obs.perturbation = [x.replace('_INTERGENIC', '').replace('INTERGENIC_', '') for x in obs.perturbation]  # remove intergenic from combis
obs.perturbation[obs.perturbation == 'INTERGENIC'] = 'control'

# obs['grna_lenient']=[x.replace('_','-') for x in obs['grna_lenient']]
adata.obs = obs.copy()

adata.obs['nperts'] = [1+p.count('_')-p.count('control')-p.count('None') if type(p)==str else 0 for p in adata.obs.perturbation]

annotate_qc(adata, species='human')
assert_annotations(adata)

adata.obs.target = adata.obs.target.astype(str)
adata.write(TEMPDIR / 'DixitRegev2016.h5ad', compression='gzip')
