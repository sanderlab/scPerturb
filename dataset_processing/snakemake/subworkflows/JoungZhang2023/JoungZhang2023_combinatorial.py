import pandas as pd
import scanpy as sc
import numpy as np
import sys
import re
import h5py

from scipy.sparse import csr_matrix, vstack

# Custom functions
sys.path.insert(1, '../')
from utils import *

### Combinatorial screen ###
print('Combinatorial screen')
tdata = sc.read_h5ad(snakemake.input['combinatorial'])
# their adata.X was scaled and their adata.raw.X was saved as dense...
# reconstruct unfiltered count matrix and make true sparse
adata = sc.AnnData(csr_matrix(tdata.raw.X), tdata.obs, tdata.raw.var, 
                   tdata.uns, tdata.obsm, obsp=tdata.obsp)
adata.obs['perturbation'] = ['+'.join(sorted(re.sub(r'TFORF[0-9]{4}-', ' ', s).replace(' ','').split(','))) for s in adata.obs.TF]  # cleaner names

# set metadata
adata.obs = adata.obs.rename({'n_genes': 'ngenes', 'n_counts': 'ncounts'}, axis=1)
adata.obs['disease'] = "healthy"
adata.obs['cancer'] = False
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "hESCs"
adata.obs["celltype"] = 'hESCs'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'ORF overexpression'
adata.obs['nperts'] = [p.count('+') + 1 - p.count('GFP') if type(p)==str else 0 for p in adata.obs.perturbation]
annotate_qc(adata, species='human')
adata.obs.index.name = 'cell_barcode'

assert_annotations(adata)
adata.write(snakemake.output['combinatorial'], compression='gzip')
