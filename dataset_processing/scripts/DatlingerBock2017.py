import scanpy as sc
import pandas as pd
import os
import sys
import numpy as np
import gzip
import matplotlib.pyplot as pl
import re

from scipy.io import mmread
from scipy.sparse import csr_matrix
from tqdm import tqdm

from process_supp import *
sys.path.append(os.path.abspath(os.path.join(os.getcwd(), '../')))
from utils import write_as_singles, read_from_singles, annotate_qc, assert_annotations

import yaml
config = yaml.safe_load(open("../../config.yaml", "r"))
DIR = config['DIR']
WDIR = config['WDIR']

# read
tab = pd.read_csv(f"{DIR}GSE92872/supp/GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz",
                 skiprows=6, header=None, index_col=0).T
X = csr_matrix(tab.values, dtype=int)
obs = pd.read_csv(f"{DIR}GSE92872/supp/GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz",
                 header=None, index_col=0, nrows=5).T.set_index('cell')

# create
adata = sc.AnnData(X, obs=obs)
adata.var_names = list(tab.columns)
adata.obs_names_make_unique()

# rename and reorder
adata.obs = adata.obs.rename({'condition': 'perturbation_2', 'gene': 'target', 'grna': 'perturbation'}, axis=1)
adata.obs = adata.obs[np.sort(adata.obs.keys())]

# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
# reform obs
adata.obs.index=adata.obs.index.rename('cell_barcode')
adata.obs['perturbation'] = ['control' if 'CTRL' in x else x for x in adata.obs['perturbation']]
adata.obs['target'][adata.obs.target=='CTRL']=None
adata.obs['celltype'] = 'T cells'
adata.obs['cell_line'] = 'Jurkat cells'
adata.obs['cancer'] = True
adata.obs['disease'] = 'acute T cell leukemia'
adata.obs['tissue_type'] = 'cell_line'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'CRISPR'
adata.obs['perturbation_type_2'] = 'TCR stimulation'

# write
adata.write(f"{WDIR}DatlingerBock2017.h5")

# inspect result
print(adata)
