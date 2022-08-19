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

# I guess only the two screens are important. (screen 3 is mouse)
# WTX means whole transcriptome seq

path = DIR + 'GSE135497/supp/'
folders = get_subfolders(path, False)
for experiment in ['TAP_SCREEN__chromosome_11_screen', 'TAP_SCREEN__chromosome_8_screen']:
    subfolders = [x for x in folders if experiment in x]
    adatas = {}
    bdatas = {}
    for folder in subfolders:
        # load data
        files = get_files(path+folder)
        X = pd.read_csv([f for f in files if '.counts.' in f][0], index_col=0).T
        adata = sc.AnnData(csr_matrix(X.values, dtype=int))
        adata.var_names = X.columns
        adata.obs_names = X.index
        # add info what was perturbed
        pert_status = pd.read_csv([f for f in files if '.pertStatus.' in f][0], index_col=0).T
        bdata = sc.AnnData(csr_matrix(pert_status.values, dtype=int))
        bdata.var_names = pert_status.columns
        bdata.obs_names = pert_status.index
        # adata.obsm['pert_status'] = csr_matrix(pert_status)
        # adata.uns['pert_status_vars'] = pert_status.columns  # Problem: this is not saved...

        # add with sample name
        adatas[folder.split('__')[-1]]=adata
        bdatas[folder.split('__')[-1]]=bdata
    adata = sc.concat(adatas, index_unique='-', label='sample')
    bdata = sc.concat(bdatas, index_unique='-', label='sample')
    break
    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')
    # reform obs
    adata.obs.index=adata.obs.index.rename('cell_barcode')
    adata.obs=adata.obs.rename({'sample': 'replicate'}, axis=1)
    adata.obs['tissue_type']='cell_line'
    adata.obs['cell_line']='K562'
    adata.obs['cancer'] = True
    adata.obs['disease'] = 'chronic myelogenous leukemia'
    adata.obs['celltype'] = 'lymphoblasts'
    adata.obs['organism'] = 'human'
    guides_per_cell = np.array(np.sum(bdata.X, axis=1)).T[0]
    best_guide = list(bdata.var_names[np.array(np.argmax(bdata.X, axis=1)).T[0]])
    adata.obs['perturbation'] = best_guide
    adata.obs['perturbation'][guides_per_cell>1]='multiplet'
    adata.obs['perturbation'][guides_per_cell==0]='control'
    adata.obs['perturbation_type'] = 'CRISPR'

    adata.write(WDIR+'SchraivogelSteinmetz2020_'+experiment+'.h5')
