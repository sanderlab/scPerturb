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

path = f'{DIR}GSE140802/supp/'
folders = get_subfolders(path)

folders = get_subfolders(path)
# only "GSM4185644_stateFate_cytokinePerturbation" has perturbations
for folder in folders:
    name = folder.split('/Supp_')[-1]
    if name == 'GSM4185644_stateFate_cytokinePerturbation':
        print('Processing ', folder)
        files = get_files(folder)

        # count matrix
        mat = csr_matrix(mmread([f for f in files if '_normed_counts.mtx' in f][0]))
        # gene names
        var = pd.read_csv([f for f in files if '_gene_names' in f][0], index_col=0, header=None, sep='\t', names=['gene_symbols'])
        # barcodes and metadata
        obs = pd.read_csv([f for f in files if '_cell_barcodes' in f][0], index_col=0, header=None, sep='\t', names=['cell_barcodes'])
        obs1 = pd.read_csv([f for f in files if '_library_names' in f][0], index_col=0, header=None, sep='\t', names=['library_names'])
        obs['library_names'] = list(obs1.index)
        metafiles = [f for f in files if '_metadata' in f]
        if len(metafiles) > 0:
            metafile = metafiles[0]
            obs2 = pd.read_csv(metafile, index_col='Cell barcode', sep='\t')
            obs = pd.concat([obs, obs2], axis=1)
        if ('library_names' in obs.columns) and ('Library' in obs.columns):
            obs = obs.drop('Library', axis=1)

        adata = sc.AnnData(mat, obs, var)

        # add clonal info
        clones = [f for f in files if '_clone_matrix.mtx' in f]
        if len(clones) > 0:
            cmat = csr_matrix(mmread(clones[0])).T
            adata.obsm['clone_matrix'] = cmat

        # reform var
        adata.var.index = adata.var.index.rename('gene_symbol')
        # reform obs
        adata.obs['cancer'] = False
        adata.obs.index=adata.obs.index.rename('cell_barcode')
        adata.obs=adata.obs.rename({'Time point': 'age', 'Cytokine condition': 'perturbation', 'Cell type annotation': 'celltype'}, axis=1)
        adata.obs.age = ['day '+str(x) for x in adata.obs.age]
        adata.obs = adata.obs.drop(['SPRING-x', 'SPRING-y'], axis=1)
        adata.obs['organism'] = 'mouse'
        adata.obs['disease'] = 'healthy'
        adata.obs['tissue_type']='primary'
        adata.obs['perturbation_type'] = 'cytokine'
        adata.obsm['clone_matrix'] = adata.obsm['clone_matrix'].A  # we can only save dense ones in obsm atm
        adata.obs_names = [x+'_'+y for x, y in zip(adata.obs_names, adata.obs.library_names)]
        write_as_singles(adata, '/fast/work/users/peidlis_c/data/perturbation_resource_paper/', 'WeinrebKlein2020', add_h5=True)
