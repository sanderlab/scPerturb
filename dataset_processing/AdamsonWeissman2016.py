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

path = DIR + 'GSE90546/supp/'
import shutil
shutil.unpack_archive(path+"GSE90546_RAW.tar", path)

for experiment in ['GSM2406675_10X001', 'GSM2406677_10X005', 'GSM2406681_10X010']:
    X = csr_matrix(mmread(path+experiment+"_matrix.mtx.txt.gz"))
    var = pd.read_csv(path+experiment+"_genes.tsv.gz", sep='\t', header=None, names=['gene_id', 'gene_name']).set_index('gene_name')
    obs = pd.read_csv(path+experiment+"_barcodes.tsv.gz", sep='\t', header=None, index_col=0)
    obs.index.names = ['cell_barcode']
    adata = sc.AnnData(X.T, var=var, obs=obs)
    adata.obs.index = [x.split('-')[0] for x in adata.obs.index]
    adata.var_names_make_unique()
    metadata = pd.read_csv(path+experiment+"_cell_identities.csv.gz", index_col=0)
    metadata.index = [x.split('-')[0] for x in metadata.index]
    adata.obs_names_make_unique()
    metadata=metadata[~metadata.index.duplicated()]
    adata.obs=pd.concat([adata.obs, metadata], axis=1, join='outer').loc[adata.obs_names]
    if 'good coverage' in adata.obs.columns:
        adata.obs['good coverage'] = np.array(adata.obs['good coverage'], dtype=str)
    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')
    adata.var = adata.var.rename({'gene_id': 'ensembl_id'}, axis=1)
    # reform obs
    adata.obs.index=adata.obs.index.rename('cell_barcode')
    adata.obs=adata.obs.rename({'guide identity': 'perturbation'}, axis=1)
    adata.obs = adata.obs.drop(['coverage', 'good coverage', 'number of cells'], axis=1)
    adata.obs['tissue_type']='cell_line'
    adata.obs['cell_line']='K562'
    adata.obs['cancer'] = True
    adata.obs['disease'] = 'chronic myelogenous leukemia'
    adata.obs['perturbation_type'] = 'CRISPR'
    adata.obs['celltype'] = 'lymphoblasts'
    adata.obs['organism'] = 'human'
    adata.write(f'{WDIR}AdamsonWeissman2016_{experiment}.h5')
