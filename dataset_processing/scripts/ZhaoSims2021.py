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

path = DIR + 'GSE148842/supp/'
folders = get_subfolders(path, True)
annotation = pd.read_excel('../tables/GSE148842_manual_annotation.xlsx')

adatas = {}
for folder in folders:
    print(folder)
    file = get_files(folder, False)[0]
    tab = pd.read_csv(folder + '/' + file, sep='\t')
    var = tab[tab.columns[:2]].set_index(tab.columns[1], drop=True)
    obs_names = tab.columns[2:]
    X = tab.drop(tab.columns[:2], axis=1).values.T
    from scipy.sparse import csr_matrix
    X = csr_matrix(X)

    adata = sc.AnnData(X, var=var)
    adata.obs_names = obs_names
    adata.var_names_make_unique('_')

    ann = annotation[annotation['GEO'] == file.split('_')[0]]
    patient = ann['Sample'].values[0][:5]
    adata.obs['patient'] = patient
    for c in ann.columns:
        adata.obs[c] = ann[c].values[0]
    adatas[file.split('.')[0].split('_')[-1]] = adata

adata = sc.concat(adatas, index_unique='-', label='library')
# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
# reform obs
adata.obs['dose_value']=[None if ('DMSO' in x or x=='none') else x.split(' ')[0] for x in adata.obs.treatment]
adata.obs['dose_unit']=[None if ('DMSO' in x or x=='none') else x.split(' ')[1] for x in adata.obs.treatment]
adata.obs['perturbation']=['control' if ('DMSO' in x or x=='none') else x.split(' ')[-1] for x in adata.obs.treatment]
adata.obs = adata.obs.drop('treatment', axis=1)
adata.obs.index=adata.obs.index.rename('cell_barcode')
adata.obs=adata.obs.rename({'patient': 'sample', 'gender': 'sex'}, axis=1)
adata.obs.perturbation[adata.obs.perturbation=='vehicle (DMSO)']='control'
adata.obs['tissue_type']='primary'
adata.obs['cancer'] = True
adata.obs['disease'] = 'glioblastoma'
adata.obs['celltype'] = 'unknown'
adata.obs['sex'] = adata.obs['sex'].replace('M', 'm').replace('F', 'f')
adata.obs['organism'] = 'human'

adata.write(WDIR+'ZhaoSims2021.h5')
