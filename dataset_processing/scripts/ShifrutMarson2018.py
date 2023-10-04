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

import tarfile
adatas = {}
folders = get_subfolders(DIR+'GSE119450/supp/')
for d in ['D1', 'D2']:
    for cond in ['Stim', 'NoStim']:
        sample_name = d+'_'+cond
        gex_path=[x for x in folders if sample_name+'_10x' in x][0]
        obs_path=[x for x in folders if sample_name+'_ReAmp' in x][0]

        # extract GEX and read
        file = [f for f in get_files(gex_path) if '.gz' in f][0]
        tar = tarfile.open(file, "r:gz")
        tar.extractall(path=gex_path)
        tar.close()
        adata = sc.read_10x_mtx(gex_path)
        adata.obs.index = [x.split('-')[0] for x in adata.obs.index]

        # read and add metadata
        tab=pd.read_csv(get_files(obs_path)[0], index_col=0)
        tab=tab.rename({'UMI.count': 'gRNA_UMI_count', 'gRNA.ID': 'gRNA_ID'}, axis=1)[['gRNA_ID', 'gRNA_UMI_count']]
        adata.obs = pd.merge(adata.obs, tab, how='left', left_index=True, right_index=True)
        adata.obs['donor'] = d
        adata.obs['condition'] = cond
        adata.obs['sample_name'] = sample_name

        adata.write(DIR+'GSE119450/supp/GSE119450' + sample_name + '_processed_supp.h5')
        adatas[sample_name] = adata
adata = sc.concat(adatas, index_unique='-')

# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
# reform obs
adata.obs['cancer'] = False
adata.obs.index=adata.obs.index.rename('cell_barcode')
adata.obs=adata.obs.rename({'gRNA_ID': 'perturbation', 'donor': 'replicate', 'condition': 'perturbation_2', 'sample_name': 'library'}, axis=1)
adata.obs['celltype'] = 'T cells'
adata.obs['organism'] = 'human'
adata.obs['disease'] = 'healthy'
adata.obs['tissue_type']='primary'
adata.obs['perturbation_type'] = 'CRISPR'
adata.obs['perturbation_type_2'] = 'TCR stimulation'
adata.obs.perturbation = adata.obs.perturbation.astype(str)
adata.obs.perturbation[adata.obs.perturbation=='nan'] = 'control'

adata.write(WDIR+'/ShifrutMarson2018.h5')
