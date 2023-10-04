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

# GSE106340
path = DIR + 'GSE122662/supp/'
folders = get_subfolders(path, False)
folders = [x for x in folders if int(x[8:15]) in np.arange(2836267, 2836288+1)]
adatas = {}
for folder in folders:
    # create adata
    files = get_files(path+folder)
    var = pd.read_csv([x for x in files if 'genes' in x][0], sep='\t', index_col=1, header=None, names=['gene_id', 'gene_name'])
    obs = pd.read_csv([x for x in files if 'barcodes' in x][0], sep='\t', index_col=0, header=None)
    try:
        X = csr_matrix(mmread([x for x in files if 'matrix.mtx' in x][0])).T
    except:
        print(f'No count matrix found for {folder}. Is this corrupted?')
        continue
    adata = sc.AnnData(X, obs, var)

    # annotate
    GSM = folder.split('_', 2)[1]
    adata.obs['GSM'] = GSM
    string = folder.split('_', 2)[2]
    import re
    adata.obs['condition'] = 'serum' if 'serum' in string else '2i' if '2i' in string else None
    adata.obs['time'] = 'iPSCs' if 'iPSCs' in string else re.search('D\d+', string).group(0)

    # clear barcode names
    adata.obs_names = [x.split('-')[0] for x in adata.obs_names]
    adata.var_names_make_unique()
    adatas[GSM] = adata
adata = sc.concat(adatas, index_unique='-')

# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
# reform obs
adata.obs.index=adata.obs.index.rename('cell_barcode')
adata.obs['cancer'] = False
adata.obs['organism'] = 'mouse'
adata.obs['disease'] = 'healthy'
adata.obs['tissue_type'] = 'stem'
adata.obs['perturbation_type'] = 'drug'
adata.obs['celltype']='mESCs'
adata.obs=adata.obs.rename({'time': 'age', 'condition': 'perturbation'}, axis=1)
adata.obs['perturbation'] = adata.obs['perturbation'].astype(str)
adata.obs['perturbation'][adata.obs.perturbation=='2i'] = '2i (LIF+MEKi+GSK3i)'

adata.write(WDIR+'SchiebingerLander2019_GSE106340.h5')


# GSE115943
path = DIR + 'GSE122662/supp/'
folders = get_subfolders(path, False)
folders = [x for x in folders if int(x[8:15]) in np.arange(3195648, 3195773+1)]
adatas = {}
for folder in folders:
    file = get_files(path+folder, True)[0]
    adata = sc.read_10x_h5(file)
    adata.obs_names = [x.split('-')[0] for x in adata.obs_names]
    adata.obs['time'] = 'iPSC' if 'iPSC' in folder else re.search('D[\d\.]+(_5)?', folder).group(0).replace('_', '.')
    adata.obs['replicate'] = re.search('C\d', folder).group(0)[-1]
    adata.obs['condition'] = folder.split('_')[-2]
    GSM = folder.split('_', 2)[1]
    adata.obs['GSM'] = GSM
    adata.var_names_make_unique()
    adatas[GSM] = adata
adata = sc.concat(adatas, index_unique='-')
# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
# reform obs
adata.obs.index=adata.obs.index.rename('cell_barcode')
adata.obs['cancer'] = False
adata.obs['organism'] = 'mouse'
adata.obs['disease'] = 'healthy'
adata.obs['tissue_type']='primary'
adata.obs['perturbation_type'] = 'drug'
adata.obs['celltype']='mESCs'
adata.obs=adata.obs.rename({'time': 'age', 'condition': 'perturbation'}, axis=1)
adata.obs['perturbation'] = adata.obs['perturbation'].astype(str)
adata.obs['perturbation'][adata.obs.perturbation=='2i'] = '2i (LIF+MEKi+GSK3i)'
adata.write(WDIR+'SchiebingerLander2019_GSE115943.h5')
