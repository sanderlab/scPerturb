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

import yaml
config = yaml.safe_load(open("config.yaml", "r"))
DIR = config['DIR']
WDIR = config['WDIR']

path = DIR+'GSE90063/'
adatas = {}


# GSM2396858_K562_TFs__7_days
path = DIR+'GSE90063/'
folder = 'GSM2396858_k562_tfs_7'


files = get_files(path+folder, full_path=False)
name = os.path.commonprefix(files)
spath = path+folder+'/'+name

X = csr_matrix(mmread(spath+'.mtx.txt.gz'))
adata = sc.AnnData(X.T)

var = pd.read_csv(spath+'_genenames.csv.gz', index_col=0)
splitted = np.array([x.split('_', 1) for x in var.values[:,0]])
adata.var_names = splitted[:,1]
adata.var['gene_id'] = splitted[:,0]
adata.var_names_make_unique()

obs = pd.read_csv(spath+'_cellnames.csv.gz', index_col=0)
splitted = np.array([x.split('_', 2) for x in obs.values[:,0]])
adata.obs_names = obs.values[:,0]
adata.obs['identifier_0'] = splitted[:,1]
adata.obs['identifier_1'] = splitted[:,2]

# annotation:
files = [spath+'_cbc_gbc_dict.csv.gz', spath+'_cbc_gbc_dict_strict.csv.gz', spath+'_cbc_gbc_dict_lenient.csv.gz']
keys = ['grna', 'grna_strict', 'grna_lenient']
for file, key in zip(files, keys):
    if os.path.isfile(file):
        cbc_gbc_dict = pd.read_csv(file, index_col=0, header=None)
        adata.obs[key]=None
        for grna in list(cbc_gbc_dict.index):
            for barcode in cbc_gbc_dict.loc[grna][1].replace(' ','').split(','):
                if barcode in adata.obs_names:
                    val = adata.obs.loc[barcode][key]
                    adata.obs.loc[barcode,key] = grna if val is None else val+' + '+grna

adata.obs['target'] = [x.replace('p_sg', '').replace('p_', '').split('_')[0] if type(x)==str else None for x in adata.obs.grna]
adata.obs = adata.obs.rename({'grna': 'perturbation'}, axis=1).drop(['identifier_1', 'identifier_0'], axis=1)
adata.obs['moi'] = 'normal'
adata.obs['time'] = 7*24  # 7 days
adata.obs['cell_line'] = 'K562'
adata.obs['celltype'] = 'lymphoblasts'
adata.obs['perturbation_type'] = 'CRISPR'
# adata.obs.perturbation[pd.isna(adata.obs.perturbation)] = 'control'
adata.obs['cancer'] = True
adata.obs['disease'] = 'myelogenous leukemia'
adatas['K562_TFs__7_days'] = adata
adata.write_h5ad("DixitRegev2016_K562_TFs_7days.h5ad",compression="gzip")

# GSM2396859_K562_TFs__13_days
path = DIR+'GSE90063/'
folder = 'GSM2396859_k562_tfs_13'

files = get_files(path+folder, full_path=False)
name = os.path.commonprefix(files)
spath = path+folder+'/'+name

X = csr_matrix(mmread(spath+'.mtx.txt.gz'))
adata = sc.AnnData(X.T)

var = pd.read_csv(spath+'_genenames.csv.gz', index_col=0)
splitted = np.array([x.split('_', 1) for x in var.values[:,0]])
adata.var_names = splitted[:,1]
adata.var['gene_id'] = splitted[:,0]
adata.var_names_make_unique()

obs = pd.read_csv(spath+'_cellnames.csv.gz', index_col=0)
splitted = np.array([x.split('_', 2) for x in obs.values[:,0]])
adata.obs_names = obs.values[:,0]
adata.obs['identifier_0'] = splitted[:,1]
adata.obs['identifier_1'] = splitted[:,2]

# annotation:
files = [spath+'_cbc_gbc_dict.csv.gz', spath+'_cbc_gbc_dict_strict.csv.gz', spath+'_cbc_gbc_dict_lenient.csv.gz']
keys = ['grna', 'grna_strict', 'grna_lenient']
for file, key in zip(files, keys):
    if os.path.isfile(file):
        cbc_gbc_dict = pd.read_csv(file, index_col=0, header=None)
        adata.obs[key]=None
        for grna in list(cbc_gbc_dict.index):
            for barcode in cbc_gbc_dict.loc[grna][1].replace(' ','').split(','):
                if barcode in adata.obs_names:
                    val = adata.obs.loc[barcode][key]
                    adata.obs.loc[barcode,key] = grna if val is None else val+' + '+grna
adata.obs['target'] = [x.replace('p_sg', '').replace('p_', '').split('_')[0] if type(x)==str else None for x in adata.obs.grna]
adata.obs = adata.obs.rename({'grna': 'perturbation'}, axis=1).drop(['identifier_1', 'identifier_0'], axis=1)
adata.obs['moi'] = 'normal'
adata.obs['time'] = 13*24  # 13 days
adata.obs['cell_line'] = 'K562'
adata.obs['celltype'] = 'lymphoblasts'
adata.obs['perturbation_type'] = 'CRISPR'
# adata.obs.perturbation[pd.isna(adata.obs.perturbation)] = 'control'
adata.obs['cancer'] = True
adata.obs['disease'] = 'myelogenous leukemia'
adatas['K562_TFs__13_days'] = adata
adata.write_h5ad("DixitRegev2016_K562_TFs_13days.h5ad",compression="gzip")


# GSM2396860_K562_TFs__High_MOI
path = DIR+'GSE90063/'
folder = 'GSM2396860_k562_tfs_highmoi'

files = get_files(path+folder, full_path=False)
name = os.path.commonprefix(files)
spath = path+folder+'/'+name

X = csr_matrix(mmread(spath+'.mtx.txt.gz'))
adata = sc.AnnData(X.T)

var = pd.read_csv(spath+'_genenames.csv.gz', index_col=0)
splitted = np.array([x.split('_', 1) for x in var.values[:,0]])
adata.var_names = splitted[:,1]
adata.var['gene_id'] = splitted[:,0]
adata.var_names_make_unique()

obs = pd.read_csv(spath+'_cellnames.csv.gz', index_col=0)
splitted = np.array([x.split('_', 2) for x in obs.values[:,0]])
adata.obs_names = obs.values[:,0]
adata.obs['identifier_0'] = splitted[:,1]
adata.obs['identifier_1'] = splitted[:,2]

# annotation:
files = [spath+'_cbc_gbc_dict.csv.gz', spath+'_cbc_gbc_dict_strict.csv.gz', spath+'_cbc_gbc_dict_lenient.csv.gz']
keys = ['grna', 'grna_strict', 'grna_lenient']
for file, key in zip(files, keys):
    if os.path.isfile(file):
        cbc_gbc_dict = pd.read_csv(file, index_col=0, header=None)
        adata.obs[key]=None
        for grna in list(cbc_gbc_dict.index):
            for barcode in cbc_gbc_dict.loc[grna][1].replace(' ','').split(','):
                if barcode in adata.obs_names:
                    val = adata.obs.loc[barcode][key]
                    adata.obs.loc[barcode,key] = grna if val is None else val+' + '+grna
adata.obs['target'] = [x.replace('p_sg', '').replace('p_', '').split('_')[0] if type(x)==str else None for x in adata.obs.grna_strict]
adata.obs = adata.obs.rename({'grna_strict': 'perturbation'}, axis=1).drop(['identifier_1', 'identifier_0'], axis=1)
adata.obs['moi'] = 'high'
adata.obs['cell_line'] = 'K562'
adata.obs['celltype'] = 'lymphoblasts'
adata.obs['perturbation_type'] = 'CRISPR'
# shift to treating non-targeting guides as control
# adata.obs.perturbation[pd.isna(adata.obs.perturbation)] = 'control'
adata.obs['cancer'] = True
adata.obs['disease'] = 'myelogenous leukemia'
adatas['K562_TFs__High_MOI'] = adata
adata.write_h5ad("DixitRegev2016_K562_TFs_highMOI.h5ad",compression="gzip")

for key in adatas:
    print(key)
    adata = adatas[key]
    adata.obs['tissue_type']='cell_line'
    adata.obs['organism'] = 'human'
    obs = adata.obs
    obs.to_csv('test_obs'+key+'.csv')
    if 'grna_lenient' in adata.obs.columns:
        adata.obs['grna_lenient']=adata.obs['grna_lenient'].str.replace(' + ',';', regex=False)
    adata.obs['guide_id']= adata.obs['perturbation']
    adata.obs['perturbation'] = adata.obs['perturbation'].astype(str)
    adata.obs.loc[adata.obs['perturbation'].str.contains('INTERGENIC') &~ adata.obs['perturbation'].str.contains('+', regex=False).astype('bool'),'perturbation']  = 'control'

    if 'grna_lenient' in adata.obs.columns:
        adata.obs['grna_lenient']=adata.obs['grna_lenient'].str.replace('_','-', regex=False)

    adata.obs['perturbation']=adata.obs['perturbation'].str.replace('_','-', regex=False)
    adata.obs['perturbation']=adata.obs['perturbation'].str.replace(' + ','_', regex=False)
    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))
    qc = sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=False)           

    adata.obs['ncounts'] = qc[0]['total_counts']
    adata.obs['ngenes'] = qc[0]['n_genes_by_counts']
    adata.obs['percent_mito'] = qc[0]['pct_counts_mt']
    adata.obs['percent_ribo'] = qc[0]['pct_counts_ribo']
    adata.var['ncounts'] = qc[1]['total_counts']
    adata.var['ncells'] = qc[1]['n_cells_by_counts']
    adata.obs['nperts'] = adata.obs['perturbation'].str.count('_') + 1
    adata.obs['tissue_type']='cell_line'
    adata.obs['organism'] = 'human'
    adata.write_h5ad('DixitRegev2016_'+key +'.h5ad', compression='gzip')





