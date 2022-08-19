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

# description under https://static-content.springer.com/esm/art%3A10.1038%2Fs41592-021-01153-z/MediaObjects/41592_2021_1153_MOESM4_ESM.xlsx

path = DIR + 'GSE168620/supp/Supp_GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells/'
lib = 'GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells'  # I guess the only relevant one
file = path+lib+'.h5ad'

# load data
if os.path.isfile(file+'.gz'):
    os.system(f"gunzip {file+'.gz'}")
adata = sc.read(file)
adata = adata[np.sum(adata.X, axis=1)>=100].copy()  # else too big

# add metadata
meta_file = path+lib+'.csv' if os.path.isfile(path+lib+'.csv') else path+lib+'.csv.gz'
metadata = pd.read_csv(meta_file)
obs = pd.merge(adata.obs, metadata, right_on='plate_well', left_on='plate_well')
obs.index = adata.obs.index
adata.obs = obs

# convert gene names
import mygene
mg = mygene.MyGeneInfo()
query = mg.querymany(list(adata.var_names) , scopes='ensembl.gene', fields='symbol', species='human', returnall=True, verbose=False)
df = pd.DataFrame(query['out'])
df = df.set_index('query')
df = df.drop_duplicates(keep='first')
var=pd.merge(adata.var, df, how='left', left_index=True, right_index=True)
var.loc[pd.isna(var['symbol']),'symbol'] = list(var.index[pd.isna(var['symbol'])])
var=var[['symbol']]
var = var.drop_duplicates(keep='first')
adata.var['ensembl_id'] = adata.var.index
adata.var_names=[var.loc[x]['symbol'] if x in var.index else x for x in adata.var_names]
adata.var_names = [str(x) for x in adata.var.index.values]  # make strings, else not writeable

# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
# reform obs
adata.obs = adata.obs.drop(['combinatorial_barcode', 'sample_name', 'material', 'gRNA_seq', 'gRNA_well'], axis=1)
adata.obs.index=adata.obs.index.rename('cell_barcode')
adata.obs=adata.obs.rename({'gRNA_ID': 'perturbation', 'TCR_status': 'perturbation_2', 'plate_well': 'sample'}, axis=1)
adata.obs['tissue_type']='cell_line'
adata.obs['cancer'] = True
adata.obs['celltype'] = 'T cells'
adata.obs['cell_line'] = 'Jurkat cells'
adata.obs['disease'] = 'acute T cell leukemia'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'CRISPR'
adata.obs['perturbation_type_2'] = 'TCR stimulation'

annotate_qc(adata)
assert_annotations(adata)

# save file
write_as_singles(adata, '/fast/work/users/peidlis_c/data/perturbation_resource_paper/', 'DatlingerBock2021', add_h5=True)
