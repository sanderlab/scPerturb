import pandas as pd
import scanpy as sc
import numpy as np
import sys
import os

from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path

TEMPDIR = Path(snakemake.config['TEMPDIR']) 
name = snakemake.wildcards['name']

time = 7 if name.endswith('7') else 3 if name.endswith('3') else None
pert_key = 'grna' if time!=None else 'grna_strict'
moi = 'normal' if time!=None else 'high'

X = csr_matrix(mmread(TEMPDIR / 'DixitRegev2016' / f'{name}.mtx.txt'))
adata = sc.AnnData(X.T)

var = pd.read_csv(TEMPDIR / 'DixitRegev2016' / f'{name}_genenames.csv', index_col=0)
splitted = np.array([x.split('_', 1) for x in var.values[:,0]])
adata.var_names = splitted[:,1]
adata.var['gene_id'] = splitted[:,0]
adata.var_names_make_unique()

obs = pd.read_csv(TEMPDIR / 'DixitRegev2016' / f'{name}_cellnames.csv', index_col=0)
splitted = np.array([x.split('_', 2) for x in obs.values[:,0]])
adata.obs_names = obs.values[:,0]
adata.obs['identifier_0'] = splitted[:,1]
adata.obs['identifier_1'] = splitted[:,2]

# annotation:
files = [TEMPDIR / 'DixitRegev2016' / f'{name}_cbc_gbc_dict.csv', TEMPDIR / 'DixitRegev2016' / f'{name}_cbc_gbc_dict_strict.csv', TEMPDIR / 'DixitRegev2016' / f'{name}_cbc_gbc_dict_lenient.csv']
keys = ['grna', 'grna_strict', 'grna_lenient']
for file, key in zip(files, keys):
    if os.path.isfile(file):
        cbc_gbc_dict = pd.read_csv(file, index_col=0, header=None)
        adata.obs[key]=None
        for grna in list(cbc_gbc_dict.index):
            for barcode in cbc_gbc_dict.loc[grna][1].replace(' ','').split(','):
                if barcode in adata.obs_names:
                    val = adata.obs.loc[barcode][key]
                    adata.obs.loc[barcode][key] = grna if val is None else val+' + '+grna
        adata.obs[key] = adata.obs[key].astype(str) # NaNs are not allowed in the obs

adata.obs['target'] = [x.replace('p_sg', '').replace('p_', '').split('_')[0] if type(x)==str else None for x in adata.obs[pert_key]]
adata.obs = adata.obs.rename({pert_key: 'perturbation'}, axis=1).drop(['identifier_1', 'identifier_0'], axis=1)
adata.obs['moi'] = moi
adata.obs['time'] = time*24  if time!=None else "None"
adata.obs['cell_line'] = 'K562'
adata.obs['celltype'] = 'lymphoblasts'
adata.obs['perturbation_type'] = 'CRISPR'
adata.obs['cancer'] = True
adata.obs['disease'] = 'myelogenous leukemia'
adata.obs['library'] = name

# we have NaNs
adata.obs.perturbation = adata.obs.perturbation.astype(str)  
adata.obs.target = adata.obs.target.astype(str)

adata.write(snakemake.output[0])
