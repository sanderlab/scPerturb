import pandas as pd
import scanpy as sc
import numpy as np
import sys
import os

from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

### OLD ###
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

path = DIR + 'GSE139944/supp/'
folders = get_subfolders(path, False)
def prepare(folder):
    files = get_files(path+folder, False)
    gene_annotations = [x for x in files if 'gene.annotations' in x][0]
    var = pd.read_csv(path+folder+'/'+gene_annotations, sep='\t', header=None, names=['gene_id', 'gene_name']).set_index('gene_name')
    # obs = pd.read_csv(path+folder+'/'+folder[5:]+'_cell.annotations.txt.gz', sep='\t', header=None)  # no additional info here?
    metadata = [x for x in files if 'pData' in x][0]
    obs2 = pd.read_csv(path+folder+'/'+metadata, sep=' ')
    counts = [x for x in files if 'UMI.count.matrix' in x][0]
    UMI_counts = pd.read_csv(path+folder+'/'+counts, sep='\t', header=None)
    X=csr_matrix((UMI_counts[2], (UMI_counts[1]-1, UMI_counts[0]-1)), shape=(len(obs2), len(var)))
    adata = sc.AnnData(X, obs=obs2, var=var)
    adata.var_names_make_unique(join=':')

    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')
    adata.var['gene_id'] = [x.split('.')[0] for x in adata.var['gene_id']]
    adata.var = adata.var.rename({'gene_id': 'ensembl_id'}, axis=1)
    return adata

# sciplex1 is just a species mixing experiment and hence irrelevant

# sciplex2
folder= 'Supp_GSM4150377_sciPlex2_A549_Transcription_Modulators'
dataset = folder.replace('Supp_', '')
adata = prepare(folder)

adata.obs = adata.obs.drop(['Cell', 'sample', 'Size_Factor'], axis=1)
df=pd.DataFrame({index: row.top_oligo.split('_') for index, row in adata.obs.iterrows() if type(row.top_oligo)==str}).T
df.columns = ['perturbation', 'dose_value', 'well']
adata.obs = pd.concat([adata.obs, df], axis=1)

# NOTE: In GEO it's mentioned there should be DMSO treated control cells. This is not explicitly named, so I guess they are the "None" ones?
adata.obs['perturbation'][pd.isna(adata.obs['perturbation'])]='control'
adata.obs['celltype'] = 'alveolar basal epithelial cells'
adata.obs['cell_line'] = 'A549'
adata.obs['cancer'] = True
adata.obs['disease'] = 'lung adenocarcinoma'
adata.obs['tissue_type'] = 'cell_line'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'drug'

adata.obs = adata.obs.rename({'n.umi': 'ncounts', 'pval': 'pval_demultiplexing', 'qval': 'qval_demultiplexing'}, axis=1)

adata.var_names_make_unique()
adata.write(WDIR+'SrivatsanTrapnell2020_sciplex2.h5')

# sciplex3
folder= 'Supp_GSM4150378_sciPlex3_A549_MCF7_K562_screen'
dataset = folder.replace('Supp_', '')
# adata = prepare(folder)
files = get_files(path+folder, False)
gene_annotations = [x for x in files if 'gene.annotations' in x][0]
var = pd.read_csv(path+folder+'/'+gene_annotations, sep='\t', header=None, names=['gene_id', 'gene_name']).set_index('gene_name')
# obs = pd.read_csv(path+folder+'/'+folder[5:]+'_cell.annotations.txt.gz', sep='\t', header=None)  # no additional info here?
metadata = [x for x in files if 'pData' in x][0]
obs2 = pd.read_csv(path+folder+'/'+metadata, sep=' ')
counts = [x for x in files if 'UMI.count.matrix' in x][0]

UMI_counts = pd.read_csv(path+folder+'/'+counts, sep='\t', header=None)
X=csr_matrix((UMI_counts[2], (UMI_counts[1]-1, UMI_counts[0]-1)), shape=(len(obs2), len(var)))  # this may crash your kernel
adata = sc.AnnData(X, obs=obs2, var=var)
adata.var_names_make_unique(join=':')

# reform var
adata.var.index = adata.var.index.rename('gene_symbol')
adata.var['gene_id'] = [x.split('.')[0] for x in adata.var['gene_id']]
adata.var = adata.var.rename({'gene_id': 'ensembl_id'}, axis=1)
adata.obs=adata.obs.drop(['cell', 'sample', 'Size_Factor', 'top_to_second_best_ratio_W', 'hash_umis_W', 'pval_W', 'qval_W', 'top_to_second_best_ratio_P', 'hash_umis_P', 'pval_P', 'qval_P'], axis=1)
adata.obs = adata.obs.drop(['top_oligo_W', 'top_oligo_P', 'rt_well', 'lig_well', 'pcr_well', 'pcr_plate', 'culture_plate',
       'rt_plate', 'lig_plate', 'Combo', 'drug_dose', 'catalog_number', 'treatment', 'dose_character', 'dose_pattern', 'vehicle'], axis=1)
adata.obs = adata.obs.rename({'n.umi': 'ncounts', 'well_oligo': 'well', 'plate_oligo': 'plate',
                  'cell_type': 'cell_line', 'time_point': 'time', 'dose': 'dose_value',
                  'product_name' : 'perturbation'}, axis=1)
adata.obs['dose_unit']='nM'  # I guess this is in nanomolar since doses are 10 nM, 100 nM, 1 μM, and 10 μM
adata.obs.perturbation[adata.obs.perturbation == 'Vehicle']='control'
adata.obs['celltype'] = ['alveolar basal epithelial cells' if line=='A549'
                       else 'mammary epithelial cells' if line=='MCF7'
                       else 'lymphoblasts' if line=='K562'
                       else 'None' for line in adata.obs.cell_line]
adata.obs['disease'] = ['lung adenocarcinoma' if line=='A549'
                      else 'breast adenocarcinoma' if line=='MCF7'
                      else 'chronic myelogenous leukemia' if line=='K562'
                      else 'None' for line in adata.obs.cell_line]
adata.obs['cancer'] = True
adata.obs['tissue_type'] = 'cell_line'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'drug'
adata.write(WDIR+'SrivatsanTrapnell2020_sciplex3.h5')

# sciplex4
folder= 'Supp_GSM4150379_sciPlex4_A549_MCF7_HDACi'
dataset = folder.replace('Supp_', '')
adata = prepare(folder)

dataset = folder.replace('Supp_', '')
hashTable = pd.read_csv(f'{path}{folder}/{dataset}_hashTable.out.txt.gz', sep='\t', header=None, index_col=1,
                        names=['sample_id', 'cell_id', 'delstring', 'axis_id', 'hash_UMI_counts'])
# split the delstring
delstring = pd.DataFrame([row.delstring.split('_') for index, row in tqdm(hashTable.iterrows())], index=hashTable.index,
                         columns=['dose_value', 'perturbation', 'dose_value_2', 'perturbation_2', 'cell_line', 'plate_id', 'well_id'])
# keep index with most hash_UMI_counts
annotation = pd.concat([hashTable.drop('delstring',axis=1), delstring], axis=1).sort_values('hash_UMI_counts', ascending=False)
annotation = annotation[~annotation.index.duplicated(keep='first')]
# overwrite obs
overloaded_obs = pd.concat([adata.obs.drop(adata.obs.columns, axis=1), annotation], axis=1)
adata.obs = overloaded_obs.loc[adata.obs_names]

adata.obs['perturbation'][adata.obs['perturbation']=='DMSO']='control'
adata.obs['perturbation_2'][adata.obs['perturbation_2']=='DMSO']='control'
adata.obs['celltype'] = ['alveolar basal epithelial cells' if line=='A549' else 'mammary epithelial cells' if line=='MCF7' else 'None' for line in adata.obs.cell_line]
adata.obs['cancer'] = True
adata.obs['disease'] = ['lung adenocarcinoma' if line=='A549' else 'breast adenocarcinoma' if line=='MCF7' else 'None' for line in adata.obs.cell_line]
adata.obs['tissue_type'] = 'cell_line'
adata.obs['organism'] = 'human'
adata.obs['perturbation_type'] = 'drug'
adata.obs['perturbation_type_2'] = 'drug'
adata.obs = adata.obs.drop(['sample_id', 'axis_id'], axis=1)
x1 = np.array([0 if x==None else 0 if x=='control' else 1 for x in adata.obs.perturbation])
x2 = np.array([0 if x==None else 0 if x=='control' else 1 for x in adata.obs.perturbation_2])
adata.obs['nperts'] = x1+x2
adata.var_names_make_unique()
adata.write(WDIR+'SrivatsanTrapnell2020_sciplex4.h5')
