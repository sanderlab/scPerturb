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

path = DIR + 'GSE149383/supp/'
# format: ['GEO', 'days', 'origin', 'condition', 'bundle', 'subseries']

X=[
# GSE134836 (Erl PC9) is weird
# ['GSM3972651', 0, 'PC9', 'control', 'PC9D0vsD3Erl', 'GSE134836'],  #PC9D0
# ['GSM3972652', 3, 'PC9', 'Erl', 'PC9D0vsD3Erl', 'GSE134836'],  #PC9D3Erl
# GSE134838 (Vem M14)
['GSM3972655', 0, 'M14', 'control', 'M14D0vsD3Vem', 'GSE134838'], #M14Day0
['GSM3972656', 3, 'M14', 'Vem', 'M14D0vsD3Vem', 'GSE134838'], #M14Day3_Vem
# GSE134839 (Erl PC9)
['GSM3972657', 0, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 0 PC9 erlotinib
['GSM3972658', 1, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 1 PC9 erlotinib
['GSM3972659', 2, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 2 PC9 erlotinib
['GSM3972660', 4, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 4 PC9 erlotinib
['GSM3972661', 9, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 9 PC9 erlotinib
['GSM3972662', 11, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 11 PC9 erlotinib
# GSE134841 (Drug holiday run Erl PC9)
['GSM3972669', 0, 'PC9', 'Erl', 'Drug_holiday', 'GSE134841'], #PC9 Day 0
['GSM3972670', 2, 'PC9', 'Erl', 'Drug_holiday', 'GSE134841'], #PC9 Day 2
['GSM3972671', 11, 'PC9', 'Erl', 'Drug_holiday', 'GSE134841'], #PC9 Day 11
['GSM3972672', 19, 'PC9', 'Erl_DMSO', 'Drug_holiday', 'GSE134841'], #PC9 D19 DMSO
['GSM3972673', 19, 'PC9', 'Erl_Erl', 'Drug_holiday', 'GSE134841'], #PC9 D19 ERL
# GSE134842 (Beads on mix)
# ['GSM3972674', 0, 'PC9+U937', 'control', 'bead_mix', 'GSE134842'], #PC9 U937 mix
# ['GSM3972675', 0, 'PC9+U937', 'EpCAM+', 'bead_mix', 'GSE134842'], #EpCAM-positive by beads
# ['GSM3972676', 0, 'PC9+U937', 'CD45-', 'bead_mix', 'GSE134842'], #CD45-negative by beads
# GSE148466 (primary tissues)
# ['GSM4472055', 0, 'P8126999_primary', 'control', 'primary', 'GSE148466'], #P8126999
# ['GSM4472056', 0, 'P8127005_primary', 'control', 'primary', 'GSE148466'], #P8127005
# ['GSM4472057', 0, 'P8127011_primary', 'control', 'primary', 'GSE148466'], #P8127011
# GSE149214 (PC9 time course Mono-drug)
['GSM4494347', 0, 'PC9_rep1', 'control', 'multirep_PC9_mono', 'GSE149214'], #Day00b scRNA-seq
['GSM4494348', 11, 'PC9_rep1', 'Erl', 'multirep_PC9_mono', 'GSE149214'], #Day11b scRNA-seq
['GSM4494349', 11, 'PC9_rep2', 'Erl', 'multirep_PC9_mono', 'GSE149214'], #Day11c scRNA-seq
# GSE149215 (PC9 time course Multi-drug)
['GSM4494350', 3, 'PC9_rep1', 'Erl', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-ERL1 scRNA-seq
['GSM4494351', 3, 'PC9_rep2', 'Erl', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-ERL2 scRNA-seq
['GSM4494352', 3, 'PC9_rep1', 'Erl+Cri', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-CRI-ERL1 scRNA-seq
['GSM4494353', 3, 'PC9_rep2', 'Erl+Cri', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-CRI-ERL2 scRNA-seq
['GSM4494354', 3, 'PC9', 'Eto', 'multirep_PC9_multi', 'GSE149215'], #PC9Day3-Eto scRNA-seq
# GSE160244 (xenografts mice)
['GSM4869650', 3, 'PC9_xeno', 'control', 'Xenograft', 'GSE160244'], #Day 3 control
['GSM4869651', 3, 'PC9_xeno', 'Cri', 'Xenograft', 'GSE160244'], #Day 3 crizotinib
['GSM4869652', 3, 'PC9_xeno', 'Osi', 'Xenograft', 'GSE160244'], #Day 3 osimertinib
['GSM4869653', 3, 'PC9_xeno', 'Osi+Cri', 'Xenograft', 'GSE160244'] #Day 3 osimertinib+crizotinib
]
GSMs = [x.split('_')[1] for x in get_subfolders(path, False)]
tab = pd.DataFrame(X, columns=['GEO', 'time', 'cell_line', 'perturbation', 'batch', 'subseries'])
tab.time = tab.time *24  # make hours
tab['replicate'] = [x.split('_')[-1] if 'rep' in x else None for x in tab.cell_line]
tab['cell_line'] = [x.split('_')[0]+'_xenograft' if 'xeno' in x else x.split('_')[0] for x in tab.cell_line]


folders = get_subfolders(path)
alldatas = {}
for subseries in pd.unique(tab.subseries):
    adatas = {}
    for index, row in tab[tab.subseries==subseries].iterrows():
        # determine file for GEO of subseries
        file = get_files([x for x in folders if row.GEO in x][0], True)[0]

        # load
        Y = pd.read_csv(file, sep='\t', index_col=0)
        adata = sc.AnnData(csr_matrix(Y.values.T, dtype=int))
        adata.var_names = Y.index
        adata.obs_names = Y.columns
        adata.var_names_make_unique()

        # annotate by table
        for col in tab.columns:
            adata.obs[col] = row[col]
        adatas[row.GEO] = adata
    adata = sc.concat(adatas, index_unique='-')

    # reform var, obs index names
    adata.var.index = adata.var.index.rename('gene_symbol')
    adata.obs.index = adata.obs.index.rename('cell_barcode')

    # abbreviations to full compound name
    for x,y in zip(['Vem', 'Erl', 'Cri', 'Eto', 'Osi'], ['vemurafenib', 'erlotinib', 'crizotinib', 'etoposide', 'osimertinib']):
        adata.obs.perturbation = [per.replace(x, y) for per in adata.obs.perturbation]

    # further metadata
    adata.obs['tissue_type']='cell_line'
    adata.obs['cancer'] = True
    adata.obs['perturbation_type'] = 'drug'
    adata.obs['disease'] = ['melanoma' if x=='M14' else 'lung adenocarcinoma' if x=='PC9' else None for x in adata.obs.cell_line]
    adata.obs['celltype'] = ['melanocytes' if x=='M14' else 'epithelial cells' if x=='PC9' else None for x in adata.obs.cell_line]
    adata.obs['organism'] = 'human'
    adata.obs['nperts'] = [1+('+' in x) for x in adata.obs.perturbation]
    alldatas[subseries] = adata

superdata = sc.concat(alldatas, index_unique='-')
superdata.write(WDIR+'AissaBenevolenskaya2021.h5')
write_as_singles(superdata, '/fast/work/users/peidlis_c/data/perturbation_resource_paper/', 'AissaBenevolenskaya2021', add_h5=True)
