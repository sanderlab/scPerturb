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

### FILES ###
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
# ['GSM4494347', 0, 'PC9_rep1', 'control', 'multirep_PC9_mono', 'GSE149214'], #Day00b scRNA-seq
# ['GSM4494348', 11, 'PC9_rep1', 'Erl', 'multirep_PC9_mono', 'GSE149214'], #Day11b scRNA-seq
# ['GSM4494349', 11, 'PC9_rep2', 'Erl', 'multirep_PC9_mono', 'GSE149214'], #Day11c scRNA-seq
# GSE149215 (PC9 time course Multi-drug)
# ['GSM4494350', 3, 'PC9_rep1', 'Erl', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-ERL1 scRNA-seq
# ['GSM4494351', 3, 'PC9_rep2', 'Erl', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-ERL2 scRNA-seq
# ['GSM4494352', 3, 'PC9_rep1', 'Erl+Cri', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-CRI-ERL1 scRNA-seq
# ['GSM4494353', 3, 'PC9_rep2', 'Erl+Cri', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-CRI-ERL2 scRNA-seq
# ['GSM4494354', 3, 'PC9', 'Eto', 'multirep_PC9_multi', 'GSE149215'], #PC9Day3-Eto scRNA-seq
# GSE160244 (xenografts mice)
['GSM4869650', 3, 'PC9_xeno', 'control', 'Xenograft', 'GSE160244'], #Day 3 control
['GSM4869651', 3, 'PC9_xeno', 'Cri', 'Xenograft', 'GSE160244'], #Day 3 crizotinib
['GSM4869652', 3, 'PC9_xeno', 'Osi', 'Xenograft', 'GSE160244'], #Day 3 osimertinib
['GSM4869653', 3, 'PC9_xeno', 'Osi+Cri', 'Xenograft', 'GSE160244'] #Day 3 osimertinib+crizotinib
]
tab = pd.DataFrame(X, columns=['GEO', 'time', 'cell_line', 'perturbation', 'batch', 'subseries'])
tab.time = tab.time *24  # make hours
tab['replicate'] = [x.split('_')[-1] if 'rep' in x else None for x in tab.cell_line]
tab['cell_line'] = [x.split('_')[0]+'_xenograft' if 'xeno' in x else x.split('_')[0] for x in tab.cell_line]
tab = tab.set_index('GEO')
# subseries GSE149215 and GSE149214 not found
tab['filename'] = ''
for GEO in tab.index:
    files = [f.name for f in (TEMPDIR / 'AissaBenevolenskaya2021/').glob(f'{GEO}*.dge.txt')]
    assert len(files) == 1
    tab.loc[GEO, 'filename'] = files[0].replace('.dge.txt', '')

### PROCESSING ###
adatas = {}
for index, row in tab.iterrows():
    # load
    print(index)
    Y = pd.read_csv(TEMPDIR / f'AissaBenevolenskaya2021/{row.filename}.dge.txt', sep='\t', index_col=0)
    adata = sc.AnnData(csr_matrix(Y.values.T, dtype=int))
    adata.var_names = Y.index
    adata.obs_names = Y.columns
    adata.var_names_make_unique()

    # annotate by table
    for col in tab.columns:
        adata.obs[col] = row[col]
    adatas[index] = adata

### MERGE ###
adata = sc.concat(adatas, index_unique='-')
# reform var, obs index names
adata.var.index = adata.var.index.rename('gene_symbol')
adata.obs.index = adata.obs.index.rename('cell_barcode')
# abbreviations to full compound name
for x,y in zip(['Vem', 'Erl', 'Cri', 'Eto', 'Osi'], ['vemurafenib', 'erlotinib', 'crizotinib', 'etoposide', 'osimertinib']):
    adata.obs.perturbation = [per.replace(x, y) for per in adata.obs.perturbation]
# further metadata
adata.obs['tissue_type']= ['cell_line' if x!='Xenograft' else 'cell_line' for x in adata.obs.batch]
adata.obs['cancer'] = True
adata.obs['perturbation_type'] = 'drug'
adata.obs['disease'] = ['melanoma' if x=='M14' else 'lung adenocarcinoma' if x=='PC9' else None for x in adata.obs.cell_line]
adata.obs['celltype'] = ['melanocytes' if x=='M14' else 'epithelial cells' if x=='PC9' else None for x in adata.obs.cell_line]
adata.obs['organism'] = 'human'
adata.obs.perturbation = adata.obs.perturbation.astype(str)
adata.obs['perturbation'] = [x.replace('+', '_') for x in adata.obs.perturbation]
adata.obs['nperts'] = [1+p.count('_')-p.count('control')-p.count('None') if type(p)==str else 0 for p in adata.obs.perturbation]

annotate_qc(adata, species='human')
assert_annotations(adata)

adata.obs.replicate=adata.obs.replicate.astype(str)
adata.write(snakemake.output[0], compression='gzip')