import scanpy as sc
import pandas as pd
import numpy as np
import re
from tqdm import tqdm
from pathlib import Path

# Custom functions
sys.path.insert(1, '../')
from utils import *

TEMPDIR = Path(snakemake.config['TEMPDIR']) 
files = sorted([file for file in (TEMPDIR / 'QinTape2023/Count Matrices & CellRanger Reports').glob('*') if 'DS_Store' not in file.name])

adatas = {}
for file in tqdm(files):
    adata = sc.read_10x_mtx(file)
    sample_name = file.name
    sample_number, sample_description = sample_name.split('_', 1)
    
    reg = re.match('[AKP]+|(WT)', sample_description)
    adata.obs['perturbation'] = reg.group() if reg else None
    adata.obs['perturbation_type'] = 'genotype'
    
    adata.obs['perturbation_2'] = 'Mac-Fib' if 'Mac-Fib' in sample_description else 'Mac' if 'Mac' in sample_description else 'Fib' if 'Fib' in sample_description else None
    adata.obs['perturbation_type_2'] = 'coculture'
    
    adata.obs['sample_number']=sample_number
    adatas[sample_name] = adata
    if int(sample_number) >= 19:
        # stop processing, additional controls with other conditions.
        # we focus on the core datasets (1-19)
        break

adata = sc.concat(adatas, label='batch')

# set metadata
adata.obs['disease'] = "colorectal cancer"
adata.obs['cancer'] = True
adata.obs['tissue_type']="cell_line"
adata.obs["cell_line"] = "CRC organoid"
adata.obs["celltype"] = 'colon epithelial cell'
adata.obs['organism'] = 'mouse'
adata.obs['nperts'] = (~adata.obs.perturbation.isna()*1 + ~adata.obs.perturbation_2.isna()*1)
annotate_qc(adata, species='mouse')
adata.obs.index.name = 'cell_barcode'

assert_annotations(adata)
adata.write(snakemake.output[0], compression='gzip')
