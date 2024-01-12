# !/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on 22.11.2023
Author: Stefan Peidli
Description: Reintegrates files created from extract.R (Rscript) into h5ad.
"""

import scanpy as sc
import pandas as pd

from argparse import ArgumentParser
from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path

parser = ArgumentParser(
                    prog = 'Remerger',
                    description = 'Reintegrates files created from extract.R (Rscript).',
                    epilog = 'Use like this: merge_from_extract.py -i path_input -o output_file -n name_prefix')
parser.add_argument('-i', '--input_path', help='Path to folder with extracted files', required=True)
parser.add_argument('-o', '--output_file', help='Path to h5ad output file', required=True)
parser.add_argument('-n', '--name_prefix', help='Prefix name of the files', default='seurat')

args = parser.parse_args()
path = Path(args.input_path)
prefix = args.name_prefix
output_file = Path(args.output_file) if '.h5' in args.output_file else Path(f'{args.output_file}.h5ad')

# load data
print('Reading raw matrix...')
X = csr_matrix(mmread(path / f'{prefix}_matrix.mtx')).T
print('Reading additional data...')
# obs = pd.read_csv(f'{path}barcodes.csv', index_col=0, usecols=[1], names=['cell_barcodes'], skiprows=1)  # not needed, contained in metadata
var = pd.read_csv(path / f'{prefix}_genes.csv', index_col=0, usecols=[1], names=['genes'], skiprows=1)
obs = pd.read_csv(path / f'{prefix}_metadata.csv', index_col=0)

# merge
print('Merging...')
adata = sc.AnnData(X, obs, var)

# Add optional data (embeddings)
for reduction in ['pca', 'umap', 'tsne', 'diffmap']:
    reduction_path = path / f'{prefix}_{reduction}.csv'
    if reduction_path.exists():
        print(f'Adding {reduction}...')
        reduction_data = pd.read_csv(reduction_path, index_col=0)
        adata.obsm[f'X_{reduction}'] = reduction_data.values

# write
print('Writing...')
adata.write(output_file)
print('Done!')
