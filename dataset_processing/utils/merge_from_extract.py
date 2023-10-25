import scanpy as sc
import pandas as pd

from argparse import ArgumentParser
from scipy.io import mmread
from scipy.sparse import csr_matrix
from pathlib import Path


parser = ArgumentParser(
                    prog = 'Remerger',
                    description = 'Reintegrates files created from extract.R (Rscript).',
                    epilog = 'Use like this: merge_from_extract.py -i /path/to/extracted/files/ -o /path/to/output_file -n name_prefix')
parser.add_argument('-i', '--input_path', help='Path to folder with extracted files', required=True)
parser.add_argument('-o', '--output_file', help='Path to h5ad output file', required=True)
parser.add_argument('-n', '--name_prefix', help='Prefix name of the files', default='seurat')
args = parser.parse_args()

input_path = Path(args.input_path)
name_prefix = args.name_prefix

# load data
print('Reading raw matrix...')
X = csr_matrix(mmread(input_path / f'{name_prefix}_matrix.mtx')).T
print('Reading additional data...')
# obs = pd.read_csv(f'{path}barcodes.csv', index_col=0, usecols=[1], names=['cell_barcodes'], skiprows=1)
obs = pd.read_csv(input_path / f'{name_prefix}_metadata.csv', index_col=0)
var = pd.read_csv(input_path / f'{name_prefix}_genes.csv', index_col=0, usecols=[1], names=['genes'], skiprows=1)
pca = pd.read_csv(input_path / f'{name_prefix}_pca.csv', index_col=0)
umap = pd.read_csv(input_path / f'{name_prefix}_umap.csv', index_col=0)

# merge
print('Merging...')
adata = sc.AnnData(X, obs, var)
adata.obsm['X_umap'] = umap.values
adata.obsm['X_pca'] = pca.values

# write
print('Writing...')
file = args.output_file if '.h5' in args.output_file else f'{args.output_file}.h5ad'
adata.write(file)
print('Done!')