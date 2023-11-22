import pandas as pd
import scanpy as sc
import numpy as np
import sys
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations
from scipy.sparse import csr_matrix

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

# Load data
adata = sc.read(snakemake.input['adata'])
prefix = snakemake.wildcards['prefix']
hashtable = snakemake.input['hashtable']

hashTable = pd.read_csv(hashtable, sep='\t', index_col=1, header=None,
                        names = ['orig_ident', 'cell_barcode', 'plate_id', 'all_ones', 'counts'])
# cell_barcodes --> count hashes for plate_id
df_hash = pd.pivot_table(hashTable, values='counts', index='cell_barcode', columns='plate_id', fill_value=0)
# just take highest
most_counts = np.argmax(df_hash.values, axis=1)
best_plates = pd.DataFrame(list(df_hash.columns[most_counts]), index = df_hash.index, columns=['best_plate_id'])

adata.write(snakemake.output[0], compression='gzip')