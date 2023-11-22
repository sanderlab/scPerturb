import pandas as pd
import scanpy as sc
import numpy as np
import sys

from scipy.io import mmread
from scipy.sparse import csr_matrix, vstack
from pathlib import Path
from tqdm import tqdm

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR'])
DOWNDIR = Path(snakemake.config['DOWNDIR'])

# maps experiments to file identifiers
sample_dict = {
    'leukemia': 'DM',
    'invivo': 'inVivo',
    'exvivo': 'LSK'
}
files = [x.name for x in (TEMPDIR / 'LaraAstiasoHuntly2023/').glob('*.h5')]

def merge_data(key):
    identifier = sample_dict[key]
    files_ = [x for x in files if f'_{identifier}_' in x]
    annot = pd.read_csv(TEMPDIR / 'LaraAstiasoHuntly2023/' / f'GSE213511_CellAnnotation_{key}.tsv.gz', sep='\t', index_col=0)
    adatas = {}
    for file in tqdm(files_):
        tempdata = sc.read_10x_h5(TEMPDIR / 'LaraAstiasoHuntly2023' / file)
        tempdata.var_names_make_unique()
        adatas[file.replace('GSE213511_', '').replace('.h5', '')] = tempdata
    adata = sc.concat(adatas, label='Sample', index_unique='-')
    # make indiced unique
    adata.obs.index = [x.replace('-1-', '-') for x in adata.obs.index]
    annot.index = [f'{x.split("-")[0]}-{sample}' for x, sample in zip(annot.index, annot.Sample)]
    # merge annotation
    if key=='invivo':
        annot = annot[[x in adata.obs.index for x in annot.index]]  # samples ['inVivo_OP2_Lin-_28d_1', 'inVivo_OP3_Lin-_28d_1'] are missing in the data
    adata.obs = pd.concat([adata.obs, annot], axis=1)
    return adata

def harmonize_data(adata, key):
    # harmonize
    adata.var.index.name = 'gene_symbol'
    adata.obs.index.name = 'cell_barcode'
    adata.obs['organism'] = 'Mus musculus'
    adata.obs['disease'] = 'leukemia' if key=='leukemia' else 'healthy'
    adata.obs['cancer'] = key=='leukemia'
    adata.obs['perturbation_type'] = 'CRISPR-cas9'
    adata.obs['tissue_type'] = 'primary'
    adata.obs['tissue'] = 'bone marrow transplant'
    
    adata.obsm['X_umap'] = adata.obs[['UMAP1', 'UMAP2']].values
    adata.obs = adata.obs.loc[:, ~adata.obs.columns.duplicated(keep='first')]  # remove duplicated "Sample" column
    adata.obs['perturbation'] = np.array([None if pd.isna(x) else 'control' if x[:3]=='NTC' else x.split('_')[0] for x in adata.obs.Guide])
    adata.obs.rename({
        'Sample': 'sample',
        'Phase': 'cellcycle_phase',
        'Clusters': 'celltype',
        'Mixscape': 'Mixscape_classification',
        'mixscape': 'Mixscape_classification',
        'Guide': 'guide_id',
        'Timepoint': 'time'
    }, axis=1, inplace=True)
    # reorder
    order = ['perturbation', 'guide_id', 'sample', 'cellcycle_phase', 'Mixscape_classification',
       'celltype', 'organism', 'disease', 'cancer',
       'perturbation_type', 'tissue_type', 'tissue']
    if 'time' in adata.obs.columns: order = ['time'] + order 
    adata.obs = adata.obs[order]
    annotate_qc(adata)
    return adata

for key in ['leukemia', 'invivo', 'exvivo']:
    adata = merge_data(key)
    bdata = harmonize_data(adata, key)
    assert_annotations(adata)
    bdata.write(DOWNDIR / f'LaraAstiasoHuntly2023_{key}.h5ad', compression='gzip')
