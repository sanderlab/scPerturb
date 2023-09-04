import pandas as pd
import anndata
from pathlib import Path
from scipy import sparse
from typing import Dict
import numpy as np
import sys

sys.path.append(str(
    (Path(__file__).parent.parent / 'code').resolve())
    )
import utils

scATAC_output_folder = Path(__file__).parent / 'output'
EXPORT_ROOT_DIR = scATAC_output_folder / 'export'

cellline2disease = {
    'MCF7' : 'adenocarcinoma',
    'K562' : 'myelogenous leukemia',
    'GM12878' : 'healthy',
    'CD4+_T_cells' : 'healthy',
}

cellline2iscancer = {
    'MCF7' : True,
    'K562' : True,
    'GM12878' : False,
    'CD4+_T_cells' : False,
}

cellline2organism = {
    'MCF7' : 'human',
    'K562' : 'human',
    'GM12878' : 'human',
    'CD4+_T_cells' : 'human',
}

dataset2perturbation_type={
    'Spear_ATAC' : 'CRISPR',
    'CRISPR_sciATAC' : 'CRISPR',
    'ASAP_seq' : 'CRISPR'
}

dataset2tissue_type = {
    'Spear_ATAC' : 'cell_line',
    'CRISPR_sciATAC' : 'cell_line',
    'ASAP_seq' : 'primary',
}


def process_bcs_Spear_ATAC(bc: pd.Series, sample_annotation: np.array, cell_line: str)->pd.Series:
    bc = bc.str.split('#').str[1]
    bc = bc.str[:-2]
    bc += ('-rep' + sample_annotation)
    bc += ('-' + cell_line)
    return bc


def compute_ns(X, frame, sum_axis):
    n1 = 'ncounts'
    if sum_axis == 0:
        n2 = 'ncells'
    elif sum_axis == 1:
        n2 = 'ngenes'
    else:
        raise Exception(f'Chosen axis {sum_axis} invalid.')

    frame[n1] = np.squeeze(np.array(X.sum(axis = sum_axis)))
    frame[n2] = np.squeeze(np.array((X != 0).sum(axis = sum_axis,  dtype = int)))
    return


def get_data_set_info():
    data_dirs = []
    for data_set_path in scATAC_output_folder.iterdir():
        if not data_set_path.is_dir():
            continue
        data_set = data_set_path.name
        if data_set == 'export': #that is the output folder
            continue
        for cell_line_path in data_set_path.iterdir():
            if not cell_line_path.is_dir():
                continue
            cell_line = cell_line_path.name
            export_dir = EXPORT_ROOT_DIR / data_set / cell_line
            if data_set == 'CRISPR_sciATAC':
                cell_line = "K562"
            data_dirs.append({
                'data_set' : data_set,
                'cell_line' : cell_line,
                'data_dir' : cell_line_path / 'export',
                'export_dir' : export_dir
                }
            )
    return data_dirs


def parse_peak_bc_matrix(data_set_info: Dict) ->anndata.AnnData:
    """Load and annotate peak-bc matrix and convert to AnnData.

    :param data_set_info: Infos about selected dataset.
    :type data_set_info: Dict
    :return: Anndata object
    :rtype: anndata.AnnData
    """

    peak_bc_dir =  data_set_info['data_dir']
    cell_line = data_set_info['cell_line']
    data_set = data_set_info['data_set']

    peak_bc = pd.read_csv(peak_bc_dir / "peak_bc.csv.gz", index_col = 0)
    peak_set = pd.read_csv(peak_bc_dir / "PeakSet.csv", index_col = 0)
    # In ArchR the order of peakSet is not preserved in peak Matrix, which was acknowledged here:
    # https://github.com/GreenleafLab/ArchR/issues/1148
    # The problem is that chromosomes are ordered as 1 10 11 12... instead of 1 2 3 4...
    # I can establish the right order by doing this
    peak_set = peak_set.sort_values(['seqnames', 'idx'])
    #which can be verified as follows:
    # peak_bc_row = pd.read_csv(data_set_info['data_dir'] / 'peak_bc_rowData.csv', index_col = 0)
    # print(np.all(peak_set['idx'].values == peak_bc_row['idx'].values))
    # peak_set_sorted.groupby('seqnames')['idx'].max()
    # peak_set_sorted.groupby('seqnames')['idx'].size()

    peak_set = peak_set.drop(
        columns=['strand', "width", "score", "replicateScoreQuantile", "groupScoreQuantile", "Reproducibility", "GroupReplicate", "idx", "N"]
        )
    peak_set = peak_set.rename(columns = {'seqnames' : 'chromosome'})
    peak_set = peak_set.reset_index(drop=True) #pythonic zero-based index
    peak_set.index.name = 'peak_index'
    #For cells it's easier because I have barcode annotation. This has all the information that the cell_meta_data.csv has.
    peak_bc_col = pd.read_csv(peak_bc_dir / 'peak_bc_colData.csv', index_col = 0)
    # cell_meta_data = pd.read_csv(data_dir['data_dir'] / "cell_meta_data.csv", index_col = 0)

    if data_set == 'Spear_ATAC':
        peak_bc_col = peak_bc_col.drop(
            columns=['PassQC', 'sgAssign2', 'sgAssign3', 'sgAssignClean', 'sgIndividual']
            )
        peak_bc_col = peak_bc_col.rename(
            columns={"sgAssign" : "guide", "sgAssignFinal" : "perturbation", 'Sample' : 'sample'}
            )

        peak_bc_col['sample'] = peak_bc_col['sample'].str.split('_').str[1]

        peak_bc_col.loc[~peak_bc_col.perturbation.isna(), 'perturbation'] = \
            peak_bc_col.loc[~peak_bc_col.perturbation.isna(), 'perturbation'].str[2:]
        peak_bc_col = peak_bc_col.replace({'perturbation': {'sgNT': 'control'}})
    elif data_set == 'CRISPR_sciATAC':
        peak_bc_col = peak_bc_col.drop(
            columns=['PassQC', 'sgAssign2', 'sgAssign3', 'sgAssignClean', 'sgIndividual', "nMonoFrags", "nMultiFrags", "NucleosomeRatio"]
            )
        peak_bc_col = peak_bc_col.rename(
            columns={"sgAssign" : "guide", "sgAssignFinal" : "perturbation", 'Sample' : 'sample'}
            )
        peak_bc_col = peak_bc_col.replace({'perturbation': {'sgsgNT': 'control'}})
    elif data_set == 'ASAP_seq':
        peak_bc_col = peak_bc_col.drop(
            columns=['PassQC', 'sgAssign2', 'sgAssign3', 'sgAssignClean', 'sgIndividual']
            )
        peak_bc_col = peak_bc_col.rename(
            columns={"sgAssign" : "guide", "sgAssignFinal" : "perturbation", 'Sample' : 'sample'}
            )
        peak_bc_col.loc[~peak_bc_col.perturbation.isna(), 'perturbation'] = \
            peak_bc_col.loc[~peak_bc_col.perturbation.isna(), 'perturbation'].str[2:]
        peak_bc_col = peak_bc_col.replace({'perturbation': {'sgNT': 'control'}})

    else:
        raise Exception('Dataset specifics not defined.')
    

    # -1 because R starts indexing from 1
    peak_bc = sparse.coo_matrix(
        (peak_bc['x'], (peak_bc['i'] - 1 , peak_bc['j'] - 1)), 
        shape= (peak_set.shape[0], peak_bc_col.shape[0]) )
    peak_bc = peak_bc.tocsr()

    scATAC = anndata.AnnData(
        X = peak_bc.T,
        obs = peak_bc_col,
        var = peak_set)

    scATAC.obs_names.name =  'cell_barcode'

    #add all the missing required fields for observations
    scATAC.obs['disease'] = cellline2disease[cell_line]
    scATAC.obs['cancer'] = cellline2iscancer[cell_line]
    scATAC.obs['tissue_type'] = dataset2tissue_type[data_set]
    scATAC.obs['cell_line'] = cell_line
    scATAC.obs['organism'] = cellline2organism[cell_line]
    scATAC.obs['perturbation_type'] = dataset2perturbation_type[data_set]
    scATAC.obs['nperts'] = 1
    scATAC.obs['percent_mito'] = np.nan
    scATAC.obs['percent_ribo'] = np.nan
    # scATAC.obs['ncounts'] = np.squeeze(np.array(scATAC.X.sum(axis = 1,  dtype = int)))
    # scATAC.obs['ngenes'] = np.squeeze(np.array((scATAC.X != 0).sum(axis = 1,  dtype = int)))
    compute_ns(scATAC.X, scATAC.obs, 1)

    #add all the missing required fields for variables
    # scATAC.var['ncounts'] = np.squeeze(np.array(scATAC.X.sum(axis = 0,  dtype = int)))
    # scATAC.var['ncells'] = np.squeeze(np.array((scATAC.X != 0).sum(axis = 0,  dtype = int)))
    compute_ns(scATAC.X, scATAC.var, 0)
    return scATAC


def parse_LSI_embedding(data_set_info: Dict) ->pd.DataFrame:
    data_dir =  data_set_info['data_dir']
    LSI = pd.read_csv(data_dir / 'LSI.csv', index_col=0)
    return LSI


def parse_gene_scores(data_set_info: Dict) ->pd.DataFrame:
    data_dir =  data_set_info['data_dir']

    gene_scores = pd.read_csv(data_dir / 'gene_scores.csv.gz', index_col = 0)
    gene_scores_row = pd.read_csv(data_dir / 'gene_scores_rowData.csv', index_col = 0)
    gene_scores_col = pd.read_csv(data_dir / 'gene_scores_colData.csv', index_col = 0)    

    gene_scores_col = gene_scores_col.index
    gene_scores_row = pd.Index(gene_scores_row['name'].values, name = 'gene_symbol')

    # -1 because R starts indexing from 1
    gene_scores = sparse.coo_matrix(
        (gene_scores['x'], (gene_scores['i'] - 1 , gene_scores['j'] - 1)), 
        shape= (gene_scores_row.shape[0], gene_scores_col.shape[0]) )
    # gene_scores = gene_scores.tocsr()

    gene_scores = pd.DataFrame.sparse.from_spmatrix(
        data = gene_scores,
        index = gene_scores_row,
        columns = gene_scores_col).T

    return gene_scores


def parse_markerpeak_target(data_set_info: Dict) ->pd.DataFrame:
    data_dir =  data_set_info['data_dir']
    data_set = data_set_info['data_set']

    markerpeak_target = pd.read_csv(data_dir / 'marker_peak_bc.csv', index_col = 0)
    if data_set == 'Spear_ATAC':
        markerpeak_target.index = markerpeak_target.index.str[2:]
        markerpeak_target.index = markerpeak_target.index.to_series().replace({'sgNT': 'control'})
    elif data_set == 'CRISPR_sciATAC':
        markerpeak_target.index = markerpeak_target.index.to_series().replace({'sgsgNT': 'control'})
    elif data_set == 'ASAP_seq':
        markerpeak_target.index = markerpeak_target.index.str[2:]
        markerpeak_target.index = markerpeak_target.index.to_series().replace({'sgNT': 'control'})
    else:
        raise Exception('Dataset specifics not defined.')

    return markerpeak_target



def parse_ChromVar(data_set_info: Dict) ->pd.DataFrame:
    data_dir =  data_set_info['data_dir']
    ChromVar = pd.read_csv(data_dir / 'ChromVar.csv.gz', index_col = 0)
    ChromVar_row = pd.read_csv(data_dir / 'ChromVar_rowData.csv', index_col = 0)
    ChromVar_col = pd.read_csv(data_dir / 'ChromVar_colData.csv', index_col = 0)
    ChromVar_row = ChromVar_row.index
    ChromVar_row.name = 'TF_motif'
    ChromVar_col = ChromVar_col.index
    ChromVar = sparse.coo_matrix(
        (ChromVar['x'], (ChromVar['i'] - 1 , ChromVar['j'] - 1)), 
        shape= (ChromVar_row.shape[0], ChromVar_col.shape[0]) )

    ChromVar = pd.DataFrame.sparse.from_spmatrix(
        data = ChromVar,
        index = ChromVar_row,
        columns = ChromVar_col).T
    return ChromVar


def separate_and_save_anndata(adata: anndata.AnnData, output_dir: Path, save = True)->None:

    # Because of h5 limitations it is not possible to save an anndata object with
    # obsm/varm/uns DataFrames that have many columns (which is the case for e.g. gene_scores)
    # Thus we agreed, that we would seperate each matrix and save it as separate files.

    adatas = {}

    #peak_bc
    X = adata.X
    obs = adata.obs.copy()
    var = adata.var.copy()
    adatas['peak_bc'] = anndata.AnnData(X = X, obs = obs, var = var)

    #LSI_embedding
    X = adata.obsm['LSI_embedding'].values
    obs = adata.obs.copy()
    compute_ns(X, obs, 1)
    var = pd.DataFrame(data = None, index = adata.obsm['LSI_embedding'].columns)
    compute_ns(X, var, 0)
    adatas['LSI_embedding'] = anndata.AnnData(X = X, obs = obs, var = var)

    #gene_scores
    X = adata.obsm['gene_scores'].sparse.to_coo()
    X = X.tocsr() # anndata does not support coo when writing to h5.
    obs = adata.obs.copy()
    compute_ns(X, obs, 1)
    var = pd.DataFrame(data = None, index = adata.obsm['gene_scores'].columns)
    compute_ns(X, var, 0)
    adatas['gene_scores'] = anndata.AnnData(X = X, obs = obs, var = var)

    #markerpeak_target
    X = adata.uns['markerpeak_target'].values
    obs = pd.DataFrame(data = None, index = adata.uns['markerpeak_target'].index)
    compute_ns(X, obs, 1)
    var = pd.DataFrame(data = None, index = adata.uns['markerpeak_target'].columns)
    compute_ns(X, var, 0)
    adatas['markerpeak_target'] = anndata.AnnData(X = X, obs = obs, var = var)


    #ChromVar
    X = adata.obsm['ChromVar'].values
    obs = adata.obs.copy()
    compute_ns(X, obs, 1)
    var = pd.DataFrame(data = None, index = adata.obsm['ChromVar'].columns)
    compute_ns(X, var, 0)
    adatas['ChromVar'] = anndata.AnnData(X = X, obs = obs, var = var)

    if save:
        for name, adata in adatas.items():
            utils.write_as_singles(adata, output_dir, name, add_h5=True)

    return adatas

def run():
    data_set_infos = get_data_set_info()

    scATACs = []
    for data_set_info in data_set_infos:
        print( f"Processing {data_set_info['data_set']} dataset with {data_set_info['cell_line']} cell line")
        if data_set_info['data_set'] != 'CRISPR_sciATAC':
            continue

        scATAC = parse_peak_bc_matrix(data_set_info)

        # Add LSI embedding #
        #####################
        LSI = parse_LSI_embedding(data_set_info)
        scATAC.obsm['LSI_embedding'] = LSI.reindex(scATAC.obs_names)

        # Add gene scores #
        ###################
        gene_scores = parse_gene_scores(data_set_info)
        scATAC.obsm['gene_scores'] = gene_scores.reindex(scATAC.obs_names)

        #add marker peaks#
        ##################
        markerpeak_target = parse_markerpeak_target(data_set_info)
        scATAC.uns['markerpeak_target'] = markerpeak_target

        #add Chrom Var#
        ###############
        ChromVar = parse_ChromVar(data_set_info)
        ChromVar = ChromVar.reindex(scATAC.obs_names)
        scATAC.obsm['ChromVar'] = ChromVar


        #fix cell barcode names.
        scATAC.obs_names = process_bcs_Spear_ATAC(
            scATAC.obs_names.to_series(),
            scATAC.obs['sample'].values,
            data_set_info['cell_line']
        )
        export_dir = data_set_info['export_dir']
        export_dir.mkdir(parents=True, exist_ok=True)
        adatas = separate_and_save_anndata(scATAC, export_dir, save = True)

        scATACs.append(scATAC)
    

    return scATACs


if __name__ == '__main__':
    run()

