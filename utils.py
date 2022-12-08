import gzip
import os
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl

from warnings import warn
from itertools import product
from tqdm import tqdm
from scipy.io import mmwrite
from scipy.sparse import issparse
from scipy.stats import zscore
from sklearn.metrics import pairwise_distances
from scipy.cluster.hierarchy import distance, linkage, dendrogram
from shutil import copyfileobj
from chembl_webresource_client.new_client import new_client


def random_split(adata, key='perturbation', inplace=True):
    # Creates a new column "{key}_X" where 50% of cells randomly get a "_X" appended to their value of {key}.
    adata = adata if inplace else adata.copy()
    N = int(np.min(adata.obs[key].value_counts())/2)
    adata.obs[f'{key}_X'] = adata.obs[key].astype(str)
    for group in pd.unique(adata.obs[key]):
        mask = np.random.choice(adata.obs_names[adata.obs[key]==group], N, replace=False)
        adata.obs.loc[mask, f'{key}_X'] = f'{group}_X'
    if not inplace:
        return adata


def query_chembl_drugs(names, qkeys=['pref_name', 'molecule_chembl_id', 'molecule_type', 'usan_stem_definition']):
    '''
    Returns results (pandas.DataFrame), amb_results (pandas.DataFrame), not_found (list).
    '''
    drug_client = new_client.molecule
    names = list(names) if isinstance(names, pd.Index) else [names] if isinstance(names, str) else names
    results = []
    amb_results = []
    not_found = []
    for drug in tqdm(names, leave=False):
        name = drug.split(' ')[0].split('.')[0]
        res = drug_client.filter(pref_name__istartswith=name).only(
            ['pref_name', 'molecule_chembl_id', 'molecule_type', 'usan_stem_definition'])

        if len(res)==0 and '(' in drug:
            name = drug.split('(')[-1].split(')')[0]
            if len(name)>1:
                res = drug_client.filter(pref_name__istartswith=name).only(
                        ['pref_name', 'molecule_chembl_id', 'molecule_type', 'usan_stem_definition'])

        if len(res)==0:
            # not found by chembl
            not_found.append(drug)
            continue

        df = pd.DataFrame(res)
        if drug.upper() in df.pref_name.values:
            # perfect match
            df = df.iloc[df.pref_name.values==drug.upper()]
        elif name.upper() in df.pref_name.values:
            # initial match
            df = df.iloc[df.pref_name.values==name.upper()]
        elif len(df.index)>1 and not isinstance(df, pd.Series):
            # ambiguous match
            df.index=[drug]*len(df)
            amb_results.append(df)
            continue
        df.index=[drug]*len(df)
        results.append(df)
    return results, amb_results, not_found


def write_as_singles(adata, path, name, add_h5=False, h5_compression='gzip'):
    """
    Write AnnData object to as folder with single file format.

    This function creates a folder "name" at "path". Then it saves
    the adata object as single files, that is:
    adata.X --> counts.mtx.gz
    adata.var --> var.csv
    adata.obs --> obs.csv
    adata.obsm --> obsm/
    adata.varm --> varm/
    Note that adata.uns is not recoverable.
    If add_h5==True, is also saves the object as an h5 file in the
    same folder as "name".h5.

    Parameters
    ----------
    adata : AnnData (Annotated data matrix) object
        Data object to write.
    path : str
        Path in which folder with data is written to.
    name : str
        Name of the folder to be created.

    Returns
    -------
    None
    """

    def create_folder(folder):
        # create directory if necessary
        if not os.path.isdir(folder):
            os.mkdir(folder)
        else:
            warn(f'{folder} already exists. Possibly overwriting files...')

    def gzip_file(file, keep_original=False):
        # gzip a file
        with open(file, 'rb') as file_in:
            with gzip.open(f'{file}.gz', 'wb') as file_out:
                copyfileobj(file_in, file_out)
        if not keep_original:
            os.remove(file)

    # create directory if necessary
    folder = f'{path}/{name}'
    create_folder(folder)

    # write annotations
    # adata.write_csvs(folder)  # this produces emptpy obsm and varm files...
    adata.obs.to_csv(f'{folder}/obs.csv')
    adata.var.to_csv(f'{folder}/var.csv')

    # write mtx
    if not issparse(adata.X):
        warn(f'Count matrix is not sparse.')
    mmwrite(f'{folder}/counts.mtx', adata.X)
    # gzip mtx
    gzip_file(f'{folder}/counts.mtx')

    # write obsm
    obsm_folder = f'{folder}/obsm/'
    create_folder(obsm_folder)
    for key in adata.obsm.keys():
        np.savetxt(f'{obsm_folder}{key}.csv.gz', adata.obsm[key], delimiter=',')

    # write varm
    varm_folder = f'{folder}/varm/'
    create_folder(varm_folder)
    for key in adata.varm.keys():
        np.savetxt(f'{varm_folder}{key}.csv.gz', adata.varm[key], delimiter=',')

    # write h5 if needed
    if add_h5:
        adata.write(f'{folder}/{name}.h5ad', compression=h5_compression)


def read_from_singles(folder, use_h5=True):
    """
    Read AnnData object from folder with single file format.

    This function reads single file format from folder. It creates an
    AnnData object like this:
    adata.X <-- counts.mtx.gz
    adata.var <-- var.csv
    adata.obs <-- obs.csv
    adata.obsm <-- obsm/
    adata.varm <-- varm/
    Note that adata.uns is not recoverable.

    Parameters
    ----------
    adata : AnnData (Annotated data matrix) object
        Data object to write.
    folder : str
        Path to folder with data in single file format.

    Returns
    -------
    AnnData
        An AnnData object reconstructed from single files in folder.
    """
    if use_h5:
        h5_files = [f for f in os.listdir(f'{folder}') if '.h5' in f]
        if len(h5_files)==1:
            return sc.read(f'{folder}/{h5_files[0]}')
        else:
            warn('Found multiple h5 files. Reading from singles instead...')

    # read from singles
    adata = sc.read_mtx(f'{folder}/counts.mtx.gz')
    adata.obs = pd.read_csv(f'{folder}/obs.csv', index_col=0)
    adata.var = pd.read_csv(f'{folder}/var.csv', index_col=0)
    if os.path.isdir(f'{folder}/obsm/'):
        for obsm_file in [f for f in os.listdir(f'{folder}/obsm/') if os.path.isfile(os.path.join(f'{folder}/obsm/', f))]:
            key = obsm_file.replace('.gz', '').replace('.csv', '')
            adata.obsm[key] = np.loadtxt(f'{folder}/obsm/{obsm_file}', delimiter=',')
    if os.path.isdir(f'{folder}/varm/'):
        for varm_file in [f for f in os.listdir(f'{folder}/varm/') if os.path.isfile(os.path.join(f'{folder}/varm/', f))]:
            key = varm_file.replace('.gz', '').replace('.csv', '')
            adata.varm[key] = np.loadtxt(f'{folder}/varm/{varm_file}', delimiter=',')
    return adata

def rsum(X, axis):
    # handles sparse sum, returns array instead of matrix object
    return np.ravel(np.sum(X, axis=axis))

def detect_organism(adata):
    if np.sum(adata.var_names.str.startswith('MT-')) > 0:
        return 'human'
    elif np.sum(adata.var_names.str.startswith('mt-')) > 0:
        return 'mouse'
    else:
        print('Could not reliably detect organism from mito genes. Setting as human.')
        return 'human'

def specify_genes(genes, source_species='mouse', target_species='human'):
    genes = genes if isinstance(genes, list) else list(genes) if (isinstance(genes, np.ndarray) or isinstance(genes, pd.Index)) else [genes]
    # load human-mouse orthology
    tab = pd.read_csv(os.path.dirname(os.path.realpath(__file__))+'/metadata/HMD_HumanPhenotype.rpt', sep='\t',
                  names=['human', 'human_gene_number (HGNC)', 'mouse', 'mouse_gene_number (MGI)', 'mouse protein'], usecols=[0,1,2,3,4])
    human_mouse_orthology_dict = dict(zip(tab.human, tab.mouse))
    mouse_human_orthology_dict = {key: value for (value, key) in human_mouse_orthology_dict.items()}
    # translate genes
    # TODO the mapping could lead to duplicate genes, is that a problem?
    if source_species=='mouse' and target_species=='human':
        return [mouse_human_orthology_dict[x] if x in mouse_human_orthology_dict.keys() else x.upper() for x in genes]
    elif source_species=='human' and target_species=='mouse':
        return [human_mouse_orthology_dict[x] if x in human_mouse_orthology_dict.keys() else x.capitalize() for x in genes]
    elif target_species == 'mouse':
        return [x.capitalize() for x in genes]
    elif target_species == 'human':
        return [x.upper() for x in genes]
    else:
        raise ValueError('Species '+target_species+' not known.')

def get_genefamily_percentage(adata, key='MT-', start=True, name='mito'):
    keys = key if isinstance(key, list) else [key, '____ignore____']
    if start:
        family_genes = np.logical_or(*[adata.var_names.str.startswith(k) for k in keys])
    else:
        family_genes = np.logical_or(*[adata.var_names.str.endswith(k) for k in keys])
    if issparse(adata.X):
        adata.obs['percent_'+name] = np.sum(
            adata[:, family_genes].X, axis=1).A1 * 100 / np.sum(adata.X, axis=1).A1
    else:
        adata.obs['percent_'+name] = np.sum(
            adata[:, family_genes].X, axis=1) * 100 / np.sum(adata.X, axis=1)

def get_mito_percentage(adata, species='human'):
    key = 'MT-' if species == 'human' else 'mt-'
    get_genefamily_percentage(adata, key=key, start=True, name='mito')

def get_ribo_percentage(adata, species='human'):
    key = specify_genes(['RPS', 'RPL'], target_species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='ribo')

def get_hemo_percentage(adata, species='human'):
    key = specify_genes(['HBA', 'HBB'], target_species=species)
    get_genefamily_percentage(adata, key=key, start=True, name='hemo')

def annotate_qc(adata, species='detect'):
    """
    Annotates quality controls and such to existing adata (inplace).

    This function adds (if they do not yet exist):
    adata.obs:
        - ncounts
        - ngenes
        - nperts
        - organism (inferred from `species`)
        - percent_mito
        - percent_ribo
        - percent_hemo (if tissue_type exists in obs and is `'primary'`)
    adata.obs.index:
        - updates name of index to be `'cell_barcode'`
    adata.var:
        - ncounts
        - ncells

    Parameters
    ----------
    adata : AnnData (Annotated data matrix) object
        Data object to write.
    species : `str` from `['detect', 'human', 'mouse']` (default: `'detect'`)
        Species for gene family percentage calling. If detect, heuristically
        tries to infer if data is mouse or human based on presence of
        either "MT-" or "mt-" genes.

    Returns
    -------
    None
    """

    # qc counts
    if 'ncounts' not in adata.obs.keys():
        adata.obs['ncounts'] = rsum(adata.X, axis=1)
    if 'ngenes' not in adata.obs.keys():
        adata.obs['ngenes'] = rsum(adata.X>0, axis=1)
    if 'ncounts' not in adata.var.keys():
        adata.var['ncounts'] = rsum(adata.X, axis=0)
    if 'ncells' not in adata.var.keys():
        adata.var['ncells'] = rsum(adata.X>0, axis=0)

    species = detect_organism(adata) if species == 'detect' else species
    if 'organism' not in adata.obs.keys():
        adata.obs['organism'] = species
    # gene modules
    if 'percent_mito' not in adata.obs.keys():
        get_mito_percentage(adata, species)
    if 'percent_ribo' not in adata.obs.keys():
        get_ribo_percentage(adata, species)
    if adata.obs['tissue_type'][0]=='tissue' and 'percent_hemo' not in adata.obs.keys():
        get_hemo_percentage(adata, species)

    # rename index
    adata.obs.index.names = ['cell_barcode']

    # perturbation stuff
    if 'nperts' not in adata.obs.keys():
        adata.obs['nperts'] = [len(x.split('_')) if isinstance(x, str) else 0 for x in adata.obs.perturbation]

def equal_subsampling_old(adata, obs_key, N_min=None, remove_low=False):
    '''
   set remove_low TRUE and N_min = 100? 1000? e.g.
   classes smaller than that get thrown out
   need this for CRISPR datasets especially

    Subsample to same class sizes. Classes given by obs_key pointing to categorical in adata.obs.
    If N_min is given, downsamples to maximum this number instead of the number of cells in the smallest class.
    This can lead to some classes having less cells than others if they are below N_min cells in that class (which won't be subsampled then).
    If remove_low is True, any class that has less than N_min cells will be discarded instead.
    '''
    N = np.min(adata.obs[obs_key].value_counts())
    N = N if N_min==None else np.max([N_min, N])
    selection = np.hstack(np.array([np.random.choice(adata.obs_names[adata.obs[obs_key]==p],
                         size=np.min([N, np.sum(adata.obs[obs_key]==p)-1]), replace=False) for p in pd.unique(adata.obs[obs_key]) if np.sum(adata.obs[obs_key]==p)>=N]))
    return adata[selection].copy()

def assert_annotations(adata):
    assert adata.obs.index.name == 'cell_barcode', "Please set your index name correctly using 'adata.obs.index.names = ['cell_barcode']'"
    for obs_key in ['disease', 'cancer', 'tissue_type', 'perturbation', 'organism', 'perturbation_type', 'ncounts', 'ngenes', 'nperts', 'percent_mito', 'percent_ribo']:
        assert obs_key in adata.obs.keys(), f"There is no column `{obs_key}` in adata.obs. Please annotate it (worst case: set it all to None)!"
    if (adata.obs.tissue_type[0] == 'cell_line') or (adata.obs.tissue_type[0] == 'organoid'):
        assert 'cell_line' in adata.obs.keys(), f"There is no column `cell_line` in adata.obs, even though this is a cell culture model. Please annotate it (worst case: set it all to None)!"
    for var_key in ['ncounts', 'ncells']:
        assert var_key in adata.var.keys(), f"There is no column `{var_key}` in adata.var. Please annotate it (worst case: set it all to None)!"

def equal_subsampling(adata, obs_key, N_min=None):
    '''
    Subsample to same class sizes. Classes given by obs_key pointing to categorical in adata.obs.
    If N_min is given, downsamples to at least this number instead of the number of cells in the smallest class
    and throws out classes with less than N_min cells.
    '''
    counts = adata.obs[obs_key].value_counts()
    groups = counts.index[counts>=N_min]  # ignore groups with less than N_min cells to begin with
    # We select downsampling target counts by min-max, i.e.
    # the largest N such that every group has at least N cells before downsampling.
    N = np.min(counts)
    N = N if N_min==None else np.max([N_min, N])
    # subsample indices per group
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=False) for group in groups]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()

def cluster_matrix(matrix, how='row', return_order=False, method='centroid'):
    '''
    Hierarchical clustering of a matrix/dataframe. `how` can be 'col', 'row' or 'both' (default: 'row').
    '''
    if how not in ['col', 'row', 'both']:
        raise ValueError('Value for "how" must be row or col.')
    if how!='both':
        M = matrix if how=='row' else matrix.T
        dist = distance.pdist(M)
        link = linkage(dist, method=method)
        dend = dendrogram(link, no_plot=True)
        order = np.array(dend['leaves'], dtype=int)
        if return_order:
            return order
        elif isinstance(matrix, pd.DataFrame):
            return matrix.iloc[order] if how=='row' else matrix.iloc[:, order]
        else:
            return matrix[order] if how=='row' else matrix[:, order]
    else:
        if return_order:
            warn('Returning order when clustering both row and col is not supported.')
        matrix_ = cluster_matrix(matrix, how='row', return_order=False, method=method)
        return cluster_matrix(matrix_, how='col', return_order=False, method=method)

def pairwise_pca_distances(adata, obs_key, obsm_key='X_pca', dist='sqeuclidean'):
    groups = pd.unique(adata.obs[obs_key])
    df = pd.DataFrame(index=groups, columns=groups, dtype=float)
    for i, p1 in enumerate(tqdm(groups)):
        x1 = adata[adata.obs[obs_key]==p1].obsm[obsm_key].copy()
        N = len(x1)
        for p2 in groups[i:]:
            x2 = adata[adata.obs[obs_key]==p2].obsm[obsm_key].copy()
            pwd = pairwise_distances(x1, x2, metric=dist)
            M = len(x2)
            factor = N*M if p1!=p2 else N**2 - N  # correct mean for zero diagonal if comparing to same set
            mean_pwd = np.sum(pwd) / factor
            df.loc[p1, p2] = mean_pwd
            df.loc[p2, p1] = mean_pwd
    return df

def pseudo_bulk(adata, keys, layer='counts', min_cells_per_group=10):
    X = []
    Y = []
    for gs in tqdm(product(*[pd.unique(adata.obs[key]) for key in keys])):
        mask = np.logical_and.reduce([adata.obs[key]==g for g, key in zip(gs, keys)])
        ncells = sum(mask)
        if ncells < min_cells_per_group: continue
        Y.append(list(gs)+[ncells])
        X_ = adata[mask].layers[layer] if layer!=None else adata[mask].X
        X.append(np.array(np.sum(X_, axis=0), dtype=int)[0])
    obs=pd.DataFrame(Y, columns=list(keys)+['ncells'])
    return sc.AnnData(np.array(X), obs=obs, var=adata.var)

def pairwise_mean_pca_distances(adata, obs_key, obsm_key='X_pca', sq_dist=True):
    groups = pd.unique(adata.obs[obs_key])
    df = pd.DataFrame(index=groups, columns=groups, dtype=float)
    for i, p1 in enumerate(tqdm(groups)):
        x1 = np.mean(adata[adata.obs[obs_key]==p1].obsm[obsm_key], axis=1)
        for p2 in groups[i:]:
            x2 = np.mean(adata[adata.obs[obs_key]==p2].obsm[obsm_key], axis=1)
            pwd = np.linalg.norm(x1 - x2) ** (1+int(sq_dist))
            df.loc[p1, p2] = pwd
            df.loc[p2, p1] = pwd
    return df
