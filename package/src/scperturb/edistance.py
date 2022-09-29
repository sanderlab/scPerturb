import pandas as pd
import numpy as np
import scanpy as sc

from tqdm import tqdm
from sklearn.metrics import pairwise_distances
from statsmodels.stats.multitest import multipletests

def pairwise_pca_distances(adata, obs_key, obsm_key='X_pca', dist='sqeuclidean', verbose=True):
    """Average of pairwise PCA distances between cells of each group in obs_key.
    For each pair of groups defined in adata.obs[obs_key] (e.g. perturbations)
    computes all pairwise distances between cells in adata.obsm[obsm_key] (e.g. PCA space)
    and averages them per group-pair. This results in a distance matrix between all groups.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys()
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.

    Returns
    -------
    pwd: pandas.DataFrame
        DataFrame with pairwise PCA distances between all groups.
    """

    if obs_key=='X_pca':
        print('PCA embedding not found, computing...')
        sc.pp.pca(adata)

    groups = pd.unique(adata.obs[obs_key])
    df = pd.DataFrame(index=groups, columns=groups, dtype=float)
    fct = tqdm if verbose else lambda x: x
    for i, p1 in enumerate(fct(groups)):
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

def edist(adata, obs_key='perturbation', obsm_key='X_pca', pwd=None, dist='sqeuclidean', verbose=True):
    """Computes the edistance to control. Accepts precomputed pwd.
    Computes the pairwise E-distances between all groups of cells defined in
    adata.obs[obs_key] (e.g. perturbations). Distances are computed in embedding
    space given by adata.obsm[obsm_key] (e.g. PCA space).

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.

    Returns
    -------
    estats: pandas.DataFrame
        DataFrame with pairwise E-distances between all groups.
    """

    pwd = pairwise_pca_distances(adata, obs_key=obs_key, obsm_key=obsm_key, dist=dist, verbose=verbose) if pwd is None else pwd
    # derive basic statistics
    sigmas = np.diag(pwd)
    deltas = pwd
    estats = 2 * deltas - sigmas - sigmas[:, np.newaxis]
    return estats

def equal_subsampling(adata, obs_key, N_min=None):
    """Subsample cells while retaining same class sizes.
    Classes are given by obs_key pointing to categorical in adata.obs.
    If N_min is given, downsamples to at least this number instead of the number
    of cells in the smallest class and throws out classes with less than N_min cells.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    N_min: `int` or `None` (default: `None`)
        If N_min is given, downsamples to at least this number instead of the number
        of cells in the smallest class and throws out classes with less than N_min cells.

    Returns
    -------
    subdata: :class:`~anndata.AnnData`
        Subsampled version of the original annotated data matrix.
    """

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

def etest(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean', control='control', alpha=0.05, runs=100, verbose=True):
    """Performs Monte Carlo permutation test with E-distance as test statistic.
    Tests for each group of cells defined in adata.obs[obs_key] if it is significantly
    different from control based on the E-distance in adata.obsm[obsm_key] space.
    Does multiple-testing correction using Holm-Sidak.

    Arguments
    ---------
    adata: :class:`~anndata.AnnData`
        Annotated data matrix.
    obs_key: `str` in adata.obs.keys() (default: `perturbation`)
        Key in adata.obs specifying the groups to consider.
    obsm_key: `str` in adata.obsm (default: `adata.obsm['X_pca']`)
        Key for embedding coordinates to use.
    dist: `str` for any distance in scipy.spatial.distance (default: `sqeuclidean`)
        Distance metric to use in embedding space.
    control: `str` (default: `'control'`)
        Defines the control group in adata.obs[obs_key] to test against.
    alpha: `float` between 0 and 1 (defaul: 0.05)
        significance cut-off for the test to annotate significance.
    runs:
    verbose: `bool` (default: `True`)
        Whether to show a progress bar iterating over all groups.

    Returns
    -------
    tab: pandas.DataFrame
        E-test results for each group in adata.obs[obs_key] with columns
        - edist: E-distance to control
        - pvalue: E-test p-value if group is different from control
        - significant: If p-value < alpha
        - pvalue_adj: Multiple-testing corrected E-test p-value
        - significant_adj: If p-value_adj < alpha
    """

    # Approximate sampling from null distribution (equal distributions)
    res = []
    groups = pd.unique(adata.obs[obs_key])
    fct = tqdm if verbose else lambda x: x
    for i in fct(range(runs)):
        # per perturbation, shuffle with control and compute e-distance
        df = pd.DataFrame(index=groups, columns=['edist'], dtype=float)
        for group in groups:
            if group==control:
                df.loc[group] = [0]
                continue
            mask = adata.obs[obs_key].isin([group, control])
            subdata = adata[mask].copy()
            subdata.obs['shuffled'] = np.random.permutation(subdata.obs[obs_key].values)
            df_ = edist(subdata, 'shuffled', obsm_key=obsm_key, dist=dist, verbose=False).loc[control]
            df.loc[group] = df_.loc[group]
        res.append(df.sort_index())

    # "Sampling" from original distribution without shuffling (hypothesis)
    df = edist(adata, obs_key, obsm_key=obsm_key, dist=dist, verbose=False).loc[control]
    df = pd.DataFrame(df.sort_index())
    df.columns = ['edist']

    # Evaluate test (hypothesis vs null hypothesis)
    results = pd.concat([r['edist'] - df['edist'] for r in res], axis=1) > 0  # count times shuffling resulted in larger e-distance
    pvalues = np.sum(results, axis=1) / runs

    # Apply multiple testing correction
    significant_adj, pvalue_adj, _, _ = multipletests(pvalues, alpha=alpha, method='holm-sidak')

    # Aggregate results
    tab = pd.concat([df['edist'], pvalues], axis=1)
    tab.columns = ['edist', 'pvalue']
    tab['significant'] = tab.pvalue < alpha
    tab['pvalue_adj'] = pvalue_adj
    tab['significant_adj'] = significant_adj
    return tab
