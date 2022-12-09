import pandas as pd
import numpy as np
import scanpy as sc

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
    if N_min is not None:
        groups = counts.index[counts>=N_min]  # ignore groups with less than N_min cells to begin with
    else:
        groups=counts.index
    # We select downsampling target counts by min-max, i.e.
    # the largest N such that every group has at least N cells before downsampling.
    N = np.min(counts)
    N = N if N_min==None else np.max([N_min, N])
    # subsample indices per group
    indices = [np.random.choice(adata.obs_names[adata.obs[obs_key]==group], size=N, replace=False) for group in groups]
    selection = np.hstack(np.array(indices))
    return adata[selection].copy()