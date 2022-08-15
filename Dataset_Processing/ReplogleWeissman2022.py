import pandas as pd
import scanpy as sc
import numpy as np 
from scipy import sparse
import utils
import mygene

import os

import collections
import tables
import itertools

from collections import defaultdict

import scipy.stats as stats


def gene_symbols_to_ensembl(gene_list, species='human', verbose=False):
    mg = mygene.MyGeneInfo()
    out = mg.querymany(gene_list, scopes='symbol', fields='ensembl.gene', species=species, verbose=verbose)
    df = pd.DataFrame([[o['query'], o['_id']] if '_id' in o.keys() else [o['query'], None] for o in out], columns=['gene_symbol', 'ensembl.gene']).set_index('gene_symbol')
    return df

def process_ReplogleWeissman2022(DIR= 'ReplogleWeissman2022/', WDIR="processed/"):
    def reformat_adata(adata, celltype="K562"):
        if celltype == "K562":
            adata.obs['disease'] = "chronic myeloid leukemia"
            adata.obs['cancer'] = True
            adata.obs["cell_line"] = "K562"
            adata.obs["sex"] = "Female"
            adata.obs["age"]=53
            


        else:
            adata.obs['disease'] = "healthy"
            adata.obs['cancer'] = False
            adata.obs['tissue_type']="cell_line"
            adata.obs["cell_line"] = "hTERT-RPE1"
            adata.obs["sex"]="Female"
            adata.obs["age"] = 1
            adata.obs["celltype"]

        adata.obs["perturbation"] = adata.obs['gene']
        adata.obs["organism"] = "human"
        adata.obs["perturbation_type"] = "CRISPR"
        adata.obs['tissue_type']="cell_line"
        
        adata.obs.rename(columns={"mitopercent": "percent_mito", "sgID_AB":"guide_id", "gem_group":"batch"},inplace=True)
        adata.obs['ncounts']= adata.X.sum(1)
        adata.obs['ngenes'] = np.count_nonzero(adata.X, axis =1)
        adata.obs["nperts"] = 1
        adata.obs.loc[adata.obs['gene'] == "non-targeting","nperts"] = 0 
        # percent_ribo
        ribomat = adata[:,adata.var["gene_name"].str.startswith(("RPS","RPL"))]
        ribomat = ribomat.X.sum(1)
        adata.obs['percent_ribo'] = np.divide(ribomat, adata.obs['ncounts'])
        # fix duplicates TBCE, HSPA14 

        adata.var['gene_name'] = adata.var['gene_name'].astype("object")
        adata.var.loc[adata.var['gene_name'] == "TBCE",'gene_name']  = 'TBCE_'+adata.var.loc[adata.var['gene_name'] == "TBCE",:].index.values
        adata.var.loc[adata.var['gene_name'] == "HSPA14",'gene_name']  = 'HSPA14_'+adata.var.loc[adata.var['gene_name'] == "HSPA14",:].index.values

        adata.var['ensembl_id'] = adata.var.index
        adata.var.set_index('gene_name', inplace=True)

        # ncounts 
        adata.var['ncounts'] = adata.X.sum(0)
        # ncells
        adata.var['ncells'] = np.count_nonzero(adata.X, axis =0)
        return adata.copy()

    adata = sc.read_h5ad(os.path.join(DIR,"K562_essential_raw_singlecell_01.h5ad"))
    adata = reformat_adata(adata)
    utils.write_as_singles(adata, WDIR, "ReplogleWeissman2022/K562_essential", add_h5=True)

    adata = sc.read_h5ad(os.path.join(DIR,"K562_gwps_raw_singlecell_01.h5ad"))
    adata = reformat_adata(adata)
    utils.write_as_singles(adata,WDIR,"ReplogleWeissman2022/K562_gwps", add_h5=True)

    adata = sc.read_h5ad(os.path.join(DIR,"rpe1_raw_singlecell_01.h5ad"))
    adata = reformat_adata(adata)
    utils.write_as_singles(adata, WDIR, "ReplogleWeissman2022/rpe1", add_h5=True)