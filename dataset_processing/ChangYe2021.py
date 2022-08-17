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

def process_ChangYe2021(DIR, WDIR):
    #unzip everything prior to running
    # current script keeps all metadata in original -- will probably want to cut this down to just ones we care about
    # ie don't need tsne components
    def grab_file_pair(name):
        mat1 = pd.read_csv(DIR+"/E-MTAB-10698.processed.1/"+name + ".count.data.tab", sep = "\t").transpose()
        obs = pd.read_csv(DIR+"/E-MTAB-10698.processed.2/" + name + ".meta.data.tab", sep = "\t")
        X = sparse.csr_matrix(mat1.values, dtype = int)
        var_names = mat1.columns
        adata = sc.AnnData(X, obs=obs)
        adata.var_names = list(var_names)
        adata.obs_names_make_unique()
        return adata
    
    #m1 = grab_file_pair("NGS1758_SAM2ebc6d44f6")
    #m2 = grab_file_pair("NGS1758_SAM6b79946c41")
    m3 = grab_file_pair("NGS2986")
   

    #need to doublecheck paper but I *think* we only care about the PC9 ones
    # set of perturbation experiments (1 none, three with drugs)
    manifest = pd.read_csv(DIR+"/E-MTAB-10698.sdrf.txt", sep = "\t")

    new_obs = pd.DataFrame(index = m3.obs['barcode'])

    new_obs = pd.DataFrame(index = m3.obs['barcode']+m3.obs['sample'], data = m3.obs['sample'].values)
    new_obs = new_obs.rename(columns={0:'sample'})
    new_obs['disease'] = 'lung adenocarcinoma'
    new_obs['cancer'] = True
    new_obs['sex'] = 'm' # googled PC9 cell line
    new_obs['perturbation'] = 'control'
    new_obs['dose_value'] = '0'
    new_obs['dose_unit'] = 'micromolar'
    new_obs['perturbation_type'] = 'drug'
    new_obs['organism'] = 'human'
    new_obs['nperts'] = 0 

    m3.obs = new_obs 
    m3.obs_names_make_unique()

    for i in range(2,6):
        if manifest['Factor Value[compound]'][i] == 'none':
            m3.obs.loc[manifest['Characteristics[sample_id]'][i] == m3.obs['sample'],('perturbation','dose_value', 'nperts')] = 'control', '0', '0'
        else:
            m3.obs.loc[manifest['Characteristics[sample_id]'][i] == m3.obs['sample'],('perturbation','dose_value', 'nperts')] = manifest['Factor Value[compound]'][i], manifest['Factor Value[dose]'][i], '1'
    # update hgnc symbols to latest version
    gene_names_updated = pd.read_csv("ChangYe2021/hgnc-symbol-check.csv", header=1)
    # fill NA with original gene name
    gene_names_updated['Approved symbol'] = gene_names_updated['Approved symbol'].fillna(gene_names_updated['Input'])

    annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["hgnc_symbol","ensembl_gene_id","chromosome_name"])#.set_index("ensembl_gene_id")
    annot['hgnc_symbol'] = annot['hgnc_symbol'].fillna(annot['ensembl_gene_id'])
    annot = annot.set_index("hgnc_symbol")

    # add ensembl id
    new_var = m3.var.copy()
    new_var['ensembl_gene_id'] = None
    for i in range(0,len(new_var)):
        gene_name = new_var.index[i]
        updated_gene_name = gene_names_updated.loc[gene_names_updated["Input"]==gene_name,"Approved symbol"].values[0]
        if updated_gene_name in annot.index:
            if isinstance(annot.loc[updated_gene_name]['ensembl_gene_id'], pd.Series):
                new_var.loc[gene_name,"ensembl_gene_id"] = annot.loc[updated_gene_name]['ensembl_gene_id'][0]
            else:
                new_var.loc[gene_name,"ensembl_gene_id"] = annot.loc[updated_gene_name]['ensembl_gene_id']

    
    m3.var['mt'] = m3.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    m3.var['ribo'] = m3.var_names.str.startswith(("RPS","RPL"))

    qc = sc.pp.calculate_qc_metrics(m3, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=False)           

    m3.var = new_var

    m3.obs['ncounts'] = qc[0]['total_counts']
    m3.obs['ngenes'] = qc[0]['n_genes_by_counts']
    m3.obs['percent_mito'] = qc[0]['pct_counts_mt']
    m3.obs['percent_ribo'] = qc[0]['pct_counts_ribo']
    m3.var['ncounts'] = qc[1]['total_counts']
    m3.var['ncells'] = qc[1]['n_cells_by_counts']
    m3.var.index.rename("gene_symbol")


    utils.write_as_singles(m3, WDIR, "ChangYe2021", add_h5=True)