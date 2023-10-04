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

def process_GehringPachter2019(DIR, WDIR):
    # note : I just downloaded the h5
    # site has two h5 files of the same size with the same name, I just got one of them
    # there are fastq available as well
    #https://data.caltech.edu/records/1311
    # process data using Jupyter Notebook at ClickTagDemuliplex20191009
    # https://colab.research.google.com/github/pachterlab/GPCTP_2019/blob/master/Colab%20Notebooks/96SampleMuliplexExperiment.ipynb
    # this will perform relevant download and process into clicktag matrices
    # add to their script/run at end: cDNAdata.write_h5ad("metadata_custom/cDNA_all.h5ad")
    # expID.to_csv("metadata_custom/perturbation_info.csv")
    perts = pd.read_csv('GehringPachter2019/metadata_custom/perturbation_info.csv', index_col=0)
    adata = sc.read_h5ad('GehringPachter2019/metadata_custom/cDNA_all.h5ad')
    # subset to cells with unambiguous assignments
    adata = adata[adata.obs.index.isin(perts.index),].copy()

    adata.obs = pd.concat([adata.obs,perts], axis=1)
    adata.obs_names = adata.obs_names.str.slice(0,16)

    adata.obs['disease'] = 'healthy'
    adata.obs['cancer'] = False
    adata.obs['tissue_type'] = 'primary'
    adata.obs['celltype'] = 'neural stem cells'
    adata.obs['perturbation'] = 'BMP4'
    adata.obs['perturbation_2'] = 'EGF and bFGF'
    adata.obs['perturbation_3'] = '1:5 Scriptaid:decitabine'
    adata.obs['perturbation_4'] = 'retinoic acid'
    adata.obs['dose_unit'] = 'ng/mL'
    adata.obs['dose_unit_2'] = 'ng/mL'
    adata.obs['dose_unit_3'] = 'uM'
    adata.obs['dose_unit_4'] = 'uM'

    adata.obs['organism'] = 'mouse'
    adata.obs['perturbation_type'] = 'drug'

    adata.obs['dose_value'] = 0
    adata.obs['dose_value_2'] = 0
    adata.obs['dose_value_3'] = 0
    adata.obs['dose_value_4'] = 0

    adata.obs['nperts'] = 0

    for i in range(0,adata.obs.shape[0]):
        vec = adata.obs.iloc[i,:].copy()
        if vec['BMP'] ==1:
            vec['dose_value'] = 8
            vec['nperts'] +=1
        elif vec['BMP'] == 2:
            vec['dose_value'] = 40
            vec['nperts'] +=1
        elif vec['BMP']==3:
            vec['dose_value'] == 200
            vec['nperts'] +=1
                
        if vec['EGF'] ==0:
            vec['dose_value_2'] = 1.6
            vec['nperts'] +=1
        elif vec['EGF'] ==1:
            vec['dose_value_2'] = 8
            vec['nperts'] +=1
        elif vec['EGF'] ==2:
            vec['dose_value_2'] = 40
            vec['nperts'] +=1
        elif vec['EGF'] ==3:
            vec['dose_value_2'] = 200
            vec['nperts'] +=1
            
        if vec['ScripDec'] ==0:
            vec['dose_value_3'] = 0
            vec['dose_value_4'] = 0 
        elif vec['ScripDec'] ==1:
            vec['nperts'] +=1
            if vec['RA']==0:
                vec['dose_value_3'] = 0.2
            else:
                vec['dose_value_4'] = 2
            
        elif vec['ScripDec'] ==2:
            vec['nperts'] +=1
            if vec['RA'] == 0:
                vec['dose_value_3'] = 1
            else:
                vec['dose_value_4'] = 10

        adata.obs.iloc[i,:] = vec
        #break
    adata.obs = adata.obs.drop(columns=['n_counts','n_countslog','BMP','EGF','ScripDec','RA','PlatePos'])
    annot = sc.queries.biomart_annotations(
        "mmusculus",
        ["mgi_symbol","ensembl_gene_id"])#.set_index("ensembl_gene_id")
    annot = annot.set_index('ensembl_gene_id')
    adata.var['ensembl_gene_id'] = adata.var.index
    adata.var.index = adata.var.index.str.split('.').str[0]

    adata.var = adata.var.join(annot)
    adata.var.index = adata.var['ensembl_gene_id']
    adata.var = pd.DataFrame(adata.var['mgi_symbol'])

    adata.var = adata.var.rename(columns={"mgi_symbol":"gene_symbol"})
    adata.var['ensembl_id'] = adata.var.index
    adata.var['gene_symbol'] = adata.var['gene_symbol'].fillna(adata.var['ensembl_id'])
    adata.var = adata.var.set_index('gene_symbol')
    adata.var['mt'] = adata.var_names.str.startswith('mt-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.str.startswith(("Rps","Rpl"))
    qc = sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=False)           

    adata.obs['ncounts'] = qc[0]['total_counts']
    adata.obs['ngenes'] = qc[0]['n_genes_by_counts']
    adata.obs['percent_mito'] = qc[0]['pct_counts_mt']
    adata.obs['percent_ribo'] = qc[0]['pct_counts_ribo']
    adata.var['ncounts'] = qc[1]['total_counts']
    adata.var['ncells'] = qc[1]['n_cells_by_counts']

    adata.var = adata.var[['ensembl_id','ncounts','ncells']]
#annot.loc[adata.var.index.str.split('.').str[0]]
    adata.obs_names_make_unique()
    adata.var_names_make_unique()
    utils.write_as_singles(adata, WDIR, "GehringPachter2019", add_h5=True)