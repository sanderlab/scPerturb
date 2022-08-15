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

def process_XieHon2017(DIR= 'XieHon2017/', WDIR="processed/"):
    GEODIR = DIR+'GSE81884_RAW/'
    SUPPDIR = DIR+'paper_supplement/'

    # functions modified from https://github.com/russellxie/Global-analysis-K562-enhancers/blob/master/VirtualFacs/sgrnaprocessing/_correct_sgrna.py
    # note need to unzip the barcode files, NOT the counts files
    def load_data_barcodes(dr=GEODIR, filename='GSM2544716_K562_dCas9-KRAB_Lenti-sgRNA_set_A_pooled_infection_biorep_3_batch_1'):
        input_fh  = open(dr+filename+'.cell_barcode_to_sgRNA_barcode.txt', 'r')
        next(input_fh)
        cell_bc_list   = []
        num_sgRNA_list = np.array([])
        sgRNAs         = []
        #umis           = []
        data_dict = nested_dict(2, list)
        for line in input_fh:
            cell_bc    = line.strip().split('\t')[0]
            sgRNA_list = line.strip().split('\t')[1].split(',')
            num_sgRNA  = len(sgRNA_list)
            umi_list   = line.strip().split('\t')[2].split(',')
            #print(cell_bc)
            #print(umi_list)
            #umi_list = list(map(int,umi_list))
            for i in zip(sgRNA_list,umi_list):
                if i[1] != 'NA':  
                    #print(type(i[1]))
                    data_dict[cell_bc][i[0]] = i[1]
            df = pd.DataFrame(data_dict).fillna(0)
        return df.astype('int')
    def nested_dict(n, type):
        '''
        Create an N-dimentional nested dictionary
        '''

        if n == 1:
            return defaultdict(type)
        else:
            return defaultdict(lambda: nested_dict(n-1, type))


    #filter the sgRNA UMI count based on the cutoff values
    def filter_umi (df):
        #calculate the cutoff value for each sgRNA in the dataset
        sgRNA_cutoff = [[turn_point(i, df)] for i in list(df.index)]
        sgRNA_cutoff = list(np.array(sgRNA_cutoff).flatten())
        for i in range(0, len(sgRNA_cutoff)):
            df.iloc[i, :].loc[df.iloc[i, :] <= sgRNA_cutoff[i]] = 0
        return df


    #Function to define the turning point of a cumulative curve by using the minimum derivative method
    def turn_point(sgRNA_name, df):
        sgRNA_count  = df.T.filter(items=[sgRNA_name]).sum(axis=1).sort_values(ascending=False)
        
        sgRNA_cumsum = sgRNA_count.cumsum()

        #get the total cell number of this sgRNA
        #cell_num = np.argwhere(sgRNA_count > 0).size
        cell_num = (sgRNA_count > 0).sum()

        #calculate the turning point by using the max derivative
        turning_point = sgRNA_cumsum.loc[((sgRNA_cumsum.diff()) / sgRNA_count.sum() > (1/cell_num))].shape
        
        return(sgRNA_count.iloc[turning_point])

    def grab_files(dr = "XieHon2017/GSE81884_RAW/", barcodes=SUPPDIR+"mmc2.xlsx",
               filename='GSM2175064_K562_dCas9-KRAB_Lenti-sgRNA_set_B_pooled_infection_batch_5'):
        #print(filename)
        codes = load_data_barcodes(dr=dr, filename=filename)
        codes = filter_umi(codes)
        
        counts = pd.read_csv(dr+filename+'.counts.txt.gz', sep = '\t')
        #codes = pd.read_csv(dr+filename+'.cell_barcode_to_sgRNA_barcode.txt.gz', sep = '\t')
        
        Xmat = counts.drop(columns=["Chr","Start","End","Strand","Length"])
        Xmat['Geneid'] = Xmat['Geneid'].str.split(';', n=1).str.get(0)
        Xmat = Xmat.groupby(['Geneid']).agg('sum')
        Xmat.columns = Xmat.columns.str.split('.', n =1).str.get(0)
        
        adata = sc.AnnData(Xmat.T)
        adata.obs['disease'] = "chronic myeloid leukemia"
        adata.obs['cancer'] = True
        adata.obs['sex'] = 'f'
        adata.obs['age'] = '53'
        adata.obs['tissue_type']= 'cell_line'
        adata.obs['cell_line'] = 'K562'
        # maybe add cell type of origin
        adata.obs['organism'] = 'human'
        adata.obs['perturbation_type'] = 'CRISPR'
        # treat set as sample
        adata.obs['sample'] = filename.split("set_",1)[1][0:1] 
        # grab batch from filename 
        adata.obs['batch'] = filename.split("batch_",1)[1]
        adata.obs['target'] = 'control'
        adata.obs['guide_id'] = 'control'
        adata.obs['perturbation'] = 'control'
        adata.obs['nperts'] = 0 
        
        

        if "biorep" in filename:
            adata.obs['replicate'] = filename.split("biorep_")[1][0]
        else: 
            adata.obs['replicate'] = '1'
        # perturbaion from other mamt
        barcode_info = pd.read_excel("XieHon2017/paper_supplement/mmc2.xlsx")
        barcode_info = barcode_info[barcode_info['batch'].str.contains(filename.split("set_",1)[1][0:1] )]

        
        barcode_info = pd.read_excel("XieHon2017/paper_supplement/mmc2.xlsx")
        xtra = barcode_info[barcode_info['batch']=="A-C"]
        barcode_info = barcode_info[barcode_info['batch']== filename.split("set_",1)[1][0:1]]
        barcode_info = barcode_info.append(xtra)
        
        # strip cells without barcode information
        adata = adata[adata.obs_names.isin(codes.columns.values),:]  
        adata.obsm['barcodes'] = codes.T.reindex(adata.obs_names)
        
        for cellindex, cell in enumerate(adata.obs_names):
            # # TODO : MODIFY THIS SO THE GUIDES ARE IN DECREASING ORDER BY COUNTS
            # # add guides and guide names into .obsm
            sgRNA = codes.loc[codes[cell]>0, cell].sort_values(ascending=False).index
            sgRNA_id = ''
            target = ''
            nperts = 0 
            if len(sgRNA) != 0:
                for index,bc in enumerate(sgRNA):
                    nperts += 1
                    ### TODO FIX THIS TO GO IN DECREASING ORDER OF COUNTS
                    if len(sgRNA_id) == 0:
                        sgRNA_id = barcode_info[barcode_info['sgRNA_barcode_seq'].str.contains(bc)]['sgRNA_id'].values[0].replace(";","_")
                        target = barcode_info[barcode_info['sgRNA_barcode_seq'].str.contains(bc)]['target'].values[0].replace("_",";")
                    else:
                        sgRNA_id += ';'+barcode_info[barcode_info['sgRNA_barcode_seq'].str.contains(bc)]['sgRNA_id'].values[0].replace(";","_")
                        target += '_'+barcode_info[barcode_info['sgRNA_barcode_seq'].str.contains(bc)]['target'].values[0].replace("_",";")
                adata.obs.loc[cell,'target'] = target
                adata.obs.loc[cell,'guide_id'] = sgRNA_id
                adata.obs.loc[cell,'nperts'] = nperts
                adata.obs.loc[cell,'perturbation'] = target

        #for cell in adata.obs_names:
        #    sgrna = codes.index[codes[cell] >0]
        
        return(adata)
    filelist = os.listdir(GEODIR)
    filelist = [ x for x in filelist if "counts" in x ]
    filelist = [x.split(".counts.txt.gz")[0] for x in filelist]

    adata = grab_files(dr= GEODIR, barcodes= SUPPDIR + "mmc2.xlsx", filename=filelist[0])
    for file in filelist[1:]:
        adata = adata.concatenate(grab_files(dr= GEODIR, barcodes= SUPPDIR + "mmc2.xlsx", filename=file), batch_key = "filename",index_unique=None, join="outer")

    adata.obs_names_make_unique()
    adata.obsm['barcodes'] = adata.obsm['barcodes'].fillna(0)
    annot = sc.queries.biomart_annotations(
        "hsapiens",
        ["hgnc_symbol","ensembl_gene_id","chromosome_name"])#.set_index("ensembl_gene_id")
    annot['hgnc_symbol'] = annot['hgnc_symbol'].fillna(annot['ensembl_gene_id'])
    annot = annot.set_index("hgnc_symbol")

    new_var = adata.var.copy()
    new_var['ensembl_gene_id'] = None
    for i in range(0,len(new_var)):
        gene_name = new_var.index[i]
        #updated_gene_name = gene_names_updated.loc[gene_names_updated["Input"]==gene_name,"Approved symbol"].values[0]
        if gene_name in annot.index:
            if isinstance(annot.loc[gene_name]['ensembl_gene_id'], pd.Series):
                new_var.loc[gene_name,"ensembl_gene_id"] = annot.loc[gene_name]['ensembl_gene_id'][0]
            else:
                new_var.loc[gene_name,"ensembl_gene_id"] = annot.loc[gene_name]['ensembl_gene_id']


    adata.var['mt'] = adata.var_names.str.startswith('MT-')  # annotate the group of mitochondrial genes as 'mt'
    adata.var['ribo'] = adata.var_names.str.startswith(("RPS","RPL"))

    qc = sc.pp.calculate_qc_metrics(adata, qc_vars=['mt','ribo'], percent_top=None, log1p=False, inplace=False)           
    adata.var = new_var

    adata.obs['ncounts'] = qc[0]['total_counts']
    adata.obs['ngenes'] = qc[0]['n_genes_by_counts']
    adata.obs['percent_mito'] = qc[0]['pct_counts_mt']
    adata.obs['percent_ribo'] = qc[0]['pct_counts_ribo']
    adata.var['ncounts'] = qc[1]['total_counts']
    adata.var['ncells'] = qc[1]['n_cells_by_counts']
    adata.var.index.rename("gene_symbol")
    adata.obs['nperts'] = [str(x) for x in adata.obs['nperts']]

    utils.write_as_singles(adata, WDIR, "XieHon2017", add_h5=True)