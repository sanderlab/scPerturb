import scanpy as sc
import pandas as pd
import os
import numpy as np
import gzip
import tarfile
from scipy.io import mmread
from scipy.sparse import issparse, csr_matrix
from csv import Sniffer
import shutil
DIR = "/fast/scratch/users/peidlis_c/perturbation_resource_paper/"

def get_subfolders(d, full_path=True):
    prefix=d if full_path else ''
    return [os.path.join(prefix, o) for o in os.listdir(d) if os.path.isdir(os.path.join(d,o))]

def get_files(d, full_path=True):
    prefix=d if full_path else ''
    return [os.path.join(prefix, f) for f in os.listdir(d) if os.path.isfile(os.path.join(d, f))]

def force_merge(df1, df2):
    # pd.concat([adata.obs, tab], axis=0, ignore_index=True) is not working as one would think
    # see https://stackoverflow.com/questions/32801806/pandas-concat-ignore-index-doesnt-work
    df = pd.DataFrame(np.concatenate([df1, df2], axis=1),
                             index=df1.index,
                             columns=list(df1.columns) + list(df2.columns))
    df = df.loc[:,~df.columns.duplicated()]  # remove duplicate columns
    return df

def peek(f):
    opener = open if f.split('.')[-1]!='gz' else lambda x: gzip.open(x, 'rb')
    # peek into file to find out length and separator
    file_length = sum(1 for line in opener(f))

    for line in opener(f):
        try:
            first_line = line.decode()
        except (UnicodeDecodeError, AttributeError):
            first_line = line
        break
    sniffer = Sniffer()
    dialect = sniffer.sniff(first_line)
    separator=dialect.delimiter
    return file_length, separator

def process_GSE(gse_path):
    common_ = os.path.commonprefix(get_subfolders(gse_path))
    for gsm_path in get_subfolders(gse_path):
        gsm_name = gsm_path.split('/')[-1]
        files = get_files(gsm_path)
        is_mtx=['mtx' in f for f in files]
        is_count=['count' in f for f in files]
        if np.sum(is_mtx)==1:
            print('Found a single mtx file. Loading as unique count matrix...')
            count_file=files[np.where(is_mtx)[0][0]]
            mat = mmread(count_file)
        elif np.sum(is_count)==1:
            print('Found a count file. Loading as unique count matrix...')
            count_file=files[np.where(is_count)[0][0]]
            if 'mtx' in count_file:
                mat = mmread(count_file)
            else:
                mat = pd.read_csv(count_file, sep='\t')

        # initialize shape
        n_obs, n_var = mat.shape

        # load all other files
        annot_files=[f for f in files if f!=count_file]

        # file name prefixes
        common = os.path.commonprefix(annot_files)

        adata = sc.AnnData(X=mat)

        for f in annot_files:
            name = f[len(common):].split('.')[0]
            # peek into file to find out length and separator
            file_length, sep = peek(f)  # TODO if the file has no sep, we get rubbish. Maybe filter this.

            category = 'var' if abs(file_length-n_var)<2 else 'obs' if abs(file_length-n_obs)<2 else 'unknown'
            print(name, '\tinfered as\t', category)

            for f in annot_files:
                name = f[len(common):].split('.')[0]
                # peek into file to find out length and separator
                file_length, sep = peek(f)

                category = 'var' if abs(file_length-n_var)<2 else 'obs' if abs(file_length-n_obs)<2 else 'unknown'
                print(name, '\tinfered as\t', category)

                if category != 'unknown':
                    # todo catch non-text files
                    header = None if file_length==eval('n_'+category) else 'infer'
                    tab = pd.read_csv(f, header=header, sep='\t')
                    if tab.shape[1]==1 and header == None:
                        tab.columns = [name]
                    if category == 'obs':
                        adata.obs = force_merge(adata.obs, tab)
                    if category == 'var':
                        adata.var = force_merge(adata.var, tab)
                elif 'mtx' in f:
                    # clonal matrix for example
                    print('Reading additional mtx file...')
                    annot_mat=mmread(f)
                    # convert to CSR, COO not writeable to h5
                    annot_mat=csr_matrix(annot_mat) if issparse(annot_mat) else annot_mat
                    m0, m1=annot_mat.shape
                    n0, n1=mat.shape
                    if (n0==m0) and (n1==m1):
                        adata.layers[name] = annot_mat
                        print(f'Added {name} to layers.')
                    elif (n0==m1) and (n1==m0):
                        adata.layers[name] = annot_mat.T
                        print(f'Added {name} to layers.')
                    elif m0 == n0:
                        adata.obsm[name] = annot_mat
                        print(f'Added {name} to obsm.')
                    elif m1 == n0:
                        adata.obsm[name] = annot_mat.T
                        print(f'Added {name} to obsm.')
                    elif m1 == n1:
                        adata.varm[name] = annot_mat
                        print(f'Added {name} to varm.')
                    elif m0 == n1:
                        adata.varm[name] = annot_mat.T
                        print(f'Added {name} to varm.')
                    else:
                        print(f'Found another mtx file, but no clue what it is. Name is {name}')
                        print('it has shape ', annot_mat.shape, ' while adata has shape ', adata.shape)
                else:
                    print(f'Could not interpret file {name}')

        # todo detect if obs and var is flipped

        # set obs indices (cell barcodes)
        potential_barcode_cols = [c for c in adata.obs if 'barcode' in c]
        if len(potential_barcode_cols)==1:
            adata.obs.set_index(potential_barcode_cols[0], inplace=True)
            print('Set ',potential_barcode_cols[0],' as (cell) obs_names.')
        elif len(potential_barcode_cols)>1:
            adata.obs.set_index(potential_barcode_cols[0], inplace=True)
            print('Found multiple cell barcode candidates in obs. Using first one: ', potential_barcode_cols[0])
        else:
            print('Did not find any cell barcodes in obs.')

        # set var indices (gene names)
        potential_genename_cols = [c for c in adata.var if 'gene' in c]
        if len(potential_genename_cols)==1:
            adata.var.set_index(potential_genename_cols[0], inplace=True)
            print('Set ',potential_genename_cols[0],' as (gene) var_names.')
        else:
            print('Did not find unique gene names. Found ',len(potential_genename_cols),' candidates.')

        # convert to CSR, COO not writeable to h5
        adata.X = csr_matrix(adata.X) if issparse(adata.X) else adata.X
        adata.write(gse_path+gsm_name+'_processed.h5')

def process_GSE92872(DIR, WDIR):
    # DatlingerBock2017
    # read
    tab = pd.read_csv(DIR+"GSE92872/supp/GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz",
                     skiprows=6, header=None, index_col=0).T
    X = csr_matrix(tab.values, dtype=int)
    obs = pd.read_csv(DIR+"GSE92872/supp/GSE92872_CROP-seq_Jurkat_TCR.digital_expression.csv.gz",
                     header=None, index_col=0, nrows=5).T.set_index('cell')
    # create
    adata = sc.AnnData(X, obs=obs)
    adata.var_names = list(tab.columns)
    adata.obs_names_make_unique()
    # rename and reorder
    adata.obs = adata.obs.rename({'condition': 'perturbation_2', 'gene': 'target', 'grna': 'perturbation'}, axis=1)
    adata.obs = adata.obs[np.sort(adata.obs.keys())]
    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')
    # reform obs
    adata.obs.index=adata.obs.index.rename('cell_barcode')
    adata.obs['perturbation'] = ['control' if 'CTRL' in x else x for x in adata.obs['perturbation']]
    adata.obs['target'][adata.obs.target=='CTRL']=None
    adata.obs['celltype'] = 'T cells'
    adata.obs['cell_line'] = 'Jurkat cells'
    adata.obs['cancer'] = True
    adata.obs['disease'] = 'acute T cell leukemia'
    adata.obs['tissue_type'] = 'cell_line'
    adata.obs['organism'] = 'human'
    adata.obs['perturbation_type'] = 'CRISPR'
    adata.obs['perturbation_type_2'] = 'TCR stimulation'
    # write
    adata.write(WDIR+"DatlingerBock2017.h5")

def process_GSE90486(DIR, WDIR):
    path = DIR + 'GSE90486/supp/'
    adatas = {}
    for folder in get_subfolders(path, full_path=False):
        spl = folder.split('_')
        name = spl[1] + '_' + spl[2]
        file = pd.read_csv(path+folder+'/'+name+'.txt.gz', sep='\t').T
        adata = sc.AnnData(csr_matrix(file.values))
        adata.var_names = file.columns
        adata.obs_names = file.index
        adatas[name]=adata
    adata = sc.concat(adatas, index_unique='-', label='library')
    adata.write(WDIR+'GSE90486.h5')

def process_GSE90063(DIR, WDIR):
    path = DIR + 'GSE90063/supp/'
    def load_GSE90063_data(folder, write=False):
        print('Processing ', folder)
        files = get_files(path+folder, full_path=False)
        name = os.path.commonprefix(files)
        spath = path+'/'+folder+'/'+name

        X = csr_matrix(mmread(spath+'.mtx.txt.gz'))
        adata = sc.AnnData(X.T)

        var = pd.read_csv(spath+'_genenames.csv.gz', index_col=0)
        splitted = np.array([x.split('_', 1) for x in var.values[:,0]])
        adata.var_names = splitted[:,1]
        adata.var['gene_id'] = splitted[:,0]

        obs = pd.read_csv(spath+'_cellnames.csv.gz', index_col=0)
        splitted = np.array([x.split('_', 2) for x in obs.values[:,0]])
        adata.obs_names = obs.values[:,0]
        adata.obs['identifier_0'] = splitted[:,1]
        adata.obs['identifier_1'] = splitted[:,2]

        # annotation:
        files = [spath+'_cbc_gbc_dict.csv.gz', spath+'_cbc_gbc_dict_strict.csv.gz', spath+'_cbc_gbc_dict_lenient.csv.gz']
        keys = ['grna', 'grna_strict', 'grna_lenient']
        for file, key in zip(files, keys):
            if os.path.isfile(file):
                cbc_gbc_dict = pd.read_csv(file, index_col=0, header=None)
                adata.obs[key]=None
                for grna in list(cbc_gbc_dict.index):
                    for barcode in cbc_gbc_dict.loc[grna][1].replace(' ','').split(','):
                        if barcode in adata.obs_names:
                            val = adata.obs.loc[barcode][key]
                            adata.obs.loc[barcode][key] = grna if val is None else val+' + '+grna
        adata.var_names_make_unique()
        if write: adata.write(path+name+'.h5')
        return adata

    # K562 TFs
    adata = load_GSE90063_data('Supp_GSM2396858_K562_TFs__7_days')
    adata.obs['moi']='normal'
    adata.obs['days_past_transduction'] = 7
    bdata = load_GSE90063_data('Supp_GSM2396859_K562_TFs__13_days')
    bdata.obs['moi']='normal'
    bdata.obs['days_past_transduction'] = 13
    cdata = load_GSE90063_data('Supp_GSM2396860_K562_TFs__High_MOI')
    cdata.obs['moi']='high'
    cdata.obs['days_past_transduction'] = 7
    superdata = sc.concat({'K562_TFs__7_days': adata, 'K562_TFs__13_days': bdata, 'K562_TFs__High_MOI': cdata},
                          index_unique='-', label='library')
    superdata.write(WDIR+'GSE90063_K562_TFs.h5')

    # K562 cellcycle
    adata = load_GSE90063_data('Supp_GSM2396861_K562_cell_cycle', write=False)
    adata.write(WDIR+'GSE90063_GSM2396861_K562_cell_cycle_supp.h5')

    # BMDC at 0h and 3h
    adata = load_GSE90063_data('Supp_GSM2396857_DC_0hr')
    adata.obs['stimulated']=False
    adata.obs['time_hours'] = 0
    bdata = load_GSE90063_data('Supp_GSM2396856_DC_3hr')
    bdata.obs['stimulated']=True
    bdata.obs['time_hours'] = 3
    superdata = sc.concat({'GSM2396857_DC_0hr': adata, 'GSM2396856_DC_3hr': bdata},
                          index_unique='-', label='library')
    superdata.write(WDIR+'GSE90063_BMDC_stim.h5')

def process_GSE139944(DIR, WDIR):
    path = DIR + 'GSE139944/supp/'
    for folder in get_subfolders(path, False):
        print(folder)
        files = get_files(path+folder, False)
        pfile = [x for x in files if 'gene.annotations' in x][0]
        var = pd.read_csv(path+folder+'/'+pfile, sep='\t', header=None, names=['gene_id', 'gene_name']).set_index('gene_name')
        # obs = pd.read_csv(path+folder+'/'+folder[5:]+'_cell.annotations.txt.gz', sep='\t', header=None)
        pfile = [x for x in files if 'pData' in x][0]
        obs2 = pd.read_csv(path+folder+'/'+pfile, sep=' ')
        pfile = [x for x in files if 'UMI.count.matrix' in x][0]
        UMI_counts = pd.read_csv(path+folder+'/'+pfile, sep='\t', header=None)
        X=csr_matrix((UMI_counts[2], (UMI_counts[1]-1, UMI_counts[0]-1)), shape=(len(obs2), len(var)))
        adata = sc.AnnData(X, obs=obs2, var=var)
        adata.var_names_make_unique(join=':')
        adata.write(WDIR+'GSE139944_'+folder[5:]+'.h5')

def process_GSE119450(DIR, WDIR):
    adatas = {}
    folders = get_subfolders(DIR+'GSE119450/supp/')
    for d in ['D1', 'D2']:
        for cond in ['Stim', 'NoStim']:
            sample_name = d+'_'+cond
            gex_path=[x for x in folders if sample_name+'_10x' in x][0]
            obs_path=[x for x in folders if sample_name+'_ReAmp' in x][0]

            # extract GEX and read
            file = [f for f in get_files(gex_path) if '.gz' in f][0]
            tar = tarfile.open(file, "r:gz")
            tar.extractall(path=gex_path)
            tar.close()
            adata = sc.read_10x_mtx(gex_path)
            adata.obs.index = [x.split('-')[0] for x in adata.obs.index]

            # read and add metadata
            tab=pd.read_csv(get_files(obs_path)[0], index_col=0)
            tab=tab.rename({'UMI.count': 'gRNA_UMI_count', 'gRNA.ID': 'gRNA_ID'}, axis=1)[['gRNA_ID', 'gRNA_UMI_count']]
            adata.obs = pd.merge(adata.obs, tab, how='left', left_index=True, right_index=True)
            adata.obs['donor'] = d
            adata.obs['condition'] = cond
            adata.obs['sample_name'] = sample_name

            adata.write(DIR+'GSE119450/supp/GSE119450' + sample_name + '.h5')
            adatas[sample_name] = adata
    adata = sc.concat(adatas, index_unique='-')
    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')
    # reform obs
    adata.obs['cancer'] = False
    adata.obs.index=adata.obs.index.rename('cell_barcode')
    adata.obs=adata.obs.rename({'gRNA_ID': 'perturbation', 'donor': 'replicate', 'condition': 'perturbation_2', 'sample_name': 'library'}, axis=1)
    adata.obs['celltype'] = 'T cells'
    adata.obs['organism'] = 'human'
    adata.obs['disease'] = 'healthy'
    adata.obs['tissue_type']='primary'
    adata.obs['perturbation_type'] = 'CRISPR'
    adata.obs['perturbation_type_2'] = 'TCR stimulation'
    adata.obs.perturbation = adata.obs.perturbation.astype(str)
    adata.obs.perturbation[adata.obs.perturbation=='nan'] = 'control'
    adata.write(WDIR+'/ShifrutMarson2018.h5')

def process_GSE168620(DIR, WDIR):
    path = DIR + 'GSE168620/supp/Supp_GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells/'
    lib = 'GSM5151370_PD213_scifi_2_CRISPR-TCR_77300_MeOH-cells'  # I guess the only relevant one
    file = path+lib+'.h5ad'

    # load data
    if os.path.isfile(file+'.gz'):
        os.system(f"gunzip {file+'.gz'}")
    adata = sc.read(file)
    adata = adata[np.sum(adata.X, axis=1)>=100].copy()  # else too big

    # add metadata
    meta_file = path+lib+'.csv' if os.path.isfile(path+lib+'.csv') else path+lib+'.csv.gz'
    metadata = pd.read_csv(meta_file)
    obs = pd.merge(adata.obs, metadata, right_on='plate_well', left_on='plate_well')
    obs.index = adata.obs.index
    adata.obs = obs

    # convert gene names
    import mygene
    mg = mygene.MyGeneInfo()
    query = mg.querymany(list(adata.var_names) , scopes='ensembl.gene', fields='symbol', species='human', returnall=True, verbose=False)
    df = pd.DataFrame(query['out'])
    df = df.set_index('query')
    df = df.drop_duplicates(keep='first')
    var=pd.merge(adata.var, df, how='left', left_index=True, right_index=True)
    var.loc[pd.isna(var['symbol']),'symbol'] = list(var.index[pd.isna(var['symbol'])])
    var=var[['symbol']]
    var = var.drop_duplicates(keep='first')
    adata.var['ensembl_id'] = adata.var.index
    adata.var_names=[var.loc[x]['symbol'] if x in var.index else x for x in adata.var_names]
    adata.var_names = [str(x) for x in adata.var.index.values]  # make strings, else not writeable
    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')
    # reform obs
    adata.obs = adata.obs.drop(['combinatorial_barcode', 'sample_name', 'material', 'gRNA_seq', 'gRNA_well'], axis=1)
    adata.obs.index=adata.obs.index.rename('cell_barcode')
    adata.obs=adata.obs.rename({'gRNA_ID': 'perturbation', 'TCR_status': 'perturbation_2', 'plate_well': 'sample'}, axis=1)
    adata.obs['tissue_type']='cell_line'
    adata.obs['cancer'] = True
    adata.obs['celltype'] = 'T cells'
    adata.obs['cell_line'] = 'Jurkat cells'
    adata.obs['disease'] = 'acute T cell leukemia'
    adata.obs['organism'] = 'human'
    adata.obs['perturbation_type'] = 'CRISPR'
    adata.obs['perturbation_type_2'] = 'TCR stimulation'
    # save file
    adata.write(WDIR+'DatlingerBock2021_scifi2.h5')

def process_GSE122662(DIR, WDIR):
    # WIP: only two subseries added
    path = DIR + 'GSE122662/supp/'

    # GSE106340
    folders = get_subfolders(path, False)
    folders = [x for x in folders if int(x[8:15]) in np.arange(2836267, 2836288+1)]
    adatas = {}
    for folder in folders:
        # create adata
        files = get_files(path+folder)
        var = pd.read_csv([x for x in files if 'genes' in x][0], sep='\t', index_col=1, header=None, names=['gene_id', 'gene_name'])
        obs = pd.read_csv([x for x in files if 'barcodes' in x][0], sep='\t', index_col=0, header=None)
        try:
            X = csr_matrix(mmread([x for x in files if 'matrix.mtx' in x][0])).T
        except:
            print(f'No count matrix found for {folder}. Is this corrupted?')
            continue
        adata = sc.AnnData(X, obs, var)

        # annotate
        GSM = folder.split('_', 2)[1]
        adata.obs['GSM'] = GSM
        string = folder.split('_', 2)[2]
        import re
        adata.obs['condition'] = 'serum' if 'serum' in string else '2i' if '2i' in string else None
        adata.obs['time'] = 'iPSCs' if 'iPSCs' in string else re.search('D\d+', string).group(0)

        # clear indices
        adata.obs_names = [x.split('-')[0] for x in adata.obs_names]
        adata.var_names_make_unique()
        adatas[GSM] = adata
    # merge
    superdata = sc.concat(adatas, index_unique='-')
    superdata.write(WDIR+'SchiebingerLander2019_GSE106340.h5')

    # GSE115943
    folders = get_subfolders(path, False)
    folders = [x for x in folders if int(x[8:15]) in np.arange(3195648, 3195773+1)]
    adatas = {}
    for folder in folders:
        file = get_files(path+folder, True)[0]
        # load data
        adata = sc.read_10x_h5(file)
        # annotate
        adata.obs['time'] = 'iPSC' if 'iPSC' in folder else re.search('D[\d\.]+(_5)?', folder).group(0).replace('_', '.')
        adata.obs['replicate'] = re.search('C\d', folder).group(0)[-1]
        adata.obs['condition'] = folder.split('_')[-2]
        GSM = folder.split('_', 2)[1]
        adata.obs['GSM'] = GSM
        # clear indices
        adata.obs_names = [x.split('-')[0] for x in adata.obs_names]
        adata.var_names_make_unique()
        adatas[GSM] = adata
    # merge
    superdata = sc.concat(adatas, index_unique='-')
    superdata.write(WDIR+'SchiebingerLander2019_GSE115943.h5')

def process_GSE149383(DIR, WDIR):
    path = DIR + 'GSE149383/supp/'

    # format: ['GEO', 'days', 'origin', 'condition', 'bundle', 'subseries']

    X=[
    # GSE134836 (Erl PC9) is weird, they fucked this up during upload.
    # ['GSM3972651', 0, 'PC9', 'control', 'PC9D0vsD3Erl', 'GSE134836'],  #PC9D0
    # ['GSM3972652', 3, 'PC9', 'Erl', 'PC9D0vsD3Erl', 'GSE134836'],  #PC9D3Erl
    # GSE134838 (Vem M14)
    ['GSM3972655', 0, 'M14', 'control', 'M14D0vsD3Vem', 'GSE134838'], #M14Day0
    ['GSM3972656', 3, 'M14', 'Vem', 'M14D0vsD3Vem', 'GSE134838'], #M14Day3_Vem
    # GSE134839 (Erl PC9)
    ['GSM3972657', 0, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 0 PC9 erlotinib
    ['GSM3972658', 1, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 1 PC9 erlotinib
    ['GSM3972659', 2, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 2 PC9 erlotinib
    ['GSM3972660', 4, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 4 PC9 erlotinib
    ['GSM3972661', 9, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 9 PC9 erlotinib
    ['GSM3972662', 11, 'PC9', 'Erl', 'D0_D1_D2_D4_D9_D11', 'GSE134839'], #Day 11 PC9 erlotinib
    # GSE134841 (Drug holiday run Erl PC9)
    ['GSM3972669', 0, 'PC9', 'Erl', 'Drug_holiday', 'GSE134841'], #PC9 Day 0
    ['GSM3972670', 2, 'PC9', 'Erl', 'Drug_holiday', 'GSE134841'], #PC9 Day 2
    ['GSM3972671', 11, 'PC9', 'Erl', 'Drug_holiday', 'GSE134841'], #PC9 Day 11
    ['GSM3972672', 19, 'PC9', 'Erl_DMSO', 'Drug_holiday', 'GSE134841'], #PC9 D19 DMSO
    ['GSM3972673', 19, 'PC9', 'Erl_Erl', 'Drug_holiday', 'GSE134841'], #PC9 D19 ERL
    # GSE134842 (Beads on mix)
    ['GSM3972674', 0, 'PC9+U937', 'control', 'bead_mix', 'GSE134842'], #PC9 U937 mix
    ['GSM3972675', 0, 'PC9+U937', 'EpCAM+', 'bead_mix', 'GSE134842'], #EpCAM-positive by beads
    ['GSM3972676', 0, 'PC9+U937', 'CD45-', 'bead_mix', 'GSE134842'], #CD45-negative by beads
    # GSE148466 (primary tissues)
    ['GSM4472055', 0, 'P8126999_primary', 'control', 'primary', 'GSE148466'], #P8126999
    ['GSM4472056', 0, 'P8127005_primary', 'control', 'primary', 'GSE148466'], #P8127005
    ['GSM4472057', 0, 'P8127011_primary', 'control', 'primary', 'GSE148466'], #P8127011
    # GSE149214 (PC9 time course Mono-drug)
    ['GSM4494347', 0, 'PC9_rep1', 'control', 'multirep_PC9_mono', 'GSE149214'], #Day00b scRNA-seq
    ['GSM4494348', 11, 'PC9_rep1', 'Erl', 'multirep_PC9_mono', 'GSE149214'], #Day11b scRNA-seq
    ['GSM4494349', 11, 'PC9_rep2', 'Erl', 'multirep_PC9_mono', 'GSE149214'], #Day11c scRNA-seq
    # GSE149215 (PC9 time course Multi-drug)
    ['GSM4494350', 3, 'PC9_rep1', 'Erl', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-ERL1 scRNA-seq
    ['GSM4494351', 3, 'PC9_rep2', 'Erl', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-ERL2 scRNA-seq
    ['GSM4494352', 3, 'PC9_rep1', 'Erl+Cri', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-CRI-ERL1 scRNA-seq
    ['GSM4494353', 3, 'PC9_rep2', 'Erl+Cri', 'multirep_PC9_multi', 'GSE149215'], #PC9D3-CRI-ERL2 scRNA-seq
    ['GSM4494354', 3, 'PC9', 'Eto', 'multirep_PC9_multi', 'GSE149215'], #PC9Day3-Eto scRNA-seq
    # GSE160244 (xenografts mice)
    ['GSM4869650', 3, 'PC9_xeno', 'control', 'Xenograft', 'GSE160244'], #Day 3 control
    ['GSM4869651', 3, 'PC9_xeno', 'Cri', 'Xenograft', 'GSE160244'], #Day 3 crizotinib
    ['GSM4869652', 3, 'PC9_xeno', 'Osi', 'Xenograft', 'GSE160244'], #Day 3 osimertinib
    ['GSM4869653', 3, 'PC9_xeno', 'Osi+Cri', 'Xenograft', 'GSE160244'] #Day 3 osimertinib+crizotinib
    ]
    GSMs = [x.split('_')[1] for x in get_subfolders(path, False)]
    tab = pd.DataFrame(X, columns=['GEO', 'days', 'origin', 'condition', 'bundle', 'subseries'])

    folders = get_subfolders(path)
    for subseries in pd.unique(tab.subseries):
        adatas = {}
        for index, row in tab[tab.subseries==subseries].iterrows():
            # determine file for GEO of subseries
            file = get_files([x for x in folders if row.GEO in x][0], True)[0]

            # load
            Y = pd.read_csv(file, sep='\t', index_col=0)
            adata = sc.AnnData(csr_matrix(Y.values.T, dtype=int))
            adata.var_names = Y.index
            adata.obs_names = Y.columns
            adata.var_names_make_unique()

            # annotate by table
            for col in tab.columns:
                adata.obs[col] = row[col]
            adatas[row.GEO] = adata
        superdata = sc.concat(adatas, index_unique='-')
        superdata.write(WDIR+'GSE149383_'+subseries+'.h5')

def process_GSE150949(DIR, WDIR):
    path = DIR + 'GSE150949/supp/'

    X = pd.read_csv(path+'GSE150949_pooled_watermelon.count.matrix.csv.gz', index_col=0).T
    obs = pd.read_csv(path+'GSE150949_pooled_watermelon.metadata.matrix.csv.gz', index_col=0)
    adata = sc.AnnData(csr_matrix(X.values, dtype=int), obs=obs)
    adata.var_names = X.columns
    # # 14 mins, I guess this is just normalized values???
    # tab = pd.read_csv(path+'GSE150949_pooled_watermelon.data.matrix.csv.gz')
    # np.max(tab.values) # 7.8710119999999995
    adata.write(WDIR+'GSE150949_pooled_watermelon_supp.h5')

    # 12 mins
    df = pd.read_csv(path+'GSE150949_pc9_count_matrix.csv.gz', index_col=0)
    adata = sc.AnnData(csr_matrix(df.values.T, dtype=int))
    adata.var_names = df.index
    adata.obs_names = df.columns
    adata.write(WDIR+'GSE150949_pc9n_supp.h5')

def process_GSE148842(DIR, WDIR):
    path = DIR + 'GSE148842/supp/'
    folders = get_subfolders(path, True)
    annotation = pd.read_excel('../tables/GSE148842_manual_annotation.xlsx')
    adatas = {}
    for folder in folders:
        file = get_files(folder, False)[0]
        tab = pd.read_csv(folder + '/' + file, sep='\t')
        var = tab[tab.columns[:2]].set_index(tab.columns[1], drop=True)
        obs_names = tab.columns[2:]
        X = tab.drop(tab.columns[:2], axis=1).values.T
        from scipy.sparse import csr_matrix
        X = csr_matrix(X)

        adata = sc.AnnData(X, var=var)
        adata.obs_names = obs_names
        adata.var_names_make_unique('_')

        ann = annotation[annotation['GEO'] == file.split('_')[0]]
        patient = ann['Sample'].values[0][:5]
        adata.obs['patient'] = patient
        for c in ann.columns:
            adata.obs[c] = ann[c].values[0]
        adatas[file.split('.')[0].split('_')[-1]] = adata
    adata = sc.concat(adatas, index_unique='-', label='library')

    # reform var
    adata.var.index = adata.var.index.rename('gene_symbol')

    # reform obs
    adata.obs['dose_value']=[None if ('DMSO' in x or x=='none') else x.split(' ')[0] for x in adata.obs.treatment]
    adata.obs['dose_unit']=[None if ('DMSO' in x or x=='none') else x.split(' ')[1] for x in adata.obs.treatment]
    adata.obs['perturbation']=['control' if ('DMSO' in x or x=='none') else x.split(' ')[-1] for x in adata.obs.treatment]
    adata.obs = adata.obs.drop('treatment', axis=1)
    adata.obs.index=adata.obs.index.rename('cell_barcode')
    adata.obs=adata.obs.rename({'patient': 'sample', 'gender': 'sex'}, axis=1)
    adata.obs.perturbation[adata.obs.perturbation=='vehicle (DMSO)']='control'
    adata.obs['tissue_type']='primary'
    adata.obs['cancer'] = True
    adata.obs['disease'] = 'glioblastoma'
    adata.obs['celltype'] = 'unknown'
    adata.obs['sex'] = adata.obs['sex'].replace('M', 'm').replace('F', 'f')
    adata.obs['organism'] = 'human'
    adata.obs['perturbation_type'] = 'drug'

    adata.write(WDIR+'ZhaoSims2021.h5')

# TODO check if clonal matrix is saved
def process_GSE140802(DIR, WDIR):
    path = f'{DIR}GSE140802/supp/'
    folders = get_subfolders(path)
    # only "GSM4185644_stateFate_cytokinePerturbation" has perturbations
    for folder in folders:
        name = folder.split('/Supp_')[-1]
        if name == 'GSM4185644_stateFate_cytokinePerturbation':
            print('Processing ', folder)
            files = get_files(folder)

            # count matrix
            mat = csr_matrix(mmread([f for f in files if '_normed_counts.mtx' in f][0]))
            # gene names
            var = pd.read_csv([f for f in files if '_gene_names' in f][0], index_col=0, header=None, sep='\t', names=['gene_symbols'])
            # barcodes and metadata
            obs = pd.read_csv([f for f in files if '_cell_barcodes' in f][0], index_col=0, header=None, sep='\t', names=['cell_barcodes'])
            obs1 = pd.read_csv([f for f in files if '_library_names' in f][0], index_col=0, header=None, sep='\t', names=['library_names'])
            obs['library_names'] = list(obs1.index)
            metafiles = [f for f in files if '_metadata' in f]
            if len(metafiles) > 0:
                metafile = metafiles[0]
                obs2 = pd.read_csv(metafile, index_col='Cell barcode', sep='\t')
                obs = pd.concat([obs, obs2], axis=1)
            if ('library_names' in obs.columns) and ('Library' in obs.columns):
                obs = obs.drop('Library', axis=1)

            adata = sc.AnnData(mat, obs, var)

            # add clonal info
            clones = [f for f in files if '_clone_matrix.mtx' in f]
            if len(clones) > 0:
                cmat = csr_matrix(mmread(clones[0])).T
                adata.obsm['clone_matrix'] = cmat

            # reform var
            adata.var.index = adata.var.index.rename('gene_symbol')
            # reform obs
            adata.obs['cancer'] = False
            adata.obs.index=adata.obs.index.rename('cell_barcode')
            adata.obs=adata.obs.rename({'Time point': 'age', 'Cytokine condition': 'perturbation', 'Cell type annotation': 'celltype'}, axis=1)
            adata.obs.age = ['day '+str(x) for x in adata.obs.age]
            adata.obs = adata.obs.drop(['SPRING-x', 'SPRING-y'], axis=1)
            adata.obs['organism'] = 'mouse'
            adata.obs['disease'] = 'healthy'
            adata.obs['tissue_type']='primary'
            adata.obs['perturbation_type'] = 'cytokine'
            adata.write(WDIR+'WeinrebKlein2020.h5')

def process_GSE135497(DIR, WDIR):
    path = DIR + 'GSE135497/supp/'
    folders = get_subfolders(path, False)
    for experiment in ['TAP_SCREEN__chromosome_11_screen', 'TAP_SCREEN__chromosome_8_screen']:
        subfolders = [x for x in folders if experiment in x]
        adatas = {}
        bdatas = {}
        for folder in subfolders:
            # load data
            files = get_files(path+folder)
            X = pd.read_csv([f for f in files if '.counts.' in f][0], index_col=0).T
            adata = sc.AnnData(csr_matrix(X.values, dtype=int))
            adata.var_names = X.columns
            adata.obs_names = X.index
            # add info what was perturbed
            pert_status = pd.read_csv([f for f in files if '.pertStatus.' in f][0], index_col=0).T
            bdata = sc.AnnData(csr_matrix(pert_status.values, dtype=int))
            bdata.var_names = pert_status.columns
            bdata.obs_names = pert_status.index
            # adata.obsm['pert_status'] = csr_matrix(pert_status)
            # adata.uns['pert_status_vars'] = pert_status.columns  # Problem: this is not saved...
            # add with sample name
            adatas[folder.split('__')[-1]]=adata
            bdatas[folder.split('__')[-1]]=bdata
        adata = sc.concat(adatas, index_unique='-', label='sample')
        bdata = sc.concat(bdatas, index_unique='-', label='sample')

        # reform var
        adata.var.index = adata.var.index.rename('gene_symbol')
        # reform obs
        adata.obs.index=adata.obs.index.rename('cell_barcode')
        adata.obs=adata.obs.rename({'sample': 'replicate'}, axis=1)
        adata.obs['tissue_type']='cell_line'
        adata.obs['cell_line']='K562'
        adata.obs['cancer'] = True
        adata.obs['disease'] = 'chronic myelogenous leukemia'
        adata.obs['celltype'] = 'lymphoblasts'
        adata.obs['organism'] = 'human'
        guides_per_cell = np.array(np.sum(bdata.X, axis=1)).T[0]
        best_guide = list(bdata.var_names[np.array(np.argmax(bdata.X, axis=1)).T[0]])
        adata.obs['perturbation'] = best_guide
        adata.obs['perturbation'][guides_per_cell>1]='multiplet'
        adata.obs['perturbation'][guides_per_cell==0]='control'
        adata.obs['perturbation_type'] = 'CRISPR'
        adata.write(DIR+'SchraivogelSteinmetz2020_'+experiment+'.h5')
        # bdata.write(DIR+'GSE135497_'+experiment+'_perturbationinfo.h5')

def process_GSE90546(DIR, WDIR):
    path = DIR + 'GSE90546/supp/'
    shutil.unpack_archive(path+"GSE90546_RAW.tar", path)
    for experiment in ['GSM2406675_10X001', 'GSM2406677_10X005', 'GSM2406681_10X010']:
        X = csr_matrix(mmread(path+experiment+"_matrix.mtx.txt.gz"))
        var = pd.read_csv(path+experiment+"_genes.tsv.gz", sep='\t', header=None, names=['gene_id', 'gene_name']).set_index('gene_name')
        obs = pd.read_csv(path+experiment+"_barcodes.tsv.gz", sep='\t', header=None, index_col=0)
        obs.index.names = ['cell_barcode']
        adata = sc.AnnData(X.T, var=var, obs=obs)
        adata.obs.index = [x.split('-')[0] for x in adata.obs.index]
        adata.var_names_make_unique()
        metadata = pd.read_csv(path+experiment+"_cell_identities.csv.gz", index_col=0)
        metadata.index = [x.split('-')[0] for x in metadata.index]
        adata.obs_names_make_unique()
        metadata=metadata[~metadata.index.duplicated()]
        adata.obs=pd.concat([adata.obs, metadata], axis=1, join='outer').loc[adata.obs_names]
        if 'good coverage' in adata.obs.columns:
            adata.obs['good coverage'] = np.array(adata.obs['good coverage'], dtype=str)
        # reform var
        adata.var.index = adata.var.index.rename('gene_symbol')
        adata.var = adata.var.rename({'gene_id': 'ensembl_id'}, axis=1)
        # reform obs
        adata.obs.index=adata.obs.index.rename('cell_barcode')
        adata.obs=adata.obs.rename({'guide identity': 'perturbation'}, axis=1)
        adata.obs = adata.obs.drop(['coverage', 'good coverage', 'number of cells'], axis=1)
        adata.obs['tissue_type']='cell_line'
        adata.obs['cell_line']='K562'
        adata.obs['cancer'] = True
        adata.obs['disease'] = 'chronic myelogenous leukemia'
        adata.obs['perturbation_type'] = 'CRISPR'
        adata.obs['celltype'] = 'lymphoblasts'
        adata.obs['organism'] = 'human'
        adata.write(f'{WDIR}AdamsonWeissman2016_{experiment}.h5')
