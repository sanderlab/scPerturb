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

