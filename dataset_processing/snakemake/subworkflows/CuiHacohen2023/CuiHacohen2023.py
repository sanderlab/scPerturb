import pandas as pd
import scanpy as sc
import sys

from tqdm import tqdm
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 

# merge
adatas = []
for f in tqdm(snakemake.input):
    adata = sc.read(f)
    adatas.append(adata)
adata = sc.concat(adatas, axis=0)

# Obs
adata.obs.rename({
    'nCount_RNA': 'ncounts', 
    'nFeature_RNA': 'ngenes',
    'nCount_HTO': 'ncounts_tags',
    'biological_replicate_number': 'bio_replicate',
    'sample': 'perturbation',
    'rep': 'bio_replicate'
}, axis=1, inplace=True)
adata.obs.drop(['orig.ident'], axis=1, inplace=True)
adata.obs.perturbation = adata.obs.perturbation.astype(str)
adata.obs['perturbation'][adata.obs['perturbation']=='PBS'] = 'control'
adata.obs['nperts'] = [1-p.count('control') if type(p)==str else 0 for p in adata.obs.perturbation]
adata.obs['perturbation_type'] = 'cytokines'
adata.obs['disease'] = "healthy"
adata.obs['cancer'] = False
adata.obs['tissue_type']="primary"
adata.obs['organism'] = 'mouse'
annotate_qc(adata, species='mouse')
assert_annotations(adata)

# Add cytokine family annotations, Thanks Tümay for the list!
cytokine_families = {"IL-1": ["IL1a","IL1b", "IL1ra","IL18", "IL33", "IL36a", "IL36RA"],
                     "Common beta chain/IL13/TSLP": ["IL2", "IL4", "IL13", "IL7", "TSLP", "IL9", "IL15", "IL21"],
                     "Common beta chain": ["IL3", "IL5", "GM-CSF"],
                     "IL-6/IL-12": ["IL6", "IL11", "IL27", "IL30", "IL31", "LIF", "OSM", "Cardiotrophin-1", "Neuropoietin", "IL12", "IL23", "IL-Y"],
                     "IL-10": ["IL10", "IL19", "IL20", "IL22", "IL24"],
                     "IL-17": ["IL17A","IL17B", "IL17C", "IL17D","IL17E", "IL17F"],
                     "Interferon": ["IFNa1", "IFNb", "IFNe","IFNk", "IFNg","IFNl2"],
                     "TNF": ["LTA1-B2", "LTA2-B1", "TNFa", "OX40L", "CD40L", "FasL", "CD27L", "CD30L", "41BBL","TRAIL","RANKL", "TWEAK", "APRIL", "BAFF", "LIGHT", "TL1A", "GITRL"],
                     "Complement": ["C3a", "C5a"],
                     "Growth factor": ["Flt3l", "IL34", "M-CSF", "G-CSF", "SCF", "EGF", "VEGF", "FGF-basic", "HGF", "IGF-I"],
                     "Other": ["TGF-beta-1", "GDNF", "Persephin", "Prolactin", "Leptin", "Adiponectin","Resistin", "Noggin", "Decorin", "TPO"],
                     "Control": ["control"]
                    }
adata.obs['cytokine_family'] = 'None'
for key, values in cytokine_families.items():
    for value in values:
        adata.obs.loc[adata.obs.perturbation==value, 'cytokine_family'] = key

adata.write(snakemake.output[0], compression='gzip')
print('Done.')
