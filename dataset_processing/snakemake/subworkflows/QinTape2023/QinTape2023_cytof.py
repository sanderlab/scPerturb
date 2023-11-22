import scanpy as sc
import pandas as pd
import numpy as np
import re
from tqdm import tqdm
from pathlib import Path

# Custom functions
sys.path.insert(1, '../../')
from utils import annotate_qc, assert_annotations

TEMPDIR = Path(snakemake.config['TEMPDIR']) 
f
#TODO: Do this. See notebook.