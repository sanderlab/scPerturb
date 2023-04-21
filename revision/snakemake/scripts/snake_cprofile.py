import scanpy as sc
print('Scanpy version:', sc.__version__)
import cProfile

import sys
# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import etest

# Setup CLI parser
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-i')      # option that takes a value
parser.add_argument('-o')      # option that takes a value
# access the arguments
args = parser.parse_args()

def exec_fun():
    adata = sc.read(args.i)
    etest(adata, obs_key='perturbation', control='control', runs=1000)

if __name__ == '__main__':
    print('starting profiling...')
    exec_fun()
