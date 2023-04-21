# Setup CLI parser
from argparse import ArgumentParser
parser = ArgumentParser()
parser.add_argument('-i')      # option that takes a value
parser.add_argument('-o')      # option that takes a value
# access the arguments
args = parser.parse_args()

fp=open(args.o,'w+')
import scanpy as sc
print('Scanpy version:', sc.__version__)
from memory_profiler import profile

import sys
# scperturb utils
sys.path.insert(1, '../../package/src/')
from scperturb import etest

@profile(stream=fp)
def exec_fun():
    etest(adata, obs_key='perturbation', control='control', runs=1000)

if __name__ == '__main__':
    adata = sc.read(args.i)
    print('starting profiling...')
    exec_fun()
