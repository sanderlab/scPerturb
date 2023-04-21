fp=open('memory_profiler.log','w+')
# run as
# mprof run profile_mem.py
# mprof plot
# or
# python -m memory_profiler profile_mem.py
import scanpy as sc
print('Scanpy version:', sc.__version__)
import cProfile
from memory_profiler import profile

import sys
sys.path.append("..")
from scperturb import edist

@profile(stream=fp)
def exec_fun():
    edist(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean')

if __name__ == '__main__':
    adata = sc.read('/fast/scratch/users/peidlis_c/perturbation_resource_paper/tmp_data_PapalexiSatija2021_eccite_RNA.h5')
    print('starting profiling...')
    #cProfile.run("estats = edist(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean')", 'cpstats.prof')
    exec_fun()