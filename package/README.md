# scperturb
A python package to compute E-distances in single-cell perturbation data and
perform E-tests.

# Install
Just install via pip:

```
pip install scperturb
```

# Usage example

```
# E-distances
estats = edist(adata, obs_key='perturbation')
# E-distances to a specific group (e.g. 'control')
estats_control = estats.loc['control']
# E-test for difference to control
df = etest(adata, obs_key='perturbation', obsm_key='X_pca', dist='sqeuclidean', control='control', alpha=0.05, runs=100)
```

Check out notebooks --> e-distance for tutorial.
