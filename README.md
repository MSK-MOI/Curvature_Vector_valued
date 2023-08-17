# Curvature_Vector_valued
This is implementation for [Multi-omic integrated curvature study on pan-cancer genomic data](https://link.springer.com/article/10.1007/s00498-023-00360-7)

## Input 
- adj: Adjecency matrix
- CNA, RNA, Methyl: omics data
- Methyl: optional 
- Implemented for 2 or 3 channels

## Curvature definitions
### Edge weights (within each layer)
- edge k := gene i ~ gene j

- Markov transition $\rho_{ij}$

- $w_k=1/\sqrt(\rho_{ij}+\rho_{ji})$
### Edge weights (cross layer)

- gamma>0: fix valued weights
- gamma<0: corr based, p-values between layer (each gene), p*50+1

### Transport 

- Nc layers transition markov distribution at i and j
- Nc layers to Nc layers

### Curvature normalization

- curv=1-Wi/cc
- $cc=\sum w_ij$ at Nc layers


## Output
vectors of length K = num of edges

- curv_i : curvature values
- curv_c_i: curvature*cc
- W_i transport distances

## Functions
- vec_curv_fun: main function
  from adj construct big adj (multilayer)
  compute transport
  compute curvature
- test_function:
  run on cluster
  check complete
  resumit jobs
- kantorovich, Kantorovich_linprog:
  Compute Wasserstein distance
  linprog version is faster
