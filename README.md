# Curvature_Vector_valued
This is implementation for [Multi-omic integrated curvature study on pan-cancer genomic data](https://link.springer.com/article/10.1007/s00498-023-00360-7)

## Input 
adj: Adjecency matrix
CNA, RNA, Methyl: omics data
Methyl: optional 
Implemented for 2 or 3 channels

## Curvature definitions
### Edge weights (within each layer)
edge k := gene i ~ gene j

Markov transition $\rho_{ij}$

$$w_k=1/sqrt(\rho_{ij}+\rho_{ji})$$

## Output
vectors of length K = num of edges

curv_i : 
