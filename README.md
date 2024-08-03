# R package B-Lightning repository
<img src="https://github.com/yirenshao/B-Lightning/blob/master/blightning.jpg?raw=true" width="100">
This repository contains the codes for the R package "B-Lightning", a novel and robust method designed to identify
heterogeneity-source-specific marker genes and the corresponding cell subpopulations that
are differentiated by a particular source of heterogeneity (e.g., cell activation state), isolating
it from other sources of heterogeneity (e.g., cell type, cell cycle phase).

<img src="https://github.com/yirenshao/B-Lightning/blob/master/blightning_workflow.jpeg?raw=true" width="300">

B-Lightning uses an iterative approach that repeatedly expands an initial set of verified feature genes; each iteration
uses differential gene expression (DGE) analysis to select candidate genes and checks
for connectivity with known markers in a gene co-expression network to exclude false discoveries.

The inputs are just a Seurat object and verified "bait" genes (split into 2 vectors: up-regulated or down-regulated). If all "bait" genes are up-regulated (or down-regulated), then leave another vector as empty vector.


## To install the package:
  
```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yirenshao/B-Lightning")
```
