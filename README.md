<img src="https://github.com/yirenshao/B-Lightning/blob/master/blightning.jpg?raw=true" width="200">
  
# R package BLightning repository
  
  
This repository contains the codes and a tutorial for the R package "BLightning" v0.1.0. 
  
B-Lightning (R package "BLightning") is a novel and robust method designed to identify
heterogeneity-source-specific marker genes and the corresponding cell subpopulations that
are differentiated by a particular source of heterogeneity (e.g., cell activation state), isolating
it from other sources of heterogeneity (e.g., cell type, cell cycle phase).
   
<img src="https://github.com/yirenshao/B-Lightning/blob/master/blightning_workflow.jpeg?raw=true">

  
B-Lightning uses an iterative approach that repeatedly expands an initial set of verified feature genes; each iteration
uses differential gene expression (DGE) analysis to select candidate genes and checks
for connectivity with known markers in a gene co-expression network to exclude false discoveries.
  
The inputs are just a Seurat object and verified "bait" genes (split into 2 vectors: up-regulated or down-regulated). If all "bait" genes are up-regulated (or down-regulated), then leave another vector as empty vector.
  

## Easy Installation:
  
```R
if (!require("devtools", quietly = TRUE))
    install.packages("devtools")
devtools::install_github("yirenshao/B-Lightning")
```

## Easy Usage:

```R
so = readRDS('toy_s1_so.rds')
upregulated.markers = c("gene 202", "gene 203" ,"gene 204", "gene 205", "gene 206", "gene 207", "gene 210")
downregulated.markers = c("gene 201" ,"gene 208", "gene 209")
ret = runBLightning(so,upregulated.markers,
                                      downregulated.markers,
                                      score = "CFS",
                                      estimated.nonfeatured.proportion = 0.92,
                                      connectivity.cutoff = 4,
                                      num.variablefeatures = 2000,
                                      alpha.genes = 0.05,
                                      max.iter = 10)
                                      
```
