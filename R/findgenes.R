
#' Find gene markers
#'
#' This function identifies differentially expressed genes in the two groups of cells.
#' @param so Seurat Object
#' @param alpha typeI error
#' @return a list consists of identified upregulated gene markers, identified
#' down regulated gene markers and their adjusted pvalues
#' @export


findgenes = function(so,alpha = 0.05){




  Idents(so) <- so@meta.data$group

  ret <- FindMarkers(so,ident.1 = "featured", ident.2 = "non-featured",test.use = "MAST",min.pct = 0.01)

  index_gene_up = ret %>% filter(p_val_adj < alpha) %>%
    filter(avg_log2FC > 0) %>% rownames()
  index_gene_down = ret %>% filter(p_val_adj < alpha) %>%
    filter(avg_log2FC < 0 ) %>% rownames()
  pvalue = ret[c(index_gene_up,index_gene_down),]
  return(list(index_gene_up,index_gene_down,pvalue))
}
