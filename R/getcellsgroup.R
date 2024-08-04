#' @title getcellsgroup
#'
#' @description
#' Regroup cells in the Seurat object
#'
#' @details This function divides cells into two groups according to their scores
#' calculated by gene markers' expression levels. The score threshold is set
#' based on the estimation of featured cells proportion in the population.
#'
#' @param so Seurat Object
#' @param index_gene_up a char vector of up-regulated markers (as the rownames of seurat object)
#' @param index_gene_down  a char vector of down-regulated markers (as the rownames of seurat object)
#' @param score cellular score to use GSVA/CFS
#' @param estimated.nonfeatured.proportion estimated featured cells proportion
#' @return Updated Seurat Object with new cell group in the metadata
#' @export
#' @import GSVA
#'

getcellsgroup <-function(so,
                          index_gene_up,index_gene_down,
                          score,estimated.nonfeatured.proportion){
  subcounts = so@assays$RNA@counts[c(index_gene_up,index_gene_down),]

  if (score == "GSVA"){

    GSVAscore <- GSVA::gsva(subcounts, list(index_gene_up,index_gene_down), verbose = FALSE)

    z <- GSVAscore[1,] - GSVAscore[2,]
    z2 = (z - mean(z))/sd(z)

  }
  if (score == "CFS"){

    z = getscore(score,index_gene_up, index_gene_down,subcounts)

    z2 = (z - mean(z))/sd(z)

  }

  threshold = as.numeric(quantile(z2,estimated.nonfeatured.proportion))
  #length(z2[z2>threshold])
  index_senecells = which(z2 >= threshold)
  #index_senecells = which(z > 0)


  cellsgroup = rep("non-featured",dim(so)[2])
  for (i in index_senecells){
    cellsgroup[i] = "featured"
  }
  so@meta.data$group = cellsgroup
  return (so)
}
