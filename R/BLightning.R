#' @title B_Lightning
#'
#' @description An iterative fishing method for gene markers
#'
#' @param so Preprocessed Seurat Object
#' @param upregulated.markers a char/index vector of established upregulated gene markers
#' @param downregulated.markers a char/index vector of established downregulated gene markers
#' @param score the type of cellular score used to group cells: CFS/GSVA
#' @param estimated.nonfeatured.proportion estimated proportion of nonfeatured cells
#' in the whole cell population
#' @param connectivity.cutoff number of established genes that candidates have to be correlated with
#' @param num.variablefeatures number of highly variable genes used to calculate the null distribution of connectivity
#' @param alpha.genes TypeI error of identifying differentally expressed genes
#'
#' @return A list that contains the number of iterations, new upregulated gene markers,
#' new downregulated gene markers, pvalues of new upregulated gene markers,
#' pvaluds of new downregualted markers
#'
#' @importFrom Seurat "FindMarkers", Seurat "AddMetaData", Seurat "FindVariableFeatures",
#' @importFrom GSVA "gsva"
#'
#'

B_Lightning <- function(so,upregulated.markers,
                        downregulated.markers,
                        score = "CFS1",
                        estimated.nonfeatured.proportion = 0.9,
                        connectivity.cutoff = 4,
                        num.variablefeatures = 2000,
                        alpha.genes = 0.05
){

  iter = 0
  
  flag = T
  memo = list(upregulated.markers,downregulated.markers,upregulated.markers,downregulated.markers)
  so <- AddMetaData(so, rep("unknown",dim(so)[2]),col.name = "group")
  so <- FindVariableFeatures(so, selection.method = "vst", nfeatures = num.variablefeatures)
  topfeatures <-  head(VariableFeatures(so), num.variablefeatures)


  while(flag){
    start_time = Sys.time()
    iter  = iter + 1
    print("iter = ")
    print(iter)
    so = getcellsgroup(so,unique(c(memo[[3]],upregulated.markers)),
                          unique(c(memo[[4]],downregulated.markers)),
                       score, estimated.nonfeatured.proportion)
    print("getcellsgroup done, now findgenes.")

    temp = findgenes(so,alpha.genes)

    #Check Connectivity!
    print("findgenes done, now check quantile connectivity.")
    candidate1 = quantile_conn(c(upregulated.markers,downregulated.markers),temp[[1]],so,connectivity.cutoff,num.variablefeatures)
    candidate2 = quantile_conn(c(upregulated.markers,downregulated.markers),temp[[2]],so,connectivity.cutoff,num.variablefeatures)
    if( length(setdiff(sort(candidate1),sort(unique(c(memo[[1]],memo[[3]]))))) == 0 && length(setdiff(sort(candidate2),sort(unique(c(memo[[2]],memo[[4]]))))) == 0){
      flag = F
      memo[[1]] = memo[[3]]
      memo[[2]] = memo[[4]]
      memo[[3]] = candidate1
      memo[[4]] = candidate2
    }else{
      memo[[1]] = memo[[3]]
      memo[[2]] = memo[[4]]
      memo[[3]] = candidate1
      memo[[4]] = candidate2
    }
    end_time = Sys.time()
    if (as.numeric(end_time - start_time) > 60 * 60 * 4){
      flag = F
    }

  }

  so = getcellsgroup(so,unique(c(memo[[3]],upregulated.markers)),
                     unique(c(memo[[4]],downregulated.markers)),
                     score, estimated.nonfeatured.proportion)

  temp = findgenes(so,alpha.genes)

  pvalue_up = temp[[3]][c(memo[[4]],memo[[3]]),] %>% filter(avg_log2FC > 0 ) %>% filter(p_val_adj < alpha.genes) %>% dplyr::select(p_val_adj)
  pvalue_down = temp[[3]][c(memo[[4]],memo[[3]]),]%>% filter(avg_log2FC < 0 ) %>% filter(p_val_adj < alpha.genes) %>% dplyr::select(p_val_adj)

  candidate1 = intersect(rownames(pvalue_up),memo[[3]])
  candidate2 = intersect(rownames(pvalue_down),memo[[4]])

  return(list("iterations" = iter,
              "new upregulated markers" = candidate1,
              "new downregulated markers" = candidate2,
              "input upregulated markers" = upregulated.markers,
              "input downregulated markers" = downregulated.markers,
              "pvalues of new upregulated markers" = pvalue_up,
              "pvalues of new downregulated markers" = pvalue_down))
  }



