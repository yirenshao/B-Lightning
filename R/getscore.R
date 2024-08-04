#' Calculate scores of cells
#'
#' This function calculates scores of each cell based on gene markers' expression levels.
#'
#' @param index_gene_up names of upregulated markers
#' @param index_gene_down names of downregulated markers
#' @param subcounts count matrix (gene markers by cells)
#' #' @return score vector
#' @export


getscore = function(score = "CFS",index_gene_up, index_gene_down, subcounts){

  n = dim(subcounts)[2]

if (score == "CFS"){
  if(length(index_gene_up) == 0){
    temp.m.u = rep(0,n)}else{
      temp.m.u = subcounts[index_gene_up,]
      if (length(index_gene_up) == 1){
        temp.m.u = (rank(temp.m.u,ties.method = "min") - 1) / n

      }else{
        for ( i in 1:length(index_gene_up))
        {
          temp.m.u[i,] = (rank(temp.m.u[i,],ties.method = "min") - 1) / n

        }
      }
      temp.m.u[temp.m.u <= 0.5] = 0
      temp.m.u = temp.m.u / length(index_gene_up)
    }
  if (length(index_gene_down) == 0){
    temp.m.d = rep(0,n)}else{
      temp.m.d = subcounts[index_gene_down,]
      if(length(index_gene_down) == 1){
        temp.m.d= (rank(-temp.m.d,ties.method = "min") - 1) / n

      }else{
        for ( i in 1:length(index_gene_down))
        {
          temp.m.d[i,] = (rank(-temp.m.d[i,],ties.method = "min") - 1) / n

        }}

      temp.m.d[temp.m.d <= 0.5] = 0
      temp.m.d = temp.m.d / length(index_gene_down)
    }

}


  temp.m = rbind(temp.m.u,temp.m.d)
  cscore = colSums(temp.m)
  return (cscore)

}
