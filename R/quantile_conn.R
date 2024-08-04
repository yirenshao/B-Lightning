#' Filter out candidates with low connectivity
#'
#'This function calculates each candidate's connectivity with known gene markers
#'and only keep the ones with high connectivity based on set threshold.
#'
#' @param index_marker_known input/known gene markers
#' @param index_marker_candidate candidate gene markers to be tested for their
#' connectivity with input/known gene markers
#' @param so Seurat Object
#' @param cutoff threshold of connectivity
#' @return a char vector of genes with high connectivity
#' @export

quantile_conn = function(index_marker_known,index_marker_candidate,so,cutoff,num.variablefeatures){
  ret = vector()
  VF = Seurat::VariableFeatures(so)
  countdat = so@assays$RNA@counts
  if(length(index_marker_candidate) == 0){
    return(ret)
  }else{
    index_marker_candidate = setdiff(index_marker_candidate,index_marker_known)
    if(length(index_marker_candidate) == 0){
      return(ret)
    }else{
      Y = countdat[unique(append(index_marker_known,index_marker_candidate)),]

      Y1 = countdat[append(index_marker_known,VF),]
      Y <- t(Y) #countdat should be gene by cell
      Y1 <- t(Y1)
      Y_c <- quantile_thres(Y) #A matrix,binary,whether the observation is above the 50% quantile of the non-zero part.
      #sourceCpp("coexp_arma.cpp")
      coex <- coexp_arma(Y_c)
      diag(coex) <- 0
      #coex_vec <- coex[upper.tri(coex, diag = F)]

      Y1_c <- quantile_thres(Y1) #A matrix,binary,whether the observation is above the 50% quantile of the non-zero part.
      coex1 <- coexp_arma(Y1_c)
      diag(coex1) <- 0
      #coex1_vec <- coex1[upper.tri(coex1, diag = F)]

      for (i in 1:length(index_marker_candidate)){
        w_gag0 = coex[i + length(index_marker_known),1:length(index_marker_known)]
        count = 0
        for (j in 1:length(w_gag0)){ #length(marker_known)
          w_g0g = tail(coex1[j,],num.variablefeatures)
          threshold1 = quantile(w_g0g,0.025,na.rm = TRUE)
          threshold2 = quantile(w_g0g,0.975, na.rm = TRUE)

          if((!is.na(threshold1) & !is.infinite(threshold1)) |(!is.na(threshold2) & !is.infinite(threshold2))){
            if(w_gag0[j] >= threshold2 | w_gag0[j] <= threshold1){
              count = count + 1
            }
          }}
        if (count > cutoff )  {
          ret = append(ret,index_marker_candidate[i])
        }
      }
      return(ret)}}
}
