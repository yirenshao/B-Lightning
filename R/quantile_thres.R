#' @title quantile_thres
#'
#' @description The function classifies count data into binary low (0) - high (1) data, based on whether the count number is greater than a threshold.
#' @param Y A binary matrix, indicating whether the observation is above the 50% quantile of the non-zero part.
#' @return A categorized count matrix
#' @details The threshold used for classification is defined by quantile regression on each gene using Frisch–Newton interior point method ("fn" option for method variable in quantreg package, rq function).
#' By default if no meta data is provided, the quantile regression would be applied on the mean expression level of each cell.
#' The quantile to be estimated in the quantile regression is set to be the estimated 50% quantile of the non-zero part of the expression level for each gene.
#' If the expression level of a gene is low (with median 0), then the threshold is set to be 0.
#' @references Koenker, R. and S. Portnoy (1997) The Gaussian Hare and the Laplacean Tortoise: Computability of Squared-error vs Absolute Error Estimators, (with discussion). Statistical Science, 12, 279-300.
#' @references Roger Koenker (2022). quantreg: Quantile Regression. R package version 5.94. https://CRAN.R-project.org/package=quantreg
#' @import quantreg
#' @import Matrix
#' @export
#'
#'
quantile_thres <- function(Y){
  n <- nrow(Y)
  Z <- rowMeans(Y)
  Zmat <- cbind(rep(1,n),Z)
  Y_c <- matrix(0, n, ncol(Y))

  for (j in 1:ncol(Y)){
    v5 <- quantile(Y[,j], 0.5)
    if (v5 == 0){
      Y_c[,j] <- (Y[,j] == 0)*1
    } else {
      quant <- sum(Y[,j] <= v5) / n
      fit <- quantreg::rq(Y[,j] ~ Z, quant, method = "fn")
      Q <- fit$fitted.values
      Y_c[,j] <- (Y[,j]<=Q)*1
    }
  }
  return(Y_c)
}
