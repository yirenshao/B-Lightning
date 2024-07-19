#library(quantreg)
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
