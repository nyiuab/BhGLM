
waldtest.bh <- function(coefs, Cov, vars=1:length(coefs), weights=NULL)
{
  b <- coefs[vars]
  V <- Cov[names(b), names(b), drop = FALSE] 
  for (i in 1:nrow(V)) V[i, i] <- ifelse(V[i,i] <= 1e-03, 1e-03, V[i, i]) 
  
  # weights can be vector or matrix
  if (is.null(weights)) {
    W <- solve(V)
    stat <- as.numeric( t(b) %*% W %*% b ) 
    p <- pchisq(stat, df=length(b), lower.tail=FALSE)
    test <- list(chisq=stat, df=length(b), pvalue=p)
    eigens <- rep(1, length(b))
  }
  else{
    W <- weights
    if (is.vector(W)) {
      if (length(W) < length(b)) W <- c(W, rep(W[length(W)], length(b) - length(W)) )
      if (length(W) > length(b)) W <- W[1:length(b)]
      W <- diag(W)
    }  
    Q <- as.numeric( t(b) %*% W %*% b )
    eigens <- as.numeric(eigen(V %*% W, only.values=TRUE)[[1]]) 
    test <- quad.chisq(Q, eigens)
  }
  test$variables <- names(b)
  
  return(test)
}


quad.chisq <- function (Q, eigens)
{
  eigens2 <- sum(eigens^2)  
  eigens3 <- sum(eigens^3) 
  chisq <- (Q - (sum(eigens) - eigens2^2/eigens3)) / (eigens3/eigens2)  
  d <- eigens2^3 / eigens3^2  
  p <- pchisq(max(chisq, 0), df=max(d, 0.1), lower.tail=FALSE)
  res <- list(chisq=max(chisq, 0), df=d, pvalue=p) 
  res  
}

#*******************************************************************************
