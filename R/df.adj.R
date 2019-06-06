
df.adj <- function(object, vars = 1:length(object$coefficients))
{
  x <- model.matrix(object)[, vars, drop = FALSE]
  w <- sqrt(object$weights)
  if(ncol(x) <= nrow(x)){ 
    x.star <- rbind(x, diag(NCOL(x))) * c(w, 1/object$prior.sd[vars])
    V <- chol2inv(qr(x.star)[[1]]) 
  }
  else{
    d <- object$prior.sd[vars]^2
    d <- ifelse(d > 1000, 1000, d)
    d <- diag(d)
    x0 <- x * w
    I <- diag(NROW(x))
    xd <- x0 %*% d %*% t(x0)
    x1 <- solve(I + xd)
    x2 <- t(x0) %*% x1 %*% x0
    x3 <- (x2 * diag(d)) * diag(d)
    V <- d - x3
  }
  df.adj <- ncol(V) - sum((1/object$prior.sd[vars]^2) * diag(V))
  names(df.adj) <- "Effective number of parameters:"
  return(round(df.adj, digits = 2))
}
