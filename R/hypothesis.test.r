

#*******************************************************************************

waldtest.bh <- function(object, vars = 2:length(object$coefficients), weights = NULL, show.vars = FALSE)
{
# weights can be vector or matrix
  Cov <- vcov.bh(object)
  b <- object$coefficients[vars]
#  if (length(b) == nrow(Cov)) V <- Cov
#  else {
#    allvars <- names(object$coefficients)
#    vars.in <- which(allvars %in% names(b))
#    vars.not <- which(! allvars %in% names(b))
#    Cov.inv <- solve(Cov)
#    I11 <- Cov.inv[vars.in, vars.in, drop = FALSE]
#    I22 <- Cov.inv[vars.not, vars.not, drop = FALSE]
#    I12 <- Cov.inv[vars.in, vars.not, drop = FALSE]
#    V <- solve(I11 - I12 %*% solve(I22) %*% t(I12))
#  }
  V <- Cov[names(b), names(b), drop = FALSE] 
  for (i in 1:nrow(V)) V[i, i] <- ifelse(V[i,i] <= 1e-03, 1e-03, V[i, i]) 
  
  if (is.null(weights)) {
    W <- solve(V)
    stat <- t(b) %*% W %*% b 
#    p <- 1 - pchisq(stat, df = length(b))
    p <- pchisq(stat, df = length(b), lower.tail = FALSE)
    test <- c(Q = stat, chisq = stat, df = length(b), pvalue = p)
    eigens <- rep(1, length(b))
  }
  else{
    W <- weights
    if (is.vector(W)) {
      if (length(W) < length(b)) W <- c(W, rep(W[length(W)], length(b) - length(W)) )
      if (length(W) > length(b)) W <- W[1:length(b)]
      W <- diag(W)
    }  
    Q <- t(b) %*% W %*% b
    eigens <- as.numeric(eigen(V %*% W, only.values = TRUE)[[1]]) 
    test <- quad.chisq(Q, eigens)
  }
  variables <- names(b)

# calculate and test average effect  
#  if (!is.null(weights)) W <- solve(V)
#  L <- rep(1, nrow(W))
#  d <- L %*% W 
#  V.ave <- 1/sum(d * L)
#  b.ave <- V.ave * (d %*% b)
  L <- rep(1, nrow(W))
  b.ave <- mean(b)
  V.ave <- t(L) %*% V %*% L / length(L)^2
  pvalue <- 2 * pnorm(-abs(b.ave/sqrt(V.ave)))
  average <- c(Estimate = b.ave, Std.Error = sqrt(V.ave), pvalue = pvalue)
  
  if (show.vars)
    res <- list(joint.test = test, average.effect = average, variables = variables)
  else 
    res <- list(joint.test = test, average.effect = average)

  return(res)
}


quad.chisq <- function (Q, eigens)
{
  eigens2 <- sum(eigens^2)  
  eigens3 <- sum(eigens^3) 
  chisq <- (Q - (sum(eigens) - eigens2^2/eigens3)) / (eigens3/eigens2)  
  d <- eigens2^3 / eigens3^2  
#  p <- 1 - pchisq(max(chisq, 0), df = max(d, 0.1), lower.tail = TRUE) 
  p <- pchisq(max(chisq, 0), df = max(d, 0.1), lower.tail = FALSE)
  res <- c(Q = Q, chisq = max(chisq, 0), df = d, pvalue = p) 
  res  
}

#*******************************************************************************

scoretest.bglm <- function (obj.null, x, check.data = TRUE, weights = NULL, joint = TRUE, verbose = TRUE) 
{
  if (any(class(obj.null) %in% "lme")) require(mgcv)
  start.time <- Sys.time()
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) = paste("x", 1:NCOL(x), sep = "")
  
  if (check.data) {
  
    func <- function(d){
      na.prop <- length(which(is.na(d)))/length(d)
      return(na.prop)
    }
    na.prop <- apply(x, 2, func)
    if(any(na.prop > 30/100)) {
      d <- which(na.prop > 30/100)
      if (verbose){
        warning(length(d), " variables have more than 30% missing values: ", call. = FALSE)  
        for(j in 1:length(d)) warning(names(d)[j], call. = FALSE)
      }
    }
  
    v.x <- apply(x, 2, var, na.rm = TRUE)
    if (length(which(v.x == 0)) != 0) {
      x <- x[, which(v.x != 0), drop = FALSE]
      d <- which(v.x == 0)
      if (verbose)
        warning(length(d), " variables with no-variation are removed!", call. = FALSE )
    }
  
    func <- function(d){
      na.index <- which(is.na(d))
      d[na.index] <- mean(d, na.rm = TRUE)
      return(d)
    }
    if (any(is.na(x))) x <- apply(x, 2, func)
  
    x0 <- obj.null$model
    if (any(class(obj.null) %in% "lme")) x0 <- obj.null$data       #model.frame(obj.null, na.action = na.pass)
    sample.inc <- as.numeric(rownames(x0)) #which(!is.na(rowSums(x0)))
    if (length(sample.inc) < nrow(x)){ 
      n <- nrow(x) - length(sample.inc)
      x <- x[sample.inc, , drop = FALSE] 
      if (verbose)
        warning(n, " samples having either missing phenotypes or covariates are excluded!", call. = FALSE )
    }
  
  }
  
  res <- score.bglm(obj.null, x, weights, joint)
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 4)
  if (verbose) 
    cat("Computational time:", minutes, "minutes \n")
  
  return(res)
}

score.bglm <- function (obj.null, x, weights, joint)
{
  if (any(class(obj.null) %in% "glm")) {
    r <- obj.null$residuals
    x0 <- model.matrix(obj.null)
    dispersion <- obj.null$dispersion
    if (is.null(dispersion)) {
      if (obj.null$family$family %in% c("poisson", "binomial")) dispersion <- 1
      else dispersion <- sum(obj.null$weights * r^2)/obj.null$df.residual
    }

    sigma.inv <- obj.null$weights/dispersion
    sigmainv.r <- sigma.inv * r 
    score <- t(x) %*% sigmainv.r 
    cov.null <- vcov(obj.null)
    d <- sigma.inv * x
    v11 <- t(x) %*% d
    v01 <- t(x0) %*% d
    cov.score <- v11 - t(v01) %*% cov.null %*% v01
  }
  
  if (any(class(obj.null) %in% "lme")) { 
    r <- obj.null$residuals[, 1] 
    x0 <- model.matrix(obj.null, data = obj.null$data)
    sigma.inv <- solve( extract.lme.cov(obj.null, obj.null$data) )  # extract.lme.cov is in package 'mgcv'
    sigmainv.r <- sigma.inv %*% r 
    score <- t(x) %*% sigmainv.r 
    d <- sigma.inv %*% x
    v11 <- t(x) %*% d
    v01 <- t(x0) %*% d
    v00 <- t(x0) %*% sigma.inv %*% x0
    cov.score <- v11 - t(v01) %*% solve(v00) %*% v01
  }
  
  if (NCOL(x) == 1) joint <- FALSE
  
  if (joint) {
    if (is.null(weights)) {
      W <- solve(cov.score)
      stat <- t(score) %*% W %*% score 
#     p <- 1 - pchisq(stat, df = NCOL(x))
      p <- pchisq(stat, df = NCOL(x), lower.tail = FALSE)
      res <- c(Q = stat, chisq = stat, df = NCOL(x), pvalue = p)
      eigens <- rep(1, NCOL(x))
    }
    else {
      if (is.vector(weights)){
        if (length(weights) < NCOL(x))
          weights <- c(weights, rep(weights[length(weights)],NCOL(x)-length(weights)))
        if (length(weights) > NCOL(x))
          weights <- weights[1:NCOL(x)] 
        W <- diag(weights)
        Q <- sum(weights * score^2)
        eigens <- as.numeric(eigen(cov.score %*% W, symmetric = TRUE, only.values = TRUE)[[1]])
        res <- quad.chisq(Q, eigens) 
      }
      if (is.matrix(weights)) {
        W <- weights
        Q <- t(score) %*% W %*% score 
        eigens <- as.numeric(eigen(cov.score %*% W, only.values = TRUE)[[1]]) 
        res <- quad.chisq(Q, eigens)
      }
    }
    rownames(W) <- colnames(W) <- colnames(x)
  }
  
  chisq <- score^2/diag(cov.score)
# p <- 1 - pchisq(chisq, df = 1)
  p <- pchisq(chisq, df = 1, lower.tail = FALSE)
  single <- cbind(chisq, p)
  colnames(single) <- c("chisq", "pvalue")
  rownames(single) <- colnames(x)
  
#  list(joint.test = res, single.test = single, score = score, cov.score = cov.score, weights = W, eigenvalues = eigens)
  if (joint) out <- list(joint.test = res, single.test = single)
  else out <- single
  
  return(out)
}

#*******************************************************************************

