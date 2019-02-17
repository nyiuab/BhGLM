
# simulation functions

sim.x <- function(n, m, group = NULL, corr = 0.6, v = rep(1, m), p = 0.5, genotype = NULL, 
                  method = c("svd", "chol", "eigen"), joint = TRUE, verbose = FALSE)
{
  start.time <- Sys.time()
  if (is.list(group)) {
    vars <- unlist(group)
    m <- length(unique(vars))
  }
  else vars <- paste("x", 1:m, sep = "")
  g <- Grouping(all.var = vars, group = group)
  group.vars <- g$group.vars
  linked.vars <- link.vars(group.vars = group.vars)
  V <- diag(m)
  rownames(V) <- colnames(V) <- unique(vars)
  for (j in 1:m) V[j, linked.vars[[j]]] <- 1
  corr1 <- corr[1]
  corr2 <- 0
  if (length(corr) > 1) corr2 <- corr[2]
  V <- ifelse(V == 1, corr1, corr2)
  diag(V) <- v
   
  method <- method[1]
  x <- matrix(0, n, m)
  colnames(x) <- unique(vars)
  if (joint) x <- rmvnorm(n = n, mean = rep(0, m), sigma = V, method = method)
  else {
    for (j in 1:length(group.vars)) {
      names <- group.vars[[j]]
      x[, names] <- rmvnorm(n = n, mean = rep(0, length(names)), sigma = V[names, names, drop = FALSE], method = method)  
    }
  }
  
  if (!is.null(genotype)){
    for(j in genotype){
      probs <- c((1 - p)^2, 2 * (1 - p) * p, p^2)
      quantiles <- quantile(x[, j], c(cumsum(probs)[1:2], 1))
      x[, j] <- as.numeric( factor(cut(x[, j], breaks = c(-Inf, quantiles, Inf))) ) - 1   
    }
  }
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose) cat("simulation time:", minutes, "minutes \n")
  
  return(data.frame(x))
}


sim.eta <- function(x, mu = 0, coefs = NULL, herit = 0.1, sigma = 1, p.neg = 0.5)
{
# p.neg = proportion of negative effects 
# herit = herit of variants
  x <- as.matrix(x)
  x.var <- apply(x, 2, var, na.rm = TRUE)
  x <- x[, which(x.var != 0), drop = FALSE]
  m <- ncol(x)
  if(sum(herit) == 0) coefs <- 0
  if(!is.null(coefs)) {
    if(length(coefs) < m) coefs <- c(coefs, rep(0, m - length(coefs)))
    coefs <- coefs[1:m]
    eta <- x %*% coefs
  }
  if(is.null(coefs)) {
    if(length(herit) == 1){
      b <- ( herit * sigma^2 / (m * mean(x.var) * (1 - herit)) )^0.5
      or1 <- 1 
      or2 <- 2 * exp(b) - or1 
      coefs <- log(runif(ncol(x), or1, or2))
      if(ncol(x) == 1) coefs <- b
    }
    if(length(herit) > 1){
      if(length(herit) < m) herit <- c(herit, rep(0, m - length(herit)))
      if(length(herit) > m) herit <- herit[1:m]
      var.y <- sigma^2 / (1 - sum(herit))
      coefs <- (herit * var.y / x.var)^0.5
    }
    if(p.neg != 0){
      m.non0 <- length(coefs[coefs != 0])
      neg <- sample(x = 1:m.non0, size = m.non0 * p.neg, replace = FALSE) 
      coefs[coefs != 0][neg] <- - coefs[coefs != 0][neg]
    }  
    eta <- x %*% coefs
  } 
  eta <- eta + mu 
  var.sum <- sum(x.var * coefs^2) + sigma^2
  herit <- x.var * coefs^2 / var.sum
  out <- list(eta = eta, coefs = coefs, sigma = sigma, herit = herit)
  return(out)             
}


sim.y <- function(x, mu = 0, coefs = NULL, herit = 0.1, p.neg = 0.5, sigma = 1, 
                  quantiles = 0.5, theta = 3, df = 3)
{
  out0 <- sim.eta(x = x, mu = mu, coefs = coefs, herit = herit, sigma = sigma, p.neg = p.neg)
  eta <- out0$eta
  sigma <- out0$sigma
  n <- length(eta)
  y.normal <- rnorm(n, eta, sigma)
  quantiles <- quantile(y.normal, quantiles)
  y.ordinal <- as.numeric( factor(cut(y.normal, breaks = c(-Inf, quantiles, Inf))) ) - 1  
  y.poisson <- rpois(n, exp(eta))
  y.nb <- rnbinom(n, mu = exp(eta), size = theta)
  
  mu.beta <- exp(eta)/(1 + exp(eta))
  y.beta <- rbeta(n, mu.beta*theta, (1-mu.beta)*theta)
  y.t <- as.vector(eta + rt(n, df)*sigma)
  
#  y <- rexp(n, exp(y.normal))
  y <- rexp(n, exp(eta))
  c <- rexp(n, exp(y.normal - eta))
  tt <- ifelse(y > c, c, y) # min(y, c)
  d <- ifelse(y > c, 0, 1) # 1: uncensored; 0: censored
#  tt <- rexp(n, exp(eta))
#  d <- sample(x = c(0, 1), size = length(tt), replace = T, prob = c(p.censored, 1 - p.censored)) # 1: uncensored; 0: censored
#  tt <- ifelse(d == 0, runif(length(tt), 0, tt), tt)
  library(survival)
  y.surv <- Surv(tt, d)
  
  out <- c(list(y.normal = y.normal, y.ordinal = y.ordinal, y.poisson = y.poisson, y.nb = y.nb, 
                y.beta = y.beta, y.t = y.t, y.surv = y.surv), out0) 

  return(out)  
}

#*******************************************************************************

sim.out <- function(coefs.p, coefs.est, alpha = c(0.05, 0.01))
{
# coefs.p: a matrix of p-values of coefficients 
# coefs.est: a matrix of coefficients estimates
# row: coefficients. column: simulations
  p.if <- coefs.p
  b <-  coefs.est
  power <- NULL
  for(j in 1:length(alpha)) {
    for(s in 1:ncol(coefs.p)) p.if[, s] <- (coefs.p[, s] <= alpha[j])
    power <- cbind(power, apply(p.if, 1, mean))
  }
  colnames(power) <- alpha 
  
  b.mean <- apply(b, 1, mean)
  b.median <- apply(b, 1, median)
  b.l <- apply(b, 1 ,quantile, prob = 0.025)
  b.h <- apply(b, 1, quantile, prob = 0.975)
  est <- cbind(b.mean, b.median, b.l, b.h)
  colnames(est) <- c("mean", "median", "2.5%", "97.5%")
   
  out <- list(power = power, est = est)
  out
}

#*******************************************************************************
# this function is package mvtnorm

rmvnorm <- function (n, mean = rep(0, nrow(sigma)), sigma = diag(length(mean)), 
                     method = c("svd", "eigen", "chol"), pre0.9_9994 = FALSE) 
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    if (length(mean) != nrow(sigma)) 
        stop("mean and sigma have non-conforming size")
    method <- match.arg(method)
    R <- if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    }
    else if (method == "svd") {
        s. <- svd(sigma)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
    }
    else if (method == "chol") {
        R <- chol(sigma, pivot = TRUE)
        R[, order(attr(R, "pivot"))]
    }
    retval <- matrix(rnorm(n * ncol(sigma)), nrow = n, byrow = !pre0.9_9994) %*% R
    retval <- sweep(retval, 2, mean, "+")
    nm <- names(mean)
    if (is.null(nm) && !is.null(colnames(sigma))) nm <- colnames(sigma)
    colnames(retval) <- nm
    if (n == 1) drop(retval)
    else retval
}


decom <- function (sigma, method = c("svd", "eigen", "chol")) 
{
    if (!isSymmetric(sigma, tol = sqrt(.Machine$double.eps), 
        check.attributes = FALSE)) {
        stop("sigma must be a symmetric matrix")
    }
    method <- match.arg(method)
    R <- if (method == "eigen") {
        ev <- eigen(sigma, symmetric = TRUE)
        if (!all(ev$values >= -sqrt(.Machine$double.eps) * abs(ev$values[1]))) {
            warning("sigma is numerically not positive definite")
        }
        t(ev$vectors %*% (t(ev$vectors) * sqrt(ev$values)))
    }
    else if (method == "svd") {
        s. <- svd(sigma)
        if (!all(s.$d >= -sqrt(.Machine$double.eps) * abs(s.$d[1]))) {
            warning("sigma is numerically not positive definite")
        }
        t(s.$v %*% (t(s.$u) * sqrt(s.$d)))
    }
    else if (method == "chol") {
        R <- chol(sigma, pivot = TRUE)
        R[, order(attr(R, "pivot"))]
    }

    R
}

#*******************************************************************************
