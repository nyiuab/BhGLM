
bpolr <- function (formula, data, weights, start, subset, na.action, 
                   method = c("logistic", "probit", "loglog", "cloglog", "cauchit"), 
                   contrasts = NULL, Hess = TRUE, prior = Student(0, 0.5, 1),  
                   verbose = FALSE, ...) 
{
  if (!requireNamespace("MASS")) install.packages("MASS")
    library(MASS)
    start.time <- Sys.time()
    
    prior.mean <- prior$mean
    prior.scale <- prior$scale
    prior.df <- prior$df
    
    model <- TRUE 
    logit <- function(p) log(p/(1 - p))
    
    dt.deriv <- function(x, mean, scale, df, log = TRUE, delta = 0.001) {
      (dt((x + delta - mean)/scale, df, log = log) - dt((x - 
        delta - mean)/scale, df, log = log))/(2 * delta)
    }
    fmin <- function(beta) {
      theta <- beta[pc + 1:q]
      gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
      eta <- offset
      if (pc > 0) eta <- eta + drop(x %*% beta[1:pc])
      pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
      if (all(pr > 0)) f <- -sum(wt * log(pr))
      else f <- Inf
      if (pc > 0) f <- f - sum(dt((beta[1:pc] - prior.mean)/prior.scale, prior.df, log = TRUE))
      return(f)
    }
    gmin <- function(beta) {
      jacobian <- function(theta) {
        k <- length(theta)
        etheta <- exp(theta)
        mat <- matrix(0, k, k)
        mat[, 1] <- rep(1, k)
        for (i in 2:k) mat[i:k, i] <- etheta[i]
        mat
      }
      theta <- beta[pc + 1:q]
      gamm <- c(-100, cumsum(c(theta[1], exp(theta[-1]))), 100)
      eta <- offset
      if (pc > 0) eta <- eta + drop(x %*% beta[1:pc])
      pr <- pfun(gamm[y + 1] - eta) - pfun(gamm[y] - eta)
      p1 <- dfun(gamm[y + 1] - eta)
      p2 <- dfun(gamm[y] - eta)
      g1 <- if (pc > 0) 
        t(x) %*% (wt * (p1 - p2)/pr)
      else numeric(0)
      xx <- .polrY1 * p1 - .polrY2 * p2
      g2 <- -t(xx) %*% (wt/pr)
      g2 <- t(g2) %*% jacobian(theta)
      if (pc > 0) g1 <- g1 - dt.deriv(beta[1:pc], prior.mean, prior.scale, prior.df, log = TRUE)
      if (all(pr > 0)) c(g1, g2)
      else rep(NA, pc + q)
    }
    
    m <- match.call(expand.dots = FALSE)
    mf <- match(c("formula", "data", "subset", "weights", "na.action", 
                  "etastart", "mustart", "offset"), names(m), 0)
    m <- m[c(1, mf)]
    method <- match.arg(method)
    pfun <- switch(method, logistic = plogis, probit = pnorm, 
                   loglog = pgumbel, cloglog = pGumbel, cauchit = pcauchy)        
    dfun <- switch(method, logistic = dlogis, probit = dnorm, 
                   loglog = dgumbel, cloglog = dGumbel, cauchit = dcauchy)        
    if (is.matrix(eval.parent(m$data))) m$data <- as.data.frame(data)
    m$start <- m$Hess <- m$method <- m$... <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval.parent(m)
    Terms <- attr(m, "terms")
    x <- model.matrix(Terms, m, contrasts)
    xint <- match("(Intercept)", colnames(x), nomatch = 0)
    n <- nrow(x)
    pc <- ncol(x)
    cons <- attr(x, "contrasts")
    if (xint > 0) {
        x <- x[, -xint, drop = FALSE]
        pc <- pc - 1
    }
    else warning("an intercept is needed and assumed")
    wt <- model.weights(m)
    if (!length(wt)) wt <- rep(1, n)
    offset <- model.offset(m)
    if (length(offset) <= 1) offset <- rep(0, n)
    y <- model.response(m)
    if (!is.factor(y)) stop("response must be a factor")
    lev <- levels(y)
    if (length(lev) <= 2) stop("response must have 3 or more levels")
    y <- unclass(y)
    q <- length(lev) - 1
    Y <- matrix(0, n, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    
    if (missing(start)) coefs <- rep(0, pc)
    else {
      coefs <- start
      if (length(start) != pc) 
        stop("'start' is not of the correct length")
    }
    xb <- as.numeric(x %*% coefs) + offset   
    thetas <- polr(as.factor(y) ~ offset(xb), weights = wt, method = method)$zeta
    start <- c(coefs, thetas)
      
    J <- NCOL(x)
    if (length(prior.mean) < J) 
      prior.mean <- c(prior.mean, rep(prior.mean[length(prior.mean)], J - length(prior.mean)) )
    if (length(prior.scale) < J) 
      prior.scale <- c(prior.scale, rep(prior.scale[length(prior.scale)], J - length(prior.scale)) )
    if (length(prior.df) < J) 
      prior.df <- c(prior.df, rep(prior.df[length(prior.df)], J - length(prior.df)) )

    prior.scale <- prior.scale / autoscale(x, min.x.sd=1e-04)
    
    prior.counts.for.bins <- 1/(q + 1)
    if (length(prior.counts.for.bins) == 1) prior.counts.for.bins <- rep(prior.counts.for.bins, q + 1)
    
    y.0 <- y
    Y.0 <- Y
    x.0 <- x
    wt.0 <- wt
    offset.0 <- offset
    .polrY1.0 <- .polrY1
    .polrY2.0 <- .polrY2
    y <- c(y.0, 1:(q + 1))
    Y <- matrix(0, n + q + 1, q)
    .polrY1 <- col(Y) == y
    .polrY2 <- col(Y) == y - 1
    x <- rbind(x.0, matrix(colMeans(x.0), nrow = (q + 1), ncol = J, byrow = TRUE))
    wt <- c(wt.0, prior.counts.for.bins)
    offset <- c(offset, rep(0, q + 1))

    res <- optim(start, fmin, gmin, method = "BFGS", hessian = Hess)

    y <- y.0
    Y <- Y.0
    x <- x.0
    wt <- wt.0
    offset <- offset.0
    .polrY1 <- .polrY1.0
    .polrY2 <- .polrY2.0

    beta <- res$par[seq_len(pc)]
    theta <- res$par[pc + 1:q]
    zeta <- cumsum(c(theta[1], exp(theta[-1])))
    deviance <- 2 * res$value
    if (pc > 0) deviance <- deviance + 2 * sum(dt((beta - prior.mean)/prior.scale, prior.df, log = TRUE)) # I add
    niter <- c(f.evals = res$counts[1], g.evals = res$counts[2])
    names(zeta) <- paste(lev[-length(lev)], lev[-1], sep = "|")
    if (pc > 0) {
        names(beta) <- colnames(x)
        eta <- drop(x %*% beta)
    }
    else {
        eta <- rep(0, n)
    }
    cumpr <- matrix(pfun(matrix(zeta, n, q, byrow = TRUE) - eta), , q)
    fitted <- t(apply(cumpr, 1, function(x) diff(c(0, x, 1))))
    dimnames(fitted) <- list(row.names(m), lev)
    fit <- list(coefficients = beta, zeta = zeta, deviance = deviance,  
                fitted.values = fitted, lev = lev, terms = Terms, df.residual = sum(wt) - pc - q, 
                edf = pc + q, n = sum(wt), nobs = sum(wt), 
                call = match.call(), method = method, convergence = res$convergence, niter = niter, lp = eta)
    if (Hess) {
        dn <- c(names(beta), names(zeta))
        H <- res$hessian
        dimnames(H) <- list(dn, dn)
        fit$Hessian <- H
    }
    if (model) fit$model <- m
    fit$na.action <- attr(m, "na.action")
    fit$contrasts <- cons
    fit$xlevels <- .getXlevels(Terms, m)
    fit$method <- method
    fit$prior <- list(prior="Stendent-t", mean=prior.mean, scale=prior.scale, df=prior.df)
    
    class(fit) <- c("polr")
    stop.time <- Sys.time()
    minutes <- round(difftime(stop.time, start.time, units="min"), 3)
    if (verbose) cat("Computational time:", minutes, "minutes \n")

    fit
}

#*********************************************************************
# from MASS
pgumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE) 
{
  q <- (q - loc)/scale
  p <- exp(-exp(-q))
  if (!lower.tail) 
    1 - p
  else p
}
dgumbel <- function (x, loc = 0, scale = 1, log = FALSE) 
{
  x <- (x - loc)/scale
  d <- log(1/scale) - x - exp(-x)
  if (!log) 
    exp(d)
  else d
}
pGumbel <- function (q, loc = 0, scale = 1, lower.tail = TRUE) 
{
  q <- (q - loc)/scale
  p <- exp(-exp(q))
  if (lower.tail) 
    1 - p
  else p
}
dGumbel <- function (x, loc = 0, scale = 1, log = FALSE) 
{
  x <- -(x - loc)/scale
  d <- log(1/scale) - x - exp(-x)
  if (!log) 
    exp(d)
  else d
}
#**************************************************************************
