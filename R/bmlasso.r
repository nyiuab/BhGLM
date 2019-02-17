

bmlasso <- function(x, y, family = c("gaussian", "binomial", "poisson", "cox"), weights = rep(1, nrow(x)), offset = NULL,
                    epsilon = 1e-04, maxit = 50, init = NULL, group = NULL, 
                    prior = c("mde", "mt"), prior.df = 1, ss = c(0.04, 0.5),
                    Warning = FALSE, verbose = TRUE) 
{
  require(glmnet)
  start.time <- Sys.time()
  call <- match.call()
  if (any(is.na(x)) | any(is.na(y))) {
    a <- apply(cbind(y,x), 1, function(z) !any(is.na(z)))
    y <- y[a]
    x <- x[a,]
  }
  x <- as.matrix(x)
  if (NROW(x) != NROW(y)) stop("nobs of 'x' and 'y' are different")
  if (is.null(colnames(x))) colnames(x) <- paste("x", 1:ncol(x), sep = "")
  family <- family[1]
  if (family == "cox")  
    if (!is.Surv(y)) stop("'y' should be a 'Surv' object")
  if (family == "gaussian") y <- (y - mean(y))/sd(y)
  if (!is.null(init) & length(init) != ncol(x)) stop("give an initial value to each coefficient (not intercept)")

  f <- bmlasso.fit(x = x, y = y, family = family, weights = weights, offset = offset, epsilon = epsilon, maxit = maxit, init = init,
                   group = group, prior = prior, prior.df = prior.df, ss = ss, Warning = Warning)
  
  f$call <- call
  class(f) <- c("glmnet", "bmlasso", "GLM")
  if (family == "cox") class(f) <- c("glmnet", "bmlasso", "COXPH") 
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose){
    cat("EM Coordinate Decent Iterations:", f$iter, "\n")
    cat("Computational time:", minutes, "minutes \n")
  }

  return(f)
}

# ******************************************************************************

bmlasso.fit <- function(x, y, family = "gaussian", weights = rep(1, nrow(x)), offset = NULL, epsilon = 1e-04, maxit = 50, 
                        init = rep(0, ncol(x)), group = NULL, prior = "mde", prior.df = 1, ss = c(0.04, 0.5), 
                        Warning = FALSE)
{  
  ss <- sort(ss)
  ss <- ifelse(ss <= 0, 0.001, ss)
  
  prior.sd <- prior.scale <- ss[length(ss)]  # used for ungrouped coefficients 
  prior <- prior[1]

  if (family == "cox") intercept <- FALSE
  else intercept <- TRUE
  x0 <- x
  if (intercept) x0 <- cbind(1, x)
  d <- prepare(x = x0, intercept = intercept, prior.mean = 0, prior.sd = prior.sd, prior.scale = prior.scale, 
               prior.df = prior.df, group = group)
  x <- d$x
  prior.sd <- d$prior.sd
  prior.scale <- d$prior.scale 
  prior.df <- d$prior.df
  sd.x <- d$sd.x
  min.x.sd <- d$min.x.sd
  group <- d$group
  group.vars <- d$group.vars
  ungroup.vars <- d$ungroup.vars
  if (intercept){
    x <- x[, -1]
    prior.sd <- prior.sd[-1]
    prior.scale <- prior.scale[-1]
    prior.df <- prior.df[-1]
  }
  
  if (prior == "mde" | prior == "mt") {
    if (length(ss) != 2) stop("ss should have two positive values")
    theta <- rep(0.5, length(group.vars))
    p <- rep(0.5, length(prior.sd))
    names(p) <- names(prior.sd)  
  }
  
  if (is.null(init)) {
    for (k in 1:5) {
      if (prior == "mde") {
        ps <- ss[1] + (k - 1) * 0.01
        if (family == "cox") ps <- min(ss[1] + (k - 1) * 0.01, 0.08)
        f <- glmnet(x = x, y = y, family = family, weights = weights, offset = offset, alpha = 0.95, lambda = 1/(nrow(x) * ps), standardize = TRUE)
      }
      if (prior == "mt"){
        ps <- ss[1]^2
        f <- glmnet(x = x, y = y, family = family, weights = weights, offset = offset, alpha = 0, lambda = 1/(nrow(x) * ps), standardize = TRUE)
      }
      b <- as.numeric(f$beta)
      if (any(b != 0)) break
    }
  }
  else b <- as.numeric(init)
  
  names(b) <- colnames(x)
  b <- ifelse(b == 0, 0.001, b)
  init <- b
  
  devold <- 0
  conv <- FALSE
  for (iter in 1:maxit){
    
    if (prior == "mde" | prior == "mt") {
      out <- mix(prior = prior, prior.df = prior.df, prior.scale = prior.scale, 
                 group.vars = group.vars, beta0 = b, ss = ss, theta = theta, p = p)
      prior.scale <- out[[1]]   
      p <- out[[2]]
      theta <- out[[3]]
    }
    if (prior == "mt") {
      prior.sd <- update.prior.sd(prior = prior, beta0 = b, prior.scale = prior.scale, 
                                  prior.df = prior.df, sd.x = sd.x, min.x.sd = min.x.sd) 
    }
    if (prior == "mde"){ 
      Pf <- 1/(prior.scale + 1e-10)
      f <- glmnet(x = x, y = y, family = family, weights = weights, offset = offset, alpha = 1, penalty.factor = Pf, 
                  lambda = sum(Pf)/(nrow(x) * ncol(x)), standardize = FALSE)
    }
    if (prior == "mt"){ 
      Pf <- 1/(prior.sd^2 + 1e-10)
      f <- glmnet(x = x, y = y, family = family, weights = weights, offset = offset, alpha = 0, penalty.factor = Pf, 
                  lambda = sum(Pf)/(nrow(x) * ncol(x)), standardize = FALSE)
    }
    
    b <- as.numeric(f$beta) #/sqrt(dispersion)
    names(b) <- colnames(x)
    dev <- deviance(f)
  
    if(abs(dev - devold)/(0.1 + abs(dev)) < epsilon & iter > 5) {
      conv <- TRUE
      break
    }
    else devold <- dev
  }
  if (Warning & !conv) warning("algorithm did not converge", call. = FALSE)
  
  f$x <- x
  f$y <- y
  f$family <- family
  f$ss <- ss
  
  f$coefficients <- as.numeric(coef(f))
  names(f$coefficients) <- rownames(coef(f))
  f$linear.predictors <- predict(f, newx = x, type = "link", offset = offset)
  if (family == "gaussian")
    f$dispersion <- bglm(y ~ f$linear.predictors - 1, weights = weights, start = 1, prior = "de", prior.mean = 1, prior.scale = 0, verbose = FALSE)$dispersion 
  
  f$iter <- iter
  f$prior <- prior
  if (prior == "mt") {
    f$prior.df <- prior.df
    f$prior.sd <- prior.sd
  }
  f$prior.scale <- prior.scale
  f$penalty.factor <- Pf
  f$group <- group
  f$group.vars <- group.vars 
  f$ungroup.vars <- ungroup.vars
  if (prior == "mde" | prior == "mt") {
    f$p <- p[unlist(group.vars)]
    f$theta <- theta
  }
  f$init <- init
  f$aic <- deviance(f) + 2 * f$df
  f$weights <- weights
  
  return(f)
}

#*******************************************************************************

