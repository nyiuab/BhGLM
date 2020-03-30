

bmlasso <- function(x, y, family=c("gaussian", "binomial", "poisson", "cox"), offset=NULL,
                    epsilon=1e-04, maxit=50, init=NULL, 
                    ss=c(0.04, 0.5), group=NULL, w.theta=NULL, 
                    Warning=FALSE, verbose=FALSE) 
{
  if (!requireNamespace("glmnet")) install.packages("glmnet")
  require(glmnet)
  start.time <- Sys.time()
  call <- match.call()
  x <- as.matrix(x)
  if (is.null(colnames(x))) colnames(x) <- paste("x", 1:ncol(x), sep = "")
  nobs <- nrow(x)
  if (NROW(y) != nobs) stop("nobs of 'x' and 'y' are different")
  inc <- apply(cbind(y, x), 1, function(z) !any(is.na(z)))
  if (!is.null(offset)) {
    if (length(offset) != nobs) stop("nobs of 'x' and 'offset' are different")
    inc <- apply(cbind(y, x, offset), 1, function(z) !any(is.na(z)))
  }
  y <- y[inc]
  x <- x[inc,]
  offset <- offset[inc]
  family <- family[1]
  if (family == "cox")  
    if (!is.Surv(y)) stop("'y' should be a 'Surv' object")
  if (family == "gaussian") y <- (y - mean(y))/sd(y)
  if (!is.null(init) & length(init) != ncol(x)) stop("give an initial value to each coefficient (not intercept)")

  f <- bmlasso.fit(x=x, y=y, family=family, offset=offset, epsilon=epsilon, maxit=maxit, init=init,
                   group=group, ss=ss, w.theta=w.theta, Warning=Warning)
  
  f$call <- call
  if (family == "cox") class(f) <- c(class(f), "bmlasso", "COXPH") 
  else class(f) <- c(class(f), "bmlasso", "GLM")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose){
    cat("EM Coordinate Decent Iterations:", f$iter, "\n")
    cat("Computational time:", minutes, "minutes \n")
  }

  return(f)
}

# ******************************************************************************

bmlasso.fit <- function(x, y, family="gaussian", offset=NULL, epsilon=1e-04, maxit=50, 
                        init=rep(0, ncol(x)), group=NULL, ss=c(0.04, 0.5), w.theta=NULL, 
                        Warning=FALSE)
{  
  ss <- sort(ss)
  ss <- ifelse(ss <= 0, 0.001, ss)
  prior.scale <- ss[length(ss)]  # used for ungrouped coefficients 

  if (family == "cox") intercept <- FALSE
  else intercept <- TRUE
  x0 <- x
  if (intercept) x0 <- cbind(1, x)
  d <- prepare(x = x0, intercept = intercept, prior.mean = 0, prior.sd = 1, prior.scale = prior.scale, 
               prior.df = 1, group = group)
  x <- d$x
  prior.scale <- d$prior.scale 
  group <- d$group
  group.vars <- d$group.vars
  ungroup.vars <- d$ungroup.vars
  if (intercept){
    x <- x[, -1]
    prior.scale <- prior.scale[-1]
  }
  
  if (length(ss) != 2) stop("ss should have two positive values")
  gvars <- unlist(group.vars)
  theta <- p <- rep(0.5, length(gvars))
  names(theta) <- names(p) <- gvars 
  if (is.null(w.theta)) w.theta <- rep(1, length(gvars))
  if (length(w.theta)!=length(gvars)) stop("all grouped variables should have w.theta")
  if (any(w.theta > 1 | w.theta < 0)) stop("w.theta should be in [0,1]")
  names(w.theta) <- gvars
  
  if (is.null(init)) {
    for (k in 1:5) {
      ps <- ss[1] + (k - 1) * 0.01
      if (family == "cox") ps <- min(ss[1] + (k - 1) * 0.01, 0.08)
      f <- glmnet(x=x, y=y, family=family, offset=offset, alpha=0.95, 
                  lambda=1/(nrow(x) * ps), standardize=TRUE)
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
    
    out <- update.scale.p(b0=b[gvars], ss=ss, theta=theta)
    prior.scale[gvars] <- out[[1]]   
    p <- out[[2]]
    if (!is.matrix(group))
      theta <- update.ptheta.group(group.vars=group.vars, p=p, w.theta=w.theta)
    else theta <- update.ptheta.network(theta=theta, p=p, w=group)
    
    Pf <- 1/(prior.scale + 1e-10)
    f <- glmnet(x = x, y = y, family = family, offset = offset, alpha = 1, 
                penalty.factor = Pf, lambda = sum(Pf)/(nrow(x) * ncol(x)), standardize = FALSE)
    
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
    f$dispersion <- bglm(y ~ f$linear.predictors-1, start=1, prior=De(1,0), verbose=FALSE)$dispersion 
  
  f$iter <- iter
  f$prior.scale <- prior.scale
  f$penalty.factor <- Pf
  f$group <- group
  f$group.vars <- group.vars 
  f$ungroup.vars <- ungroup.vars
  f$p <- p
  f$ptheta <- theta
  f$init <- init
  f$aic <- deviance(f) + 2 * f$df
  f$offset <- offset
  
  return(f)
}

#*******************************************************************************

