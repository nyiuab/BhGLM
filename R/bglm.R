
bglm <- function (formula, family=gaussian, data, offset, weights, subset, na.action, 
           start=NULL, etastart, mustart, control=glm.control(epsilon=1e-04, maxit=50), 
           prior=Student(), group=NULL, method.coef, 
           theta.weights=NULL, inter.hierarchy=NULL, inter.parents=NULL,
           prior.sd=0.5, dispersion=1, Warning=FALSE, verbose=FALSE)  
{
  start.time <- Sys.time()
  
  autoscale <- prior$autoscale
  if (is.null(autoscale)) autoscale <- FALSE
  prior.mean <- prior$mean
  prior.scale <- prior$scale
  if (is.null(prior.scale)) prior.scale <- 0.5
  prior.df <- prior$df
  if (is.null(prior.df)) prior.df <- 1
  ss <- prior$ss
  if (is.null(ss)) ss <- c(0.04, 0.5)
  b <- prior$b
  prior <- prior[[1]]
  if (missing(method.coef)) method.coef <- NULL 
  
  contrasts <- NULL
  call <- match.call()
  if (is.character(family))
    family <- get(family, mode = "function", envir = parent.frame())
  if (is.function(family))
    family <- family()
  if (is.null(family$family)) {
    print(family)
    stop("'family' not recognized")
  }
  if (missing(data))
    data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "weights", "na.action",
        "etastart", "mustart", "offset"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  mf$na.action <- NULL
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  mt <- attr(mf, "terms")
  Y <- model.response(mf, "any")
  if (length(dim(Y)) == 1) {
    nm <- rownames(Y)
    dim(Y) <- NULL
    if (!is.null(nm))
      names(Y) <- nm
  }
  X <- if (!is.empty.model(mt)) 
    model.matrix(mt, mf, contrasts)
  else matrix(, NROW(Y), 0L)
  weights <- as.vector(model.weights(mf))
  if (!is.null(weights) && !is.numeric(weights)) 
    stop("'weights' must be a numeric vector")
  if (!is.null(weights) && any(weights < 0))
    stop("negative weights not allowed")
  offset <- as.vector(model.offset(mf))
  if (!is.null(offset)) {
    if (length(offset) != NROW(Y)) 
      stop(gettextf("number of offsets is %d should equal %d (number of observations)", 
            length(offset), NROW(Y)), domain = NA)
  }
  mustart <- model.extract(mf, "mustart")
  etastart <- model.extract(mf, "etastart")
  
  fit <- bglm.fit(x=X, y=Y, weights=weights, start=start,
                  etastart=etastart, mustart=mustart, offset=offset, 
                  family=family, control=control, intercept=attr(mt, "intercept") > 0,
                  prior=prior, group=group, method.coef=method.coef, 
                  dispersion=dispersion, prior.mean=prior.mean, prior.sd=prior.sd, 
                  prior.scale=prior.scale, prior.df=prior.df, autoscale=autoscale, ss=ss, b=b,
                  theta.weights=theta.weights, inter.hierarchy=inter.hierarchy, inter.parents=inter.parents,
                  Warning=Warning)
  
  fit$model <- mf
  fit$na.action <- attr(mf, "na.action")
  fit <- c(fit, list(call = call, formula = formula, terms = mt,
           data = data, offset = offset, control = control, 
           contrasts = attr(X, "contrasts"), xlevels = .getXlevels(mt, mf)) )
    
  class(fit) <- c("glm", "lm")
  if (family[[1]]==NegBin()[[1]]) class(fit) <- c("negbin", "glm", "lm")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose) {
    cat("EM-IWLS iterations:", fit$iter, "\n")
    cat("Computational time:", minutes, "minutes \n")
  }
  
  fit
}

#*******************************************************************************

bglm.fit <- function (x, y, weights=rep(1, nobs), start=NULL, etastart=NULL, mustart=NULL, 
               offset=rep(0, nobs), family=gaussian(), control=glm.control(), intercept=TRUE,  
               prior="de", group=NULL, method.coef=1, 
               dispersion=1, prior.mean=0, prior.sd=0.5, prior.scale=1, prior.df=1, autoscale=TRUE,
               ss=c(0.05, 0.1), b=1, theta.weights=NULL, inter.hierarchy=NULL, inter.parents=NULL,
               Warning=FALSE)   
{
    ss <- sort(ss)
    ss <- ifelse(ss <= 0, 0.001, ss)
    if (prior == "mde" | prior == "mt")
      prior.sd <- prior.scale <- ss[length(ss)]  # used for ungrouped coefficients
    
    if (is.null(dispersion)) dispersion <- 1
    
    d <- prepare(x=x, intercept=intercept, prior.mean=prior.mean, prior.sd=prior.sd, prior.scale=prior.scale, 
                 prior.df=prior.df, group=group)
    x <- d$x
    prior.mean <- d$prior.mean 
    prior.sd <- d$prior.sd 
    prior.scale <- d$prior.scale 
    prior.df <- d$prior.df 
    sd.x <- d$sd.x
    min.x.sd <- d$min.x.sd
    group <- d$group
    group.vars <- d$group.vars
    ungroup.vars <- d$ungroup.vars
    
    if (autoscale){
      prior.scale <- prior.scale / auto_scale(x, min.x.sd)
      if (family[[1]]=="gaussian") prior.scale <- prior.scale * sd(y)
    }
    
    x0 <- x
    if (intercept) x0 <- x[, -1, drop = FALSE] 
    g0 <- Grouping(all.var = colnames(x0), group = method.coef)
    group0 <- g0$group.vars
    covars0 <- g0$ungroup.vars  
    if (intercept) covars0 <- c(colnames(x)[1], covars0)
    method.coef <- "joint"
    if (length(group0) > 1) method.coef <- "group"
    
    # for mixture prior
    if (prior == "mde" | prior == "mt") {
      if (length(ss) != 2) stop("ss should have two positive values")
      gvars <- unlist(group.vars)
      theta <- p <- rep(0.5, length(gvars))
      names(theta) <- names(p) <- gvars
    
      if (is.null(theta.weights)) theta.weights <- rep(1, length(gvars))
      if (length(theta.weights)!=length(gvars)) stop("all grouped variables should have theta.weights")
      if (any(theta.weights > 1 | theta.weights < 0)) stop("theta.weights should be in [0,1]")
      names(theta.weights) <- gvars
      
      if (length(b) < length(group.vars)) 
        b <- c(b, rep(b[length(b)], length(group.vars) - length(b)) )
      b <- b[1:length(group.vars)]
    }
    
    # for negative binomial model
    nb <- FALSE
    if (family[[1]] == NegBin()[[1]]) nb <- TRUE
    if (nb){
      if (!requireNamespace("MASS")) install.packages("MASS")
      library(MASS)
    }
       
   # *************************      
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    conv <- FALSE
    nobs <- NROW(y)
    nvars <- ncol(x)
    EMPTY <- nvars == 0
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)
    variance <- family$variance
    dev.resids <- family$dev.resids
    aic <- family$aic
    linkinv <- family$linkinv
    mu.eta <- family$mu.eta
    if (!is.function(variance) || !is.function(linkinv))
        stop("'family' argument seems not to be a valid family object")
    valideta <- family$valideta
    if (is.null(valideta))
        valideta <- function(eta) TRUE
    validmu <- family$validmu
    if (is.null(validmu))
        validmu <- function(mu) TRUE
    if (is.null(mustart)) {
        eval(family$initialize)
    }
    else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }
    if (EMPTY) {
        eta <- rep.int(0, nobs) + offset
        if (!valideta(eta))
            stop("invalid linear predictor values in empty model")
        mu <- linkinv(eta)
        if (!validmu(mu))
            stop("invalid fitted means in empty model")
        dev <- sum(dev.resids(y, mu, weights))
        w <- ((weights * mu.eta(eta)^2)/variance(mu))^0.5
        residuals <- (y - mu)/mu.eta(eta)
        good <- rep(TRUE, length(residuals))
        boundary <- conv <- TRUE
        coef <- numeric(0)
        iter <- 0
    }
    else {
        coefold <- NULL
        if (!is.null(etastart)) {
            eta <- etastart
        }
        else if (!is.null(start)) {
            if (length(start) != nvars)
                stop(gettextf("length of 'start' should equal %d and correspond to initial coefs for %s",
                  nvars, paste(deparse(xnames), collapse = ", ")),
                  domain = NA)
            else {
                eta <- offset + as.vector(if (NCOL(x) == 1)
                  x * start
                else x %*% start)
                coefold <- start
            }
        }
        else {
            eta <- family$linkfun(mustart)
        }
        mu <- linkinv(eta)
        if (!(validmu(mu) && valideta(eta)))
            stop("cannot find valid starting values: please specify some")
        devold <- sum(dev.resids(y, mu, weights))
        boundary <- conv <- FALSE
        
        
        if (!is.null(start)){
          coefs.hat <- start
          names(coefs.hat) <- names(prior.mean)
        }
        else coefs.hat <- prior.mean
        dispersionold <- dispersion
        for (iter in 1:control$maxit) {
        
            good <- weights > 0
            varmu <- variance(mu)[good]
            varmu <- ifelse(varmu == 0, 1e-04, varmu)
            if (any(is.na(varmu)))
                stop("NAs in V(mu)")
            if (any(varmu == 0))
                stop("0s in V(mu)")
            mu.eta.val <- mu.eta(eta)
            mu.eta.val <- ifelse(mu.eta.val == 0, 1e-04, mu.eta.val)  
            if (any(is.na(mu.eta.val[good])))
                stop("NAs in d(mu)/d(eta)")
            good <- (weights > 0) & (mu.eta.val != 0)
            if (all(!good)) {
                conv <- FALSE
                warning("no observations informative at iteration ",
                  iter)
                break
            }
            z <- (eta - offset)[good] + (y - mu)[good]/mu.eta.val[good]
            w <- sqrt((weights[good] * mu.eta.val[good]^2)/varmu[good])
            ngoodobs <- as.integer(nobs - sum(!good))
            
            w <- ifelse(w == 0, 1e-04, w) # I add
            
            if (iter > 1) {
              beta0 <- (coefs.hat - prior.mean)/sqrt(dispersion)
              
              if (prior == "mde" | prior == "mt") {
                out <- update.scale.p(prior=prior, df=prior.df[gvars], b0=beta0[gvars], ss=ss, theta=theta)
                prior.scale[gvars] <- out[[1]]   
                p <- out[[2]]
                if (!is.matrix(group))
                  theta <- update.ptheta.group(group.vars=group.vars, p=p, w=theta.weights, b=b)
                else theta <- update.ptheta.network(theta=theta, p=p, w=group) 
                
                if (!is.null(inter.hierarchy))
                  theta.weights <- update.theta.weights(gvars=gvars, 
                                                        theta.weights=theta.weights, 
                                                        inter.hierarchy=inter.hierarchy, 
                                                        inter.parents=inter.parents, 
                                                        p=p)
              }
              
              prior.sd <- update.prior.sd(prior=prior, beta0=beta0, prior.scale=prior.scale, 
                                          prior.df=prior.df, sd.x=sd.x, min.x.sd=min.x.sd) 
            }
          
          
          if (method.coef == "joint") {  
            z.star <- c(z, prior.mean)
            x.prior <- diag(NCOL(x)) 
            colnames(x.prior) <- colnames(x)
            x.star <- rbind(x, x.prior)
            w.star <- c(w, 1/(prior.sd + 1e-04) )
#            good.star <- c(good, rep(TRUE, NCOL(x.prior)))
#            fit <- qr(x.star[good.star, , drop = FALSE] * w.star, tol = min(1e-07, control$epsilon/1000))
            fit <- qr(x.star * w.star, tol = min(1e-07, control$epsilon/1000))
            fit$coefficients <- qr.coef(fit, z.star * w.star)
            coefs.hat <- fit$coefficients
            if (any(!is.finite(fit$coefficients))) {
              conv <- FALSE
              warning("non-finite coefficients at iteration ", iter)
              break
            }
          }  
           
          if (method.coef != "joint") {  
            for (j in 1:length(group0)) {
              vars <- c(covars0, group0[[j]])
              if (iter <= 5 | any((abs(coefs.hat[vars] - prior.mean[vars])) > 1e-03)) { 
                if (iter > 5) vars <- vars[abs(coefs.hat[vars] - prior.mean[vars]) > 1e-03]
                x0 <- x[, vars, drop = FALSE]
                eta0 <- x0 %*% coefs.hat[vars]
                z0 <- z - (eta - eta0) + offset  
                z0.star <- c(z0, prior.mean[vars])
                x0.prior <- diag(NCOL(x0)) 
                colnames(x0.prior) <- vars
                x0.star <- rbind(x0, x0.prior)
                w0.star <- c(w, 1/(prior.sd[vars] + 1e-04)) 
#                 good.star <- c(good, rep(TRUE, NCOL(x0.prior)))
#                 fit <- qr(x0.star[good.star, , drop = FALSE] * w0.star, tol = min(1e-07, control$epsilon/1000))
                fit <- qr(x0.star * w0.star, tol = min(1e-07, control$epsilon/1000))
                coefs.hat[vars] <- qr.coef(fit, z0.star * w0.star)
                eta <- eta - eta0 + x0 %*% coefs.hat[vars]
              } 
            }
            fit$coefficients <- coefs.hat
          }
            
#            start[fit$pivot] <- fit$coefficients
            start <- fit$coefficients 
            eta <- drop(x %*% start)
            mu <- linkinv(eta <- eta + offset)
            dev <- sum(dev.resids(y, mu, weights))
            if (control$trace)
                cat("Deviance =", dev, "Iterations -", iter,
                  "\n")
            boundary <- FALSE
            if (!is.finite(dev)) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated due to divergence",
                  call. = FALSE)
                ii <- 1
                while (!is.finite(dev)) {
                  if (ii > control$maxit)
                    stop("inner loop 1; cannot correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                  dev <- sum(dev.resids(y, mu, weights))
                }
                boundary <- TRUE
                if (control$trace)
                  cat("Step halved: new deviance =", dev, "\n")
            }
            if (!(valideta(eta) && validmu(mu))) {
                if (is.null(coefold))
                  stop("no valid set of coefficients has been found: please supply starting values",
                    call. = FALSE)
                warning("step size truncated: out of bounds",
                  call. = FALSE)
                ii <- 1
                while (!(valideta(eta) && validmu(mu))) {
                  if (ii > control$maxit)
                    stop("inner loop 2; cannot correct step size")
                  ii <- ii + 1
                  start <- (start + coefold)/2
                  eta <- drop(x %*% start)
                  mu <- linkinv(eta <- eta + offset)
                }
                boundary <- TRUE
                dev <- sum(dev.resids(y, mu, weights))
                if (control$trace)
                  cat("Step halved: new deviance =", dev, "\n")
            }
            
            if ( !(family[[1]] %in% c("binomial", "poisson")) & !nb ){
              Sum <- sum((w * (z - (eta - offset)[good]))^2) + sum((coefs.hat - prior.mean)^2/(prior.sd^2 + 1e-04)) 
              n.df0 <- nobs
              n.df <- n.df0 - length(coefs.hat[prior.sd >= 1e+04])
              if (n.df <= 0) n.df <- n.df0
              dispersion <- Sum/n.df 
            } 
            dispersion <- ifelse(dispersion > 1e+04,1e+04, dispersion) 
            dispersion <- ifelse(dispersion < 1e-04, 1e-04, dispersion)
            
            if(nb)  # for negative binomial model
            {
              th <- suppressWarnings( theta.ml(y=y, mu=mu, n=sum(weights), weights=weights, limit=10, trace=FALSE) )
              if (is.null(th)) th <- family$theta 
              family <- NegBin(theta = th)
              
              variance <- family$variance
              dev.resids <- family$dev.resids
              aic <- family$aic
              linkinv <- family$linkinv
              mu.eta <- family$mu.eta
              valideta <- family$valideta
              if (is.null(valideta))
                valideta <- function(eta) TRUE
              validmu <- family$validmu
              if (is.null(validmu))
                validmu <- function(mu) TRUE
            }

            if (iter > 2 & abs(dev - devold)/(0.1 + abs(dev)) < control$epsilon) {
                conv <- TRUE
                coef <- start
                break
            }
            else {                
                devold <- dev
                dispersionold <- dispersion
                coef <- coefold <- start
            }
            
        }  # iter end
        
        nvars <- ncol(fit$qr)
        
        if (Warning) {
          if (!conv)
              warning("algorithm did not converge", call. = FALSE)
          if (boundary)
              warning("algorithm stopped at boundary value", call. = FALSE)
          eps <- 10 * .Machine$double.eps
          if (family$family == "binomial") {
              if (any(mu > 1 - eps) || any(mu < eps))
                  warning("fitted probabilities numerically 0 or 1 occurred", call. = FALSE)
          }
          if (family$family == "poisson" | nb) {
              if (any(mu < eps))
                  warning("fitted rates numerically 0 occurred", call. = FALSE)
          }
        }
        if (fit$rank < nvars)
            coef[fit$pivot][seq(fit$rank + 1, nvars)] <- NA
        xxnames <- xnames[fit$pivot]
        residuals <- rep.int(NA, nobs)
        residuals[good] <- z - (eta - offset)[good]
        fit$qr <- as.matrix(fit$qr)
        nr <- min(sum(good), nvars)
        if (nr < nvars) {
            Rmat <- diag(nvars)
            Rmat[1:nr, 1:nvars] <- fit$qr[1:nr, 1:nvars]
        }
        else Rmat <- fit$qr[1:nvars, 1:nvars]
        Rmat <- as.matrix(Rmat)
        Rmat[row(Rmat) > col(Rmat)] <- 0
        names(coef) <- xnames
        colnames(fit$qr) <- xxnames
        dimnames(Rmat) <- list(xxnames, xxnames)
    } # end
    names(residuals) <- ynames
    names(mu) <- ynames
    names(eta) <- ynames
    wt <- rep.int(0, nobs)
    wt[good] <- w^2
    names(wt) <- ynames
    names(weights) <- ynames
    names(y) <- ynames
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else linkinv(offset)
    nulldev <- sum(dev.resids(y, wtdmu, weights))
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok
    if (all(prior.sd >= 1e+04)) nulldf <- n.ok - as.integer(intercept)
    rank <- if (EMPTY)
        0
    else fit$rank
    if (method.coef != "joint") rank <- ncol(x)
    resdf <- n.ok
    if (all(prior.sd >= 1e+04)) resdf <- n.ok - rank    
    aic.model <- aic(y, n, mu, weights, dev) + 2 * rank
    
    loglik <- -(aic.model - 2 * rank)/2

    if (intercept) {
      prior.mean <- prior.mean[-1]
      prior.scale <- prior.scale[-1]
      prior.df <- prior.df[-1]
    }
    
    out <- list(coefficients = coef, residuals = residuals, fitted.values = mu,
                effects = if (!EMPTY) fit$effects, R = if (!EMPTY) Rmat,
                rank = rank, qr = if (!EMPTY) structure(fit[c("qr", "rank", "qraux", "pivot")], class = "qr"), 
                linear.predictors = eta, deviance = dev, aic = aic.model, loglik = loglik,
                null.deviance = nulldev, iter = iter, weights = wt, prior.weights = weights,
                df.residual = resdf, df.null = nulldf, y = y, z = z, converged = conv, boundary = boundary, 
                intercept = intercept, 
                prior.sd = prior.sd, dispersion = dispersion, group = group, group.vars = group.vars, 
                ungroup.vars = ungroup.vars, method.coef = method.coef, family = family )
    
    if (prior == "t") 
      out$prior <- list(prior=prior, mean=prior.mean, scale=prior.scale, df=prior.df)
    if (prior == "de") 
      out$prior <- list(prior=prior, mean=prior.mean, scale=prior.scale)
    if (prior == "mde" | prior == "mt") {
      out$prior.scale <- prior.scale
      out$p <- p
      out$ptheta <- theta
      if (prior == "mde")
        out$prior <- list(prior=prior, mean=prior.mean, s0=ss[1], s1=ss[2], b=b)
      if (prior == "mt") 
        out$prior <- list(prior=prior, mean=prior.mean, s0=ss[1], s1=ss[2], df=prior.df, b=b)
      out$theta.weights <- theta.weights
    }
    if (nb){ 
      out$theta <- as.vector(th)
      out$SE.theta <- attr(th, "SE")
      out$twologlik <- 2 * loglik
    }
    
    return(out)
}

#*******************************************************************************

update.prior.sd <- function (prior, beta0, prior.scale, prior.df, sd.x, min.x.sd) 
{
  prior.scale <- prior.scale + 1e-04
  J <- length(beta0)
  if (prior == "t" | prior == "mt")   
    prior.sd <- sqrt((beta0^2 + prior.df * prior.scale^2)/(1 + prior.df)) 
  if (prior == "de" | prior == "mde")     # prior.scale = lamda in Exp(1/(2*lamda^2) )   
    prior.sd <- sqrt(abs(beta0) * prior.scale)
  
  prior.sd <- ifelse(prior.sd > 1e+04, 1e+04, prior.sd) 
  prior.sd <- ifelse(prior.sd < 1e-04, 1e-04, prior.sd)
  prior.sd <- ifelse(sd.x < min.x.sd, 1e-04, prior.sd)
  if (names(beta0)[1] == "(Intercept)") prior.sd[1] <- 1e+10
  prior.sd            
}

update.scale.p <- function(prior="mde", df=1, b0, ss, theta) 
{
  if (prior == "mde"){
    den0 <- (2 * ss[1])^(-1) * exp(-abs(b0)/ss[1]) # de density
    den1 <- (2 * ss[2])^(-1) * exp(-abs(b0)/ss[2]) 
  }
  if (prior == "mt"){
    den0 <- (ss[1])^(-1) * (1 + b0^2/(df * ss[1]^2))^(-(df + 1)/2) # t density
    den1 <- (ss[2])^(-1) * (1 + b0^2/(df * ss[2]^2))^(-(df + 1)/2)
  }
  p <- theta * den1 / (theta * den1 + (1 - theta) * den0 + 1e-10)
  scale <- 1/((1 - p)/ss[1] + p/ss[2] + 1e-10)
  
  list(scale=scale, p=p)
}

update.ptheta.group <- function(group.vars, p, w, b) # group-specific probability
{
  f <- function(theta, w, p, bb) { # theta ~ beta(1,b)  
    sum(p*log(w*theta) + (1-p)*log(1-w*theta)) + mean((bb-1)*log(1-theta))
  }
  theta <- p
  for (j in 1:length(group.vars)) {  
    vars <- group.vars[[j]]
#    theta[vars] <- mean(p[vars])  # posterior mode with theta~beta(1,1)
    theta[vars] <- optimize(f, interval=c(0, 1), 
                            w=w[vars], p=p[vars], bb=b[j], maximum=T)$maximum
  } 
  theta <- ifelse(theta < 0.01, 0.01, theta)
  theta <- ifelse(theta > 0.99, 0.99, theta)
  theta <- w * theta
  
  theta
}

update.ptheta.network <- function(theta, p, w) 
{
  phi <- 2
  for (j in 1:length(theta)) {  
    mu <- w %*% theta
    m <- mu[j] - w[j,j]*theta[j]
    a <- m*phi
    b <- (1-m)*phi
    theta[j] <- (p[j] + a)/(1 + a + b) # posterior mean
  } 
  theta <- ifelse(theta < 0.01, 0.01, theta)
  theta <- ifelse(theta > 0.99, 0.99, theta)
  
  theta
}

update.theta.weights <- function (gvars, theta.weights, inter.hierarchy, inter.parents, p)
{
  if (is.null(inter.parents)) 
    stop("'inter.parents' should be given")
  if (!is.list(inter.parents))
    stop("'inter.parents' should be a list") 
  xnames <- strsplit(gvars, split=":", fixed=T)
  inter <- unlist(lapply(xnames, function(x){length(x)}))
  if (length(inter.parents)!=length(inter[inter==2]))
    stop("interactions are not correctly specified in formula or inter.parents")
  
  p.main <- p[inter==1]
  if (inter.hierarchy=="strong")
    ww <- lapply(inter.parents, 
                function(x, p.main){ p.main[x[1]] * p.main[x[2]] }, 
                p.main)
  if (inter.hierarchy=="weak")
    ww <- lapply(inter.parents, 
                 function(x, p.main){ (p.main[x[1]] + p.main[x[2]])/2 }, 
                 p.main)
  theta.weights[inter==2] <- unlist(ww)
  theta.weights
}
  

# from MASS
NegBin <- function (theta=3, link="log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt"))
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",
                       linktemp))
  }
  .Theta <- theta
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", theta, envir = env)
  variance <- function(mu) mu + mu^2/.Theta
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1,
                                                           y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
  aic <- function(y, n, mu, wt, dev) {
    term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
      lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
      lgamma(.Theta + y)
    2 * sum(term * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    rnegbin(nsim * length(ftd), ftd, .Theta)
  }
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  famname <- paste("NegBin(", format(round(theta,
                                           4)), ")", sep = "")
  structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, theta = .Theta),
                 class = "family")
}

Student <- function(mean=0, scale=0.5, df=1, autoscale=TRUE)
{
  if (any(scale < 0)) stop("'scale' cannot be negative")
  if (any(df < 0)) stop("'df' cannot be negative")
  list(prior="t", mean=mean, scale=scale, df=df, autoscale=autoscale)
}
De <- function(mean=0, scale=0.5, autoscale=TRUE)
{
  if (any(scale < 0)) stop("'scale' cannot be negative")
  list(prior="de", mean=mean, scale=scale, autoscale=autoscale)
}
mde <- function(mean=0, s0=0.04, s1=0.5, b=1)
{
  if (s0 < 0 | s1 < 0) stop("scale cannot be negative")
  if (s0 > s1) stop("s0 should be smaller than s1")
  list(prior="mde", mean=mean, ss=c(s0,s1), b=b)
}
mt <- function(mean=0, s0=0.04, s1=0.5, df=1, b=1)
{
  if (s0 < 0 | s1 < 0) stop("scale cannot be negative")
  if (s0 > s1) stop("s0 should be smaller than s1")
  if (any(df < 0)) stop("'df' cannot be negative")
  list(prior="mt", mean=mean, ss=c(s0,s1), df=df, b=b)
}

auto_scale <- function(x, min.x.sd=1e-04)
{
  scale <- apply(x, 2, sd)
  scale <- ifelse(scale<=min.x.sd, 1, scale)
  two <- which(apply(x, 2, function(u) {length(unique(u))==2}))
  scale[two] <- apply(x[, two, drop=F], 2, function(u){max(u)-min(u)})
  scale
}


#************************************************************************************

prepare <- function(x, intercept, prior.mean, prior.sd, prior.scale, prior.df, group)
{
    x0 <- x
    if (intercept) x0 <- x[, -1, drop = FALSE] 
    g <- Grouping(all.var = colnames(x0), group = group)
    group <- g$group
    group.vars <- g$group.vars
    ungroup.vars <- g$ungroup.vars
    covars <- g$ungroup.vars  
    
    if (is.list(group)) { # for overlap groups
      if (length(unlist(group)) > length(unique(unlist(group)))) {
        x1 <- as.data.frame(x0)
        x1 <- x1[, c(covars, unlist(group))]
        g <- c(length(ungroup.vars), length(ungroup.vars) + cumsum(lapply(group, length)))
        for (j in 1:(length(group)-1))
          group.vars[[j]] <- colnames(x1[, (g[j]+1):g[j+1]])
        x1 <- as.matrix(x1)
        x <- x1 
        if (intercept) {
          x <- cbind(1, x)
          colnames(x)[1] <- "(Intercept)"
        }
      }
    }
    
    J <- NCOL(x)
    
    if (intercept & J > 1) {
      prior.mean <- c(0, prior.mean)
      prior.scale <- c(prior.scale[1], prior.scale)
      prior.df <- c(prior.df[1], prior.df)
    }
    
    if (length(prior.mean) < J) 
      prior.mean <- c(prior.mean, rep(prior.mean[length(prior.mean)], J - length(prior.mean)) )
    if (length(prior.scale) < J) 
      prior.scale <- c(prior.scale, rep(prior.scale[length(prior.scale)], J - length(prior.scale)) )
    if (length(prior.df) < J) 
      prior.df <- c(prior.df, rep(prior.df[length(prior.df)], J - length(prior.df)) )
    prior.mean <- prior.mean[1:J]
    prior.scale <- prior.scale[1:J]
    prior.df <- prior.df[1:J]
    prior.df <- ifelse(prior.df==Inf, 1e+10, prior.df)
    
    if (is.null(prior.sd)) prior.sd <- prior.scale + 0.2   ## + 0.2 to avoid prior.sd=0
    if (length(prior.sd) < J)  
      prior.sd <- c(prior.sd, rep(prior.sd[length(prior.sd)], J - length(prior.sd)) )
    prior.sd <- prior.sd[1:J]
    sd.x <- apply(x, 2, sd, na.rm=TRUE)
    min.x.sd <- 1e-04
    prior.sd <- ifelse(sd.x < min.x.sd, 1e-04, prior.sd)
    if (intercept) prior.sd[1] <- 1e+10 

    names(prior.mean) <- names(prior.scale) <- names(prior.df) <- names(prior.sd) <- colnames(x)

    if (intercept) covars <- c(colnames(x)[1], covars)
    if (!is.null(covars)) prior.mean[covars] <- 0
    
    list(x=x, prior.mean=prior.mean, prior.sd=prior.sd, prior.scale=prior.scale, prior.df=prior.df, 
         sd.x=sd.x, min.x.sd=min.x.sd,
         group=group, group.vars=group.vars, ungroup.vars=ungroup.vars)
}  

Grouping <- function(all.var, group) 
{ 
  n.vars <- length(all.var)
  group.vars <- list()
  
  if (is.matrix(group))
  {
    if (nrow(group)!=ncol(group) | ncol(group)>n.vars) 
      stop("wrong dimension for 'group'")
    if (any(rownames(group)!=colnames(group)))
      stop("rownames should be the same as colnames")
    if (any(!colnames(group)%in%all.var))
      stop("variabe names in 'group' not in the model predictors")
    group.vars <- colnames(group)
    group <- abs(group)
    wcol <- rowSums(group) - diag(group)
    group <- group/wcol
  }
  else{
    if (is.list(group)) group.vars <- group
    else
    {
      if (is.numeric(group) & length(group)>1) { 
        group <- sort(group)  
        if (group[length(group)] > n.vars) stop("wrong grouping")
      }
      if (is.numeric(group) & length(group)==1)
        group <- as.integer(seq(0, n.vars, length.out = n.vars/group + 1))
      if (is.null(group)) group <- c(0, n.vars)
      group <- unique(group)
      for (j in 1:(length(group) - 1))
        group.vars[[j]] <- all.var[(group[j] + 1):group[j + 1]]
    }
  }
  all.group.vars <- unique(unlist(group.vars))
  
  if (length(all.group.vars) == n.vars) ungroup.vars <- NULL
  else ungroup.vars <- all.var[which(!all.var %in% all.group.vars)]
  
  group.new <- c(length(ungroup.vars), length(ungroup.vars) + cumsum(lapply(group.vars, length)))
  var.new <- c(ungroup.vars, unlist(group.vars))
  
  list(group=group, group.vars=group.vars, ungroup.vars=ungroup.vars, 
       group.new=group.new, var.new=var.new) 
}



# only used in simulation
link.vars <- function(group.vars) {
  all.group.vars <- unique(unlist(group.vars))
  n.vars <- length(all.group.vars)
  n.groups <- length(group.vars)
  linked.vars <- vector(mode = "list", length = n.vars)
  names(linked.vars) <- all.group.vars
  for (i in 1:n.vars) {
    for (j in 1:n.groups)
      if (all.group.vars[i] %in% group.vars[[j]])
        linked.vars[[i]] <- unique(c(linked.vars[[i]], group.vars[[j]])) 
      d <- which(linked.vars[[i]] %in% all.group.vars[i])
      linked.vars[[i]] <- linked.vars[[i]][-d]
  }
  linked.vars
}

#*************************************************************************
# not used
truncated <- function(y, s.x, eta, Sd) #for extreme phenotype sampling
{
  s.x <- ifelse (s.x == -Inf, -1e+10, s.x)
  s.x <- ifelse (s.x == Inf, 1e+10, s.x)
  t1 <- (s.x[, 1] - eta)/Sd
  t2 <- (s.x[, 2] - eta)/Sd
  t3 <- (s.x[, 3] - eta)/Sd
  t4 <- (s.x[, 4] - eta)/Sd
  y.add <- (dnorm(t2) - dnorm(t1) + dnorm(t4) - dnorm(t3)) / (pnorm(t2) - pnorm(t1) + pnorm(t4) - pnorm(t3) + 1e-10) 
  v.add <- (t2 * dnorm(t2) - t1 * dnorm(t1) + t4 * dnorm(t4) - t3 * dnorm(t3))/(pnorm(t2) - pnorm(t1) + pnorm(t4) - pnorm(t3) + 1e-10)
  w <- sqrt(abs(1 - v.add - y.add^2))
  z <- eta + (y - eta + Sd * y.add)/(w^2 + 1e-10)
  out <- list(z = z, w = w)
  out
}  

#*************************************************************************
