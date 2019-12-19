
bcoxph <- function (formula, data, weights, subset, na.action, init, 
                    control = coxph.control(eps = 1e-04, iter.max = 50), ties = c("breslow", "efron"), tt,  
                    prior = Student(0, 0.5, 1), group = NULL, method.coef, 
                    Warning = FALSE, verbose = FALSE, ...) 
{
  if (!requireNamespace("survival")) install.packages("survival")
    require(survival)
    start.time <- Sys.time()
    prior.mean <- prior$mean
    prior.scale <- prior$scale
    if (is.null(prior.scale)) prior.scale <- 0.5
    prior.df <- prior$df
    if (is.null(prior.df)) prior.df <- 1
    ss <- prior$ss
    if (is.null(ss)) ss <- c(0.04, 0.5)
    prior <- prior[[1]]
    prior.sd <- 0.5 
    if (missing(method.coef)) method.coef <- NULL 
    
    robust <- FALSE
    ties <- ties[1]
#    ties <- match.arg(ties)
    Call <- match.call()
    extraArgs <- list(...)
    if (length(extraArgs)) {
      controlargs <- names(formals(coxph.control))
      indx <- pmatch(names(extraArgs), controlargs, nomatch = 0L)
      if (any(indx == 0L)) 
        stop(gettextf("Argument %s not matched", names(extraArgs)[indx == 
                                                                    0L]), domain = NA)
    }
    if (missing(control)) control <- coxph.control(...)
    indx <- match(c("formula", "data", "weights", "subset", "na.action"), 
        names(Call), nomatch = 0)
    if (indx[1] == 0) 
        stop("A formula argument is required")
    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")
    special <- c("strata", "cluster", "tt")
    temp$formula <- if (missing(data)) 
        terms(formula, special)
    else terms(formula, special, data = data)
    if (!is.null(attr(temp$formula, "specials")$tt)) {
        coxenv <- new.env(parent = environment(formula))
        assign("tt", function(x) x, env = coxenv)
        environment(temp$formula) <- coxenv
    }
    mf <- eval(temp, parent.frame())
    if (nrow(mf) == 0) 
        stop("No (non-missing) observations")
    Terms <- terms(mf)
    Y <- model.extract(mf, "response")
    if (!inherits(Y, "Surv")) 
        stop("Response must be a survival object")
    type <- attr(Y, "type")
    if (type != "right" && type != "counting") 
        stop(paste("Cox model doesn't support \"", type, "\" survival data", 
            sep = ""))
    data.n <- nrow(Y)
    if (control$timefix) Y <- aeqSurv(Y)

    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
        stemp <- untangle.specials(Terms, "strata", 1)
        if (length(stemp$vars) == 1) 
            strata.keep <- mf[[stemp$vars]]
        else strata.keep <- strata(mf[, stemp$vars], shortlabel = TRUE)
        strats <- as.numeric(strata.keep)
    }
    timetrans <- attr(Terms, "specials")$tt
    if (missing(tt)) 
        tt <- NULL
    if (length(timetrans)) {
        timetrans <- untangle.specials(Terms, "tt")
        ntrans <- length(timetrans$terms)
        if (is.null(tt)) {
            tt <- function(x, time, riskset, weights) {
                obrien <- function(x) {
                  r <- rank(x)
                  (r - 0.5)/(0.5 + length(r) - r)
                }
                unlist(tapply(x, riskset, obrien))
            }
        }
        if (is.function(tt)) 
            tt <- list(tt)
        if (is.list(tt)) {
            if (any(!sapply(tt, is.function))) 
                stop("The tt argument must contain function or list of functions")
            if (length(tt) != ntrans) {
                if (length(tt) == 1) {
                  temp <- vector("list", ntrans)
                  for (i in 1:ntrans) temp[[i]] <- tt[[1]]
                  tt <- temp
                }
                else stop("Wrong length for tt argument")
            }
        }
        else stop("The tt argument must contain a function or list of functions")
        if (ncol(Y) == 2) {
            if (length(strats) == 0) {
                sorted <- order(-Y[, 1], Y[, 2])
                newstrat <- rep.int(0L, nrow(Y))
                newstrat[1] <- 1L
            }
            else {
                sorted <- order(strats, -Y[, 1], Y[, 2])
                newstrat <- as.integer(c(1, 1 * (diff(strats[sorted]) != 
                  0)))
            }
            if (storage.mode(Y) != "double") 
                storage.mode(Y) <- "double"
            counts <- .Call(Ccoxcount1, Y[sorted, ], as.integer(newstrat))
            tindex <- sorted[counts$index]
        }
        else {
            if (length(strats) == 0) {
                sort.end <- order(-Y[, 2], Y[, 3])
                sort.start <- order(-Y[, 1])
                newstrat <- c(1L, rep(0, nrow(Y) - 1))
            }
            else {
                sort.end <- order(strats, -Y[, 2], Y[, 3])
                sort.start <- order(strats, -Y[, 1])
                newstrat <- c(1L, as.integer(diff(strats[sort.end]) != 
                  0))
            }
            if (storage.mode(Y) != "double") 
                storage.mode(Y) <- "double"
            counts <- .Call(Ccoxcount2, Y, as.integer(sort.start - 
                1L), as.integer(sort.end - 1L), as.integer(newstrat))
            tindex <- counts$index
        }
        Y <- Surv(rep(counts$time, counts$nrisk), counts$status)
        type <- "right"
        mf <- mf[tindex, ]
        strats <- rep(1:length(counts$nrisk), counts$nrisk)
        weights <- model.weights(mf)
        if (!is.null(weights) && any(!is.finite(weights))) 
            stop("weights must be finite")
        tcall <- attr(Terms, "variables")[timetrans$terms + 2]
        pvars <- attr(Terms, "predvars")
        pmethod <- sub("makepredictcall.", "", as.vector(methods("makepredictcall")))
        for (i in 1:ntrans) {
          newtt <- (tt[[i]])(mf[[timetrans$var[i]]], Y[, 1], 
                             strats, weights)
          mf[[timetrans$var[i]]] <- newtt
          nclass <- class(newtt)
          if (any(nclass %in% pmethod)) {
            dummy <- as.call(list(as.name(class(newtt)[1]), 
                                  tcall[[i]][[2]]))
            ptemp <- makepredictcall(newtt, dummy)
            pvars[[timetrans$terms[i] + 2]] <- ptemp
          }
        }
        attr(Terms, "predvars") <- pvars
    }
    cluster <- attr(Terms, "specials")$cluster
    if (length(cluster)) {
        robust <- TRUE
        tempc <- untangle.specials(Terms, "cluster", 1:10)
        ord <- attr(Terms, "order")[tempc$terms]
        if (any(ord > 1)) 
            stop("Cluster can not be used in an interaction")
        cluster <- strata(mf[, tempc$vars], shortlabel = TRUE)
        Terms <- Terms[-tempc$terms]
        xlevels <- .getXlevels(Terms[-tempc$terms], mf)
    }
    else {
        if (missing(robust)) 
            robust <- FALSE
        xlevels <- .getXlevels(Terms, mf)
    }
    contrast.arg <- NULL
    attr(Terms, "intercept") <- 1
    adrop <- 0
    dropterms <- NULL
    stemp <- untangle.specials(Terms, "strata", 1)
    if (length(stemp$vars) > 0) {
        hasinteractions <- FALSE
        for (i in stemp$vars) {
            if (any(attr(Terms, "order")[attr(Terms, "factors")[i, 
                ] > 0] > 1)) 
                hasinteractions <- TRUE
        }
        if (!hasinteractions) 
            dropterms <- c(dropterms, stemp$terms)
        else adrop <- c(0, match(stemp$var, colnames(attr(Terms, 
            "factors"))))
    }
    if (length(dropterms)) {
        temppred <- attr(terms, "predvars")
        Terms2 <- Terms[-dropterms]
        if (!is.null(temppred)) {
            attr(Terms2, "predvars") <- temppred[-(1 + dropterms)]
        }
        X <- model.matrix(Terms2, mf, constrasts = contrast.arg)
        renumber <- match(colnames(attr(Terms2, "factors")), 
            colnames(attr(Terms, "factors")))
        attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
    }
    else X <- model.matrix(Terms, mf, contrasts = contrast.arg)
    Xatt <- attributes(X)
    xdrop <- Xatt$assign %in% adrop
    X <- X[, !xdrop, drop = FALSE]
    attr(X, "assign") <- Xatt$assign[!xdrop]
    attr(X, "contrasts") <- Xatt$contrasts
    offset <- model.offset(mf)
    if (is.null(offset) | all(offset == 0)) 
        offset <- rep(0, nrow(mf))
    else if (any(!is.finite(offset))) 
        stop("offsets must be finite")
    weights <- model.weights(mf)
    if (!is.null(weights) && any(!is.finite(weights))) 
        stop("weights must be finite")
    assign <- attrassign(X, Terms)
    contr.save <- attr(X, "contrasts")

    if (missing(init))
      init <- rep(0, NCOL(X))
    else {
      if (length(init) != NCOL(X)) 
        stop("wrong length for init argument")
      temp <- X %*% init - sum(colMeans(X) * init)
      if (any(temp < .Machine$double.min.exp | temp > .Machine$double.max.exp)) 
#        stop("initial values lead to overflow or underflow of the exp function")
        warning("initial values lead to overflow or underflow of the exp function")
    }

    fit <- bcoxph.fit(X, Y, offset = offset, weights = weights, init = init, 
                      control = control, strats = factor(strats),
                      ties = ties, prior = prior, group = group, method.coef = method.coef, 
                      prior.mean = prior.mean, prior.sd = prior.sd, 
                      prior.scale = prior.scale, prior.df = prior.df, ss = ss, Warning = Warning) 
    
    na.action <- attr(mf, "na.action")
    if (length(na.action)) 
      fit$na.action <- na.action
    if (length(timetrans)) {
      mf[[".surv."]] <- Y
      mf[[".strata."]] <- strats
      stop("Time transform + model frame: code incomplete")
    }
    fit$model <- mf
    fit$x <- X
    if (length(strats)) {
      if (length(timetrans)) 
        fit$strata <- strats
      else fit$strata <- strata.keep
    }

    fit$terms <- Terms
    fit$assign <- assign
    if (!is.null(weights) && any(weights != 1)) fit$weights <- weights
    fit$formula <- formula(Terms)
    if (length(xlevels) > 0) fit$xlevels <- xlevels
    fit$contrasts <- contr.save
    if (any(offset != 0)) fit$offset <- offset
    fit$call <- Call

    stop.time <- Sys.time()
    minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
    if (verbose) {
      cat("EM Newton-Raphson iterations:", fit$iter[2], "\n")
      cat("Computational time:", minutes, "minutes \n")
    }
    
    fit
}



bcoxph.fit <- function(x, y, offset = rep(0, nobs), weights = rep(1, nobs), init = 0, strats = NULL,
                       control = coxph.control(eps = 1e-04, iter.max = 50), ties = "efron", 
                       prior = "t", group = NULL, method.coef = NULL, 
                       prior.mean = 0, prior.sd = 1, prior.scale = 1, prior.df = 1, ss = c(0.05, 0.1), 
                       Warning = FALSE) 
{
  ss <- sort(ss)
  ss <- ifelse(ss <= 0, 0.001, ss)
  if (prior == "mde")
    prior.sd <- prior.scale <- ss[length(ss)]  # used for ungrouped coefficients
  
  d <- prepare(x = x, intercept = FALSE, prior.mean = prior.mean, prior.sd = prior.sd, prior.scale = prior.scale, 
               prior.df = prior.df, group = group)
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

  prior.scale <- prior.scale / autoscale(x, min.x.sd)
    
  g0 <- Grouping(all.var = colnames(x), group = method.coef)
  group0 <- g0$group.vars
  covars0 <- g0$ungroup.vars  
  method.coef <- "joint"
  if (length(group0) > 1) method.coef <- "group"
  
  # for mixture prior
  if (prior == "mde") {
    if (length(ss) != 2) stop("ss should have two positive values")
    gvars <- unlist(group.vars)
    theta <- p <- rep(0.5, length(gvars))
    names(theta) <- names(p) <- gvars
  }
  
  names(init) <- colnames(x)
  coefs.hat <- init
  eta <- x %*% coefs.hat
  
  means <- init
  var <- array(0, c(ncol(x), ncol(x)))
  rownames(var) <- colnames(var) <- colnames(x)
  var2 <- var
  
  offset0 <- offset
    
  devold <- -100
  conv <- FALSE
  for (iter in 1:control$iter.max){
    
    if (iter > 1) {
      beta0 <- coefs.hat - prior.mean
      
      if (prior == "mde") {
        out <- update.scale.p(b0=beta0[gvars], ss=ss, theta=theta)
        prior.scale[gvars] <- out[[1]]   
        p <- out[[2]]
        theta <- update.ptheta.group(group.vars=group.vars, p=p)
      }
      
      prior.sd <- update.prior.sd(prior = prior, beta0 = beta0, prior.scale = prior.scale, 
                                  prior.df = prior.df, sd.x = sd.x, min.x.sd = min.x.sd) 
    } 
       
    if (method.coef == "joint") {
      formula <- y ~ b.ridge(x, theta = 1/prior.sd^2, prior.mean = prior.mean)
      if (any(offset0 != 0))
        formula <- y ~ offset(offset0) + b.ridge(x, theta = 1/prior.sd^2, prior.mean = prior.mean)
      if (length(strats) > 1) { 
        formula <- y ~ b.ridge(x, theta = 1/prior.sd^2, prior.mean = prior.mean) + strata(strats)
        if (any(offset0 != 0))
          formula <- y ~ offset(offset0) + b.ridge(x, theta = 1/prior.sd^2, prior.mean = prior.mean) + strata(strats)
      }
      fit <- coxph(formula = formula, init = init, weights = weights, 
                   control = coxph.control(iter.max=1, outer.max=1), ties = ties, method = ties)
      init <- coefs.hat <- fit$coefficients
      names(init) <- names(coefs.hat) <- colnames(x)
      if (iter == 1) loglik0 <- fit$loglik[1]
    }

    if (method.coef != "joint") {  
      for (j in 1:length(group0)) {
        vars <- c(covars0, group0[[j]])
        if (iter <= 5 | any((abs(coefs.hat[vars] - prior.mean[vars])) > 1e-03)) { 
          if (iter > 5) vars <- vars[abs(coefs.hat[vars] - prior.mean[vars]) > 1e-03]
          x0 <- x[, vars, drop = FALSE]
          eta0 <- eta - x0 %*% coefs.hat[vars]
          formula <- y ~ offset(eta0) + b.ridge(x0, theta = 1/prior.sd[vars]^2, prior.mean = prior.mean[vars])
          if (any(offset0 != 0))
            formula <- y ~ offset(offset0 + eta0) + b.ridge(x0, theta = 1/prior.sd[vars]^2, prior.mean = prior.mean[vars])
          if (length(strats) > 1){ 
            formula <- y ~ offset(eta0) + b.ridge(x0, theta = 1/prior.sd[vars]^2, prior.mean = prior.mean[vars]) + strata(strats)
            if (any(offset0 != 0))
              formula <- y ~ offset(offset0 + eta0) + b.ridge(x0, theta = 1/prior.sd[vars]^2, prior.mean = prior.mean[vars]) + strata(strats)
          }
          fit <- coxph(formula = formula, init = init[vars], weights = weights,
                       control = coxph.control(iter.max=1, outer.max=1), ties = ties, method = ties)
          init[vars] <- coefs.hat[vars] <- fit$coefficients 
          eta <- eta0 + x0 %*% coefs.hat[vars]
          means[vars] <- fit$means
          var[vars, vars] <- fit$var
          var2[vars, vars] <- fit$var2
        } 
        if (iter == 1 & j == 1) loglik0 <- fit$loglik[1]
      }
      fit$coefficients <- coefs.hat
      fit$means <- means
      fit$var <- var
      fit$var2 <- var2
    }
          
    dev <- -2 * fit$loglik[2]
    if(abs(dev - devold)/(0.1 + abs(dev)) < control$eps) {
      conv <- TRUE
      break
    }
    else devold <- dev
  } # iteration end
  if (Warning & !conv) warning("algorithm did not converge", call. = FALSE)
    
  names(fit$coefficients) <- colnames(x)
    
  fit$iter[1] <- fit$iter[2] <- iter
  fit$loglik[1] <- loglik0
  fit$deviance <- dev
  fit$prior.sd <- prior.sd 
  fit$group <- group
  fit$group.vars <- group.vars
  fit$ungroup.vars <- ungroup.vars
  fit$method.coef <- method.coef
  
  if (prior == "t") 
    fit$prior <- list(prior="Stendent-t", mean=prior.mean, scale=prior.scale, df=prior.df)
  if (prior == "de") 
    fit$prior <- list(prior="Double-exponential", mean=prior.mean, scale=prior.scale)
  if (prior == "mde") {
    fit$p <- p
    fit$ptheta <- theta
    fit$prior <- list(prior="mixture double-exponential", mean=prior.mean, s0=ss[1], s1=ss[2])
  }
  
  fit
}


b.ridge <- function (..., theta, prior.mean)
{
    x <- cbind(...)
    nvar <- ncol(x)
    xname <- as.character(parse(text = substitute(cbind(...))))[-1]
    vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
    class(x) <- "coxph.penalty"
    scale <- FALSE
    prior.mean <- prior.mean
    
    pfun <- function(coef, theta, ndead, scale) {
              list(penalty = sum((coef - prior.mean)^2 * theta/2), first = theta *
                  (coef - prior.mean), second = theta, flag = FALSE)
            }
    
    temp <- list( pfun = pfun, diag = TRUE, 
              cfun = function(parms, iter, history) {
                        list(theta = parms$theta, done = TRUE)}, 
              cparm = list(theta = theta), pparm = vars, varname = paste("ridge(",
                           xname, ")", sep = "") )
    
    attributes(x) <- c(attributes(x), temp)
    
    x
}

#*******************************************************************************

b.ridge0 <- function (..., theta, df, eps = 0.1, scale = FALSE)
{
    x <- cbind(...)
    nvar <- ncol(x)
    xname <- as.character(parse(text = substitute(cbind(...))))[-1]
    vars <- apply(x, 2, function(z) var(z[!is.na(z)]))
    class(x) <- "coxph.penalty"
    if (!missing(theta) && !missing(df))
        stop("Only one of df or theta can be specified")
    if (scale)
        pfun <- function(coef, theta, ndead, scale) {
            list(penalty = sum(coef^2 * scale * theta/2), first = theta *
                coef * scale, second = theta * scale, flag = FALSE)
        }
    else pfun <- function(coef, theta, ndead, scale) {
        list(penalty = sum(coef^2 * theta/2), first = theta *
            coef, second = theta, flag = FALSE)
    }
    if (!missing(theta)) {
        temp <- list(pfun = pfun, diag = TRUE, cfun = function(parms,
            iter, history) {
            list(theta = parms$theta, done = TRUE)
        }, cparm = list(theta = theta), pparm = vars, varname = paste("ridge(",
            xname, ")", sep = ""))
    }
    else {
        temp <- list(pfun = pfun, diag = TRUE, cfun = frailty.controldf,
            cargs = "df", cparm = list(df = df, eps = eps, thetas = 0,
                dfs = nvar, guess = 1), pparm = vars, varname = paste("ridge(",
                xname, ")", sep = ""))
    }
    attributes(x) <- c(attributes(x), temp)
    x
}

#*******************************************************************************
