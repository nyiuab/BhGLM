
#*******************************************************************************

cv.bh <- function(object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE)
{
  start.time <- Sys.time()
  
  if (any(class(object) %in% "glm")) 
    out <- cv.bh.glm(object = object, nfolds = nfolds, foldid = foldid, ncv = ncv, verbose = verbose)
  
  if (any(class(object) %in% "coxph")) 
    out <- cv.bh.coxph(object = object, nfolds = nfolds, foldid = foldid, ncv = ncv, verbose = verbose)
  
  if (any(class(object) %in% "glmNet") | any(class(object) %in% "bmlasso")) 
    out <- cv.bh.lasso(object = object, nfolds = nfolds, foldid = foldid, ncv = ncv, verbose = verbose)
  
  if (any(class(object) %in% "polr")) 
    out <- cv.bh.polr(object = object, nfolds = nfolds, foldid = foldid, ncv = ncv, verbose = verbose)
  
  stop.time <- Sys.time()
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  if(verbose) 
    cat("\n Cross-validation time:", Time, "minutes \n")
    
  out
}


### for bglm, glm
cv.bh.glm <- function(object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE)
{
  data.obj <- model.frame(object)
  x.obj <- model.matrix(object)
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  offset <- model.offset(data.obj)
  
  measures0 <- lp0 <- y.fitted0 <- foldid0 <- NULL
  fold <- foldid
  if (!is.null(foldid)) {
    fold <- as.matrix(foldid)
    nfolds <- max(foldid)
    ncv <- ncol(fold)
  }
  if (nfolds > n) nfolds <- n
  if (nfolds == n) ncv <- 1
  j <- 0
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv * nfolds, "models: \n")
  for (k in 1:ncv) {
    
    y.fitted <- lp <- rep(NA, n)
    deviance <- NULL
    
    if (!is.null(fold)) foldid <- fold[, k]
    else foldid <- sample(rep(seq(nfolds), length = n)) #sample(1:nfolds, size = n, replace = TRUE)
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid == i)
      subset1[omit] <- FALSE
      if (!is.null(object$prior.sd)) fit <- update(object, subset = subset1, verbose = FALSE)
      else fit <- update(object, subset = subset1) 
      lp[omit] <- x.obj[omit, , drop = FALSE] %*% fit$coefficients
      if (!is.null(offset)) lp[omit] <- lp[omit] + offset[omit]
      if (any(class(object) %in% "negbin")) 
        dd <- measure.nb(lp=lp[omit], y=y.obj[omit], theta=fit$nb.theta, linkinv=object$family$linkinv)
      else  
        dd <- measure.glm(lp=lp[omit], y=y.obj[omit], family=object$family, dispersion=fit$dispersion) 
      y.fitted[omit] <- dd$y.fitted
      deviance <- c(deviance, dd$measures["deviance"])
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
#        J <- nfolds * ncv
#        j <- j + 1
#        pre <- rep("\b", J)
#        cat(pre, "\n", j, "/", J, "\n", sep = "")
#        flush.console()
      }
    }
    
    if (any(class(object) %in% "negbin")) 
      measures <- measure.nb(lp=lp, y=y.obj, linkinv=object$family$linkinv)$measures
    else  
      measures <- measure.glm(lp=lp, y=y.obj, family=object$family)$measures 
    measures["deviance"] <- sum(deviance)
    
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
    y.fitted0 <- cbind(y.fitted0, y.fitted)
    foldid0 <- cbind(foldid0, foldid)
    
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$y.obs <- y.obj
  out$lp <- rowMeans(lp0, na.rm = TRUE)
  out$y.fitted <- rowMeans(y.fitted0, na.rm = TRUE)
  out$foldid <- foldid0
  if (ncv > 1){
    rownames(measures0) <- NULL
    out$detail <- list(measures = measures0, lp = lp0)
  }
  
  out
}

# for bcoxph, coxph
cv.bh.coxph <- function(object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE)
{
  data.obj <- model.frame(object)
  x.obj <- model.matrix(object)
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  offset <- model.offset(data.obj)
  
  measures0 <- lp0 <- foldid0 <- NULL
  fold <- foldid
  if (!is.null(foldid)) {
    fold <- as.matrix(foldid)
    nfolds <- max(foldid)
    ncv <- ncol(fold)
  }
  if (nfolds > n) nfolds <- n
  if (nfolds == n) ncv <- 1
  j <- 0
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv * nfolds, "models: \n")
  for (k in 1:ncv) {
    
    lp <- rep(NA, n)
    pl <- NULL
    
    if (!is.null(fold)) foldid <- fold[, k]
    else foldid <- sample(rep(seq(nfolds), length = n)) 
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid == i)
      subset1[omit] <- FALSE
      if (!is.null(object$prior.sd)) fit <- update(object, subset = subset1, verbose = FALSE)
      else fit <- update(object, subset = subset1)
      xb <- x.obj %*% fit$coefficients
      if (!is.null(offset)) xb <- xb + offset
      dd1 <- coxph(y.obj ~ xb, init = 1, control = coxph.control(iter.max=1), method = object$method)
      dd2 <- coxph(y.obj ~ xb, init = 1, control = coxph.control(iter.max=1), subset = subset1, method = object$method)
      lp[omit] <- xb[omit]        
      pl <- c(pl, dd1$loglik[1] - dd2$loglik[1])
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
  
    measures <- c(sum(pl), measure.cox(lp=lp, y=y.obj))
    names(measures) <- c("CVPL", "pl", "Cindex")
    
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
    foldid0 <- cbind(foldid0, foldid)
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$y.obs <- y.obj
  out$lp <- rowMeans(lp0, na.rm = TRUE)
  out$foldid <- foldid0
  if (ncv > 1){
    rownames(measures0) <- NULL
    out$detail <- list(measures = measures0, lp = lp0)
  }
  
  out
}

# for lasso, mlasso
cv.bh.lasso <- function(object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE)
{ 
  family <- object$family
  if (family=="gaussian") fa <- gaussian()
  if (family=="binomial") fa <- binomial()
  if (family=="poisson") fa <- poisson()
  
  x.obj <- object$x
  y.obj <- object$y
  n <- NROW(y.obj)
  offset <- object$offset
  if (!offset) offset <- NULL
  init <- object$coefficients
  init <- init[!names(init)%in%"(Intercept)"]
  
  measures0 <- lp0 <- y.fitted0 <- foldid0 <- NULL
  fold <- foldid
  if (!is.null(foldid)) {
    fold <- as.matrix(foldid)
    nfolds <- max(foldid)
    ncv <- ncol(fold)
  }
  if (nfolds > n) nfolds <- n
  if (nfolds == n) ncv <- 1
  j <- 0
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv * nfolds, "models: \n")
  for (k in 1:ncv) {
    
    y.fitted <- lp <- rep(NA, n)
    deviance <- pl <- NULL
    
    if (!is.null(fold)) foldid <- fold[, k]
    else foldid <- sample(rep(seq(nfolds), length = n)) #sample(1:nfolds, size = n, replace = TRUE)
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid == i)
      subset1[omit] <- FALSE
      if (any(class(object) %in% "glmNet"))
        fit <- update(object, x = x.obj[-omit, ], y = y.obj[-omit], weights = object$weights[-omit], offset = offset[-omit],
                      lambda = object$lambda, verbose = FALSE)
      if (any(class(object) %in% "bmlasso"))
        fit <- update(object, x = x.obj[-omit, ], y = y.obj[-omit], weights = object$weights[-omit], offset = offset[-omit], 
                      init = init, verbose = FALSE)
      if (any(class(object) %in% "GLM")) {
        lp[omit] <- cbind(1, x.obj[omit, , drop = FALSE]) %*% fit$coefficients
        if (!is.null(offset)) lp[omit] <- lp[omit] + offset[omit]
        dd <- measure.glm(lp=lp[omit], y=y.obj[omit], family=fa, dispersion=fit$dispersion) 
        y.fitted[omit] <- dd$y.fitted
        deviance <- c(deviance, dd$measures["deviance"])
      }
      if (any(class(object) %in% "COXPH")) {
        xb <- x.obj %*% fit$coefficients 
        if (!is.null(offset)) xb <- xb + offset
        dd1 <- coxph(y.obj ~ xb, init = 1, control = coxph.control(iter.max=1), method = "breslow")
        dd2 <- coxph(y.obj ~ xb, init = 1, control = coxph.control(iter.max=1), subset = subset1, method = "breslow")
        lp[omit] <- xb[omit]
        pl <- c(pl, dd1$loglik[1] - dd2$loglik[1])
      }
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    if (any(class(object) %in% "GLM")) {
      measures <- measure.glm(lp=lp, y=y.obj, family=fa)$measures 
      measures["deviance"] <- sum(deviance)
    }
    if (any(class(object) %in% "COXPH")) {
      measures <- c(sum(pl), measure.cox(lp=lp, y=y.obj))
      names(measures) <- c("CVPL", "pl", "Cindex")
    }
    
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
    y.fitted0 <- cbind(y.fitted0, y.fitted)
    foldid0 <- cbind(foldid0, foldid)
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$y.obs <- y.obj
  out$lp <- rowMeans(lp0, na.rm = TRUE)
  out$y.fitted <- rowMeans(y.fitted0, na.rm = TRUE)
  out$foldid <- foldid0
  
  if (ncv > 1){
    rownames(measures0) <- NULL
    out$detail <- list(measures = measures0, lp = lp0)
  }
  
  out
}


### for bpolr, polr
cv.bh.polr <- function(object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE)
{ 
  data.obj <- model.frame(object)
  x.obj <- data.obj[, -1, drop = FALSE]
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  
  measures0 <- lp0 <- foldid0 <- NULL
  y.fitted0 <- list()
  fold <- foldid
  if (!is.null(foldid)) {
    fold <- as.matrix(foldid)
    nfolds <- max(foldid)
    ncv <- ncol(fold)
  }
  if (nfolds > n) nfolds <- n
  if (nfolds == n) ncv <- 1
  j <- 0
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv * nfolds, "models: \n")
  for (k in 1:ncv) {
    
    y.fitted <- array(0, c(n, length(levels(y.obj))))
    
    if (!is.null(fold)) foldid <- fold[, k]
    else foldid <- sample(rep(seq(nfolds), length = n)) #sample(1:nfolds, size = n, replace = TRUE)
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid == i)
      subset1[omit] <- FALSE
      if (!is.null(object$prior.scale)) fit <- update(object, subset=subset1, Hess=FALSE, verbose=FALSE)
      else fit <- update(object, subset=subset1, Hess=FALSE) 
      dd <- predict.bh(fit, new.x=x.obj[omit, , drop=FALSE], new.y=y.obj[omit])
      y.fitted[omit, ] <- dd$y.fitted
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    measures <- measure.polr(pred=y.fitted, y=y.obj)
    
    measures0 <- rbind(measures0, measures)
    y.fitted0[[k]] <- y.fitted
    foldid0 <- cbind(foldid0, foldid)
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$y.obs <- y.obj
  out$y.fitted <- array(0, c(n, length(levels(y.obj))))
  for (k in 1:ncv) out$y.fitted <- out$y.fitted + y.fitted0[[k]]/ncv 
  out$foldid <- foldid0
  if (ncv > 1) out$detail <- measures0
  
  out
}

#***********************************************************************************

