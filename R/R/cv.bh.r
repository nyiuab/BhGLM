
#*******************************************************************************

cv.bh <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{
  start.time <- Sys.time()
  
  if (!"gam" %in% class(object))
  {
    if (any(class(object) %in% "glm")) 
      out <- cv.bh.glm(object=object, nfolds=nfolds, foldid=foldid, ncv=ncv, verbose=verbose)
  
    if (any(class(object) %in% "coxph")) 
      out <- cv.bh.coxph(object=object, nfolds=nfolds, foldid=foldid, ncv=ncv, verbose=verbose)
  
    if (any(class(object) %in% "glmNet") | any(class(object) %in% "bmlasso")) 
      out <- cv.bh.lasso(object=object, nfolds=nfolds, foldid=foldid, ncv=ncv, verbose=verbose)
  
    if (any(class(object) %in% "polr")) 
      out <- cv.bh.polr(object=object, nfolds=nfolds, foldid=foldid, ncv=ncv, verbose=verbose)
  }
  else
  {
    fam <- object$family$family 
    if (substr(fam, 1, 17) == "Negative Binomial") fam <- "NegBin"
    gam.fam <- c("gaussian", "binomial", "poisson", "quasibinomial", "quasipoisson", "NegBin", 
                 "Cox PH")
    if (! fam %in% gam.fam)
      stop("Cross-validation for this family has not been implemented yet")
    if (fam %in% gam.fam[1:6])
      out <- cv.gam.glm(object=object, nfolds=nfolds, foldid=foldid, ncv=ncv, verbose=verbose)
    if (fam == "Cox PH")
      out <- cv.gam.coxph(object=object, nfolds=nfolds, foldid=foldid, ncv=ncv, verbose=verbose)
  }
  
  stop.time <- Sys.time()
  Time <- round(difftime(stop.time, start.time, units = "min"), 3)
  if(verbose) 
    cat("\n Cross-validation time:", Time, "minutes \n")
    
  out
}

generate.foldid <- function(nobs, nfolds=10, foldid=NULL, ncv=1)
{
  if (nfolds > nobs) nfolds <- nobs
  if (nfolds == nobs) ncv <- 1
  if (is.null(foldid)) {
   foldid <- array(NA, c(nobs, ncv)) 
   for(j in 1:ncv) 
     foldid[, j] <- sample(rep(seq(nfolds), length=nobs))
  }
  foldid <- as.matrix(foldid)
  nfolds <- max(foldid)
  ncv <- ncol(foldid)
  
  list(foldid=foldid, nfolds=nfolds, ncv=ncv)
}
  

### for bglm, glm
cv.bh.glm <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{
  data.obj <- model.frame(object)
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  fol <- generate.foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- lp0 <- y.fitted0 <- NULL
  j <- 0
  if (!is.null(object$offset)) {
    data.obj <- object$data
    if (is.null(object$data)) stop("'data' not given in object")
  }
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv*nfolds, "models: \n")
  for (k in 1:ncv) {
    y.fitted <- lp <- rep(NA, n)
    deviance <- NULL
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      fit <- update(object, subset = subset1)
      lp[omit] <- predict(fit, newdata=data.obj[omit, , drop=FALSE])
      y.fitted[omit] <- object$family$linkinv(lp[omit])
      if (any(class(object) %in% "negbin")) fit$dispersion <- fit$theta
      dd <- suppressWarnings( measure.glm(y.obj[omit], y.fitted[omit], family=object$family$family, dispersion=fit$dispersion) ) 
      deviance <- c(deviance, dd["deviance"])
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    measures <- measure.glm(y.obj, y.fitted, family=object$family$family) 
    measures["deviance"] <- sum(deviance)
    
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
    y.fitted0 <- cbind(y.fitted0, y.fitted)
    
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$measures <- round(out$measures, digits=3)
  out$y.obs <- y.obj
  out$lp <- lp0
  out$y.fitted <- y.fitted0
  out$foldid <- foldid
  
  out
}

# for bcoxph, coxph
cv.bh.coxph <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{
  data.obj <- model.frame(object)
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  fol <- generate.foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- lp0 <- NULL
  j <- 0
  if (!is.null(object$offset)) {
    data.obj <- object$data
    if (is.null(object$data)) stop("'data' not given in object")
  }
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv*nfolds, "models: \n")
  for (k in 1:ncv) {
    lp <- rep(NA, n)
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      fit <- update(object, subset = subset1)
      lp[omit] <- predict(fit, newdata=data.obj[omit, , drop=FALSE])
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
  
    measures <- measure.cox(y.obj, lp)
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$measures <- round(out$measures, digits=3)
  out$y.obs <- y.obj
  out$lp <- lp0
  out$foldid <- foldid
  
  out
}

# for lasso, mlasso
cv.bh.lasso <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{ 
  family <- object$family
  x.obj <- object$x
  y.obj <- object$y
  n <- NROW(y.obj)
  offset <- object$offset
  init <- object$coefficients
  init <- init[!names(init)%in%"(Intercept)"]
  
  fol <- generate.foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- lp0 <- y.fitted0 <- NULL
  j <- 0
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv*nfolds, "models: \n")
  for (k in 1:ncv) {
    y.fitted <- lp <- rep(NA, n)
    deviance <- NULL
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      if (any(class(object) %in% "glmNet"))
        fit <- update(object, x=x.obj[-omit, ], y=y.obj[-omit], offset=offset[-omit],
                      lambda=object$lambda, verbose=FALSE)
      if (any(class(object) %in% "bmlasso"))
        fit <- update(object, x=x.obj[-omit, ], y=y.obj[-omit], offset=offset[-omit], 
                      init=init, verbose=FALSE)
      if (is.null(fit$offset)) fit$offset <- FALSE
      else fit$offset <- TRUE
      xx <- x.obj[omit, , drop=FALSE]
      off <- offset[omit]
      lp[omit] <- as.vector(predict(fit, newx=xx, newoffset=off))
      if (any(class(object) %in% "GLM")) {
        y.fitted[omit] <- as.vector(predict(fit, newx=xx, type="response", newoffset=off))
        dd <- suppressWarnings( measure.glm(y.obj[omit], y.fitted[omit], family=family, dispersion=fit$dispersion) )
        deviance <- c(deviance, dd["deviance"])
      }
       
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    if (any(class(object) %in% "GLM")) {
      measures <- measure.glm(y.obj, y.fitted, family=family)
      measures["deviance"] <- sum(deviance)
      y.fitted0 <- cbind(y.fitted0, y.fitted)
    }
    if (any(class(object) %in% "COXPH")) 
      measures <- measure.cox(y.obj, lp)
    
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$measures <- round(out$measures, digits=3)
  out$y.obs <- y.obj
  out$lp <- lp0
  if (any(class(object) %in% "GLM")) out$y.fitted <- y.fitted0
  out$foldid <- foldid

  out
}


### for bpolr, polr
cv.bh.polr <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{ 
  data.obj <- model.frame(object)
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  fol <- generate.foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- NULL
  y.fitted0 <- list()
  j <- 0
  if (!is.null(object$offset)) {
    data.obj <- object$data
    if (is.null(object$data)) stop("'data' not given in object")
  }

  if (verbose) cat("Fitting", "ncv*nfolds =", ncv * nfolds, "models: \n")
  for (k in 1:ncv) {
    y.fitted <- array(0, c(n, length(levels(y.obj))))
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      fit <- update(object, subset=subset1, Hess=FALSE) 
      y.fitted[omit, ] <- predict(fit, newdata=data.obj[omit, , drop=FALSE], type="probs")
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    measures <- measure.polr(y.obj, y.fitted)
    
    measures0 <- rbind(measures0, measures)
    y.fitted0[[k]] <- y.fitted
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$measures <- round(out$measures, digits=3)
  out$y.obs <- y.obj
#  out$y.fitted <- array(0, c(n, length(levels(y.obj))))
#  for (k in 1:ncv) out$y.fitted <- out$y.fitted + y.fitted0[[k]]/ncv 
  out$y.fitted <- y.fitted0
  out$foldid <- foldid
  
  out
}

#***********************************************************************************
# for gam from mgcv

cv.gam.glm <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{
  data.obj <- model.frame(object)
  y.obj <- model.response(data.obj)
  n <- NROW(y.obj)
  fol <- generate.foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- lp0 <- y.fitted0 <- NULL
  j <- 0
  if (!is.null(object$offset)) {
    data.obj <- object$data
    if (is.null(object$data)) stop("'data' not given in object")
  }
  
  fam <- object$family$family
  if (substr(object$family$family, 1, 17) == "Negative Binomial")
    fam <- "NegBin"
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv*nfolds, "models: \n")
  for (k in 1:ncv) {
    y.fitted <- lp <- rep(NA, n)
    deviance <- NULL
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      fit <- update(object, subset = subset1)
      lp[omit] <- predict(fit, newdata=data.obj[omit, , drop=FALSE])
      y.fitted[omit] <- object$family$linkinv(lp[omit])
      disp <- fit$sig2
      if (fam[[1]] == "NegBin") disp <- fit$family$getTheta(TRUE)
      dd <- suppressWarnings( measure.glm(y.obj[omit], y.fitted[omit], family=fam, dispersion=disp) ) 
      deviance <- c(deviance, dd["deviance"])
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    measures <- measure.glm(y.obj, y.fitted, family=fam) 
    measures["deviance"] <- sum(deviance)
    
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
    y.fitted0 <- cbind(y.fitted0, y.fitted)
    
  }
  
  out <- list()
  if (nrow(measures0)==1) out$measures <- colMeans(measures0, na.rm=TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm=TRUE), apply(measures0, 2, sd, na.rm=TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$measures <- round(out$measures, digits=3)
  out$y.obs <- y.obj
  out$lp <- lp0
  out$y.fitted <- y.fitted0
  out$foldid <- foldid
  
  out
}

cv.gam.coxph <- function(object, nfolds=10, foldid=NULL, ncv=1, verbose=TRUE)
{
  require(survival)
  data.obj <- model.frame(object)
  y.obj <- Surv(model.response(data.obj), object$prior.weights)
  n <- NROW(y.obj)
  fol <- generate.foldid(nobs=n, nfolds=nfolds, foldid=foldid, ncv=ncv)
  foldid <- fol$foldid
  nfolds <- fol$nfolds
  ncv <- fol$ncv
  measures0 <- lp0 <- NULL
  j <- 0
  if (!is.null(object$offset)) {
    data.obj <- object$data
    if (is.null(object$data)) stop("'data' not given in object")
  }
  
  if (verbose) cat("Fitting", "ncv*nfolds =", ncv*nfolds, "models: \n")
  for (k in 1:ncv) {
    lp <- rep(NA, n)
    
    for (i in 1:nfolds) {
      subset1 <- rep(TRUE, n)
      omit <- which(foldid[, k] == i)
      subset1[omit] <- FALSE
      fit <- update(object, subset=subset1)
      lp[omit] <- predict(fit, newdata=data.obj[omit, , drop=FALSE])
      
      if (verbose) {
        j <- j + 1
        cat(j, "")
      }
    }
    
    measures <- measure.cox(y.obj, lp)
    measures0 <- rbind(measures0, measures)
    lp0 <- cbind(lp0, lp)
  }
  
  out <- list()
  if (nrow(measures0) == 1) out$measures <- colMeans(measures0, na.rm = TRUE)
  else {
    out$measures <- rbind(colMeans(measures0, na.rm = TRUE), apply(measures0, 2, sd, na.rm = TRUE))
    rownames(out$measures) <- c("mean", "sd")
  }
  out$measures <- round(out$measures, digits=3)
  out$y.obs <- y.obj
  out$lp <- lp0
  out$foldid <- foldid
  
  out
}

