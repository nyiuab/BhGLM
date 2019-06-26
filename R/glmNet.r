
#*******************************************************************************

glmNet <- function (x, y, family = c("gaussian", "binomial", "poisson", "cox"), offset = NULL,
                    alpha = c(1, 0.5, 0), lambda, penalty.factor = rep(1, ncol(x)), 
                    nfolds = 10, ncv = 10, verbose = FALSE)
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
  alpha <- alpha[1]
  if (length(penalty.factor) != ncol(x)) stop("give each predictor a penalty")
  names(penalty.factor) <- colnames(x)
  penalty.factor <- ncol(x) * penalty.factor/sum(penalty.factor) 
  
  if (nfolds == NROW(x)) ncv <- 1 

  if (missing(lambda)) {
    lam <- rep(NA, ncv)
    for (i in 1:ncv) { 
      cv.f <- cv.glmnet(x = x, y = y, family = family, offset = offset, alpha = alpha, nfolds = nfolds, grouped = TRUE,
                        penalty.factor = penalty.factor, standardize = FALSE)
      lam[i] <- cv.f$lambda.min
      if (verbose & ncv > 1) cat(i, "/", ncv, "\n")
    }
    lambda <- mean(lam)
  }
  f <- glmnet(x = x, y = y, family = family, offset = offset, alpha = alpha, lambda = lambda,
              penalty.factor = penalty.factor, standardize = FALSE)
    
  f$coefficients <- as.numeric(coef(f))
  names(f$coefficients) <- rownames(coef(f))
  f$linear.predictors <- predict(f, newx = x, type ="link", offset = offset)
  if (family == "gaussian")
   f$dispersion <- bglm(y ~ f$linear.predictors - 1, start = 1, prior = "de", prior.mean = 1, prior.scale = 0, verbose = FALSE)$dispersion
  
  f$x <- x
  f$y <- y
  f$family <- family 
  f$alpha <- alpha
  f$penalty.factor <- penalty.factor 
  if (alpha == 1) f$prior.scale <- mean(1/(f$lambda * penalty.factor * nrow(x)))
  if (alpha == 0) f$prior.sd <- mean(sqrt(1/(f$lambda * penalty.factor * nrow(x))))
  f$offset <- offset
  f$aic <- deviance(f) + 2 * f$df
  
  f$call <- call
  
  if (family == "cox") class(f) <- c(class(f), "glmNet", "COXPH") 
  else class(f) <- c(class(f), "glmNet", "GLM")
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)
  if (verbose)
    cat("Computational time:", minutes, "minutes \n")

  return(f)
}

#*******************************************************************************

