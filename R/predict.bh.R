
predict.bh <- function (object, new.x, new.y, new.offset) 
{
  if (missing(new.y)){ 
    y <- object$y
    if (any(class(object) %in% "polr")) y <- object$model[, 1]
  }  
  else y <- new.y
  if (!missing(new.x)) 
    if (NROW(new.x)!=NROW(y))
      stop("'new.x' and 'new.y' should be the same length.", call. = FALSE)
  if (missing(new.offset)) new.offset <- NULL
  if (!is.null(new.offset)) 
    if(NROW(new.offset)!=NROW(y))
      stop("'new.offset' and 'new.y' should be the same length.", call. = FALSE)
  if (!is.null(new.offset) & !missing(new.x) & any(class(object) %in% c("glm", "coxph", "polr"))) 
    new.x <- data.frame(new.x, new.offset)
  
  if (any(class(object) %in% "glm")) 
  {
    if (!is.numeric(y)) stop("'new.y' must be numeric")
    lp <- predict(object, newdata=new.x) 
    mu <- object$family$linkinv(lp)
    if (any(class(object) %in% "negbin"))
      measures <- measure.nb(pred=mu, obs=y, theta=object$theta)
    else  
      measures <- measure.glm(pred=mu, obs=y, family=object$family, dispersion=object$dispersion) 
    res <- list(lp=lp, y.fitted=mu, measures=measures)
  }
  
  if (any(class(object) %in% "glmNet") | any(class(object) %in% "bmlasso"))
  {
    family <- object$family
    if (family=="cox")
      if (!is.Surv(y)) stop("'new.y' must be a Surv object")
    if (missing(new.x)) lp <- object$linear.predictors 
    else{
      if (is.vector(new.x)) new.x <- t(as.matrix(new.x))
      coefs <- object$coefficients
      if (names(coefs)[1]=="(Intercept)") new.x <- cbind(1, new.x)
      lp <- as.matrix(new.x) %*% coefs
      if (!is.null(new.offset)) lp <- lp + new.offset
    }
    lp <- as.numeric(lp)
    names(lp) <- 1:length(lp)
    mu <- lp
    if (family=="binomial") mu <- exp(lp)/(1 + exp(lp))
    if (family=="poisson") mu <- exp(lp)
    if (family=="cox")
      res <- list(lp=lp, measures=measure.cox(pred=lp, obs=y))
    else {
      if (family=="gaussian") fa <- gaussian()
      if (family=="binomial") fa <- binomial()
      if (family=="poisson") fa <- poisson()
      measures <- measure.glm(pred=mu, obs=y, family=fa, dispersion=object$dispersion)
      res <- list(lp=lp, y.fitted=mu, measures=measures)
    }
  }
  
  if (any(class(object) %in% "coxph")) 
  {
    if (!is.Surv(y)) stop("'new.y' must be a Surv object")
    lp <- predict(object, newdata=new.x)  
    res <- list(lp=lp, measures=measure.cox(pred=lp, obs=y))
  }
  
  if (any(class(object) %in% "polr")) 
  {
    if (!is.factor(y)) stop("'new.y' must be a factor")
    pred <- predict(object, newdata=new.x, type="probs")
    res <- list(y.fitted=pred, measures=measure.polr(pred=pred, obs=y))
  }
  
  res    
}

#********************************************************************

measure.glm <- function (pred, obs, family, dispersion = 1) 
{
  y <- obs
  mu <- pred
  family <- family$family
  if (family == "gaussian") logL <- dnorm(y, mu, sqrt(dispersion), log = TRUE)
  if (family == "binomial") logL <- dbinom(y, 1, mu, log = TRUE)
  if (family == "poisson") logL <- dpois(y, mu, log = TRUE)
    
  logL <- sum(logL, na.rm = TRUE)
  deviance <- -2 * logL
    
  mse <- mean((y - mu)^2, na.rm = TRUE)
  mae <- mean(abs(y - mu), na.rm = TRUE)
  measures <- list(deviance = deviance, mse = mse, mae = mae)
  if (family == "gaussian") {
      R2 <- (var(y, na.rm = TRUE) - mse)/var(y, na.rm = TRUE)
      measures <- list(deviance = deviance, mse = mse, R2 = R2)
  }
  if (family == "binomial") {
    nna <- !is.na(y)&!is.na(mu)
    auc <- roc.auc(y[nna], mu[nna], plot = FALSE)$AUC
    misclassification <- mean(abs(y - mu) >= 0.5, na.rm = TRUE)
    measures <- list(deviance = deviance, auc = auc, mse = mse, 
                     misclassification = misclassification)
  }
  
  round(unlist(measures), digits=3)
}

measure.nb <- function (pred, obs, theta = 1) 
{
  y <- obs
  mu <- pred
  logL <- dnbinom(y, size = theta, mu = mu, log = TRUE)
  logL <- sum(logL, na.rm = TRUE)
  deviance <- -2 * logL
  
  mse <- mean((y - mu)^2, na.rm = TRUE)
  mae <- mean(abs(y - mu), na.rm = TRUE)
  measures <- list(deviance = deviance, mse = mse, mae = mae)
  
  round(unlist(measures), digits=3)
}


measure.cox <- function (pred, obs) 
{
  y <- obs
  lp <- pred
  pl <- coxph(y ~ lp, init = 1, control = coxph.control(iter.max=1), method = "breslow")$loglik[1]
  nna <- !is.na(y)&!is.na(lp)
  cindex <- Cindex(y[nna], lp[nna])$cindex
  measures <- list(loglik = pl, Cindex = cindex)
  round(unlist(measures), digits=3)
}


measure.polr <- function (pred, obs) 
{
  y <- obs
  if (is.vector(pred)) pred <- t(as.matrix(pred))
  auc <- mse <- misclassification <- 0
  y.level <- levels(y)
  for (k in 1:NCOL(pred)) {
    y1 <- ifelse(y == y.level[k], 1, 0)
    auc <- auc + roc.auc(y1, pred[, k], plot = FALSE)$AUC
    misclassification <- misclassification + mean(abs(y1 - pred[, k]) > 0.5, na.rm = TRUE)
    mse <- mse + mean((y1 - pred[, k])^2, na.rm = TRUE)
  }
  auc <- auc/NCOL(pred)
  mse <- mse/NCOL(pred)
  misclassification <- misclassification/NCOL(pred)
  L <- rep(NA, NROW(pred))
  for (i in 1:NROW(pred)){
    y2 <- rep(0, NCOL(pred))
    for (k in 1:NCOL(pred)) y2[k] <- ifelse(y[i]==y.level[k], 1, 0)
    L[i] <- sum(y2*pred[i,])
  }
  L <- ifelse(L==0, 1e-04, L)
  deviance <- -2 * sum(log(L))
  
  measures <- list(deviance = deviance, auc = auc, mse = mse, misclassification = misclassification)
  round(unlist(measures), digits=3)
}

#***************************************************************************
