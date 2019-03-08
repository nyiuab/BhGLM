
summary.bh <- function (object, digits = 3)
{
  if (any(class(object) == "glm")) res <- summary.bglm(object, digits=digits)
  if (any(class(object) == "coxph")) res <- summary.bcoxph(object, digits=digits)
  if (any(class(object) == "polr")) res <- summary.bpolr(object, digits=digits)
  res
}

summary.bglm <- function (object, digits)
{
  if (is.null(object$method.coef)) object$method.coef <- "joint"
  if (object$method.coef == "joint") 
    res <- summary(object, dispersion = object$dispersion)$coefficients[, c(1, 2, 4), drop = FALSE]
  else {
    coefs <- object$coefficients
    covmat <- vcov.bh(object)
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coefs/s.err
    pvalue <- 2 * pnorm(-abs(tvalue))
    res <- cbind(coefs, s.err, pvalue)
  }
  
  if (object$family[[1]] == "binomial" & object$family[[2]] == "logit") {
    OR <- exp(res[, 1])
    lower.95 <- exp(res[, 1] - 2 * res[, 2])
    upper.95 <- exp(res[, 1] + 2 * res[, 2])
    OR <- cbind(OR, lower.95, upper.95)
    res <- cbind(res[, 1:2], OR, res[,3])
  }
  
  if (ncol(res) == 3)
    colnames(res) <- c("coef", "se(coef)", "pvalue")
  if (ncol(res) == 6)
    colnames(res) <- c("coef", "se(coef)", "exp(coef)", "lower.95", "upper.95", "pvalue")
  res <- cbind(round(res[, 1:(ncol(res)-1), drop=FALSE], digits=digits),
               signif(res[, ncol(res), drop=FALSE], digits=digits))
  
  return(res)
}

summary.bcoxph <- function (object, digits)
{
  if (is.null(object$method.coef)) object$method.coef <- "joint"
  if (object$method.coef == "joint") {
    out <- summary(object)
    coefs <- out$coefficients
    if (ncol(coefs) == 5) {
      coefs <- coefs[, c(1, 3, ncol(coefs)), drop = FALSE]
    }
    if (ncol(coefs) == 6) {
      coefs[, ncol(coefs)] <- pchisq(coefs[, 4], df = 1, lower.tail = FALSE)
      coefs <- coefs[, c(1, 2, ncol(coefs)), drop = FALSE]
    }
    conf.int <- out$conf.int[, c(1, 3, 4), drop = FALSE]
    res <- cbind(coefs, conf.int) 
  }
  else {
    coefs <- object$coefficients
    covmat <- vcov.bh(object)
    var.cf <- diag(covmat)
    s.err <- sqrt(var.cf)
    tvalue <- coefs/s.err
    pvalue <- 2 * pnorm(-abs(tvalue))
    res <- cbind(coefs, s.err, pvalue)
    OR <- exp(res[, 1])
    lower.95 <- exp(res[, 1] - 2 * res[, 2])
    upper.95 <- exp(res[, 1] + 2 * res[, 2])
    OR <- cbind(OR, lower.95, upper.95)
    res <- cbind(res, OR)
  }
  
  colnames(res) <- c("coef", "se(coef)", "pvalue", "exp(coef)", "lower.95", "upper.95")
  res <- res[,c("coef", "se(coef)", "exp(coef)", "lower.95", "upper.95", "pvalue"), drop = FALSE]
  res <- cbind(round(res[, 1:(ncol(res)-1), drop=FALSE], digits=digits),
               signif(res[, ncol(res), drop=FALSE], digits=digits))
  return(res)
}

summary.bpolr <- function (object, digits)
{
  require(MASS)
  out <- summary(object)$coefficients[names(object$coefficients), , drop = FALSE]
  tvalue <- out[, 3]
  pvalue <- 2 * pnorm(-abs(tvalue))
  pvalue <- as.matrix(pvalue)
  out <- cbind(out[, -3, drop = FALSE], pvalue)
  OR <- exp(out[, 1])
  lower.95 <- exp(out[, 1] - 2 * out[, 2])
  upper.95 <- exp(out[, 1] + 2 * out[, 2])
  OR <- cbind(OR, lower.95, upper.95)
  res <- cbind(out, OR)

  colnames(res) <- c("coef", "se(coef)", "pvalue", "exp(coef)", "lower.95", "upper.95")
  res <- res[,c("coef", "se(coef)", "exp(coef)", "lower.95", "upper.95", "pvalue"), drop = FALSE]
  res <- cbind(round(res[, 1:(ncol(res)-1), drop=FALSE], digits=digits),
               signif(res[, ncol(res), drop=FALSE], digits=digits))
  return(res)
}

vcov.bh <- function (object)
{
  if (is.null(object$method.coef)) object$method.coef <- "joint"
  if (object$method.coef == "joint" | any(class(object) == "coxph")) 
  {
    vn <- names(object$coefficients)
    v <- vcov(object)[vn, vn]
  }
  else {
    X <- model.matrix(object)
    w <- sqrt(object$weights)
    object$prior.sd <- ifelse(object$prior.sd > 100, 100, object$prior.sd)
    d <- diag(object$prior.sd^2)
    x0 <- X * w
  
    if (ncol(X) <= nrow(X))
      v <- solve(t(x0) %*% x0 + diag(1/object$prior.sd^2))
  
    if (ncol(X) > nrow(X)) {
      I <- diag(NROW(X))
      xd <- x0 %*% d %*% t(x0)
      x1 <- solve(I + xd)
      x2 <- t(x0) %*% x1 %*% x0
      x3 <- (x2 * diag(d)) * diag(d)
      v <- d - x3
    }
    v <- v * object$dispersion
  }

  v
}


#****************************************************************************************



