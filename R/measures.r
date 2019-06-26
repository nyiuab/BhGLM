
measure.bh <- function(object, new.x, new.y, new.offset) 
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
    mu <- predict(object, newdata=new.x, type="response") 
    if (any(class(object) %in% "negbin")) object$dispersion <- object$theta
    measures <- measure.glm(y, mu, family=object$family$family, dispersion=object$dispersion) 
  }
  
  if (any(class(object) %in% "glmNet") | any(class(object) %in% "bmlasso"))
  {
    family <- object$family
    if (family=="cox")
      if (!is.Surv(y)) stop("'new.y' must be a Surv object")
    if (is.null(object$offset)) object$offset <- FALSE
    else object$offset <- TRUE
    if (missing(new.x)) new.x <- object$x
    newx <- as.matrix(new.x)
    if (family=="cox"){
      lp <- predict(object, newx=newx, newoffset=new.offset)
      lp <- as.vector(lp)
      measures <- measure.cox(y, lp)
    }
    else{ 
      mu <- predict(object, newx=newx, newoffset=new.offset, type="response")
      mu <- as.vector(mu)
      measures <- measure.glm(y, mu, family=family, dispersion=object$dispersion)
    }
  }
  
  if (any(class(object) %in% "coxph")) 
  {
    if (!is.Surv(y)) stop("'new.y' must be a Surv object")
    lp <- predict(object, newdata=new.x)  
    measures <- measure.cox(y, lp)
  }
  
  if (any(class(object) %in% "polr")) 
  {
    if (!is.factor(y)) stop("'new.y' must be a factor")
    probs <- predict(object, newdata=new.x, type="probs")
    measures <- measure.polr(y, probs)
  }
  
  measures    
}

#********************************************************************

measure.glm <- function(y, y.fitted, family, dispersion = 1) 
{
  if (NROW(y)!=NROW(y.fitted))
    stop("y and y.fitted should be of the same length", call. = FALSE)
  if (is.null(dispersion)) dispersion <- 1
  
  mu <- y.fitted
  if (substr(family, 1, 6)=="NegBin" | substr(family, 1, 17)=="Negative Binomial"
      | family=="nb") 
    family <- "NegBin"
  fam <- c("gaussian", "binomial", "poisson", "quasibinomial", "quasipoisson", "NegBin")
  if (! family %in% fam)
    stop("Measures for this family have not been implemented yet")
  
  if (family=="gaussian") 
    logL <- dnorm(y, mean=mu, sd=sqrt(dispersion), log=TRUE)
  if (family=="binomial" | family=="quasibinomial") 
    logL <- dbinom(y, size=1, prob=mu, log=TRUE)
  if (family=="poisson" | family=="quasipoisson") 
    logL <- dpois(y, lambda=mu, log=TRUE)
  if (family == "NegBin")
    logL <- dnbinom(y, size=dispersion, mu=mu, log=TRUE) 
   
  logL <- sum(logL, na.rm=TRUE)
  deviance <- -2 * logL
    
  mse <- mean((y - mu)^2, na.rm = TRUE)
  mae <- mean(abs(y - mu), na.rm = TRUE)
  measures <- list(deviance=deviance, mse=mse, mae=mae)
  if (family=="gaussian") {
      R2 <- (var(y, na.rm = TRUE) - mse)/var(y, na.rm = TRUE)
      measures <- list(deviance=deviance, R2=R2, mse=mse, mae=mae)
  }
  if (family=="binomial" | family=="quasibinomial") {
    if (!requireNamespace("pROC")) install.packages("pROC")
    require(pROC)
    AUC <- suppressMessages( pROC::auc(y, mu) )
    AUC <- as.numeric(AUC)
    misclassification <- mean(abs(y - mu) >= 0.5, na.rm = TRUE)
    measures <- list(deviance=deviance, auc=AUC, mse=mse, mae=mae, 
                     misclassification = misclassification)
  }
  
  round(unlist(measures), digits=3)
}

measure.polr <- function(y, y.fitted) 
{
  if (!requireNamespace("pROC")) install.packages("pROC")
  require(pROC)
  pred <- y.fitted
  if (is.vector(pred)) pred <- t(as.matrix(pred))
  if (NROW(pred)!=NROW(y))
    stop("y and y.fitted should be the same length.", call. = FALSE)
  AUC <- mse <- misclassification <- 0
  y.level <- levels(y)
  for (k in 1:NCOL(pred)) {
    y1 <- ifelse(y == y.level[k], 1, 0)
    AUC <- AUC + as.numeric( suppressMessages(pROC::auc(y1, pred[, k])) )
    misclassification <- misclassification + mean(abs(y1 - pred[, k]) > 0.5, na.rm = TRUE)
    mse <- mse + mean((y1 - pred[, k])^2, na.rm = TRUE)
  }
  AUC <- AUC/NCOL(pred)
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
  
  measures <- list(deviance = deviance, auc = AUC, mse = mse, misclassification = misclassification)
  round(unlist(measures), digits=3)
}

measure.cox <- function(y, lp) 
{
  if (NROW(y)!=NROW(lp))
    stop("y and lp should be of the same length", call. = FALSE)
  ff <- bcoxph(y ~ lp, prior.scale=0, prior.mean=1, verbose=FALSE)
  deviance <- -2 * ff$loglik[2]
  cindex <- summary(ff)$concordance[[1]]
  measures <- list(deviance = deviance, Cindex = cindex)
  round(unlist(measures), digits=3)
}

#***************************************************************************

surv.curves <- function (y, lp, probs=0.50, mark.time=TRUE, main=" ", 
                         lwd=2, lty=1, col=c(3, 4), add=FALSE) 
{
  library(survival)
  group = as.numeric(cut(lp, c(-Inf, quantile(lp, probs = probs), Inf)))   
  sf = survfit(y ~ group)
  if (!add)
    plot(sf, conf.int=FALSE, mark.time=mark.time, main=main, lwd=lwd, lty=lty, col=col, 
         xlab="Time", ylab="Survival probability")
  else
    lines(sf, conf.int=FALSE, mark.time=mark.time, lwd=lwd, lty=lty, col=col)
  logrank = survdiff(y ~ group)
  p = signif(pchisq(logrank$chisq, df=length(unique(group))-1, lower.tail=F), digits=3)
  legend("topright", paste("p =", p), cex = 0.8, bty = "n")
  return(logrank)
}


aucCox <- function (y, lp, main=" ", lwd=2, lty=1, col="black", add=FALSE) 
{
  time <- y[, 1]
  status <- y[, 2]
  tt <- sort(unique(time[status == 1]))
  nt <- length(tt)
  x <- lp
  AUCt <- rep(NA, nt)
  numsum <- denomsum <- 0
  for (i in 1:nt) {
    ti <- tt[i]
    Y <- sum(time >= ti)
    R <- which(time > ti)
    xi <- x[time == ti]
    num <- sum(x[R] < xi) + 0.5 * sum(x[R] == xi)
    AUCt[i] <- num/(Y - 1)
    numsum <- numsum + num
    denomsum <- denomsum + Y - 1
  }
  AUC <- numsum/denomsum
  if (!add) 
    plot(tt, AUCt, ylim = c(min(AUCt), max(AUCt)), main = main, xlab = "Time", ylab = "AUC(t)", type = "n")
  lines(lowess(data.frame(tt, AUCt)), lwd = lwd, lty = lty, col = col)
  abline(h = 0.5, lty = 3)
  
  return(list(AUCt = data.frame(time = tt, AUC = AUCt), AUC = AUC))
}


# the following functions need the package dynpred
peCox <- function (y, lp, FUN=c("Brier", "KL"), main="", lwd=2, lty=1, col="black", add=FALSE) 
{
  if (!requireNamespace("dynpred")) install.packages("dynpred")
  library(dynpred)
  FUN <- FUN[1]
  time <- y[, 1]
  status <- y[, 2]
  n <- length(time)
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  tt <- sort(unique(time[status == 1]))
  nt <- length(tt)
  
  x <- lp[ord]
  cox1 <- bcoxph(Surv(time, status) ~ x, init = 1, prior = "de",
                 prior.mean = 1, prior.scale = 0, verbose = FALSE)
  if (sum(x^2) == 0) 
    sf <- survfit(cox1, newdata = data.frame(x = x), type = "kalbfl")
  else sf <- survfit(cox1, newdata = data.frame(x = x))
  tt <- sf$time
  survmat <- sf$surv
  
  if (tt[1] > 0) {
    tsurv <- c(0, tt)
    survmat <- rbind(rep(1, n), survmat)
  }
  else tsurv <- tt
  nsurv <- length(tsurv)
  
  coxcens <- coxph(Surv(time, 1 - status) ~ 1)
  ycens <- coxcens[["y"]]
  p <- ncol(ycens)
  tcens <- ycens[, p - 1]
  dcens <- ycens[, p]
  xcens <- coxcens$linear.predictors
  coxcens <- coxph(Surv(tcens, dcens) ~ xcens)
  if (sum(xcens^2) == 0) 
    sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens), type = "kalbfl")
  else sfcens <- survfit(coxcens, newdata = data.frame(xcens = xcens))
  tcens <- sfcens$time
  censmat <- sfcens$surv
  if (tcens[1] > 0) {
    tcens <- c(0, tcens)
    censmat <- rbind(rep(1, n), censmat)
  }
  ncens <- length(tcens)
  
  tout <- unique(c(tsurv, tcens))
  tout <- sort(tout)
  nout <- length(tout)
  res <- pe(time, status, tsurv, survmat, tcens, censmat, FUN, tout)
  d <- which(res[, 2] == "NaN")
  if (length(d) > 0) res <- res[-d, ]
  
  if (!add)
    plot(res$time[res$Err > 0], res$Err[res$Err > 0], type = "s", ylim = c(0, max(res$Err)),
         lwd = lwd, lty = lty, col = col, main = main, xlab = "Time", ylab = "Prediction error")
  else 
    lines(res$time[res$Err > 0], res$Err[res$Err > 0], type = "s", ylim = c(0, max(res$Err)),
          lwd = lwd, lty = lty, col = col)
  
  return(res)
}

#*****************************************************

