
ROC <- function(y, y.fitted, main = " ", xlab="False positive rate", ylab="True positive rate",
                lwd=2, lty=1, col="black", type="s", pch=20,
                m=100, plot=TRUE, add=FALSE)
{
  dn <- y
  pre <- y.fitted
  # draw roc curve
  n <- length(dn)                   # sample size
  d <- length(dn[dn != 0])          # number of disease (dn=1)
  sr <- (max(pre) - min(pre))/m     # size of reclassifying interval
  tpf <- rep(NA, m)                 # true positive fraction
  fpf <- rep(NA, m)                 # false positive fraction
  for (i in 1:(m + 1)){
    tpf[i] <- sum(pre >= max(pre) - sr*(i - 1) & dn == 1)/d
    fpf[i] <- sum(pre >= max(pre) - sr*(i - 1) & dn == 0)/(n-d)
  }
  alpha <- max(pre) - sr*c(0:m)
  
  if (add) plot <- TRUE
  if (plot){
    if(!add) plot(c(0,1), c(0,1), xlim = c(0, 1), ylim = c(0, 1), main = main, xlab = xlab, ylab = ylab,
                  type = "l", lty = 3, col = "gray")    # draw a diagonal line
    lines(fpf, tpf, type = type, lwd = lwd, lty = lty, col = col, pch = pch)
  }
  
  # calculate area under curve based on non-parametric statistic (Mann-Whitney U)
  tp2 <- rep(0, m)                  # number of true positive for each interval of cut-off values
  fp2 <- rep(0, m)                  # number of false positive for each interval of cut-off values
  for (i in 1:m){
    tp2[i] <- sum((pre >= min(pre) + sr*(i - 1) & pre < min(pre) + sr * (i)) & dn == 1)
    fp2[i] <- sum((pre >= min(pre) + sr*(i - 1) & pre < min(pre) + sr * (i)) & dn == 0)
  }
  
  dtp2 <- rep(0, m)                 # degression of tp2
  afp2 <- rep(0, m)                 # acumulation of fp2
  dtp2[1] <- (d - tp2[1])
  for (i in 1:(m - 1)){
    dtp2[i + 1] <- dtp2[i] - tp2[i + 1]
    afp2[i + 1] <- afp2[i] + fp2[i]
  }
  
  a <- rep(0, m)
  for (i in 1:m) {a[i] <- fp2[i] * dtp2[i] + (fp2[i] * tp2[i]) / 2}
  auc <- sum(a) / (d * (n - d) + 1e-06)     # area under curve
  if (auc <= 0) auc <- 0.5
  
  # test hypothesis of auc=0.5 (Mann-Whitney U test)
  q1 <- rep(0, m)                   # degression (leijian) of tp2
  q2 <- rep(0, m)                   # acumulation (leijia) of fp2
  for (i in 1:m){
    q1[i] <- fp2[i] * ((dtp2[i])^2 + tp2[i] * dtp2[i] + ((tp2[i])^2) / 3)
    q2[i] <- tp2[i] * ((afp2[i])^2 + fp2[i] * afp2[i] + ((fp2[i])^2) / 3)
  }
  sq1 <- sum(q1) / (d^2 * (n - d) + 1e-06)
  sq2 <- sum(q2) / (d * (n - d)^2 + 1e-06)
  se <- z <- p <- lci <- rci <- NULL
  se <- (auc * (1 - auc) + (d - 1) * (sq1 - auc^2) + ((n - d) - 1) * (sq2 - auc^2)) / (d * (n - d) + 1e-06)
  if(se > 0){
    se <- sqrt(se)
    z <- (auc - 0.5) / se             # z score
    p <- 2 * pnorm(-abs(z))           # p-value for two-tailed test
    lci <- auc - 1.96 * se            # left boundary of 95% confidence interval
    rci <- auc + 1.96 * se            # right boundary of 95% confidence interval
  }
  w = list(alpha = alpha, TPF = tpf, FPF = fpf, AUC = auc, SE = se, 
           Test.statistics = z, P.value = p, CI = c(lci, rci))
  w
}

#*******************************************************************************

surv.curves <- function (y, lp, probs = 0.50, main = " ", lwd = 2, lty = 1, 
                         col = c(3, 4), add = FALSE) 
{
  library(survival)
  group = as.numeric(cut(lp, c(-Inf, quantile(lp, probs = probs), Inf)))   
  sf = survfit(y ~ group)
  if (!add)
    plot(sf, conf.int = FALSE, mark.time = TRUE, main = main, lwd = lwd, lty = lty, col = col, 
         xlab = "Time", ylab = "Survival probability")
  else
    lines(sf, conf.int = FALSE, mark.time = TRUE, lwd = lwd, lty = lty, col = col)
  logrank = survdiff(y ~ group)
  p = signif(pchisq(logrank$chisq, df = length(unique(group))-1, lower.tail = F), digits = 3)
  legend("topright", paste("p =", p), bty = "n")
  return(logrank)
}

Cindex <- function (y, lp) 
{
  time <- y[, 1]
  status <- y[, 2]
  x <- lp
  n <- length(time)
  ord <- order(time, -status)
  time <- time[ord]
  status <- status[ord]
  x <- x[ord]
  wh <- which(status == 1)
  total <- concordant <- 0
  for (i in wh) {
    for (j in ((i + 1):n)) {
      tt <- (time[j] > time[i])
      if (is.na(tt)) tt <- FALSE
      if (tt) {
        total <- total + 1
        if (x[j] < x[i]) concordant <- concordant + 1
        if (x[j] == x[i]) concordant <- concordant + 0.5
      }
    }
  }
  return(list(concordant = concordant, total = total, cindex = concordant/total))
}


aucCox <- function (y, lp, main = " ", lwd = 2, lty = 1, col = "black", add = FALSE) 
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
#    num <- 0
#    for (k in 1:length(xi))
#      num <- num + sum(x[R] < xi[k]) + 0.5 * sum(x[R] == xi[k])
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
peCox <- function (y, lp, FUN = c("Brier", "KL"), main = "", lwd = 2, lty = 1, col = "black", add = FALSE) 
{
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


#*******************************************************************************

#*******************************************************************************

aucwCox <- function (y, lp, width = 2, plot = TRUE) 
{
  time <- y[, 1]
  status <- y[, 2]
  tt <- sort(unique(time[status == 1]))
  x <- lp
  ttw <- c(tt, tt - width)
  ttw <- ttw[ttw > 0]
  ttw <- sort(unique(ttw))
  ntw <- length(ttw)
  AUCw <- rep(NA, ntw)
  for (j in 1:ntw) {
    twj <- ttw[j]
    ttj <- tt[((tt >= twj) & (tt <= twj + width))]
    ntj <- length(ttj)
    AUCt <- rep(NA, ntj)
    numsum <- denomsum <- 0
    for (i in 1:ntj) {
      ti <- ttj[i]
      Y <- sum(time >= ti)
      R <- which(time > ti)
      xi <- x[time == ti]
      num <- sum(x[R] < xi) + 0.5 * sum(x[R] == xi)
      numsum <- numsum + num
      denomsum <- denomsum + Y - 1
    }
    AUCw[j] <- numsum/denomsum
  }
  if (plot){
    plot(ttw, AUCw, type = "s", lwd = 2, ylim = c(min(AUCw), max(AUCw)),
         xlab = "Time", ylab = "Dynamic C")
  }
  
  res <- data.frame(time = ttw, AUCw = AUCw)
  attr(res, "width") <- width
  return(res)
}

pewCox <- function (y, lp, width = 2, FUN = c("Brier", "KL"), plot = TRUE) 
{
  library(BhGLM)
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
  xcens <- coxcens$linear.predictors[ord]
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
  tout <- c(tout, tout - width)
  tout <- tout[tout >= 0]
  if (min(tout) > 0) tout <- c(0, tout)
  tout <- sort(unique(tout))
  nout <- length(tout)
  if (nout > 0) {
    tout <- tout[-nout]
    nout <- nout - 1
  }
  else {
    tout <- sort(tout)
    nout <- length(tout)
  }
  res <- pew(time, status, tsurv, survmat, tcens, censmat, width, FUN, tout)
  d <- which(res[, 2] == "NaN")
  if (length(d) > 0) res <- res[-d, ]
  
  #    interval.err <- evalstep(res$time, res$Err, 0:max(res$time))
  #    f0 <- survfit(Surv(time, status) ~ 1)
  #    KMstart <- evalstep(f0$time, f0$surv, 0:max(res$time), subst = 0)
  #    total.err <- sum(KMstart * interval.err)
  
  if (plot)
    plot(res$time[res$Err > 0], res$Err[res$Err > 0], type = "s", ylim = c(0, max(res$Err)),
         lwd = 2, xlab = "Time", ylab = "Prediction error")
  
  return(res)
}

#*******************************************************************************************
