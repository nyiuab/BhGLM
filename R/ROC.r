
ROC <- function(y, y.fitted, main = " ", xlab="False positive rate", ylab="True positive rate",
                lwd=2, lty=1, col="black", type="s", pch=20, m=100,
                plot=TRUE, add=FALSE)
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
    if(!add) plot(c(0,1), c(0,1), xlim=c(0, 1), ylim=c(0, 1), main=main, xlab=xlab, ylab=ylab,
                  type="l", lty=3, col="gray")    # draw a diagonal line
    lines(fpf, tpf, type=type, lwd=lwd, lty=lty, col=col, pch=pch)
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
  w = list(AUC=auc, se=se, ci=c(lci, rci), pvalue=p)
  w
}
