
mc.adj <- function (p, df.adj = length(p), digits = 10)
{
  adjust <- c("none", "bonferroni", "holm", "hochberg", "hommel", "BH", "BY")
  adjust.p <- NULL
  for (j in 1:length(adjust)) adjust.p <- cbind(adjust.p, p.adjust(p = p, adjust[j]))
  colnames(adjust.p) <- adjust
  if (!is.null(df.adj)){
    if (df.adj > length(p)) df.adj <- length(p) 
    adjust.p <- cbind(adjust.p, pmin(1, p * df.adj))
    colnames(adjust.p) <- c(adjust, "bonferroni.adj") 
  }
  out <- adjust.p <- round(adjust.p, digits = digits)
  return(out)
}