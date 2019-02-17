
plot.bh <- function(object = NULL, coefs, sds = NULL, pvalues = NULL, vars.rm = NULL, 
                    threshold = 0.05, show.all.vars = FALSE, show.pvalues = TRUE, gap = 0, 
                    main = " ", cex.main = 0.9, xlim = NULL, cex.axis = 0.8, cex.var = 0.8, cex.pts = 1, 
                    pch.pts = 20, type = "p", lwd = 1, lty = 1, line = 0, col.pts = "black", 
                    OR = FALSE, add = FALSE)
{
  if (threshold <= 0) stop("'threshold' should be greater than 0")
  if (!is.null(object)) {
    out <- summary.bh(object) 
    if(!is.null(vars.rm)) out <- out[-vars.rm, , drop = FALSE]
    coefs <- out[, 1]
    sds <- out[, 2]
    pvalues <- out[, 3]
  }
  if (is.null(names(coefs))) {
    names(coefs) <- paste("v", 1:length(coefs), sep = "")
    warning("variables are not given!")
  }
  if (is.null(sds)) sds <- rep(0, length(coefs))
  if (threshold >= 1) {  # plot top coefs
    threshold <- as.integer(threshold)
    top <- sort(abs(coefs), decreasing = TRUE)[1:threshold]
    if (!is.null(pvalues)) threshold <- max(pvalues[names(top)])
    if (is.null(pvalues)) {
      pvalues <- rep(1, length(coefs)) 
      names(pvalues) <- names(coefs)
      pvalues[names(top)] <- 0
      show.pvalues <- FALSE
      threshold <- 0.05
    }
  }
  n <- length(coefs)
  if (is.numeric(sds)) {
    coef.l <- coefs - 2*sds
    coef.h <- coefs + 2*sds
  }
  if (is.matrix(sds)) { # can be posterior quantiles
    coef.l <- sds[, 1]
    coef.h <- sds[, 2]
  }
  coefs <- coefs[n:1]
  coef.l <- coef.l[n:1] 
  coef.h <- coef.h[n:1]
  if (OR) {
    coefs <- exp(coefs) 
    coef.l <- exp(coef.l) 
    coef.h <- exp(coef.h)
  }
  if (!is.null(pvalues)) {
    pvalues <- pvalues[n:1]
    pvalues <- ifelse(pvalues == "NaN", 1, pvalues)
  }
  p <- rep(-1, n)
  if (!is.null(pvalues)) p <- signif(pvalues, 2)
  varnames <- names(coefs)
  if (is.null(varnames)) varnames <- paste("v", n:1, sep = "")
  if (!show.all.vars) varnames <- ifelse(p <= threshold, varnames, "")
  if (length(col.pts) == 1) col.pts[2] <- "gray"

  z <- c(1:n)
  if (n > 1)
    for(i in 2:n) z[i] <- ifelse(p[i] <= threshold, z[i-1] + 1 + gap, z[i-1] + 1)
  if (is.null(xlim)) xlim = c(min(coef.l), max(coef.h))
  if (!add) plot(coefs, z, main = main, cex.main = cex.main, xlab = "", xlim = xlim,
                 ylab = "", yaxt = "n", type = "n", frame.plot = FALSE, cex.axis = cex.axis)
  axis(side = 3, cex.axis = cex.axis)
  if (!OR) lines(c(0, 0), c(-10, max(z) + 10), lty = 2, lwd = 1)
  else lines(c(1, 1), c(-10, max(z) + 10), lty = 2, lwd = 1)
  col <- rep(col.pts[1], n)
  for (i in 1:n) {
    col[i] <- ifelse(p[i] <= threshold, col.pts[1], col.pts[2])
    lines(c(coef.l[i], coef.h[i]), c(z[i], z[i]), col = col[i], lwd = lwd, lty = lty)
    cex <- ifelse(p[i] <= threshold, cex.var, cex.var - 0.1)
    if (!add) mtext(paste(varnames[i]), side = 2, at = z[i], cex = cex, las = 1, col = col[i], line = line)
    cex <- ifelse(p[i] <= threshold, cex.var - 0.1, cex.var - 0.2)
    if(show.pvalues & p[1] > 0 & varnames[i] != "")
      mtext(p[i], side = 4, at = z[i], cex = cex, las = 1, col = col[i], line = line/2)
  }
  points(coefs, z, pch = pch.pts, cex = cex.pts, col = col, type = type, lwd = lwd, lty = lty)

}



