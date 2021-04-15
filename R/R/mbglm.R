

mbglm <- function(y, formula, data, family=NegBin(), prior=Student(0, 1),
                  min.p=0, verbose=TRUE)
{
  start.time <- Sys.time()
  
  call <- match.call()
  y <- nonzero(y=y, min.p=min.p, sort=FALSE)$y.filter
  if (missing(data)) data <- NULL
  
  if (verbose) cat("Analyzing", ncol(y), "responses: \n")
  fm <- y.one ~ .
  fm[[3]] <- formula[[2]]
  fit <- vector(mode="list", length=NCOL(y))
  names(fit) <- colnames(y)
  for (j in 1:ncol(y)){
    y.one <- y[, j]
    data1 <- data.frame(cbind(y.one, data))
    tryCatch( {
      fit[[j]] <- bglm(formula=fm, data=data1, family=family, prior=prior, 
                       verbose=FALSE ) 
      if (verbose) cat(j, "")
    }, error = function(e) {message("\n", "y", j, " error: ", conditionMessage(e), sep="")} )
  } 
  fit <- fit[!sapply(fit, is.null)]
  responses <- names(fit)
  variables <- names(fit[[1]]$coefficients)

  res <- list(fit=fit, responses=responses, variables=variables, call=call)
  
  class(res) <- "mbglm"
  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units="min"), 3)
  if (verbose) 
    cat("\n Computational time:", minutes, "minutes \n")
  
  res
}

summary.mbglm <- function(object) 
{
  obj.fit <- object$fit
  obj.resp <- object$responses
  obj.vars <- object$variables
  out <- lapply(obj.fit, summary.bh)
  
  res <- responses <- NULL
  for (j in 1:length(out)){
    res <- rbind(res, out[[j]])
    responses <- c(responses, rep(obj.resp[j], length(obj.vars)))
  }
  variables <- rownames(res)
  res <- data.frame(responses, variables, res)
  rownames(res) <- paste(res[,1], "--", res[,2], sep="")
  res$padj <- res$pvalue
  for(j in 1:length(obj.vars))
  {
    p <- res[res[,"variables"]==obj.vars[j], "pvalue"]
    nam <- rownames(res[res[,2]==obj.vars[j], ])
    res[nam, "padj"] <- signif(p.adjust(p, method="fdr"), 2)
  }
  colnames(res)[3:4] <- c("Estimate", "Std.Error")
  out <- res[res[,2]!="(Intercept)", ] 
  
  return(out)
}

nonzero <- function(y, total, min.p=0, sort=TRUE, plot=FALSE)
{
  y <- as(y, "matrix")
  if (is.null(colnames(y))) colnames(y) <- paste("y", 1:ncol(y), sep = "")
  if(missing(total)) total <- rowSums(y)
  total <- unlist(total)
  nonzero.p <- apply(y, 2, function(z) {length(z[z != 0])/length(z)} )
  if (sort) nonzero.p <- sort(nonzero.p, decreasing=T)
  sub <- names(nonzero.p[nonzero.p > min.p])
  if (length(sub) == 0) stop("min.p is too large") 
  else y0 <- y[, sub, drop=FALSE]
  total.mean <- round(mean(total), 2)
  total.sd <- round(sd(total), 2)
  if (plot){
    par(mfrow = c(1, 2), mar = c(5, 4, 4, 4))
    hist(total, nclass=1000, xlab="Total reads", main="")
    zero <- sort(1 - nonzero.p, decreasing=F)
    plot(x=1:length(zero), y=zero, xlab="Taxa", ylab="Zero Proportion")
  }
  list(nonzero.p=nonzero.p, total=total, total_mean_sd=c(total.mean, total.sd), y.filter=y0)
}


