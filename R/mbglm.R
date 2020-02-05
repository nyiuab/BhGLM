

mbglm <- function(y, formula, data, family=NegBin(), prior=Student(0, 1),
                  verbose=TRUE)
{
  start.time <- Sys.time()
  
  call <- match.call()
  y <- as.matrix(y)
  if (NCOL(y)==1) stop("'y' should include multiple columns (responses)")
  if (is.null(colnames(y))) colnames(y) <- paste("y", 1:ncol(y), sep="")
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

summary.mbglm <- function(object, adj.p=FALSE, vr.name=NULL) 
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
  rownames(res) <- NULL
  
  if (adj.p)
  {
    for(j in 1:length(obj.vars))
    {
      p <- res[res[,"variables"]==obj.vars[j], "pvalue"]
      nam <- rownames(res[res[,2]==obj.vars[j], ])
      res[nam, "pvalue"] <- signif(p.adjust(p, method="fdr"), 2)
    }
  }
  out <- res
  
  if (!is.null(vr.name))
  {
    if (!vr.name %in% c(obj.resp, obj.vars)) stop("wrong name given")
    if (vr.name %in% obj.resp) res0 <- res[res[, 1]==vr.name, ]
    if (vr.name %in% obj.vars) res0 <- res[res[, 2]==vr.name, ]
    out <- list(res=res, res0=res0)
  }
  
  return(out)
}

