
# construct interaction matrix

make.inter <- function(x1, x2, back = 1, na.action = na.pass)
{
#back = 0, 1
  if(class(x1) == "factor"){
    f <- ~ x1
    mf <- model.frame(f, na.action = na.action)
    x1 <- model.matrix(f, mf)[, -1, drop = FALSE]
  }
  if(class(x2) == "factor"){
    f <- ~ x2
    mf <- model.frame(f, na.action = na.action)
    x2 <- model.matrix(f, mf)[, -1, drop = FALSE]
  }
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  if( is.null(colnames(x1)) | is.null(colnames(x2)) ) stop("no colnames for the inputs")
  f <- ~ x1:x2 - 1
  mf <- model.frame(f, na.action = na.action)
  z <- model.matrix(f, mf)
  if(ncol(x1) == 1 | ncol(x2) == 1) colnames(z) <- paste( paste("x1", colnames(x1), sep = ""),":",
                                                   paste("x2", colnames(x2), sep = ""), sep = "" )
  znames <- strsplit(colnames(z), split = ":", fixed = TRUE)
  znames0 <- lapply(znames, function(x) substr(x, 3, nchar(x) - back))
  index <- unlist(lapply(znames0, function(x) (x[1] != x[2])))
  z1 <- z[, index, drop = FALSE]
  z1names <- strsplit(colnames(z1), split = ":", fixed = TRUE)
  z1names0 <- lapply(z1names, function(x) substr(x, 3, nchar(x)))
  z1names1 <- lapply(z1names0, function(x) sort(c(x[1], x[2])))
  z1names2 <- unique(z1names1)
  z2 <- z1[, match(z1names2, z1names1), drop = FALSE]
  colnames(z2) <- lapply(z1names2, function(x) paste(x[1], ":" , x[2], sep = ""))
  data.frame(z2)
}