
# construct interaction matrix

make.inter <- function(x1, x2, back=0)   #back = 0, 1
{
  x1 <- as.matrix(x1)
  x2 <- as.matrix(x2)
  if( is.null(colnames(x1)) | is.null(colnames(x2)) ) 
    stop("no colnames for the inputs")
  f <- ~ x1:x2 - 1
  mf <- model.frame(f, na.action=na.pass)
  z <- model.matrix(f, mf)
  if(ncol(x1)==1 | ncol(x2)==1) 
    colnames(z) <- paste( paste("x1", colnames(x1), sep=""), ":",
                          paste("x2", colnames(x2), sep=""), sep="" )
  znames <- strsplit(colnames(z), split=":", fixed=TRUE)
  znames <- lapply(znames, function(x) substr(x, 3, nchar(x)-back))
  inc <- unlist(lapply(znames, function(x) (x[1]!=x[2])))
  z1 <- z[, inc, drop=FALSE]
  z1names <- strsplit(colnames(z1), split=":", fixed=TRUE)
  z1names <- lapply(z1names, function(x) substr(x, 3, nchar(x)))
  z1names1 <- lapply(z1names, function(x) sort(c(x[1], x[2])))
  z1names2 <- unique(z1names1)
  z2 <- z1[, match(z1names2, z1names1), drop=FALSE]
  z2names <- strsplit(colnames(z2), split=":", fixed=TRUE)
  z2names <- lapply(z2names, function(x) substr(x, 3, nchar(x)))
  colnames(z2) <- lapply(z2names, function(x) paste(x[2], ":" , x[1], sep=""))
  as.data.frame(z2)
}


