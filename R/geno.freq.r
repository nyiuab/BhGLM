
# geno is a matrix of all genotypes
geno.freq <- function(geno, freq = TRUE, digits = 4, verbose = FALSE)
{
  geno <- make.main(geno, model = "additive", fill.missing = FALSE, verbose = verbose)
  counts <- apply(geno, 2, table)
  if(class(counts) == "matrix") counts <- t(counts)
  if(class(counts) == "list"){
    counts <- lapply(counts, function(x) {length(x) = 3; x = x})
    counts <- do.call(rbind, counts)
    counts <- apply(counts, 2, function(x) {x[is.na(x)] = 0; x = x})
  }
  colnames(counts) <- c("common", "heterzygote", "rare")
  MAF <- apply(counts, 1, function(x) pa = 1 - (x[1] + x[2]/2)/sum(x))
  E <- cbind((1 - MAF)^2, 2 * (1 - MAF) * MAF, MAF^2)
  E <- E * rowSums(counts)
  chisq <- rowSums((counts - E)^2/E)
  p.HWE <- 1 - pchisq(chisq, df = 1)
  if(freq) counts <- t(apply(counts, 1, function(x) x/sum(x))) 
  missing <- apply(geno, 2, function(x) length(x[is.na(x)]))
  results <- round(cbind(missing, counts, MAF, p.HWE), digits = digits)
  cat("\n")
  results
}


#*******************************************************************************
geno.pheno <- function(pheno, geno, sd = TRUE, digits = 2)
{
  geno <- make.main(geno = geno, fill.missing = F, geno.order = TRUE, model = "additive") + 1
  geno.list <- as.list(geno)
  pheno.mean <- lapply(geno.list, function(w) tapply(pheno, w, mean, na.rm = T))
  pheno.mean <- lapply(pheno.mean, function(x) {length(x) = 3; x = x})
  pheno.mean <- do.call(rbind, pheno.mean)
  colnames(pheno.mean) = c("common", "heterzygote", "rare")
  pheno.sd = NULL
  if(sd){
    pheno.sd = lapply(geno.list, function(w){
                         z = table(w)
                         w[w==names(z[z==1])] = NA
                         tapply(pheno, w, sd, na.rm=T) }
                      )
    pheno.sd = lapply(pheno.sd, function(x) {length(x)=3; x=x})
    pheno.sd = do.call(rbind, pheno.sd)
    colnames(pheno.sd) = c("common", "heterzygote", "rare")
  }
  results = round(cbind(pheno.mean, pheno.sd), digits=digits)
  results
}