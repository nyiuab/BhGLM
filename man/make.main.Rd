
\name{make.main}
\Rdversion{1.1}
\alias{make.main}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Making Design Matrix of Main Effects from Genotypic Data
}

\description{
    This function is to construct a design matrix of main effects from genotypic data of genetic markers. The genotypic data can include missing values.
  Different genetic models can be used, which transform the three-level (or two-level for a backcross) genotypic data
  to main-effect predictors.    
}

\usage{
make.main(geno, model = c("Cockerham", "codominant", "additive", "dominant", "recessive", "overdominant"),
          fill.missing = TRUE, ind.group = NULL, geno.order = TRUE, 
          loci.names = c("marker", "position"), imprint = TRUE, verbose = FALSE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{ 
    For human association data, it is a matrix or data frame of genotypes with dimension \code{n*J}, where \code{n} is the number of individuals and \code{J} is the number of markers. 
    The genotype data can be any three numbers or characters that indicate three genotypes (but the heterozygote should be the second by sort).
    For experimental crosses, it is an object of class \code{cross}. See \code{\link[qtl]{read.cross}} in the package \code{qtl} for details.
}
  \item{model}{ 
  a genetic model to construct main-effect predictors. 
}
  \item{fill.missing}{ 
   logical. If \code{TRUE}, fill in missing genotypic data. The default is \code{TRUE}. If \code{FALSE}, individuals with missing data will be removed from the analysis.   
}
  \item{ind.group}{ 
  a vector of length \code{n}, indicating groups of individuals (e.g., case-control status, race). 
  If provided, we fill in missing data based on observed genotype data for each group separately.
  The default is \code{NULL}; we don't use group information.  
}
  \item{geno.order}{ 
   logical. If \code{TRUE}, re-code genotypes as 1: common homozygote, 2: heterozygote, 3: rare homozygote.  
}
  \item{loci.names}{ 
  the way to name main-effect predictors; use marker names or chromosome positions.  
}
  \item{imprint}{
  logical. Indicates whether consider imprinting effects. 
}
\item{verbose}{
  logical. If \code{TRUE}, print out markers removed or with only two genotypes.
}
}

\details{
This function provides different genetic models to code the main-effect predictors. 
    
Denote common homozygote (i.e., the homozygote with higher frequency), heterozygote, and rare homozygote for each SNP by c, h, and r, respectively.

The \code{Cockerham} model defines two main effects for each SNP (with suffix 'a' and 'd'): an additive predictor as -1, 0, and 1 for c, h, and r, and a dominance predictor as -0.5 for c and r and 0.5 for h.

The \code{codominant} model also introduces two main effects for each SNP (with suffix 'r' and 'h'), with the two main-effect predictors being two indicator variables with the common homozygote c chosen as the reference group:
'r' and 'h' represent indicators for rare homozygote and heterozygote, respectively.

The \code{additive} model defines a main-effect predictor for each SNP, equal to 0, 1, 2 for c, h, r, respectively.

The \code{dominant} model defines a main-effect predictor for each SNP, equal to 1 for r and h, and 0 for c.

The \code{recessive} model defines a main-effect predictor for each SNP, equal to 1 for r, and 0 for h and c.

The \code{overdominant} model defines a main-effect predictor for each SNP, equal to 1 for h, and 0 for c and r.
  
For missing genotypes, we first calculate the genotypic probabilities of missing genotypes conditioning on the observed marker data, and then use these conditional probabilities to construct the main-effect predictors as above.
For QTL mapping in experimental crosses, we use the multipoint method as implemented in \code{R/qtl} (see \code{\link[qtl]{calc.genoprob}}) and \code{R/qtlbim} (see \code{\link[qtlbim]{qb.genoprob}}).    
For human association data, we simply replace missing genotypes by their expected values (i.e., dosages) based only on the observed genotypes for that marker.

This function removes markers with only one genotype or more than three genotypes, and for markers with only two genotypes, always uses genotype indicator variables.
}

\value{
This function returns a data frame consisting of values of all main-effect predictors.         
}

\author{
Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link[qtl]{read.cross}}, \code{\link[qtl]{calc.genoprob}}, \code{\link[qtlbim]{qb.genoprob}} 
}

\examples{
library(BhGLM)

x = sim.x(n=100, m=10, genotype=6:10)
geno = x[, 6:10] #get genotype data
x.g = make.main(geno=geno, model="additive", fill.missing=T)
x.g = make.main(geno=geno, model="Cockerham", fill.missing=T)

}

