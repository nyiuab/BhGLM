
\name{geno.freq}
\Rdversion{1.1}
\alias{geno.freq}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Summarizing Genotypic Data
}

\description{
This function is to calculate frequecies of missing data, three genotypes, and MAF (minor allelic frequecies),
and to test Hardy-Weinberg Equilibrium (HWE) for each marker.
}

\usage{
geno.freq(geno, freq = TRUE, digits = 4, verbose = FALSE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geno}{ 
   a matrix of genotypic data. 
  }
  \item{freq}{
  logical. Shows frequencies or counts.
  }
  \item{digits}{
  digits to round with print.
  }
  \item{verbose}{
  logical. If \code{TRUE}, print out markers removed or with only two genotypes.
  }
}

\value{
    This function returns a matrix consisting of frequecies of missing data, three genotypes, MAF (minor allelic frequecies), and p-value for testing HWE.         
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\examples{
library(BhGLM)

x = sim.x(n=100, m=10, genotype=6:10)
geno = x[, 6:10] #get genotype data
geno.freq(geno, freq=TRUE)
geno.freq(geno, freq=FALSE)

}

