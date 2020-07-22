
\name{make.inter}
\Rdversion{1.1}
\alias{make.inter}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Making Design Matrix of Interactions (GxG and GxE)
}

\description{
    This function is to construct a design matrix of interactions from two main-effect design matrices.
}

\usage{
make.inter(x1, x2, back = 0)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{ 
  design matrices of main-effect predictors. 
}
  \item{back}{
  For GxE interactions, \code{back = 0}. For the Cockerham and codominant GxG models, \code{back = 1}, and for other GxG models, \code{back = 0}.     
  if \code{back = 1}, interaction between rs123a and rs123d is not allowed. This is only need for the Cockerham and codominant models, because they defines two main effects for each marker. 
}
}

\value{
  This function returns a data.frame consisting of values of all interaction predictors.         
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link[BhGLM]{make.main}}
}

\examples{
library(BhGLM)

x = sim.x(n=100, m=10, genotype=6:10)
geno = x[, 6:10] #get genotype data
cov1 = x[, 1, drop=F] #get the first covariate

x.main = make.main(geno=geno, model="Cockerham", fill.missing=T)
x.GE = make.inter(x.main, cov1, back=0) #interactions between cov1 and all main effects
x.GG1 = make.inter(x.main[, 1:4], x.main[, 5:10], back=1) # some GxG interactions 
x.GG2 = make.inter(x.main, x.main, back=1) # all GxG interactions 

}
