
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
make.inter(x1, x2, back = 1, na.action = na.pass)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x1, x2}{ 
  design matrices of main-effect predictors. Can be a factor. For a factor, its design matrix is first constructed internally.     
}
  \item{back}{
  For GxE interactions, \code{back = 0}. For the Cockerham and codominant GxG models, \code{back = 1}, and for other GxG models, \code{back = 0}.     
  if \code{back = 1}, interaction between rs123a and rs123d is not allowed. This is only need for the Cockerham and codominant models, 
  because they defines two main effects for each marker. 
}
  \item{na.action}{
  how NAs are treated. The default is \code{na.pass}, i.e., including all NAs. This is a general augument. See \code{\link{glm}}.
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

data(fake.cv) #load data from BhGLM
cov1 = fake.cv[, 3, drop = F] #get the first covariate
geno = fake.cv[, -c(1:4)] #get genotype data
x.main = make.main(geno = geno, model = "Cockerham", fill.missing = T)
x.GxE = make.inter(cov1, x.main, back = 0) #interactions between cov1 and all main effects
x.GxG1 = make.inter(x.main[, 1:10], x.main[, 11:20], back = 1) # some GxG interactions 
x.GxG2 = make.inter(x.main, x.main, back = 1) # all GxG interactions 

}
