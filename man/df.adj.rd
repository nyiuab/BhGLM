
\name{df.adj}
\Rdversion{1.1}
\alias{df.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating Adjusted Degree of Freedom
}

\description{
This function calculates adjusted degrees of freedom (i.e., effective number of parameters) 
contributed by specified variables based on an object from \code{\link{bglm}}. 
}

\usage{
df.adj(object, vars = 1:length(object$coefficients))
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  a fitted object from \code{\link{bglm}}. 
}
  \item{vars}{
  a vector of variable index or names; default is 1:length(object$coefficients), i.e., all coefficients.
  }
}

\details{
   In classical models, the degree of freedom equals the actual number of parameters. However, the degree of freedom (the effective number of parameters) in a hierarchical model can be much smaller than the actual number of parameters. In a hierarchical model, the effective number of parameters is defined as the difference between  the posterior mean of the deviance and the deviance at the posterior means of parameters of interest. This function calculates the effective number of parameters based on this definition. 
}

\value{
This function returns the adjusted degree of freedom (effective number of parameters) contributed by specified variables.

}
\references{
Yi, N., Xu, S., Lou, X.Y., and Mallick, H. (2013) Multiple Comparisons in Genetic Association Studies: A Hierarchical Modeling Approach. Statistical Applications in Genetics and Molecular Biology.

Spiegelhalter, D.J., Best, N.G., Carlin, B.P., and Linde, A.v.d. (2002) Bayesian measures of model complexity and fit (with discussion). Journal of the Royal Statistical Society Series B 64: 583-639.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  
}

\examples{

}

