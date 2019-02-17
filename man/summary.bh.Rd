
\name{summary.bh}
\Rdversion{1.1}
\alias{summary.bh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Numerically Summarizing Bayesian Hierarchical Model Fits
}

\description{
This function is to numerically summarize the results of the function \code{\link{bglm}} or \code{\link{bcoxph}} or \code{\link{bpolr}}.
}

\usage{
summary.bh(object)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  an object from \code{\link{bglm}}, \code{\link{bpolr}}, or \code{\link{bcoxph}}. 
}
  
}

\details{

}

\value{
This function returns: 

\item{coef}{coefficient estimates;} 

\item{se(coef)}{standard errors of coefficient estimates;} 

\item{pvalue}{p values for testing coefficients;} 

For logistic regression and Cox model, this function also returns:

\item{exp(coef)}{exponential coefficient estimates, i.e., odds ratio or hazards ratio;} 

\item{lower.95, upper.95}{exp(coef - 2 * se(coef)), exp(coef + 2 * se(coef));} 
       
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{

}

\examples{

}

