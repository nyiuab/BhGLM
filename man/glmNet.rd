
\name{glmNet}
\Rdversion{1.1}
\alias{glmNet}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Fit GLMs or Cox Models with Lasso or Elastic Net
}

\description{
 Fit a generalized linear model or Cox model via the cyclic coordinate descent algorithm using the functions 
 \code{\link{glmnet}} and \code{\link{cv.glmnet}} in the package \bold{glmnet}. 
}

\usage{    
glmNet(x, y, family = c("gaussian", "binomial", "poisson", "cox"), weights = rep(1, nrow(x)), offset = NULL,
       alpha = c(1, 0.5, 0), lambda,
       penalty.factor = rep(1, ncol(x)), nfolds = 10, ncv = 10, 
       verbose = TRUE)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x, y, family, weights, offset, alpha, lambda, penalty.factor, nfolds}{ 
  These arguments are the same as in the functions \code{\link{glmnet}} and \code{\link{cv.glmnet}} in the package \bold{glmnet}.
  }
  \item{ncv}{
  repeated number of cross-validation.  
  }
  \item{verbose}{
  logical. If \code{TRUE}, print out the computational time and progress.
  }
}

\details{
  The function \code{\link{cv.glmnet}} performs cross-validation to determine an optimal penalty \code{lambda}.  
  Since the folds are selected at random, the estimate of the optimal penalty is not stable and depends on the folds. 
  This function does K-fold cross-validation \code{ncv} times and uses the mean of the \code{ncv} penalty values as the estimate of the optimal penalty \code{lambda}, 
  and then fits the elastic net model (including lasso and ridge) using the optimal penalty \code{lambda}.   
}

\value{
  This function returns all outputs from the function \code{\link{glmnet}}, and also \code{prior.scale}, which can be used in Bayesian hierarchical models.
}

\references{
 Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw 33, 1-22.
 
 Simon, N., Friedman, J., Hastie, T. & Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software 39, 1-13.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{glmnet}}, \code{\link{bmlasso}}, \code{\link{cv.glmnet}}
}

\examples{

 see examples in \code{\link{bmlasso}} and \code{\link{cv.bh}}. 

}