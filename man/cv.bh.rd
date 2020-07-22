
\name{cv.bh}
\Rdversion{1.1}
\alias{cv.bh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-Validation for Bayesian Models or Elastic Net
}

\description{
The function \code{\link{cv.bh}} performs K-fold cross-validation and calculates cross-validated predictive measures
for Bayesian hierarchical GLMs and Cox survival model, or for elastic net from the package \bold{glmnet}.  
}

\usage{
cv.bh(object, nfolds = 10, foldid = NULL, ncv = 1, verbose = TRUE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  a fitted object. 
  }
  \item{nfolds}{
  number of folds(groups) into which the data should be split to estimate the cross-validation prediction error. 
  default is 10. Although nfolds can be as large as the sample size (leave-one-out CV), it is not recommended for large datasets. 
  }
  \item{foldid}{
  an optional vector (if ncv = 1) or matrix (if ncv > 1) of values between 1 and nfolds identifying what fold each observation is in. 
  If supplied, nfolds can be missing.If \code{foldid = NULL}, \code{nfolds} subsets will be generated randomly.  
  }
  \item{ncv}{
  repeated number of cross-validation.  
  }
  \item{verbose}{
  logical. If \code{TRUE}, print out computational time and progress.
  }
}

\details{
The data is divided randomly into \code{nfolds} subsets with equal (or approximately equal) numbers of indivudals. 
For each subset, the model is fit to data omitting that subset, 
and then predict the omitted responses and calculate various prediction errors in that subset. 
Since the folds are selected at random, the cross-validation results are random. Users can reduce this randomness by running cross-validation several times, and averaging the predictive values. 
Cross-validation is repeated \code{ncv} times. 
}

\value{

The returned values include:
\item{y.obs}{The observed responses.}
\item{lp}{linear predictors of all observations.}
\item{foldid}{a vector (if ncv = 1) or matrix (if ncv > 1) indicating folds.} 
\item{measures}{various predictive values.} 

For GLMs and polr, also include:
\item{y.fitted}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}

For all GLMs and polr, \code{measures} includes:
\item{deviance}{estimate of deviance.}
\item{mse}{estimate of mean squared error.}

For binomial and polr models, \code{measures} also includes:
\item{auc}{area under ROC curve.}
\item{misclassification}{estimate of misclassification error.}

For Cox models, \code{measures} includes:
\item{deviance}{deviance using cross-validated prognostic index.}
\item{Cindex}{concordance index.}

}

\references{
  Steyerberg, E. W., 2009 Clinical Prediction Models: A Practical Approch to Development, Validation, and Updates. Springer, New York.
  
  van Houwelinggen, H.G. & Putter, H. Dynamic Prediction in Clinical Survival Analysis, (CRC Press, 2012).
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{bglm}}, \code{\link{bcoxph}}, \code{\link{glmNet}}, \code{\link{bmlasso}}
}

\examples{
library(BhGLM)
library(survival)
library(glmnet)

N = 1000
K = 30
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu = 0, herit=h, p.neg=0.5, sigma=1.6) # simulate responses
yy$coefs


#y = yy$y.normal; fam = "gaussian"
y = yy$y.ordinal; fam = "binomial"

f = glmNet(x, y, family = fam, alpha = 1, ncv = 2)
cv = cv.bh(f, ncv = 2)
cv$measures 

f1 = bglm(y ~ ., data = x, family = fam, prior = De(scale=f$prior.scale))
cv1 = cv.bh(f1, foldid = cv$foldid)
cv1$measures

par(mfrow = c(1, 2), cex.axis = 1, mar = c(3, 4, 4, 4))
plot.bh(coefs = f$coef[-1], threshold = 10, gap = 10) 
plot.bh(f1, vars.rm = 1, gap = 10, col.pts = c("red", "black"), threshold = 0.01) 


# censored survival data
y = yy$y.surv

f = glmNet(x, y, family = "cox", alpha = 1, ncv = 2)
cv = cv.bh(f, ncv = 2)
cv$measures 

f1 = bcoxph(y ~ ., data = x, prior = De(scale=f$prior.scale))
cv1 = cv.bh(f1, foldid = cv$foldid)
cv1$measures

par(mfrow = c(1, 2), cex.axis = 1, mar = c(3, 4, 4, 4))
plot.bh(coefs = f$coef, threshold = 10, gap = 10) 
plot.bh(f1, gap = 10, col.pts = c("red", "black"), threshold = 0.01) 

}

