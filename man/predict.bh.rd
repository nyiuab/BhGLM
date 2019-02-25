
\name{predict.bh}
\Rdversion{1.1}
\alias{predict.bh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Predictions for a Bayesian hierarchical model
}

\description{
This function obtains various predictive measures for new (or existing) data from a fitted model. 
}

\usage{
predict.bh(object, new.x, new.y, new.offset) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  a fitted object.   
  }
  \item{new.x}{
  data frame or matrix of new values for variables used in \code{object}. for object from \code{\link{bglm}} or \code{\link{bpolr}} or \code{\link{bcoxph}}, it is data frame; for object from \code{\link{glmNet}} or \code{\link{bmlasso}}, it is matrix.
}
  \item{new.y}{  
  vector of new response values corresponding to \code{new.x}. 
  If \code{new.x} or \code{new.y} are omitted, the fitted linear predictors are used for prediction. 
}
  \item{new.offset}{ 
  data frame or vector of offset values for new data points.  
  }
}

\details{
  This function calculates various predictive values, including linear predictor, predictive responses or probabilities,
predictive likelihood, and more. It can be used to evaluate the predictive performance of the model.    
}

\value{

The returned values include:
\item{lp}{linear predictors of all observations.}
\item{measures}{various predictive values.} 

For GLMs and polr, also include:
\item{y.fitted}{the fitted mean values, obtained by transforming the linear predictors by the inverse of the link function.}

For all GLMs and polr, \code{measures} includes:
\item{deviance}{estimate of deviance.}
\item{mse}{estimate of mean squared error.}

For binomial and polr models, \code{measures} also includes:
\item{auc}{area under ROC curve.}
\item{misclassification}{estimate of misclassification.}

For Cox models, \code{measures} includes:
\item{loglik}{partial log-likelihood.}
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
  \code{\link{cv.bh}}, \code{\link{predict.glm}}, \code{\link{predict.coxph}}
}

\examples{
 
library(BhGLM)

N = 1000
K = 50
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5, sigma=1.6) # simulate responses
yy$coefs

y = yy$y.ordinal; fam = binomial

# partition the data into two parts; 
# fit a model using the first part and evaluate the prediction using the second part 
x1 = x[1:(N/2),]; y1 = y[1:(N/2)]
x2 = x[(N/2):N,]; y2 = y[(N/2):N]

f1 = bglm(y1 ~ ., data=x1, family=fam, prior="de", prior.scale=0.5)   
plot.bh(f1, vars.rm = 1, threshold = 0.01, gap = 10)  

e1 = predict.bh(f1, x2, y2)
e1$measures


y = yy$y.surv
y1 = y[1:(N/2)]
y2 = y[(N/2):N]
f1 = bcoxph(y1 ~ ., data=x1, prior="de", prior.scale=0.5)   
plot.bh(f1, threshold = 0.01, gap = 10)  

e1 = predict.bh(f1, x2, y2)
e1$measures
 
}

