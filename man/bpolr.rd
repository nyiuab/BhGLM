
\name{bpolr}
\Rdversion{1.1}
\alias{bpolr}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayesian Ordinal Models
}

\description{
  This function is to set up and fit Bayesian hierarchical ordered logistic or probit regressions for ordinal response (e.g., disease severity)
  with Student-t prior on the coefficients. The default model is Bayesian proportional odds logistic regression, after which the function is named. 
  The Bayesian hierarchical ordered logistic or probit models include classical ordered logistic or probit regression as special case.
   
}

\usage{
bpolr(formula, data, weights, start, subset, na.action, 
      method = c("logistic", "probit", "loglog", "cloglog", "cauchit"), 
      contrasts = NULL, Hess = TRUE, 
      prior.mean = 0, prior.scale = 0.5, prior.df = 1, 
      verbose = FALSE, ...) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula, data, weights, subset, na.action, method, contrasts, Hess}{
  These arguments are the same as in \code{\link[MASS]{polr}}.
}
  \item{start}{
  Initial values for the coefficients (not for zeta).
}
  \item{prior.mean, prior.scale, prior.df, verbose }{
  These arguments are the same as in \code{\link{bglm}}.
}
\item{\dots}{
  further arguments for \code{\link[MASS]{polr}}.
}
  
}

\details{
  This function is an alteration of \code{\link[MASS]{polr}} for fitting classical ordered logistic or probit regressions. It uses the quasi-Newton algorithm (BFGS) as implemented in \code{\link{optim}} to fit the model by maximizing the posterior distribution. Bayesian ordinal models allow us to jointly analyze many correlated predictors. The function includes all the \code{\link[MASS]{polr}} arguments and also some new arguments for the hierarchical modeling. 
}
\value{
  This function returns an object of class "polr", including all outputs from the function \code{\link[MASS]{polr}}.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link[MASS]{polr}}, \code{\link{bglm}} 
}

\examples{

library(BhGLM)


N = 1000
K = 100
x = sim.x(n = N, m = K, corr = 0.6) # simulate correlated variables  
x = as.matrix(x)
h = rep(0.1, 4) # assign non-zero effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu = 10, herit=h, p.neg=0.5, sigma=1.6, quantiles = c(0.3, 0.6)) # simulate responses
yy$coefs
y = as.factor(yy$y.ordinal)
table(y)


par(mfrow = c(1, 3), cex.axis = 1, mar = c(3, 4, 4, 4))
# classical ordered logistic regression
library(MASS)
f1 = polr(y ~ ., data = x, Hess = T)
plot.bh(f1, gap = 5)
summary.bh(f1)

# equivalent to classical ordered logistic regression
f2 = bpolr(y ~ ., data = x, prior.scale = Inf)
plot.bh(f2, gap = 5)
summary.bh(f2)

# hierarchical ordered logistic regression
f3 = bpolr(y ~ ., data = x, prior.scale = 0.05)
plot.bh(f3, gap = 5)
summary.bh(f3)

}

