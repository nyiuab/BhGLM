
\name{bmlasso}
\Rdversion{1.1}
\alias{bmlasso}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayesian Spike-and-Slab Lasso or Ridge GLMs and Cox Models
}

\description{
  This function is to set up Bayesian GLMs or Cox models with spike-and-slab mixture double-exponential (mde) or t prior (mt) (called the Bayesian spike-and-slab mixture lasso or ridge), 
  and to fit the model by incorporating EM steps into the fast coordinate descent algorithm. 
}

\usage{    
  bmlasso(x, y, family = c("gaussian", "binomial", "poisson", "cox"), weights = rep(1, nrow(x)), offset = NULL, epsilon = 1e-04, maxit = 50, 
          init = NULL, group = NULL, prior = c("mde", "mt"), prior.df = 1, ss = c(0.04, 0.5),
          Warning = FALSE, verbose = TRUE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{ 
  input matrix, of dimension nobs x nvars; each row is an observation vector.
  both x and y cannot include missing values.
  }
  \item{y}{
   response variable. Quantitative for family="gaussian", or family="poisson" (non-negative counts). 
   For family="gaussian", y is always been standardized. 
   For family="binomial", y should be either a factor with two levels, or a two-column matrix of counts or proportions (the second column is treated as the target class; for a factor, the last level in alphabetical order is the target class). 
   For family="cox", y should be a two-column matrix with columns named 'time' and 'status'. The latter is a binary variable, with '1' indicating death, and '0' indicating right censored. The function Surv() in package \bold{survival} produces such a matrix.  
  }
  \item{family}{
   Response type (see above).
  }
  \item{epsilon}{
   positive convergence tolerance e; the iterations converge when |dev - dev_{old}|/(|dev| + 0.1) < e.
  }
  \item{maxit}{
   integer giving the maximal number of EM iterations.
  }
  \item{init}{
   vector of initial values for all coefficients (not for intercept). If not given, it will be internally produced. 
  }
  \item{group}{
  a numeric vector, or an integer, or a list indicating the groups of predictors. 
  If \code{group = NULL}, all the predictors form a single group.
  If \code{group = K}, the predictors are evenly divided into groups each with \code{K} predictors.
  If \code{group} is a numberic vector, it defines groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2: \code{(group[2]+1):group[3]}, Group 3: \code{(group[3]+1):group[4]}, .....  
  If \code{group} is a list of variable names, \code{group[[k]]} includes variables in the k-th group. 
  The mixture double-exponential prior is only used for grouped predictors. 
  For ungrouped predictors, the prior is double-exponential or t. 
}
  \item{prior}{
  Two types of priors for the coefficents; spike-and-slab mixture double-exponential (\code{mde}),
  and spike-and-slab mixture t (\code{mt}). 
}
  \item{prior.df}{
  prior degrees of freedom in the \code{mt} prior: default is \code{1} (leading to Cauchy prior).  
}
  \item{ss}{
   a vector of two positive scale values for the spike-and-slab mixture double-exponential or t prior, allowing for different scales for different predictors, 
   leading to different amount of shrinkage. smaller scale values give stronger shrinkage.
}
  \item{Warning}{
  logical. If \code{TRUE}, show the error messages of not convergence and identifiability.
}
  \item{verbose}{
  logical. If \code{TRUE}, print out number of iterations and computational time.
}
}

\details{
  This function sets up Bayesian GLMs and Cox models with spike-and-slab mixture double-exponential or t prior (Bayesian spike-and-slab mixture lasso or ridge), and fits the model 
by incorporating EM steps into the fast coordinate descent algorithm implemented in the package \bold{glmnet}. 
It is an alteration of the function \code{\link{glmnet}} in the package \bold{glmnet}. 
  
  The mixture double-exponential or t prior is used for grouped predictors: \code{coefficients ~ DE(0, (1-d)*ss[1]+d*ss[2]), d ~ Bin(1; theta)}, 
  or \code{coefficients ~ t(prior.df, 0, (1-d)*ss[1]+d*ss[2]), d ~ Bin(1; theta)};
The mixture priors allow different shrinkage (penalty) values for different predictors.
   
}

\value{
  This function returns all outputs from the function \code{\link{glmnet}}, and some other values used in Bayesian hierarchical models.
}

\references{
  Friedman, J., Hastie, T. and Tibshirani, R. (2010) Regularization Paths for Generalized Linear Models via Coordinate Descent. J Stat Softw 33, 1-22.
 
  Simon, N., Friedman, J., Hastie, T. & Tibshirani, R. (2011) Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent. Journal of Statistical Software 39, 1-13.

  Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Generalized Linear Models for Prediction and Associated Genes Detection. Genetics 205, 77-88.
  
  Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Cox Models for Survival Prediction and Associated Genes Detection. Bioinformatics, 33(18), 2799-2807.
  
  Zaixiang Tang, et al. (2017) Spike-and-Slab Lasso Generalized Linear Models for Disease Prediction and Associated Genes Detection by Incorporating Pathway Information. Bioinformatics.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{glmnet}}, \code{\link{glmNet}}, \code{\link{bglm}}, \code{\link{bcoxph}}
}

\examples{

library(BhGLM)
library(survival)
library(glmnet)


N = 1000
K = 100
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5, sigma=1.6) # simulate responses
yy$coefs

# y = yy$y.normal; fam = "gaussian"; y = scale(y)
# y = yy$y.ordinal; fam = "binomial"
y = yy$y.surv; fam = "cox" 

group = NULL
#group = rep(0, 21)
#for(j in 1:length(group)) group[j] = (j-1) * K/(length(group)-1)

# lasso and mixture lasso

f1 = glmNet(x, y, family = fam, ncv = 1) 

ps = f1$prior.scale; ps 
ss = c(ps, 0.5)
f2 = bmlasso(x, y, family = fam, prior = "mde", ss = ss, group = group)

par(mfrow = c(1, 2), mar = c(3, 4, 4, 4))
gap = 10
plot.bh(coefs = f1$coef, threshold = f1$df, gap = gap, main = "lasso") 
plot.bh(coefs = f2$coef, threshold = f2$df, gap = gap, main = "mixture lasso") 

}