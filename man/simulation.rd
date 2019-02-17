
\name{simulation}
\Rdversion{1.1}
\alias{simulation}

\title{
Simulating Predictor Variables and Response Data and Summarizing the Analysis of Simulated Data
}

\description{
These functions simulate predictor variables and response outcomes and summarize the analysis of simulated data.  
}

\usage{
sim.x(n, m, group = NULL, corr = 0.6, v = rep(1, m), p = 0.5, genotype = NULL)

sim.eta(x, mu = 0, coefs = NULL, herit = 0.1, sigma = 1, p.neg = 0.5)

sim.y(x, mu = 0, coefs = NULL, herit = 0.1, p.neg = 0.5, sigma = 1, quantiles = 0.5, , theta = 3, df = 3)

sim.out(coefs.p, coefs.est, alpha = c(0.05, 0.01))
}

\arguments{
  \item{n}{ 
  number of simulated data points (individuals). 
  }
  \item{m}{ 
  number of simulated continuous variables or discrete genetic markers. 
  }
  \item{group}{
  a numeric vector, or an integer, or a list indicating the groups of variables. 
  If \code{group = NULL}, all the variables form a single group.
  If \code{group = K}, the predictors are evenly divided into groups each with \code{K} predictors.
  If \code{group} is a numberic vector, it defines groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2: \code{(group[2]+1):group[3]}, Group 3: \code{(group[3]+1):group[4]}, .....  
  If \code{group} is a list of variable names, \code{group[[k]]} includes variables in the k-th group. 
  }
  \item{corr}{ 
  correlation between variables. If \code{length(corr)=1}, within-group and between correlations are \code{corr} and zero. 
  If \code{length(corr)=2}, within-group and between correlations are \code{corr[1]} and \code{corr[2]}.
  }
  \item{v}{
  variances of simulated variables.
  }
  \item{p}{
  minor allelic frequencies for simulated markers. 
  }
  \item{genotype}{ 
  transform some continuous variables to three-level genotypes. 
  }
  \item{x}{ 
  a design matrix of simulated variables. 
  }
  \item{mu}{ 
  intercept. 
  }
  \item{coefs}{ 
  coefficients of variables. If \code{coefs = NULL}, the coefficients are calculated by \code{herit}.
  If \code{length(coefs) < ncol(x)}, all other coefficients are set to zero, i.e., \code{coefs} is expanded to 
  \code{c(coefs, 0, ..., 0)}. The linear predictors equal \code{mu + x * coefs}. 
  }
  \item{herit}{
  heritabilities of variables (proportions of variance explained by variables), which is used to calculate the coefficients. 
  If \code{coefs} is given, \code{herit} is not used.
  If \code{herit} is a single value (for example, \code{herit = 0.05}), it is the total heritability of all variables.
  If \code{herit} is a vector, it gives the heritabilities of the corresponding variables 
  (for example, if \code{herit = c(h1, h2, h3)}, the heritabilities are h1, h2 and h3, for the first three variables, 
  and zero for other variables). 
  }
  \item{p.neg}{ 
  proportion of negative coefficients. 
  }
  \item{sigma}{ 
  residual standard deviation for normal or t response. 
  }
  \item{quantiles}{ 
  quantiles for generating binary or ordinal responses. 
  }
  \item{theta}{ 
  shape parameter for negative binomial or beta responses. 
  }
  \item{df}{ 
  degree of freedom of t response. 
  }
  \item{coefs.p}{ 
  a matrix of p-values of coefficients. The rows and columns are coefficients and simulations, respectively. 
  }
  \item{coefs.est}{ 
  a matrix of coefficient estimates. The rows and columns are coefficients and simulations, respectively. 
  }
  \item{alpha}{ 
   significance levels for calculating power.
  }
}

\details{
The function \code{sim.x()} simulates \code{m} variables or genotypes of \code{m} markers for \code{n} individuals.

The function \code{sim.y()} calculates coefficient values if \code{coefs = NULL} 
and the linear predictors \code{eta = mu + x * coefs}, simulates \code{n} normal phenotypes with mean \code{eta} 
and variance \code{sigma^2}, categorizes the normal phenotypes to binary or ordinal phenotypes, 
simulates survival outcome, count outcomes from Poisson or negative binomial, Student-t outcomes, and beta outcome. 

The function \code{sim.out()} calculates statistical powers and estimates for all variables.
   
}

\value{
\code{sim.x()} returns a \code{n x m} data frame of continuous values or genotypes 0, 1, 2.

\code{sim.y()} returns a list of normal outcome \code{y.normal}, 
binary or ordinal outcome \code{y.ordinal}, poisson outcome \code{y.poisson}, negative binomial outcome \code{y.nb}, t outcome \code{y.t}, beta outcome \code{y.beta}, survival outcome \code{y.surv},
linear predictor values \code{eta}, coefficients \code{coefs}, 
residual standard deviation \code{sigma} and heritabilities \code{herit}.

\code{sim.out()} returns a list of power and estimate for each variable.
  
}
\references{

}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  
}

\examples{
See examples in the functions \code{\link{bglm}}, \code{\link{bpolr}}, \code{\link{bcoxph}}.

}

