
\name{df.adj}
\Rdversion{1.1}
\alias{df.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculating Adjusted Degree of Freedom (Effective Number of Parameters) in Bayesian Hierarchical GLMs
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

data(fake.cv) #load simulated common variants dataset
y = fake.cv[, 1] # binary trait
covar = fake.cv[, c(3, 4)]   # covariates
x.main = make.main(geno = fake.cv[, -c(1:4)], model = "Cockerham", fill.missing = T, geno.order = T)
x = cbind(covar, x.main)
f = bglm(y ~ ., data = x, family = binomial(link = "logit"), prior = "t", prior.scale = 0.5)
plot.bh(f, gap = 5)
df.adj(f)  # for the whole fitted model
df.adj(f, f$group.vars[[1]]) # for the predictors in the first group

}

