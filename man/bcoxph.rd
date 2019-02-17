
\name{bcoxph}
\Rdversion{1.1}
\alias{bcoxph}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayesian Hierarchical Cox Proportonal Hazards Model
}

\description{
  This function is to set up Bayesian hierarchical Cox proportonal hazards model and to fit the model using the EM Newton-Raphson algorithm. 
  As in \code{\link{bglm}}, several types of priors on the coefficients can be used.   
  The Bayesian hierarchical Cox model includes classical Cox model and ridge Cox regression as special cases. 
  It can be used for analyzing general survival data and large-scale and highly-correlated variables 
  (for example, detecting disease-associated factors and predicting survival times).
}

\usage{    
bcoxph(formula, data, weights, subset, na.action, init, 
       control = coxph.control(eps = 1e-04, iter.max = 50), 
       ties = c("breslow", "efron"), tt,  
       prior = c("de", "t", "mde", "mt"), group = NULL, method.coef, 
       prior.sd = 1, prior.scale = 0.5, prior.df = 1, prior.mean = 0, ss = c(0.04, 0.5), 
       Warning = FALSE, verbose = TRUE, ...)
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula, data, weights, subset, na.action, init, control, ties, tt}{ 
  These arguments are the same as in the function \code{\link{coxph}} in the package \bold{survival}.
}
  \item{prior, group, method.coef, prior.sd, prior.scale, prior.df, prior.mean, ss,  
        Warning, verbose}{ 
  These arguments are the same as in the function \code{\link{bglm}}.
}
\item{\dots}{
  further arguments for \code{\link{coxph}}.
}

}

\details{
  This function incorporates four types of prior distributions into the conventional Cox proportonal hazards model.
It is an alteration of the standard function \code{\link{coxph}} for fitting classical Cox proportonal hazards model, and includes all the \code{\link{coxph}} arguments and also some new arguments for the hierarchical modeling. 
The standard procedure for fitting the conventional Cox proportonal hazards model is the Newton-Raphson algorithm, 
The function incorporates an EM algorithm for updating hyper-parameters 
into the standard Newton-Raphson procedure as implemented in the function \code{\link{coxph}}.   
}

\value{
  This function returns an object of class "coxph.penal" and "coxph", including all outputs from the function \code{\link{coxph}} and also results for the additional parameters in the hierarchical models.
}

\references{
  van Houwelinggen, H.G. & Putter, H. Dynamic Prediction in Clinical Survival Analysis, (CRC Press, 2012).

  Yi, N. and Banerjee, S. (2009). Hierarchical generalized linear models for multiple quantitative trait locus mapping. Genetics 181, 1101-1113.

  Yi, N., Kaklamani, V. G. and Pasche, B. (2011). Bayesian analysis of genetic interactions in case-control studies, with application to adiponectin genes and colorectal cancer risk. Ann Hum Genet 75, 90-104. 

  Yi, N. and Ma, S. (2012). Hierarchical Shrinkage Priors and Model Fitting Algorithms for High-dimensional Generalized Linear Models. Statistical Applications in Genetics and Molecular Biology 11 (6), 1544-6115. 
 
  Gelman, A., Jakulin, A., Pittau, M. G. and Su, Y. S. (2008). A weakly informative default prior distribution for logistic and other regression models. Annals of Applied Statistics 2, 1360-1383.

  Armagan, A., Dunson, D. and Lee, J. (2010) Bayesian generalized double Pareto shrinkage. Biometrika
  
  Gelman, A. et al. (2014) Bayesian Data Analysis. Chapman & Hall/CRC Press, New York. 

  Rockova, V. and George, E. I. (2014) EMVS: The EM Approach to Bayesian Variable Selection. JASA 109: 828-846.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{coxph}}, \code{\link{bglm}} 
}

\examples{

library(BhGLM)
library(survival)

N = 1000
K = 100
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5) # simulate responses
yy$coefs

y = yy$y.surv
d = table(y[,2]); d[1]/sum(d) # cencoring proportion


# jointly update
# Compare with conventional and ridge Cox model
par(mfrow = c(2, 2), cex.axis = 1, mar = c(3, 4, 4, 4))
gap = 10

ps = 0.05
f1 = bcoxph(y ~ ., data = x, prior = "de", prior.scale = ps)
# summary.bh(f1)
plot.bh(f1, threshold = 0.01, gap = gap, main = "Cox with double-exponential")

f2 = bcoxph(y ~ ., data = x, prior = "t", prior.scale = ps/1.4)
# summary.bh(f2)
plot.bh(f2, threshold = 0.01, gap = gap, main = "Cox with t")

ss = c(0.04, 0.5) 
f3 = bcoxph(y ~ ., data = x, prior = "mde", ss = ss)
# summary.bh(f3)
plot.bh(f3, threshold = 0.01, gap = gap, main = "Cox with mixture double exponential") 

ss = c(0.04, 0.5) 
f4 = bcoxph(y ~ ., data = x, prior = "mt", ss = ss)
# summary.bh(f4)
plot.bh(f4, threshold = 0.01, gap = gap, main = "Cox with mixture t") 


# group-wise update
ps = 0.05
f1 = bcoxph(y ~ ., data = x, prior = "de", method.coef = 50, prior.scale = ps)
#summary.bh(f1)
plot.bh(coefs = f2$coefficients, threshold = 10, gap = gap)  

}