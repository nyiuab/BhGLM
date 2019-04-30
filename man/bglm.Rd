
\name{bglm}
\Rdversion{1.1}
\alias{bglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayesian Generalized Linear Models (GLMs)
}

\description{
  This function is to set up Bayesian hierarchical GLMs, and to fit the model using the EM-IWLS algorithm. 
  Three types of priors on each coefficient can be used: Student-t (t), double-exponential (de), and spike-and-slab mixture double-exponential (mde).  
  The Bayesian hierarchical GLMs include various models as special cases, e.g., classical GLMs, ridge regression, and Bayesian lasso. 
  It can be used for analyzing general data and large-scale and highly-correlated variables 
  (for example, detecting disease-associated factors and predicting phenotypes).
}

\usage{    
bglm(formula, family = gaussian, data, offset, weights, subset, na.action, 
    start = NULL, etastart, mustart, control = glm.control(epsilon = 1e-04, maxit = 50), 
    prior = c("t", "de", "mde"), group = NULL, method.coef,
    dispersion = 1, prior.sd = 0.5, prior.scale = 0.5, prior.df = 1, prior.mean = 0, ss = c(0.04, 0.5), 
    Warning = FALSE, verbose = FALSE, ...)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{formula, family, data, offset, weights, subset, na.action, start, etastart, mustart, control}{ 
  These arguments are the same as in \code{\link{glm}}.
}
  \item{family}{
  can be all the standard families defined in \code{\link{glm}}. 
  also can be Negative Binomial (\code{NegBin or "NegBin"}).
}
  \item{prior}{
  Several types of priors for the coefficents; double-exponetial (\code{"de"}), Student-t (\code{"t"}), and spike-and-slab mixture double-exponential (\code{"mde"}). 
  The mixture priors are only used for grouped predictors (defined by \code{group}). 
  For ungrouped predictors, the priors are double-exponential or t with scale = prior.scale.  
}
\item{group}{
  a numeric vector, or an integer, or a list indicating the groups of predictors. 
  If \code{group = NULL}, all the predictors form a single group.
  If \code{group = K}, the predictors are evenly divided into groups each with \code{K} predictors.
  If \code{group} is a numberic vector, it defines groups as follows: Group 1: \code{(group[1]+1):group[2]}, Group 2: \code{(group[2]+1):group[3]}, Group 3: \code{(group[3]+1):group[4]}, .....  
  If \code{group} is a list of variable names, \code{group[[k]]} includes variables in the k-th group. 
}
\item{method.coef}{
  jointly updating all coefficients or updating coefficients group by group. The default is jointly updating.
  If \code{method.coef = NULL} or \code{method.coef} is missing, jointly updating.
  If \code{method.coef = K}, update \code{K} coefficients at a time.
  \code{method.coef} can be a numeric vector or a list of variable names (as defined by \code{group}) that defines groups.
  If the number of coefficients is large, the group-by-group updating method can be much faster than the jointly updating.
}
  \item{dispersion}{
  dispersion parameter. If a value is provided, it is the starting value of \code{dispersion}. For Poisson and Binomial models, \code{dispersion} equals 1.
}
  \item{prior.sd}{
  prior standard deviations in the normal priors of the coefficients. If provided, they are starting values of the prior standard deviations for the iterative algorithms. 
}
  \item{prior.scale}{
  scale parameters in the priors \code{de} and \code{t}. \code{prior.scale = c(a1,a2,...,ak)}; if \code{k < J} (\code{J} is the total number of predictors , not counting the intercept),
  it is internally expanded to \code{c(a1,a2,...,ak, rep(ak,J-k))}. 
  For both \code{t} and \code{de} prior, smaller scale induces stronger shrinkage
  (weak shrinkage: \code{prior.scale > 1}, strong shrinkage: \code{prior.scale < 0.3}).  
}
  \item{prior.df}{
  prior degrees of freedom in the \code{t} and \code{mt} priors: default is \code{1} (leading to Cauchy prior).  
}
  \item{prior.mean}{
  prior mean for each coefficient.
}
  \item{ss}{
   a vector of two positive scale values for the mixture prior \code{mde}, allowing for different scales for different predictors, 
   leading to different amount of shrinkage. smaller scale values give stronger shrinkage.
}
  \item{Warning}{
  logical. If \code{TRUE}, show the error messages of not convergence and identifiability.
}
  \item{verbose}{
  logical. If \code{TRUE}, print out number of iterations and computational time.
}
  \item{\dots}{
  further arguments for \code{\link{glm}}.
}
}

\details{
  This function sets up Bayesian hierarchical generalized linear models and fits the model using the EM-IWLS algorithm. 
It is an alteration of the standard function \code{\link{glm}} for fitting classical GLMs, and includes all the \code{\link{glm}} arguments and also some new arguments for the hierarchical modeling. 
  
  Several types of prior distributions on the coefficients can be used:  
for t, \code{coefficients ~ t(prior.df, prior.mean, prior.sd^2)};
for de, \code{coefficients ~ DE(prior.mean, prior.scale)};
for mde, \code{coefficients ~ DE(prior.mean, (1-d)*ss[1]+d*ss[2]), d ~ Bin(1; theta)};

These priors are expressed hierarchically. Therefore, EM algorithms can be used to fit the models.

  These prior distributions include various priors as special cases, for example,
1) If \code{prior.scale = Inf}, the t and double-exponential become Uniform, leading to classical GLMs; 
2) If \code{prior.df = Inf}, the t is Normal \code{N(prior.mean, prior.scale^2)}, leading to ridge regression; 
3) If \code{prior.df = 1}, the t is Cauchy; 
4) If \code{prior.scale = prior.df = 0}, the t is Jeffrey prior.
   
  The EM-IWLS algorithm estimates posterior modes of the parameters, and incorporates an EM algorithm into the iterative weighted least squares (IWLS) as implemented in the function \code{\link{glm}}. 
The algorithm treats the variances \code{prior.sd^2} and the indicator variables (for the mixture prior \code{mde}) as missing data.

  The computational cost of the algorithm is mainly in the step of jointly updating coefficients, because this step needs to calculate the inverse of a large matrix. 
This function implements two methods to update the coefficients, i.e.,  jointly updating all coefficients, or updating a group of coefficients at a time and proceeding by cycling through all the groups at each iteration. 
If the number of coefficients is large (say > 1000), the group-by-group method can be much faster.   

  The models can include both ungrouped and grouped variables. Genetic and other scientific studies routinely generate very many predictor variables, which can be naturally grouped, with predictors in a group being biologically related or statistically correlated. 
In genetic studies, we can use biological pathways to construct groups. The argument \code{group} defines the groups. 

  For genetic association analysis, first run \code{\link{make.main}} to construct design matrix of main effects from genotypic data of genetic markers and/or \code{\link{make.inter}} to make interactions.
}

\value{
  This function returns an object of class "glm", including all outputs from the function \code{\link{glm}}, and also results for the additional parameters in the hierarchical models.
}

\references{
 Yi, N. and Banerjee, S. (2009). Hierarchical generalized linear models for multiple quantitative trait locus mapping. Genetics 181, 1101-1113.

 Yi, N., Kaklamani, V. G. and Pasche, B. (2011). Bayesian analysis of genetic interactions in case-control studies, with application to adiponectin genes and colorectal cancer risk. Ann Hum Genet 75, 90-104. 

 Yi, N. and Ma, S. (2012). Hierarchical Shrinkage Priors and Model Fitting Algorithms for High-dimensional Generalized Linear Models. Statistical Applications in Genetics and Molecular Biology 11 (6), 1544-6115. 
 
 Gelman, A., Jakulin, A., Pittau, M. G. and Su, Y. S. (2008). A weakly informative default prior distribution for logistic and other regression models. Annals of Applied Statistics 2, 1360-1383.

 Rockova, V. and George, E. I. (2014) EMVS: The EM Approach to Bayesian Variable Selection. JASA 109: 828-846.

 Zaixiang Tang, Yueping Shen, Xinyan Zhang, Nengjun Yi (2017) The Spike-and-Slab Lasso Generalized Linear Models for Prediction and Associated Genes Detection. Genetics 205, 77–88.
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{glm}}, \code{\link[MASS]{glm.nb}}
}

\examples{
library(BhGLM)

N = 1000
K = 100
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5, sigma=1.6, theta=2) # simulate responses
yy$coefs


# y = yy$y.normal; fam = gaussian; y = scale(y)
# y = yy$y.ordinal; fam = binomial
y = yy$y.nb; fam = NegBin

# jointly fit all variables (can be slow if m is large)
par(mfrow = c(2, 2), cex.axis = 1, mar = c(3, 4, 4, 4))
gap = 10

ps = 0.05
f1 = bglm(y ~ ., data = x, family = fam, prior = "de", prior.scale = ps)   
plot.bh(f1, vars.rm = 1, threshold = 0.01, gap = gap, main = "de")  

f2 = bglm(y ~ ., data = x, family = fam, prior = "t", prior.scale = ps/1.4)   
plot.bh(f2, vars.rm = 1, threshold = 0.01, gap = gap, main = "t")  

ss = c(0.04, 0.5) 
f3 = bglm(y ~ ., data = x, family = fam, prior = "mde", ss = ss)   
plot.bh(f3, vars.rm = 1, threshold = 0.01, gap = gap, main = "mde")  

# group-wise update (can be much faster if m is large)
ps = 0.05
f1 = bglm(y ~ ., data = x, family = fam, prior = "de", method.coef = 50, prior.scale = ps) # update 50 coefficients at a time
plot.bh(f1, vars.rm = 1, threshold = 0.01, gap = gap)  

}
