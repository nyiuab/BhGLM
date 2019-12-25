
\name{bglm}
\Rdversion{1.1}
\alias{bglm}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Bayesian Generalized Linear Models (GLMs)
}

\description{
  This function is to set up Bayesian hierarchical GLMs, and to fit the model using the EM-IWLS algorithm. 
  Three types of priors on each coefficient can be used: Student-t, double-exponential, and spike-and-slab mixture double-exponential.  
  The Bayesian hierarchical GLMs include various models as special cases, e.g., classical GLMs, ridge regression, and Bayesian lasso. 
  It can be used for analyzing general data and large-scale and highly-correlated variables 
  (for example, detecting disease-associated factors and predicting phenotypes).
}

\usage{    
bglm(formula, family = gaussian, data, offset, weights, subset, na.action, 
    start = NULL, etastart, mustart, control = glm.control(epsilon = 1e-04, maxit = 50), 
    prior = Student(0, 0.5, 1), group = NULL, method.coef,
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
  Prior distribustions for the coefficents. Three types of priors can be used; Student-t: \code{Student(mean, scale, df)} (default: mean=0, scale=0.5, df=1), Double-exponetial: \code{De(mean, scale)} (default: mean=0, scale=0.5), and mixture double-exponential: \code{mde(mean, s0, s1)}, s0 < s1, default: mean=0, s0=0.04, s1=0.5. 
  The \code{mean}, \code{scale} and \code{df} can be a vector. For example, scale = c(a1,a2,...,ak); if k < the total number of predictors, it is internally expanded to c(a1,a2,...,ak, rep(ak,J-k)).
}
\item{group}{
  a numeric vector, or an integer, or a list definining the groups of predictors. Only used for \code{mde}. 
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
  This function sets up Bayesian hierarchical GLMs and fits the model using the EM-IWLS algorithm. It is an alteration of the standard function \code{\link{glm}} for classical GLMs, and includes all the \code{\link{glm}} arguments and also some new arguments for Bayesian modeling. 

  For \code{Student(mean, scale)} and \code{De(mean, scale)}, \code{scale} is modified internally: for a predictor with only one value nothing is changed; for a predictor x with exactly two unique values, we divide the user-specified scale \code{scale} by the range of x; for a predictor x with more than two unique values, we divide the user-specified scale by sd(x). For gaussian models, \code{scale} is further multiplied by sd(y).  
  
  The argument \code{group} is used for the prior \code{mde(mean, s0, s1)}, allowing for group-specific inclusion probabilities and thus incorporating group information (e.g. bilogical pathways). With the prior \code{mde(mean, s0, s1)}, the predictors can be grouped or ungrouped. For ungrouped predictors, the prior is double-exponential with scale \code{s1}. 
   
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
 
 Zaixiang Tang, Yueping Shen, Yan Li, Xinyan Zhang, Jia Wen, Chen'ao Qian, Wenzhuo Zhuang, Xinghua Shi, and Nengjun Yi (2018) Group Spike-and-Slab Lasso Generalized Linear Models for Disease Prediction and Associated Genes Detection by Incorporating Pathway Information. Bioinformatics 34(6): 901-910. 
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
f1 = bglm(y ~ ., data = x, family = fam, prior = De(0, ps))   
plot.bh(f1, vars.rm = 1, threshold = 0.01, gap = gap, main = "de")  

f2 = bglm(y ~ ., data = x, family = fam, prior = Student(0, ps/1.4))   
plot.bh(f2, vars.rm = 1, threshold = 0.01, gap = gap, main = "t")  

ss = c(0.04, 0.5) 
f3 = bglm(y ~ ., data = x, family = fam, prior = mde(0, ss[1], ss[2]))   
plot.bh(f3, vars.rm = 1, threshold = 0.01, gap = gap, main = "mde")  

# group-wise update (can be much faster if m is large)
ps = 0.05
f1 = bglm(y ~ ., data = x, family = fam, prior = De(0, ps), method.coef = 50) # update 50 coefficients at a time
plot.bh(f1, vars.rm = 1, threshold = 0.01, gap = gap)  

}
