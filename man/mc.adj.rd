
\name{mc.adj}
\Rdversion{1.1}
\alias{mc.adj}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Adjust P-values for Multiple Comparisons
}

\description{
Given p-values and adjusted degrees of freedom for multiple variables, 
calculate p-values adjusted for multiple comparisons using several methods.
}

\usage{
mc.adj(p, df.adj = length(p), digits = 10)  
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{p}{ 
   a vector of p-values corresponding to a set of variables (coefficients). 
   These p-values can be obtained from \code{\link{summary.bh}}.
  }
  \item{df.adj}{ 
   adjusted degrees of freedom (effective number of parameters) corresponding to 
   the variables included in the \code{p} vector. 
   This value can be obtained from \code{\link{df.adj}}.
  }
  \item{digits}{
  digits to round with print.
  }
}

\details{
The adjustment methods include the Bonferroni correction ("bonferroni") in which the p-values are multiplied by the number of comparisons. 
Less conservative corrections are also included by Holm (1979) ("holm"), Hochberg (1988) ("hochberg"), Hommel (1988) ("hommel"), 
Benjamini & Hochberg (1995) ("BH" or its alias "fdr"), Benjamini & Yekutieli (2001) ("BY"), and the adjusted Bonferroni correction ("bonferroni.adj") using the adjusted degree of freedom, respectively. 
A pass-through option ("none") is also included. 

The first four methods are designed to give strong control of the family wise error rate. There seems no reason to use the unmodified Bonferroni 
correction because it is dominated by Holm's method, which is also valid under arbitrary assumptions. 

Hochberg's and Hommel's methods are valid when the hypothesis tests are independent or when they are non-negatively associated (Sarkar, 1998; Sarkar and Chang, 1997). 
Hommel's method is more powerful than Hochberg's, but the difference is usually small and the Hochberg p-values are faster to compute. 

The "BH" (aka "fdr") and "BY" method of Benjamini, Hochberg, and Yekutieli control the false discovery rate, the expected proportion of false discoveries amongst the rejected hypotheses. 
The false discovery rate is a less stringent condition than the family wise error rate, so these methods are more powerful than the others. 

The "bonferroni.adj" method is similar to the "bonferroni", but uses the the adjusted degree of freedom as the number of tests. 

}

\value{
  a matrix of adjusted p-values with row names copied from \code{p}.         
}

\references{
Yi, N., Xu, S., Lou, X.Y., and Mallick, H. (2013) Multiple Comparisons in Genetic Association Studies: A Hierarchical Modeling Approach. Statistical Applications in Genetics and Molecular Biology.

Benjamini, Y., and Hochberg, Y. (1995). Controlling the false discovery rate: a practical and powerful approach to multiple testing. Journal of the Royal Statistical Society Series B, 57, 289–300. 

Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics 29, 1165–1188. 

Holm, S. (1979). A simple sequentially rejective multiple test procedure. Scandinavian Journal of Statistics, 6, 65–70. 

Hommel, G. (1988). A stagewise rejective multiple test procedure based on a modified Bonferroni test. Biometrika, 75, 383–386. 

Hochberg, Y. (1988). A sharper Bonferroni procedure for multiple tests of significance. Biometrika, 75, 800–803. 

Shaffer, J. P. (1995). Multiple hypothesis testing. Annual Review of Psychology, 46, 561–576. (An excellent review of the area.) 

Sarkar, S. (1998). Some probability inequalities for ordered MTP2 random variables: a proof of Simes conjecture. Annals of Statistics, 26, 494–504. 

Sarkar, S., and Chang, C. K. (1997). Simes' method for multiple hypothesis testing with positively dependent test statistics. Journal of the American Statistical Association, 92, 1601–1608. 

Wright, S. P. (1992). Adjusted P-values for simultaneous inference. Biometrics, 48, 1005–1013. (Explains the adjusted P-value approach.)   
}
\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{p.adjust}}, \code{\link{summary.glm}}, \code{\link{df.adj}}  
}

\examples{
data(swiss)  

f = bglm(Fertility ~ ., data = swiss, prior.scale = Inf)
p = summary.bh(f)[, 3]
df.adj = df.adj(f)
mc.adj(p = p, df.adj = df.adj)

f = bglm(Fertility ~ ., data = swiss, prior.scale = 1)
p = summary.bh(f)[, 3]
df.adj = df.adj(f)
mc.adj(p = p, df.adj = df.adj)
p = summary.bh(f)[, 3][-1]
df.adj = df.adj(f, vars = 2:length(f$coef))
mc.adj(p = p, df.adj = df.adj)

}

