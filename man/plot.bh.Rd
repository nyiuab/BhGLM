
\name{plot.bh}
\Rdversion{1.1}
\alias{plot.bh}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Graphically Summarizing Bayesian Model Fits
}

\description{
This function is to graphically summarize Bayesian hierarchical model fits by plotting the coefficients (or odds ratios) and the standard deviation (2*sd).
It also can show the p-values. 
}

\usage{
plot.bh(object = NULL, coefs, sds = NULL, pvalues = NULL, vars.rm = NULL, 
        threshold = 0.05, show.all.vars = FALSE, show.pvalues = TRUE, gap = 0, 
        main = " ", cex.main = 0.9, xlim = NULL, cex.var = 0.8, cex.pts = 1, 
        pch.pts = 20, type = "p", lwd = 1, lty = 1, line = 0, col.pts = "black", 
        OR = FALSE, add = FALSE) 
}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{
   an object from \code{\link{bglm}} or \code{\link{bpolr}} or \code{\link{bcoxph}}. 
   If specified, the function get the estimates of coefficints, sd, and p-values, otherwise (i.e., \code{object = NULL}), the coefficints, sd and p-values
   should be given in the arguments \code{coefs}, \code{sds}, and \code{pvalues}. 
}
  \item{coefs, sds, pvalues}{
   vector of coefficints (\code{coefs}), vector of standard deviations (or n*2 matrix of interval estimates) (\code{sds}) and vector of p-values (\code{pvalues}) for variables. 
   The names of variables can be given via coefs.
   If \code{sds} is a vector, interval estimates of coefficients are [coefs - 2*sds, coefs + 2*sds]. 
   If \code{sds} is a n*2 matrix, interval estimates of coefficients are [sds[,1], sds[,2]]; this is particularly useful for summarizing MCMC samples.  
}
  \item{vars.rm}{
   a vector of variables to be removed; default is \code{NULL}, display all the variables (including the intercept).
}
  \item{threshold}{
   a p-value or a positive integer value to determine variables to be shown: variables with p-value < \code{threshold} or top \code{threshold} coefficients. 
}
  \item{show.all.vars}{
  logical. If \code{TRUE},show names of all variables. The default is \code{FALSE}.
}
  \item{show.pvalues}{
  logical. If \code{TRUE}, show p-values for variables. The default is \code{TRUE}. The p-values are shown on the right side of the plot.
}
   \item{gap}{
   a value for the distance between two flanking significant variables. 
}
  \item{main, cex.main, xlim, type, lwd, lty, line}{ 
  These arguments are the same as in \code{\link[graphics]{plot}} and \code{\link[graphics]{par}}.
}
  \item{cex.var}{
  the fontsize of the varible names, default = \code{0.8}.
}
  \item{cex.pts}{
  the size of points, default = \code{1}.
}
  \item{col.pts}{
  color of points and segments, default is black. It can be a vector with two elements, which will use specified colors for different 
  variables determined by the augument \code{threshold}.
}
  \item{pch.pts}{
  symbol of points, default is solid dot.
}
  \item{OR}{
  logical. If \code{TRUE}, show odds ratios for variables. The default is \code{FALSE}. 
}
  \item{add}{
  logical. if \code{TRUE}, plot over the existing plot. The default is \code{FALSE}. 
}
  
}

\details{
This function plots the estimates of coefficients, intervals and p-values from a fitted GLMs or Cox model.
It uses different colors to distinguish between siginificant and insignificant variables based on a threshold.  
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link[graphics]{plot}}, \code{\link[graphics]{par}}
}

\examples{

}

