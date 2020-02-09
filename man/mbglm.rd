
\name{mbglm}
\Rdversion{1.0}
\alias{mbglm}

\title{
  Fitting Bayesian GLMs for Many Responses 
}

\description{
  \code{mbglm} fits Bayesian GLMs separately for many responses, and \code{summary.mbglm} summarizes the fitted models.    
}

\usage{  
mbglm(y, formula, data, family=NegBin(), prior=Student(0,1),
      min.p=0, verbose=TRUE)
      
summary.mbglm(object, vr.name=NULL, sort=FALSE)
}

\arguments{
  \item{y}{
  a matrix of responses; each column is a response and rows are samples. 
}
  \item{formula}{ 
  a one-sided formula of the form \code{~ x} (i.e., the respose is omitted); the right side of \code{~} is the same as in \code{\link{bglm}}. 
}
  \item{data, family, prior, verbose}{ 
  These arguments are the same as in \code{\link{bglm}}.
}
  \item{min.p}{
  a value in [0, 1). The responses with the proportion of non-zero values > min.p are analyzed.
}

  \item{object}{ 
  an object from \code{\link{mbglm}}.
}
  \item{vr.name}{
  name of a variable or response. 
}
  \item{sort}{
  sort by the adjusted p-values into ascending order.
}

}

\details{
 This function analyzes the responses in \code{y} by repeated calls to \code{\link{bglm}}. 
 
}

\value{
  The function \code{mbglm} returns a list including \code{fit},  \code{responses}, \code{variables}, and \code{call}:

  \item{fit}{fitted models for all the responses;}
  \item{responses}{names of all the responses;}
  \item{variables}{names of all covariates;}

  If \code{vr.name=NULL}, the function \code{summary.mbglm} returns a data frame consisting of responses, variables, estimates, standard deviations, p-vlaues and adjusted p-vlaues (using FDR method) for all coefficients. 
  If \code{vr.name} is specified, \code{summary.mbglm} returns a list with two components, for all coefficients and for the specified coefficient, respectively.
}

\references{
  
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{bglm}}
}

\examples{

library(BhGLM)
library(NBZIMM)

# load data in NBZIMM
data(Romero)
names(Romero)

otu = Romero$OTU; dim(otu)
sam = Romero$SampleData; dim(sam)
colnames(sam)

Days = sam$GA_Days; Days = scale(Days)
Age = sam$Age; Age = scale(Age)
Race = sam$Race
preg = sam$pregnant; table(preg)

# calculate size factor using the cumulative sum scaling (CSS) method in metagenomeSeq
library(metagenomeSeq)
obj = newMRexperiment(counts=t(otu), phenoData=NULL, featureData=NULL,
                      libSize=NULL, normFactors=NULL)
obj = cumNorm(obj, p=cumNormStat(obj)) #set p=1: generate total reads
s = normFactors(obj)
s = unlist(s)

# analyze taxa with nonzero proportion > min.p
f = mbglm(y=otu, formula=~ Days + Age + Race + preg + offset(log(s)), 
          family=NegBin(), prior=Student(0, 1), min.p=0.1)

out = summary.mbglm(f, vr.name="preg", sort=T)[[2]]
#out = summary.mbglm(f, vr.name="Finegoldia.magna")[[2]]
out
coefs = out[,1]
sds = out[,2]
padj = out[,4]
plot.bh(coefs=coefs, sds=sds, pvalues=padj, threshold=0.001, gap=1000)

}