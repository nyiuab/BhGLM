
\name{ROC}
\Rdversion{1.1}
\alias{ROC}

\title{
Calculating AUC and ROC curves
}

\description{
This function calculates AUC and ROC curves for logistic models.  
}

\usage{
ROC(y, y.fitted, main = " ", xlab="False positive rate", ylab="True positive rate", lwd=2, lty=1, col="black", type="s", pch=20, m=100, plot=TRUE, add=FALSE)
}

\arguments{
  \item{y}{ 
  observed binary response. 
  }
  \item{y.fitted}{ 
  a vector of fitted values (probabilities). 
  }
  \item{plot}{
  logical. if \code{TRUE}, plot ROC curve.
  }
  \item{add}{
  logical. if \code{TRUE}, plot over the existing plot. 
  }
  \item{main, xlab, ylab, lwd, lty, col, type, pch}{ 
  same as in \code{\link[graphics]{plot}}.
}
}

\details{
 
}

\value{ 
It returns AUC, its standard error (se), confidence interval (ci), and p value for testing AUC=0.5, and displays ROC curve.  
}

\references{

}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
\code{\link{pROC}}  
}

\examples{

}

