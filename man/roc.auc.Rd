
\name{roc.auc}
\Rdversion{1.1}
\alias{roc.auc}

\title{
Calculating AUC and ROC Curves
}

\description{
This function calculates AUC and ROC curves for logistic models.  
}

\usage{
roc.auc(y, y.fitted, main = " ", xlab="False positive rate", ylab="True positive rate", lwd=2, lty=1, col="black", type="s", pch=20, m=100, plot=TRUE, add=FALSE, legend=TRUE)
}

\arguments{
  \item{y}{ 
  observed binary response. 
  }
  \item{y.fitted}{ 
  a vector of fitted values (probabilities). 
  }
  \item{main, xlab, ylab, lwd, lty, col, type, pch}{ 
  same as in \code{\link[graphics]{plot}}.
  }
  \item{plot}{
  logical. if \code{TRUE}, plot ROC curve.
  }
  \item{add}{
  logical. if \code{TRUE}, plot over the existing plot. 
  }
  \item{legend}{
  logical. if \code{TRUE}, show auc and pvlaue. 
  }
}

\details{
 
}

\value{ 
It returns auc, its standard error (se), and p value for testing auc=0.5, and displays ROC curve.  
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

