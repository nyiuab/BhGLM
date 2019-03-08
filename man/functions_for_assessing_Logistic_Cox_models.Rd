
\name{functions_for_assessing_Logistic_Cox_models}
\Rdversion{1.1}
\alias{functions_for_assessing_Logistic_Cox_models}

\title{
Functions for summarizing a fitted logistic or Cox survival model or cross-validation
}

\description{
These functions display survival curves, C-index,  curve of Area Under the Curve values over time, and prediction error curve, from a fitted Cox survival model or cross-validation.  
}

\usage{
ROC(y, y.fitted, main = " ", xlab = "False positive rate", ylab = "True positive rate", lwd = 2, lty = 1, col = "black", plot = TRUE, add = FALSE)

surv.curves(y, lp, probs = 0.50, main = " ", lwd = 2, lty = 1, col = c("black", "red"), add = FALSE) 

Cindex(y, lp)

aucCox(y, lp, main = " ", lwd = 2, lty = 1, col = "black", add = FALSE) 

peCox(y, lp, FUN = c("Brier", "KL"), main = "", lwd = 2, lty = 1, col = "black", add = FALSE) 
}

\arguments{
  \item{y}{ 
  For function \code{ROC}, observed binary response; For others,  observed survival response, a '\code{Surv}' object. 
  }
  \item{y.fitted}{ 
  a vector of fitted values from a fitted logistic model or cross-validation. 
  }
  \item{lp}{ 
  a vector of prognostic index (linear predictor) from a fitted Cox model or cross-validation. 
  }
  \item{probs}{ 
  numeric value or vector of probabilities with values in [0,1] for grouping the patients based on \code{lp}. Default \code{probs = 0.50}: partition samples to two groups according to the percentile of 50% of \code{lp}.   
  }
  \item{FUN}{
  The error function, either "KL" (default) for Kullback-Leibler or "Brier" for Brier score.
  }
  \item{add}{
  logical. if \code{TRUE}, plot over the existing plot. The default is \code{FALSE}.
  }
  \item{main, lwd, lty, col}{ 
  same as in \code{\link[graphics]{plot}}.
}
}

\details{
The function \code{ROC} calculates AUC (area under the ROC curve), its interval estimates, p_value for testing AUC = 0.5, etc., and plots ROC curve.

The function \code{surv.curves} plots survival curves of groups of individuals and test the difference between curves using the log-rank method in \code{\link{survdiff}}.

\code{Cindex} calculates Harrell's c-index.

\code{aucCox} calculates model-free curve of Area Under the Curve values over time.

\code{peCox} calculates prediction error curve. It is an alteration of \code{\link{dynpred}}, which should be installed.  
}

\value{ 

}
\references{
van Houwelinggen, H.G. & Putter, H. Dynamic Prediction in Clinical Survival Analysis, (CRC Press, 2012).
}

\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
\code{\link{dynpred}}  
}

\examples{

}

