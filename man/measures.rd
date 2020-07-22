
\name{measures}
\Rdversion{1.1}
\alias{measures}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Evaluating Fitted Models
}

\description{
These functions provide various measures to evaluate fitted models. 
}

\usage{
measure.bh(object, new.x, new.y, new.offset) 

measure.glm(y, y.fitted, family, dispersion=1)

measure.polr(y, y.fitted)

measure.cox(y, lp)

surv.curves(y, lp, probs=0.50, mark.time=FALSE, main=" ", lwd=2, lty=1, col=c("black", "red"), add=FALSE) 

aucCox(y, lp, main=" ", lwd=2, lty=1, col="black", add=FALSE) 

peCox(y, lp, FUN=c("Brier", "KL"), main="", lwd=2, lty=1, col="black", add=FALSE) 

}

%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{object}{ 
  a fitted object.   
  }
  \item{new.x}{
  data frame or matrix of new values for variables used in \code{object}. for object from \code{\link{bglm}} or \code{\link{bpolr}} or \code{\link{bcoxph}}, it is data frame; for object from \code{\link{glmNet}} or \code{\link{bmlasso}}, it is matrix.
}
  \item{new.y}{  
  vector of new response values corresponding to \code{new.x}. 
  If \code{new.x} or \code{new.y} are omitted, the fitted linear predictors are used for prediction. 
}
  \item{new.offset}{ 
  data frame or vector of offset values for new data points. If \code{new.x} includes offset, do not need to set \code{new.offset}.  
  }
  \item{y}{
  observed response values.
  }
  \item{y.fitted}{
  predicted (estimated) response values for GLMs or probabilties of response values for ordinal models from a fitted model or cross-validation.  
  }
  \item{family}{
  family in GLMs.
  }
  \item{dispersion}{
  dispersion in GLMs. theta in negative binomial model.
  }
  \item{lp}{ 
  a vector of prognostic index (linear predictor) from a fitted Cox model or cross-validation. 
  }
  \item{probs}{ 
  numeric value or vector of probabilities with values in [0,1] for grouping the patients based on \code{lp}. Default \code{probs = 0.50}: partition samples to two groups according to the percentile of 50% of \code{lp}.   
  }
  \item{mark.time}{ 
  controls the labeling of the curves. If set to FALSE, no labeling is done. If TRUE, then curves are marked at each censoring time.
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
The functions, \code{measure.bh, measure.glm, measure.polr, measure.cox}, provide various measures to evaluate fitted models. For all GLMs and polr models, return: \code{deviance}: estimate of deviance, \code{mse}: estimate of mean squared error. For binomial and polr models, also return: \code{auc}: area under ROC curve, \code{misclassification}: estimate of misclassification. For Cox models, return: \code{deviance}: estimate of deviance, \code{Cindex}: concordance index.

The function \code{surv.curves} plots survival curves of groups of individuals and test the difference between curves using the log-rank method in \code{\link{survdiff}}.

\code{aucCox} calculates model-free curve of Area Under the Curve values over time.

\code{peCox} calculates prediction error curve. It is an alteration of \code{\link{dynpred}}, which should be installed. 

}

\value{

}

\references{
  Steyerberg, E. W., 2009 Clinical Prediction Models: A Practical Approch to Development, Validation, and Updates. Springer, New York.
  
  van Houwelinggen, H.G. & Putter, H. Dynamic Prediction in Clinical Survival Analysis, (CRC Press, 2012).
}


\author{
  Nengjun Yi, nyi@uab.edu
}

\seealso{
  \code{\link{cv.bh}}, \code{\link{predict.glm}}, \code{\link{predict.coxph}}
}

\examples{
library(BhGLM)

N = 1000
K = 50
x = sim.x(n=N, m=K, corr=0.6) # simulate correlated continuous variables  
h = rep(0.1, 4) # assign four non-zero main effects to have the assumed heritabilty 
nz = as.integer(seq(5, K, by=K/length(h))); nz
yy = sim.y(x=x[, nz], mu=0, herit=h, p.neg=0.5, sigma=1.6) # simulate responses
yy$coefs

y = yy$y.ordinal; fam = binomial

# partition the data into two parts; 
# fit a model using the first part and evaluate the prediction using the second part 
x1 = x[1:(N/2),]; y1 = y[1:(N/2)]
x2 = x[(N/2):N,]; y2 = y[(N/2):N]

f1 = bglm(y1 ~ ., data=x1, family=fam, prior=De(0,0.5))   
plot.bh(f1, vars.rm=1, threshold=0.01, gap=10)  

measure.bh(f1, x2, y2)


y = yy$y.surv
y1 = y[1:(N/2)]
y2 = y[(N/2):N]
f1 = bcoxph(y1 ~ ., data=x1, prior=De(0,0.5))   
plot.bh(f1, threshold=0.01, gap=10)  

measure.bh(f1, x2, y2)
 
}

