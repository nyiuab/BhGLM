

#*******************************************************************************

bnegbin.fit <- function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL, mustart = NULL, 
                   offset = rep(0, nobs), control = glm.control(), intercept = TRUE,
                   prior = "de", group = NULL, method.coef = 1, dispersion = 1, prior.mean = 0, prior.sd = 0.5,
                   prior.scale = 1, prior.df = 1, ss = c(0.05, 0.1), Warning = FALSE)
{
  if (!requireNamespace("MASS")) install.packages("MASS")
  library(MASS)
  nobs <- NROW(y)
  if (is.null(weights)) weights <- rep.int(1, nobs)
  off <- offset
  if (is.null(offset)) off <- rep.int(0, nobs)
  th <- NegBin()$theta
  
  devold <- -100
  conv <- FALSE
  for (iter in 1:control$maxit){
    family <- NegBin(theta = th)
    fit <- suppressWarnings( bglm.fit(x = x, y = y, weights = weights, start = start,
                  etastart = etastart, mustart = mustart, offset = offset,
                  family = family, control = glm.control(maxit = 2), intercept = TRUE,
                  prior = prior, group = group, method.coef = method.coef,
                  dispersion = dispersion, prior.mean = prior.mean, prior.sd = prior.sd,
                  prior.scale = prior.scale, prior.df = prior.df, ss = ss, Warning = FALSE) )
    start <- fit$coefficients
    etastart <- fit$linear.predictors
    prior.sd <- fit$prior.sd
    prior.scale <- fit$prior.scale
    mu <- fit$fitted.values
    dispersion <- fit$dispersion
    
    dev <- fit$deviance
    if( abs(dev - devold)/(0.1 + abs(dev)) < control$eps ) {
      conv <- TRUE
      break
    }
    else {
      devold <- dev
      th <- suppressWarnings( theta.ml(y = y, mu = mu, n = sum(weights), weights = weights, limit = 10, trace = FALSE) )
      if (is.null(th)) th <- family$theta 
    }
  }   # theta.ml can stop sometimes!
  if (Warning & (iter == control$maxit)) warning("alternation limit reached")
  fit$theta <- family$theta 
  fit$iter <- iter
  fit$family <- family

  fit
}

#*******************************************************************************
# from MASS

NegBin <- function (theta = 3, link = "log")
{
    linktemp <- substitute(link)
    if (!is.character(linktemp))
        linktemp <- deparse(linktemp)
    if (linktemp %in% c("log", "identity", "sqrt"))
        stats <- make.link(linktemp)
    else if (is.character(link)) {
        stats <- make.link(link)
        linktemp <- link
    }
    else {
        if (inherits(link, "link-glm")) {
            stats <- link
            if (!is.null(stats$name))
                linktemp <- stats$name
        }
        else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",
            linktemp))
    }
    .Theta <- theta
    env <- new.env(parent = .GlobalEnv)
    assign(".Theta", theta, envir = env)
    variance <- function(mu) mu + mu^2/.Theta
    validmu <- function(mu) all(mu > 0)
    dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1,
        y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
    aic <- function(y, n, mu, wt, dev) {
        term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
            lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
            lgamma(.Theta + y)
        2 * sum(term * wt)
    }
    initialize <- expression({
        if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
        n <- rep(1, nobs)
        mustart <- y + (y == 0)/6
    })
    simfun <- function(object, nsim) {
        ftd <- fitted(object)
        rnegbin(nsim * length(ftd), ftd, .Theta)
    }
    environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
    famname <- paste("NegBin(", format(round(theta,
        4)), ")", sep = "")
    structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
        linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
        aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
        validmu = validmu, valideta = stats$valideta, simulate = simfun, theta = .Theta),
        class = "family")
}

#*******************************************************************************

