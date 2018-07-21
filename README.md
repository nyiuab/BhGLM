# BhGLM: Bayesian hierarchical GLMs and survival models, with applications to Genomics and Epidemiology 

This R package provides functions for setting up and fitting various Bayesian hierarchical models (generalized linear models (GLMs), Cox survival models, negative binomial models, and ordered logistic or probit regressions), for numerically and graphically summarizing the fitted models, and for evaluating the predictive performance. Several types of priors on the coefficients can be used: double-exponential, Student-t, mixture double-exponential, and mixture t. The models are fitted by using fast algorithms for estimating posterior modes rather than MCMC. The methods can be used to analyze not only general data but also large-scale molecular data (i.e., detecting disease-associated genes or variants, predictive and prognostic modeling of diseases and traits, etc).
       

License: GPL

Author: Nengjun Yi <nyi@uab.edu>

Maintainer: Nengjun Yi <nyi@uab.edu>; Xinyan Zhang <xzhang@georgiasouthern.edu>
