# BhGLM: Bayesian hierarchical GLMs and survival models, with applications to Genomics and Epidemiology 

# Overview

This R package provides functions for setting up and fitting various Bayesian hierarchical models (generalized linear models (GLMs), Cox survival models, negative binomial models, and ordered logistic or probit regressions), for numerically and graphically summarizing the fitted models, and for evaluating the predictive performance. Several types of priors on the coefficients can be used: double-exponential, Student-t, mixture double-exponential, and mixture t. The models are fitted by using fast algorithms for estimating posterior modes rather than MCMC. The methods can be used to analyze not only general data but also large-scale molecular data (i.e., detecting disease-associated genes or variants, predictive and prognostic modeling of diseases and traits, etc).

Author: Nengjun Yi <nyi@uab.edu>;  Maintainer: Nengjun Yi <nyi@uab.edu>

# Installation

Two ways to install the package in R:

1. With Vignettes (must install packages: devtools, knitr, R.rsp)
```{r}
devtools::install_github("nyiuab/BhGLM", build_opts = c("--no-resave-data", "--no-manual"), force = T)
```
2. Without Vignettes (must install package: devtools) 
```{r}
devtools::install_github("nyiuab/BhGLM", force = T)
```
