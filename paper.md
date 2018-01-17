---
title: 'registr: Registration for Exponential Family Functional Data'
date: "17 January 2018"
tags:
- statistical analysis
- R
authors:
- name: Julia Wrobel
orcid: 000-0001-6783-1421
affiliation: 1
affiliations:
name: Columbia University
index: 1
bibliography: paper.bib
---

# Summary

- what is exponential family functional data
- what is registration
- what is a warping function
- difference between t and t*
- model
- talk about model parameters
- iterative algorithm to estimate parameters from model
- given a dataset in format outlined in our vignette `register_fpca` will register exponential family functional data of a specified family. This function returns on object of class `registration`.

- implemented for speed and efficiency (Key pieces of the code are written in C++ using Rcpp and RcppArmadillo libraries, cite)
- compatible with refund.shiny

$$
\begin{eqnarray*}
E\left[Y_i\left(h_i^{-1}(t_i^*)\right) | c_i, h_i^{-1} \right] &=& \mu_i(t) \\
g\left[\mu_i(t)\right]&=& \alpha(t) + \sum_{k = 1}^K c_{ik}\psi_k(t).
\end{eqnarray*}
$$

where $Y_i(t)$ is the functional response for subject $i$, the $X_{ik}(t)$ are functional predictors, the $\beta_k(t)$ are functional coefficients of interest, and the $\delta_i(t)$ are possibly correlated errors. This package implements two statistical methods (with and without variable selection) for estimating the parameters in the functional linear concurrent model; 



these methods are described in detail [@wrobel2018].

Given tidy datasets containing functional responses and predictors for all subjects, `vb_concurrent` and `vbvs_concurrent` fit the functional linear concurrent model using variational Bayes and variational Bayes with variable selection, respectively. These functions produce objects of class `flcm` and have the same structure. Coefficients and predictions can be extracted or computed from `flmc` objects using `coef` and `predict`. Interactive visualizations of `flmc` objects are supported through the `refund.shiny` package [@refund.shiny, @wrobel2016].


# References
