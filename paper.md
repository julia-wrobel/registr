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

Functional data analysis is a set of tools for understanding patterns and variability in data where the basic unit of observation is a curve measured over time, space, or another domain. Classic functional data analysis assumes that each curve is continuous or comes from a Gaussian distribution. However, applications with exponential family functional data -- curves that arise from any exponential family distribution, and have a smooth latent mean -- are increasingly common.  

Often in a functional dataset curves have similar underlying patterns but the main features of each curve, such as the minimum and maximum, have shifts such that the data appear misaligned. This misalignment can obscure patterns shared across curves and produce messy summary statistics. Registration methods reduce variability in functional data and clarify underlying patterns by aligning curves. Our method estimates a map, called a **warping function**, which transforms the domain from so that curves are aligned. The model for registration can be written

$$
\begin{eqnarray*}
E\left[Y_i\left(h_i^{-1}(t_i^*)\right) | c_i, h_i^{-1} \right] &=& \mu_i(t) \\
g\left[\mu_i(t)\right]&=& \alpha(t) + \sum_{k = 1}^K c_{ik}\psi_k(t).
\end{eqnarray*}
$$

For subject $i$, warping function $h_i^{-1}$ maps the domain on which curves are misaligned, $t_i^*$, to aligned domain $t$ such that $h_i^{-1}(t_i^*) = t$. Then $Y_i\left(t_i^*\right)$ and $Y_i\left(h_i^{-1}(t_i^*)\right)$ are the unregistered and registered functional response curves, respectively. The $\mu_i(t)$ are subject-specific means related to the population-level mean $\alpha(t)$ and a linear combination of population-level basis functions $\psi(t)$ and subject-specific scores $c_i$ through a known link function $g$. 

The `registr` package estimates warping functions and other parameters in this model via a two-step iterative algorithm which is detailed in @wrobel2018. The main function is `register_fpca`, which registers functional data from a specified exponential family distribution. `register_fpca` reads in a long-format functional dataset and outputs an object of class `registration`.

To enhance computational efficency, key algorithm components are implemented in C++ using the R libraries `Rcpp` and `RcppArmadillo` [@rcpp, @rcppArma]. Interactive visualizations are enabled with the `refund.shiny` package [@refund.shiny, @wrobel2016]. 



# References
