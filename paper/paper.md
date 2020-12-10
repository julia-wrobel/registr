---
title: 'registr 2.0: Incomplete Curve Registration for Exponential Family Functional Data'
tags:
  - R
  - Statistical analysis
  - Functional data
  - Partially observed curves
authors:
  - name: Julia Wrobel
    orcid: 0000-0001-6783-1421
    affiliation: 1
  - name: Alexander Bauer
    orcid: 0000-0003-3495-5131
    affiliation: 2
affiliations:
 - name: Columbia University, USA
   index: 1
 - name: Department of Statistics, LMU Munich, Germany
   index: 2
date: 10 December 2020
bibliography: paper.bib
---

# Introduction

Functional data are observed in many different fields.
One classic example are longer-term panel studies where a sequence of measurements
is observed for each subject.
Functional data can however also be observed over other domains like space.
Compared to classical longitudinal analyses and by setting the focus to the one
observed _curve_ per subject, functional data analysis sets the focus
on the analysis of the observed shapes of the time-dependent processes.
E.g., one can analyze the speed of growth of children until adulthood
in the Berkeley child growth study:

![\label{fig:data}](figures/1_data.png)

Functional data have different modes of variation that underly them.
In the mentioned example, not only can growth spurts be more or less pronounced
regarding the happening growth (i.e., _amplitude variation_), but each spurt
can also be delayed for some months / years for individual subjects (i.e., _phase variation_).
Even if focus is mostly on amplitude variation, observed curves often have to be
preprocessed by a _registration method_ to eliminate phase variation.

Most registration methods can only handle continuous data or data with a Gaussian
structure. However, functional data are often non-Gaussian or even categorical.
E.g., each time point can comprise a binary indicator for physical (non-)activity
of patients (compare @wrobel2019).
Moreover, most registration approaches rely on completely observed curves where
the measurement comprises the underlying process from its very start to its very
end.

# Exponential Family-based Registration

The `registr` package is based on the methods outlined in @wrobel2019.
Registration is performed using a likelihood-based approach and estimates
_inverse warping functions_ $h_i^{-1}: t_i^* \mapsto t$ that map the observed
time domain $t_i^*$ of subject $i$ to the common time domain $t$.
The overall model results to

$$
\begin{aligned}
E\left[Y_i\left(h_i^{-1}(t_i^*)\right) | c_i, h_i^{-1} \right] &= \mu_i(t), \\
g\left[\mu_i(t)\right] &= \alpha(t) + \sum_{k = 1}^K c_{ik}\psi_k(t),
\end{aligned}
$$

with $Y_i\left(t_i^*\right)$ and $Y_i\left(h_i^{-1}(t_i^*)\right)$ the unregistered and registered curves, respectively,
and $\mu_i(t)$ the subject-specific means.
The latter are estimated with regard to a lower-dimensional representation based on
a population-level mean $\alpha(t)$ and a linear combination of population-level basis functions $\psi(t)$
and subject-specific scores $c_i$, based on a known link function $g$.
The lower-dimensional representation is estimated using a likelihood-based
approach for Generalized Functional Principal Component Analysis (GFPCA).

The overall model is estimated with function `register_fpca()`, which iterates 
between the estimation of warping
functions for the registration step (implemented in function `registr()`)
and GFPCA estimation (functions `fpca_gauss()` or `bfpca()` for Gaussian or binomial data, respectively).
The latter GFPCA functions are partly implemented in C++ using R libraries `Rcpp` and `RcppArmadillo` [@rcpp, @rcppArma]
to enhance computational efficiency.
The package also includes an implementation of the _two-step GFPCA approach_
of @gertheiss2017 to also handle further exponential family distributions.
The respective implementation is based on the existing `gfpca` package of @goldsmith2016.

The registration and GFPCA approaches in the `registr` package are implemented
for gaussian, binomial, gamma and poisson data.
When performing a registration the template function (i.e. the reference curve
to which all curves are mapped) can be flexibly defined by the user,
using the argument `Y_template` in `registr()` and `register_fpca()`.
The results can be interactively visualized with the ``refund.shiny` package [@refund.shiny, @wrobel2016]. 

# Incomplete Curve Registration

The `registr` package extends the approach of @wrobel2019 to also handle
incomplete curves where the underlying process of interest was not observed
from its very beginning (i.e., _leading incompleteness_), until its very end
(_trailing incompleteness_), or both (_full incompleteness_).

It often is a quite strict assumption in incomplete data
settings that all warping functions start and/or end on the diagonal, i.e. that the individual,
observed part of the whole time domain is not (to some extent) distorted.
Therefore, the `registr` package gives the additional option to estimate
warping functions without the constraint that their starting point and/or endpoint
lies on the diagonal.

On the other hand, if we fully remove the constraint for the starting points / endpoints
to lie on the diagonal, this can lead to very extreme and unrealistic distortions
of the time domain. This problem is further accompanied by the fact that
the assessment of some given warping to be realistic or unrealistic can heavily
vary between different applications.
As of this reason, our method includes a penalization parameter $\lambda$ that
has to be set manually to specify which kinds of distortions are deemed realistic
in the application at hand.

Mathematically speaking, we add a penalization term to the likelihood $\ell(i)$ 
for curve $i$. For a setting with **full incompleteness** (i.e., where both the starting
point and endpoint are free to vary from the diagonal) this results in

$$
\begin{aligned}
\ell_{pen}(i) &= \ell(i) - \lambda \cdot pen(i), \\
\text{with} \ \ \ 
pen(i) &= \left( \hat{h}^{-1}_i(t^*_{min,i}) - t^*_{min,i} \right)^2 +
\left( \hat{h}^{-1}_i(t^*_{max,i}) - t^*_{max,i} \right)^2,
\end{aligned}
$$

where $t^*_{min,i},t^*_{max,i}$ are the minimum / maximum of the observed time domain of curve $i$ and
$\hat{h}^{-1}_i(t^*_{min,i}), \hat{h}^{-1}_i(t^*_{max,i})$ the inverse warping function evaluated at this
minimum / maximum.

The higher the penalization parameter $\lambda$, the more the starting point and endpoint
of the warping function are forced towards an identity warp, i.e. the diagonal line.
Given a specific application, $\lambda$ should be chosen s.t.
unrealistic distortions of the time domain are prevented.
To do so, the user has to run the registration approach multiple times with
different $\lambda$'s to find an optimal value.

The results of our registration approach applied to the growth curve data
are shown in Figure \autoref{fig:registration}.

![\label{fig:registration}](figures/2_registration.png)

# Acknowledgements

We thank Fabian Scheipl for valuable methodological contributions.

# References