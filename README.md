
<!-- README.md is generated from README.Rmd. Please edit that file -->

# registr <img src="README_files/figures/registr.png" align="right" height = "150" />

<!-- badges: start -->

[![](https://www.r-pkg.org/badges/version/registr)](https://cran.r-project.org/package=registr)
[![](http://cranlogs.r-pkg.org/badges/grand-total/registr?color=green)](https://cran.r-project.org/package=registr)
[![Codecov test
coverage](https://codecov.io/gh/julia-wrobel/registr/branch/master/graph/badge.svg)](https://codecov.io/gh/julia-wrobel/registr/coverage.svg?branch=master)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02964/status.svg)](https://doi.org/10.21105/joss.02964)

<!-- badges: end -->

Registration for incomplete exponential family functional data.

- Authors: [Julia Wrobel](https://julia-wrobel.github.io/), Alexander
  Bauer, [Erin McDonnell](http://eimcdonnell.com/), and [Jeff
  Goldsmith](https://jeffgoldsmith.com/)
- License: [MIT](https://opensource.org/licenses/MIT). See the
  [LICENSE](LICENSE) file for details
- Version: 2.1

### What it does

------------------------------------------------------------------------

Functional data analysis is a set of tools for understanding patterns
and variability in data where the basic unit of observation is a curve
measured over some domain such as time or space. An example is an
accelerometer study where intensity of physical activity was measured at
each minute over 24 hours for 50 subjects. The data will contain 50
curves, where each curve is the 24-hour activity profile for a
particular subject.

Classic functional data analysis assumes that each curve is continuous
or comes from a Gaussian distribution. However, applications with
exponential family functional data – curves that arise from any
exponential family distribution, and have a smooth latent mean – are
increasingly common. For example, take the accelerometer data just
mentioned, but assume researchers are interested in *sedentary behavior*
instead of *activity intensity*. At each minute over 24 hours they
collect a binary measurement that indicates whether a subject was active
or inactive (sedentary). Now we have a *binary curve* for each subject –
a trajectory where each time point can take on a value of 0 or 1. We
assume the binary curve has a smooth latent mean, which in this case is
interpreted as the probability of being active at each minute over 24
hours. This is a example of exponential family functional data.

Often in a functional dataset curves have similar underlying patterns
but the main features of each curve, such as the minimum and maximum,
have shifts such that the data appear misaligned. This misalignment can
obscure patterns shared across curves and produce messy summary
statistics. Registration methods reduce variability in functional data
and clarify underlying patterns by aligning curves.

This package implements statistical methods for registering exponential
family functional data. The basic methods are described in more detail
in our
[paper](https://julia-wrobel.github.io/Downloads/registration_ef.pdf)
and were further adapted to (potentially) incomplete curve settings
where (some) curves are not observed from the very beginning and/or
until the very end of the common domain. For details on the incomplete
curve methodology and how to use it see the corresponding package
vignette. Instructions for installing the software and using it to
register simulated binary data are provided below.

### Installation

------------------------------------------------------------------------

To install from `CRAN`, please use:

``` r
install.packages("registr")
```

To install the latest version directly from Github, please use:

``` r
install.packages("devtools")
devtools::install_github("julia-wrobel/registr")
```

The `registr` package includes vignettes with more details on package
use and functionality. To install the latest version and pull up the
vignettes please use:

``` r
devtools::install_github("julia-wrobel/registr", build_vignettes = TRUE)
vignette(package = "registr")
```

### How to use it

------------------------------------------------------------------------

This example registers simulated binary data. More details on the use of
the package can be found in the vignettes mentioned above.

The code below uses `registr::simulate_unregistered_curves()` to
simulate curves for 100 subjects with 200 timepoints each, observed over
domain $(0, 1)$. All curves have similar structure but the location of
the peak is shifted. On the observed domain $t^*$ the curves are
unregistered (misaligned). On the domain $t$ the curves are registered
(aligned).

``` r
library(registr)

registration_data = simulate_unregistered_curves(I = 100, D = 200, seed = 2018)
```

The plot below shows the unregistered curves and registered curves.

<img src="README_files/figure-gfm/plot_sim_data-1.png" alt="" style="display: block; margin: auto;" />

Continuously observed curves are shown above in order to illustrate the
misalignment problem and our simulated data; the simulated dataset also
includes binary values which have been generated by using these
continuous curves as probabilities. The unregistered and registered
binary curves for two subjects are shown below.

<img src="README_files/figure-gfm/plot_2subjs-1.png" alt="" style="display: block; margin: auto;" />

Our software registers curves by estimating $t$. For this we use the
function `registration_fpca()`.

``` r
binary_registration = register_fpca(Y = registration_data, family = "binomial", 
                                    Kt = 6, Kh = 4, npc  = 1)
## Running initial registration step
## current iteration: 1
## Running final FPCA step
```

The plot below shows unregistered, true registered, and estimated
registered binary curves for two subjects after fitting our method.

<img src="README_files/figure-gfm/plot_fit-1.png" alt="" style="display: block; margin: auto;" />

### Citation

If you like our software, please cite it in your work! To cite the
latest `CRAN` version of the package with `BibTeX`, use

    @Manual{,
        title = {registr: Registration for Exponential Family Functional Data},
        author = {Julia Wrobel and Alexander Bauer and Erin McDonnell and Jeff Goldsmith},
        year = {2022},
        note = {R package version 2.1.0},
        url = {https://CRAN.R-project.org/package=registr},
      }

To cite the 2021 Journal of Open Source Software paper, use

    @article{wrobel2021registr,
      title={registr 2.0: Incomplete Curve Registration for Exponential Family Functional Data},
      author={Wrobel, Julia and Bauer, Alexander},
      journal={Journal of Open Source Software},
      volume={6},
      number={61},
      pages={2964},
      year={2021}
    }

To cite the 2018 Journal of Open Source Software paper, use

    @article{wrobel2018regis,
      title={registr: Registration for Exponential Family Functional Data},
      author={Wrobel, Julia},
      journal={The Journal of Open Source Software},
      volume={3},
      year={2018}
    }

### Contributions

------------------------------------------------------------------------

If you find small bugs, larger issues, or have suggestions, please file
them using the [issue
tracker](https://github.com/julia-wrobel/registr/issues) or email the
maintainer at <julia.wrobel@cuanschutz.edu>. Contributions (via pull
requests or otherwise) are welcome.
