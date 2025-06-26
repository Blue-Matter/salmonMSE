
<!-- README.md is generated from README.Rmd. Please edit that file -->

# salmonMSE

> Management Strategy Evaluation for Salmon Species

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/salmonMSE)](https://cran.r-project.org/package=salmonMSE)
[![R-CMD-check](https://github.com/Blue-Matter/salmonMSE/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Blue-Matter/salmonMSE/actions/workflows/R-CMD-check.yaml)
[![License: GPL (\>=
3)](https://img.shields.io/badge/license-GPL%20(%3E=%203)-blue.svg)](https://cran.r-project.org/web/licenses/GPL-3)

<!-- badges: end -->

salmonMSE is a decision-support tool for Pacific salmon focusing on
strategic trade-offs among harvest, hatchery and habitat management
levers.

Funding for the development of salmonMSE is provided through the
[Pacific Salmon Strategy
Initiative](https://www.dfo-mpo.gc.ca/campaign-campagne/pss-ssp/index-eng.html)
in collaboration with the Department of Fisheries and Oceans Canada.

## Installation

salmonMSE can be downloaded from [CRAN](https://CRAN.R-project.org)
with:

``` r
install.packages("salmonMSE")
```

And the development version is available on
[GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("Blue-Matter/salmonMSE")
```

$$
\textrm{Smolt} =
\begin{cases}
p N_1 &, N_1 \le \frac{N_1}{N_2}C\\
\frac{N_1}{N_2}C &, N_1 \gt \frac{N_1}{N_2}C
\end{cases}
$$
