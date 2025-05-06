
# svycoxme

<!-- badges: start -->
[![R-CMD-check](https://github.com/bdrayton/svycoxme/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bdrayton/svycoxme/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of svycoxme is to fit mixed-effects proportional hazards models to data from 
complex samples. Most of the work is done by the `coxme` package. The svycoxme package 
provides wrappers to fit models using survey designs from the `survey` package, and 
provides variances estimation by Taylor series linearisation or replicate weights. 


## Installation

You can install svycoxme from CRAN with:

``` r
install.packages("svycoxme")
```

You can install the development version of svycoxme from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bdrayton/svycoxme")
```

## Example

This is a basic example using the [samp_srcs](/man/samp_srcs.Rd) dataset provided with the package.

``` r
library(survey)
library(svycoxme)
des <- svydesign(ids = ~group_id, weights = ~weight, data = samp_srcs)

fit1 <- svycoxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | group_id), design = des)

summary(fit1)

```

