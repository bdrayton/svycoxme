
# svycoxme

<!-- badges: start -->
[![R-CMD-check](https://github.com/bdrayton/svycoxme/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/bdrayton/svycoxme/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of svycoxme is to fit mixed-effects proportional hazards models to data from 
complex samples. Most of the work is done by the `coxme` package. The svycoxme package 
provides wrappers to fit models using survey designs from the `survey` package, and provides 
a range of variance estimators and provides variances estimation by Taylor series 
linearisation or replicate weights. 


## Installation

You can install the development version of svycoxme from [GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("bdrayton/svycoxme")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(svycoxme)
## basic example code
```

