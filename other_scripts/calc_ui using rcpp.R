

# see http://adv-r.had.co.nz/Rcpp.html#rcpp-package



# convert ui_calc. to a rcpp function.

vignette(package = "Rcpp")

# calc_ui.coxme_parts
# these are the parts of parts used in calc_ui.coxme_parts.
# X
# Z
# stat
# time_start
# time_stop
# weightts
# exp_risk_score
# S0
# S1_X
# S1_Z
# ui_penalty

#include <Rcpp.h>

// Function declaration with export tag
// [[Rcpp::export]]

Rcpp::NumericMatrix
calc_ui_coxme(Rcpp::NumericMatrix X,
              Rcpp::NumericMatrix Z,
              Rcpp::NumericMatrix stat,
              Rcpp::NumericMatrix time_start,
              Rcpp::NumericMatrix time_stop,
              Rcpp::NumericMatrix weights,
              Rcpp::NumericMatrix exp_risk_score,
              Rcpp::NumericMatrix S0,
              Rcpp::NumericMatrix S1_X,
              Rcpp::NumericMatrix S1_Z,
              Rcpp::NumericMatrix ui_penalty) {

//Preallocate storage

  # need to get dimensions from inputs

  Rcpp::NumericMatrix lin_score();





}






