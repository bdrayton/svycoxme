#include <Rcpp.h>
using namespace Rcpp;

//' C_calc_ui
//' 
//' Calculate score residuals for a mixed effect cox model. Only calculates for fixed effects.
//'
//' @param time_start start time for an observation
//' @param time_stop stop time for an observation
//' @param stat indicator, 1 for event, 0 for censoring
//' @param weights observation-level sampling weights
//' @param exp_risk_score exp(betaX + bZ) 
//' @param S0 sum of exp_risk_score for observations at risk at time_stop
//' @param S1_X sum of X * exp_risk_score for  observations at risk at time_stop
//' @param X matrix of covariate data
//' @param weighted logical scalar indicating if the residuals should be weighted
//' @export
// [[Rcpp::export]]

Rcpp::NumericMatrix C_calc_ui(Rcpp::NumericMatrix& time_start,
                              Rcpp::NumericMatrix& time_stop,
                              Rcpp::NumericMatrix& stat,
                              Rcpp::NumericMatrix& weights,
                              Rcpp::NumericMatrix& exp_risk_score,
                              Rcpp::NumericMatrix& S0,
                              Rcpp::NumericMatrix& S1_X,
                              Rcpp::NumericMatrix& X,
                              bool weighted) {

    int nrow = X.nrow();
    int X_ncol = X.ncol();

    bool Yi_at_tj;
    double temp;
    Rcpp::NumericMatrix S0_inv(S0.nrow(), S0.ncol());  

    Rcpp::NumericMatrix resid(nrow, X_ncol);

    // calculate the second term
    for (int i = 0; i < nrow; i++) {
        S0_inv[i] = 1.0 / S0[i];
        for (int j = 0; j < nrow; j++) {
            Yi_at_tj = ((time_start[i] < time_stop[j]) & (time_stop[i] >= time_stop[j]));
            temp = stat[j] * weights[j] * Yi_at_tj * exp_risk_score[i] * S0_inv[j];
            for (int k = 0; k < X_ncol; k++) {
                resid(i, k) += temp * (X(i, k) - S1_X(j, k) * S0_inv[j]);
            }
        }
    }

    // calculate the first term, minus the second term, and weight if needed.
    for (int i = 0; i < nrow; i++) {
        for (int k = 0; k < X_ncol; k++) {
            resid(i, k) = stat[i] * (X(i, k) - S1_X(i, k) * S0_inv[i]) - resid(i, k);
            resid(i, k) = resid(i, k) * (weighted * weights[i] + !weighted);
        }
    }

    return resid;

}
