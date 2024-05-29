#include <Rcpp.h>
using namespace Rcpp;

//' C_calc_S0_S1X
//' 
//' Calculate S0 and S1X for a mixed effect cox model. Only calculates for fixed effects.
//'
//' @param time_start start time for an observation
//' @param time_stop stop time for an observation
//' @param X matrix of covariate data
//' @export
// [[Rcpp::export]]

Rcpp::List  C_calc_S0_S1X(  Rcpp::NumericMatrix& time_start,
                            Rcpp::NumericMatrix& time_stop,
                            Rcpp::NumericMatrix& exp_risk_score,
                            Rcpp::NumericMatrix& X) {

    int nrow = X.nrow();
    int X_ncol = X.ncol();

    bool in_risk_set;
    double temp;

    Rcpp::NumericMatrix S0(nrow, 1);
    Rcpp::NumericMatrix S1_X(nrow, X_ncol);
    
    for(int i = 0; i < nrow; i++){
        for(int j = 0; j < nrow; j++){
             // start and stop tests.
            in_risk_set = (time_stop[i] > time_start[j]) & (time_stop[i] <= time_stop[j]); 
            temp = in_risk_set * exp_risk_score[j];
            S0[i] = S0[i] + temp;

            for(int k = 0; k < X_ncol; k++){
                S1_X(i, k) = S1_X(i, k) + temp * X(j, k); 
            }
        }
    }

    Rcpp::List L = List::create(S0,
                                S1_X);
                            
    // remember to update the method return bit
    return L; 

}
