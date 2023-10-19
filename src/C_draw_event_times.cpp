#include <Rcpp.h>
using namespace Rcpp;

//' C_draw_event_times
//' 
//' Draw event times times using possibly time-varying covariate data
//'
//' @param subject_id unique identifier for subject
//' 
//' @export
// [[Rcpp::export]]

Rcpp::IntegerVector C_draw_event_times(Rcpp::IntegerVector id,
                                       Rcpp::NumericVector start_time,
                                       Rcpp::NumericVector end_time,
                                       Rcpp::IntegerVector status,
                                       Rcpp::NumericMatrix X,
                                       Rcpp::NumericVector risk_score,
                                       Rcpp::NumericVector t,
                                       double end_of_follow_up,
                                       double origin) {

    Rcpp::IntegerVector unique_id = Rcpp::unique(id);

    int n_unique_id = unique_id.length();
    int n_id = id.length();

    // working variables.

    bool temp_in;
    Rcpp::NumericVector temp_start_time;
    Rcpp::NumericVector temp_end_time;
    Rcpp::IntegreVector temp_status;
    Rcpp::NumericMatrix temp_X;
    Rcpp::NumericVector temp_risk_score;
    Rcpp::NumericVector temp_start_time;

    // loop over each unique subject
    for (int i = 0; i < n_unique_id; i++) {

        // populate the temp variables
        for (int j = 0; j < n_id; j++) {
            temp_in = id[j] == unique_id[i];

            if (temp_in) {
                temp_start_time.push_back(start_time[j]);
                temp_end_time.push_back(end_time[j]);
                temp_status.push_back(status[j]);
   //             temp_X ;  will need to figure out this bit.
                temp_risk_score.push_back(risk_score[j]);
            }

        }

     // clear temp vars. 

        
    }

    // remember to update the method return bit
    return unique_subject_id;

}


