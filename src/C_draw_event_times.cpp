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

Rcpp::NumericVector C_draw_event_times(Rcpp::IntegerVector id,
                                       Rcpp::NumericVector start_time,
                                       Rcpp::NumericVector end_time,
                                       Rcpp::IntegerVector status,
                                       Rcpp::NumericMatrix X,
                                       Rcpp::NumericVector risk_score,
                                       Rcpp::NumericVector baseline_hazard,
                                       Rcpp::NumericVector baseline_hazard_start,
                                       Rcpp::NumericVector baseline_hazard_end,
                                       double end_of_follow_up,
                                       double origin) {

    Rcpp::IntegerVector unique_id = Rcpp::unique(id);

    int n_unique_id = unique_id.length();
    int n_id = id.length();

    int p = X.ncol();

    int n_baseline_hazard = baseline_hazard.length();

    Rcpp::NumericVector all_hazards;

    bool hazard_in_interval;
    bool temp_in;
    double temp_hazard;

    // loop over each unique subject
    for (int i = 0; i < n_unique_id; i++) {
    // for (int i = 0; i < 1; i++) {
        // declare the working variables here.
        // they should die at the end of each iteration.

        Rcpp::NumericVector hazards;
        Rcpp::NumericVector changes;
        Rcpp::NumericVector ends;
        int temp_n_rows = 0;

//         bool temp_in;
//         Rcpp::NumericVector temp_start_time;
//         Rcpp::NumericVector temp_end_time;
//         Rcpp::IntegerVector temp_status;
//         Rcpp::NumericMatrix temp_X;
//         Rcpp::NumericVector temp_risk_score;

        // populate the temp variables
        // iterate over data to generate hazards and times for subject i
        for (int j = 0; j < n_id; j++) {
            // indicates if a row is associated with the current subject.
            temp_in = (id[j] == unique_id[i]);

            temp_n_rows += temp_in;

            // assign data from subject i to the temporary variables.
            if (temp_in) {
                // temp_start_time.push_back(start_time[j]);
                // temp_end_time.push_back(end_time[j]);
                // temp_status.push_back(status[j]);
                // temp_risk_score.push_back(risk_score[j]);
                // for(int k = 0; k < p; k++){
                //   temp_X(temp_n_rows - 1, k) = X(j,k);
                // }

                for (int k = 0; k < n_baseline_hazard; k++) {

                    // identify if hazard is in same time interval as the risk score.
                    hazard_in_interval = ((baseline_hazard_end[k] > start_time[j]) & (baseline_hazard_start[k] < end_time[j]));

                    if (hazard_in_interval) {
                      temp_hazard = baseline_hazard[k] * exp(risk_score[j]);
                      hazards.push_back(temp_hazard);
                      all_hazards.push_back(temp_hazard);


                      // changes.push_back(std::max(1, 2));
                      // ends.push_back(std::min(3,4));

                      // changes.push_back(std::max(temp_start_time[j], baseline_hazard_start[k]));
                      // ends.push_back(std::min(temp_end_time[j]  , baseline_hazard_end[k]));
                    }
                }
            }
        }
    }

    // remember to update the method return bit
    return all_hazards;

}


