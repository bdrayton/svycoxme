#include <Rcpp.h>
#include <svycoxme_RcppExports.h>
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
                                       double origin,
                                       int single) {

    Rcpp::IntegerVector unique_id = Rcpp::unique(id);

    // these get returned in ... maybe a list

    Rcpp::NumericVector event_times ;
    Rcpp::NumericMatrix new_X ;
    Rcpp::IntegerVector new_status;

    //

    int n_unique_id = unique_id.length();
    int n_id = id.length();

    int p = X.ncol();

    int n_baseline_hazard = baseline_hazard.length();

    bool hazard_in_interval;
    bool temp_in;
    double temp_hazard;

    double temp_start_time;
    double temp_baseline_hazard_start;

    // loop over each unique subject
    for (int i = 0; i < n_unique_id; i++) {
    // for (int i = 0; i < 1; i++) {
        // declare the working variables here.
        // they should die at the end of each iteration.

        Rcpp::NumericVector hazards;
        Rcpp::NumericVector changes;
        Rcpp::NumericVector ends;
        int temp_n_rows = 0;
        Rcpp::NumericVector subject_event_times;
        bool current_event_time_after_end_of_follow_up;
        Rcpp::NumericVector current_event_time_NV;
        double current_event_time;
        int event_count {0};

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

                for (int k = 0; k < n_baseline_hazard; k++) {

                    // identify if hazard is in same time interval as the risk score.
                    hazard_in_interval = ((baseline_hazard_end[k] > start_time[j]) & (baseline_hazard_start[k] < end_time[j]));

                    if (hazard_in_interval) {
                        temp_hazard = baseline_hazard[k] * exp(risk_score[j]);
                        hazards.push_back(temp_hazard);

                        temp_start_time = start_time[j];
                        temp_baseline_hazard_start = baseline_hazard_start[k];

                        changes.push_back(std::max(temp_start_time, temp_baseline_hazard_start));

                      // ends.push_back(std::min(temp_end_time[j]  , baseline_hazard_end[k]));
                    }


                    if(single){
                        current_event_time_NV = svycoxme::C_rpexp(1, hazards, changes, changes[0]);
                        current_event_time = current_event_time_NV[0];

                        subject_event_times.push_back(current_event_time);

                        event_count++;

                    } else {

                      current_event_time = origin;

                      while (!current_event_time_after_end_of_follow_up) {

                        current_event_time_NV = svycoxme::C_rpexp(1, hazards, changes, current_event_time);
                        current_event_time = current_event_time_NV[0];

                        subject_event_times.push_back(current_event_time);

                        current_event_time_after_end_of_follow_up = current_event_time >= end_of_follow_up;

                        event_count++;

                      }

                      // drop the last event time.
                      subject_event_times.erase(event_count-1);

                    }


                }
            }
        }
    }

    // remember to update the method return bit
    return event_times;

}


