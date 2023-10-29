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

Rcpp::List          C_draw_event_times(Rcpp::IntegerVector id,
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

    //Rcpp::IntegerVector unique_id = Rcpp::unique(id);
    // slower but preserves order of ids.
    Rcpp::IntegerVector unique_id = id[duplicated(id) == 0];

    // these get returned in ... maybe a list
    Rcpp::IntegerVector new_id;
    Rcpp::NumericVector new_start_time ;
    Rcpp::NumericVector new_end_time ;
    Rcpp::IntegerVector new_status;

    Rcpp::IntegerVector all_subjects_new_X_rows;

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
        double current_event_time {origin};

        int subject_id {unique_id[i]};
        Rcpp::NumericVector subject_start_time;
        Rcpp::NumericVector subject_end_time;
        Rcpp::IntegerVector subject_status;

        // these get the existing data, plus any events
        Rcpp::NumericVector subject_new_start_time;
        Rcpp::NumericVector subject_new_end_time;
        Rcpp::IntegerVector subject_new_status;
        Rcpp::IntegerVector subject_j;

        // populate the temp variables
        // iterate over data to generate hazards and times for subject i
        for (int j = 0; j < n_id; j++) {
            // indicates if a row is associated with the current subject.
            temp_in = (id[j] == subject_id);
          // Rcout << "The value of temp_in is:" << temp_in << "\n";

            temp_n_rows += temp_in;

            // assign data from subject i to the temporary variables,
            // generate the baseline hazard
            if (temp_in) {

              subject_start_time.push_back(start_time[j]);
              subject_end_time.push_back(end_time[j]);
              subject_status.push_back(status[j]);

              subject_j.push_back(j);

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
                }
            }
        }

        // populate the subject_X matrix
        Rcpp::NumericMatrix subject_X (subject_j.length(), p);

        for (int k = 0; k < subject_j.length(); k++){
            // for (int pi = 0; pi<p; pi++){
            //     subject_X(k, pi) = X(subject_j[k], pi);
            // }
            subject_X(k, _) = X(subject_j[k], _);
        }

        // we now have
        // a vector of hazards for subject i
        // a vector of change points for subject i

        // next we generate one or more event times, depending on user specification
        // generate one, and then more if recurrent.

        current_event_time_NV = svycoxme::C_rpexp(1, hazards, changes, changes[0]);
        current_event_time = current_event_time_NV[0];

        subject_event_times.push_back(current_event_time);

        current_event_time_after_end_of_follow_up = current_event_time >= end_of_follow_up;

        if(!single){
            int reps {0};
            while (!current_event_time_after_end_of_follow_up) {

                current_event_time_NV = svycoxme::C_rpexp(1, hazards, changes, current_event_time);
                current_event_time = current_event_time_NV[0];

                current_event_time_after_end_of_follow_up = current_event_time >= end_of_follow_up;

                if(!current_event_time_after_end_of_follow_up){
                    subject_event_times.push_back(current_event_time);
                }

                if(reps > 10) {
                  break;
                }

            }

        }

        // for each event time, need to loop subject's data and interleaf
        // the events, filling in id, X (just get the row index, build later),
        // and status too.

        bool after; // event is after this interval ends.
        bool in_interval ; // event occurs in this interval
        bool ends_at_same_time;
        int  comp_start = 0;

        for (int j = 0; j < subject_event_times.length(); j++){

            for (int k = 0; k < subject_start_time.length(); k++){

                if(k < comp_start) continue;

                after = subject_event_times[j] > subject_end_time[k];
                in_interval = ((subject_event_times[j] > subject_start_time[k]) & (subject_event_times[j] <= subject_end_time[k]));
                ends_at_same_time = subject_event_times[j] == subject_end_time[k];

                if (after){
                    new_id.push_back(subject_id);
                    new_start_time.push_back(subject_start_time[k]);
                    new_end_time.push_back(subject_end_time[k]);
                    new_status.push_back(subject_status[k]);
                    all_subjects_new_X_rows.push_back(k);

                } else if (in_interval) {
                    if (ends_at_same_time){
                        new_id.push_back(subject_id);
                        new_start_time.push_back(subject_start_time[k]);
                        new_end_time.push_back(subject_end_time[k]);
                        new_status.push_back(1);
                        all_subjects_new_X_rows.push_back(k);
                    } else {
                        new_id.push_back(subject_id);
                        new_start_time.push_back(subject_start_time[k]);
                        new_end_time.push_back(subject_event_times[j]);
                        new_status.push_back(1);
                        all_subjects_new_X_rows.push_back(k);
                        // edit start time. comparisons will start here with the next event time.
                        subject_start_time[k] = subject_event_times[k];
                        comp_start = k;
                    }
                    // hoping this pops the for loop.
                    break;
                }
            }
        }
    }

    int n_all_subjects_new_X_rows = all_subjects_new_X_rows.length();

    Rcpp::NumericMatrix new_X (n_all_subjects_new_X_rows, p);

    // build the new X
    for (int k = 0; k < n_all_subjects_new_X_rows; k++){
      new_X(k, _)  = X(all_subjects_new_X_rows[k], _);
    }

    Rcpp::List L = List::create(new_id,
                                new_start_time,
                                new_end_time,
                                new_status,
                                new_X);

    // remember to update the method return bit
    return L;

}


