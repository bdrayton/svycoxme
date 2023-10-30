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
                                       double origin,
                                       int single,
                                       int maximum_events = 1000) {

    //Rcpp::IntegerVector unique_id = Rcpp::unique(id);
    // slower but preserves order of ids.
    Rcpp::IntegerVector unique_id = id[duplicated(id) == 0];

    // these get returned in a list
    Rcpp::IntegerVector new_id;
    Rcpp::NumericVector new_start_time ;
    Rcpp::NumericVector new_end_time ;
    Rcpp::IntegerVector new_status;

    // for generating a new_X matrix later (I need to calculate number of rows first)
    Rcpp::IntegerVector all_subjects_new_X_rows;

    // for iterating over the data
    int n_unique_id = unique_id.length();
    int n_id = id.length();
    int p = X.ncol();

    // for subject specific calcs
    bool temp_in;
    double end_of_follow_up;

    // subject specific hazards
    Rcpp::NumericVector baseline_hazard_end = clone(baseline_hazard_start);
    baseline_hazard_end.erase(0);
    baseline_hazard_end.push_back(R_PosInf);

    int n_baseline_hazard = baseline_hazard.length();
    bool hazard_in_interval;
    double temp_hazard;
    double temp_start_time;
    double temp_baseline_hazard_start;

    // loop over each unique subject
    for (int i = 0; i < n_unique_id; i++) {
        // declare more working variables here.
        // they should die at the end of each iteration.

        // subject's data
        int subject_id {unique_id[i]};
        Rcpp::NumericVector subject_start_time;
        Rcpp::NumericVector subject_end_time;
        Rcpp::IntegerVector subject_status;

        // outputs
        Rcpp::NumericVector subject_event_times;

        // for time to event generation
        Rcpp::NumericVector hazards;
        Rcpp::NumericVector changes;
        bool current_event_time_before_end_of_follow_up;
        Rcpp::NumericVector current_event_time_NV;
        double current_event_time {origin};

        // populate the temp variables
        // iterate over data to generate hazards and times for subject i
        for (int j = 0; j < n_id; j++) {
            // indicates if a row is associated with the current subject.
            temp_in = (id[j] == subject_id);

            // assign data from subject i to the temporary variables,
            // generate the subject's hazard
            if (temp_in) {

              subject_start_time.push_back(start_time[j]);
              subject_end_time.push_back(end_time[j]);
              subject_status.push_back(status[j]);

                for (int k = 0; k < n_baseline_hazard; k++) {

                    // identify if hazard is in same time interval as the risk score.
                    hazard_in_interval = ((baseline_hazard_end[k] > start_time[j]) & (baseline_hazard_start[k] < end_time[j]));

                    if (hazard_in_interval) {
                        temp_hazard = baseline_hazard[k] * exp(risk_score[j]);
                        hazards.push_back(temp_hazard);

                        temp_start_time = start_time[j];
                        temp_baseline_hazard_start = baseline_hazard_start[k];

                        changes.push_back(std::max(temp_start_time, temp_baseline_hazard_start));

                    }
                }
            }
        }

        end_of_follow_up = Rcpp::max(subject_end_time);

        // we now have
        // a vector of hazards for subject i
        // a vector of change points for subject i

        // next we generate one or more event times, depending on user specification
        // generate one, and then more if recurrent.

        current_event_time_NV = svycoxme::C_rpexp(1, hazards, changes, changes[0]);
        current_event_time = current_event_time_NV[0];

        subject_event_times.push_back(current_event_time);

        current_event_time_before_end_of_follow_up = current_event_time < end_of_follow_up;

        if(single == 0){
            int reps {0};
            while (current_event_time_before_end_of_follow_up) {

                current_event_time_NV = svycoxme::C_rpexp(1, hazards, changes, current_event_time);
                current_event_time = current_event_time_NV[0];

                current_event_time_before_end_of_follow_up = current_event_time < end_of_follow_up;

                subject_event_times.push_back(current_event_time);

                reps++;

                if(reps > maximum_events) {
                  Rcout << "Events went past maximum_events for subject_id : " << subject_id << "\n";
                  break;
                }

            }

        }

        // for each event time, need to loop subject's data and interleave
        // the events, filling in id, X (just get the row index, build later),
        // and status too.

        bool after; // event is after this interval ends.
        bool in_interval ; // event occurs in this interval
        bool ends_at_same_time; // idk if this will ever actually happen, but just in case.
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
                        subject_start_time[k] = subject_event_times[j];
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


