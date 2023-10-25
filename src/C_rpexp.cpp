#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::interfaces(r, cpp)]]

//' C_rpexp
//'
//' Draw a single random number from a piecewise exponential distribution.
//' See msm::rpexp()
//'
//' @param rate vector of event rates
//' @param t vector same length as rate, giving the times the rate changes. The values of t should be ascending order.
//' @param start numeric scalar; delayed entry time. The random deviates will be left truncated.
//'
//' @export
// [[Rcpp::export]]
Rcpp::NumericVector C_rpexp(int n,
               Rcpp::NumericVector rate,
               Rcpp::NumericVector t,
               double start) {

  Rcpp::NumericVector random_draw (n);

  int n_t = t.length();
  int n_rate = rate.length();

  int n_deleted_rate = 0;
  double the_rate;
  double last_rate = rate[n_rate-1];

  Rcpp::NumericVector original_rate = rate;

  double diff_t;
  Rcpp::NumericVector H = {0};

  double H_temp;

  Rcpp::NumericVector e (n);

  if ( (n_rate == 1) | (start > t[n_t-1]) ) {

    for(int i = 0; i < n; i++){
      random_draw[i] = R::rexp(1/last_rate);
    }

    return(random_draw);

  }

  if (start > t[0]) {
    // remove rates and times before the start time.
    for (int i = (n_t-1); i > -1; i--){
      if(start > t[i]){
          t.erase(i);
          rate.erase(i);
          n_deleted_rate++;
      }
    }
    // append start to t, and add the correct rate
    the_rate = original_rate[n_deleted_rate-1];

    t.push_front(start);
    rate.push_front(the_rate);
  }

  for (int i = 0; i < (t.length()-1); i++){
    diff_t = t[i+1] - t[i];

    H_temp = H[i] + rate[i] * diff_t;

    H.push_back(H_temp);

  }

  e = Rcpp::rexp(n);

  for (int i = 0; i < n; i++) {
    // compare e to each value of H.
    for (int j = H.length()-1; j > -1; j--) {
      if(e[i]>=H[j]){
        random_draw[i] = t[j] + (e[i] - H[j])/rate[j];
        break;
      }
    }
  }

  return(random_draw);

}

