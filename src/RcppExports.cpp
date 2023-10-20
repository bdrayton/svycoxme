// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// C_calc_ui
Rcpp::NumericMatrix C_calc_ui(Rcpp::NumericMatrix& time_start, Rcpp::NumericMatrix& time_stop, Rcpp::NumericMatrix& stat, Rcpp::NumericMatrix& weights, Rcpp::NumericMatrix& exp_risk_score, Rcpp::NumericMatrix& S0, Rcpp::NumericMatrix& S1_X, Rcpp::NumericMatrix& X, bool weighted);
RcppExport SEXP _svycoxme_C_calc_ui(SEXP time_startSEXP, SEXP time_stopSEXP, SEXP statSEXP, SEXP weightsSEXP, SEXP exp_risk_scoreSEXP, SEXP S0SEXP, SEXP S1_XSEXP, SEXP XSEXP, SEXP weightedSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type time_start(time_startSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type time_stop(time_stopSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type stat(statSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type weights(weightsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type exp_risk_score(exp_risk_scoreSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type S0(S0SEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type S1_X(S1_XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix& >::type X(XSEXP);
    Rcpp::traits::input_parameter< bool >::type weighted(weightedSEXP);
    rcpp_result_gen = Rcpp::wrap(C_calc_ui(time_start, time_stop, stat, weights, exp_risk_score, S0, S1_X, X, weighted));
    return rcpp_result_gen;
END_RCPP
}
// C_draw_event_times
Rcpp::NumericVector C_draw_event_times(Rcpp::IntegerVector id, Rcpp::NumericVector start_time, Rcpp::NumericVector end_time, Rcpp::IntegerVector status, Rcpp::NumericMatrix X, Rcpp::NumericVector risk_score, Rcpp::NumericVector baseline_hazard, Rcpp::NumericVector baseline_hazard_start, Rcpp::NumericVector baseline_hazard_end, double end_of_follow_up, double origin);
RcppExport SEXP _svycoxme_C_draw_event_times(SEXP idSEXP, SEXP start_timeSEXP, SEXP end_timeSEXP, SEXP statusSEXP, SEXP XSEXP, SEXP risk_scoreSEXP, SEXP baseline_hazardSEXP, SEXP baseline_hazard_startSEXP, SEXP baseline_hazard_endSEXP, SEXP end_of_follow_upSEXP, SEXP originSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type id(idSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type start_time(start_timeSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type end_time(end_timeSEXP);
    Rcpp::traits::input_parameter< Rcpp::IntegerVector >::type status(statusSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type X(XSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type risk_score(risk_scoreSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type baseline_hazard(baseline_hazardSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type baseline_hazard_start(baseline_hazard_startSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type baseline_hazard_end(baseline_hazard_endSEXP);
    Rcpp::traits::input_parameter< double >::type end_of_follow_up(end_of_follow_upSEXP);
    Rcpp::traits::input_parameter< double >::type origin(originSEXP);
    rcpp_result_gen = Rcpp::wrap(C_draw_event_times(id, start_time, end_time, status, X, risk_score, baseline_hazard, baseline_hazard_start, baseline_hazard_end, end_of_follow_up, origin));
    return rcpp_result_gen;
END_RCPP
}
// C_rpexp
Rcpp::NumericVector C_rpexp(int n, Rcpp::NumericVector rate, Rcpp::NumericVector t, double start);
RcppExport SEXP _svycoxme_C_rpexp(SEXP nSEXP, SEXP rateSEXP, SEXP tSEXP, SEXP startSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type n(nSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type rate(rateSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type t(tSEXP);
    Rcpp::traits::input_parameter< double >::type start(startSEXP);
    rcpp_result_gen = Rcpp::wrap(C_rpexp(n, rate, t, start));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _svycoxme_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_svycoxme_C_calc_ui", (DL_FUNC) &_svycoxme_C_calc_ui, 9},
    {"_svycoxme_C_draw_event_times", (DL_FUNC) &_svycoxme_C_draw_event_times, 11},
    {"_svycoxme_C_rpexp", (DL_FUNC) &_svycoxme_C_rpexp, 4},
    {"_svycoxme_rcpp_hello_world", (DL_FUNC) &_svycoxme_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_svycoxme(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
