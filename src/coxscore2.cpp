#include <Rcpp.h>
using namespace Rcpp;



// /*
 // ** This is a port of the survival script coxscore2.c
 // ** The substantive changes are in variable initiation,
 // ** and changing indexing of matrices. There is no need to protect things
 // ** and reassign passed in objects.
 // ** SEXP matrix has variables in the
 // ** rows and observations in the columns, ordering inherited from Fortran. THis
 // ** has been changed to the R ordering, using RCPP matricies.
 // **
 // ** Compute the score residuals for a Cox model
 // **
 // ** Input
 // **      y       matrix of time and status values
 // **      strata  non-negative integer, unique value for each strata
 // **      covar2  the matrix of covariates
 // **      score   the vector of subject scores, i.e., exp(beta*X + b*M)
 // **      weights case weight
 // **      method  ==1 for efron method
 // **
 // ** Output
 // **      resid   a matrix of the same shape as x
 // **
 // ** Scratch
 // **      scratch,  from which a and a2 are carved
 // **
 // ** Data must be sorted by strata, ascending time within strata.
 // **
 // ** Updated 4/2023 to be O(np) instead of O(n^2 p):
 // **   the score is sum (x_i -xbar(t)) (dN_i(t) - Y_i(t) risk_i lambda(t))
 // **   keep cumhaz = sum lambda(t) and xhaz = sum xbar(t) lambda(t) as running
 // **     totals, the second is a vector.  xbar requires loops that go from
 // **     last time to first so that we can accumulate it.
 // ** The main loop will go in "spurts": process all obs for this ending time,
 // **   all obs for the next ending time, etc.
 // **   1, when we find a new time, for all obs at this time set initial resid to
 // **       risk * (x_i * cumhaz - xhaz)  (using the old cumhaz and xhaz)
 // **   2. find #events, hazard at this time point, mean at this time point, etc.
 // **     Update xbar, cumhaz and xhaz.
 // **   3. if obs is a death at this time, add (x_i - xbar(t))
 // **   4. at the end of a strata, subtract risk * (x_i*cumhaz - xhaz) for all
 // **     in the strata, and zero temporaries.
 // */




// replace y2 with time and status being passed in.
// covar2 is replaced with covar
// strata2 is replaced with strata
// score2 is replaced with score
// weights2 is replaced whith weights
// method2 is replaced with method

// [[Rcpp::export]]

Rcpp::NumericMatrix coxscore2(Rcpp::NumericVector time,
                              Rcpp::NumericVector status,
                              Rcpp::NumericMatrix covar,
                              Rcpp::IntegerVector strata,
                              Rcpp::NumericVector score,
                              Rcpp::NumericVector weights,
                              int method)
{
  int i,j, k, stratastart;
  int currentstrata;
  double temp;
  double deaths, newtime;
  int dd;
  double xbar;
  double denom=0.0, e_denom;
  double risk;
  double hazard, cumhaz, meanwt;
  double downwt, temp2;

  int n = covar.nrow();
  int nvar = covar.ncol();

  Rcpp::NumericMatrix resid(n, nvar);

  /* scratch space */
  Rcpp::NumericVector a(nvar);
  Rcpp::NumericVector a2(nvar);
  Rcpp::NumericVector xhaz(nvar);

  denom=0.0; cumhaz=0.0;

  for (i=0; i<nvar; i++) {
    a2[i] =0;
    a[i] =0;
    xhaz[i] =0;
  }

  stratastart = n-1;
  currentstrata = strata[n-1];

  for (i=n-1; i >=0; ) {
    newtime = time[i];
    deaths =0; e_denom=0; meanwt =0;
    for (j=0; j< nvar; j++) a2[j] =0;

    for (; i>=0 && time[i]== newtime && strata[i] == currentstrata; i--) {
      /* walk through any tied times */
      risk = score[i] * weights[i];
      denom += risk;
      for (j=0; j<nvar; j++) {
        /* future accumulated risk that new entries don't get */
        resid(i, j) = score[i] * (covar(i, j)*cumhaz - xhaz[j]);
        a[j] += risk * covar(i, j); /* running sum for covariates */
      }
      if (status[i]==1) {
        deaths++;
        e_denom += risk;
        meanwt += weights[i];
        for (j=0; j<nvar; j++)
          a2[j] += risk*covar(i, j);
      }
    }

    if (deaths > 0) { /* update cumhaz and etc */
        if (deaths <2 || method==0) {
          hazard = meanwt/denom;
          cumhaz += hazard;
          for (j=0; j<nvar; j++)  {
            xbar = (a[j]/denom);     /* xbar for this variable */
        xhaz[j] += xbar * hazard;
        for (k=1+i; k<= i+ deaths; k++)
          resid(k, j) += covar(k, j) - xbar;
          }
        }
        else {  /* the harder case, Efron approx */
        /* If there are 3 deaths, the risk set includes all of
         **  them, then 2/3 of each, then 1/3 of each: think of it as
         **  3 separate deaths.  The censored people get all the cumhaz
         **  the deaths only a portion; we 'pre charge' them for the
         **  part of cumhaz and xhaz that they should not get at the
         **  end of the strata.
         */
        meanwt /= deaths;
          for (dd=0; dd<deaths; dd++) {
            downwt = dd/deaths;
            temp = denom - downwt* e_denom;  /* working denominator */
        hazard = meanwt/temp;
        cumhaz += hazard;
        for (j=0; j<nvar; j++) {
          xbar = (a[j] - downwt*a2[j])/ temp;
          xhaz[j] += xbar*hazard;
          for (k=1+i ; k<= i+ deaths; k++) {
            temp2 = covar(k, j) - xbar;
            resid(k, j) += temp2/deaths;
            resid(k, j) += temp2 * score[k] * hazard * downwt;
          }
        }
          }
        }
    }

    if (i<0 || strata[i] != currentstrata) { /* end of a strata */
        /* final work for each obs in the stratum */
        for (k= stratastart; k> i; k--) {
          for (j=0; j < nvar; j++ )
            resid(k, j) += score[k]* (xhaz[j] - covar(k, j)* cumhaz);
        }
        /* reset */
        denom =0; cumhaz=0;
        for (j=0; j<nvar; j++) {
          a[j] =0;
          xhaz[j] =0;
        }
        stratastart = i;
        currentstrata= strata[i];
    }
  }

  return (resid);

}
