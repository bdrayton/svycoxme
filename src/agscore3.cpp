#include <Rcpp.h>
using namespace Rcpp;

//' agscore3
//'
//' Calculate score residuals for a mixed effect cox model. Only calculates for fixed effects.
//'
//' @param tstart  start time
//'	@param tstop   stop time
//'	@param event   event
//' @param strata  unique non-negative for each stratum
//' @param covar2  the matrix of covariates
//' @param score   the vector of subject scores, i.e., exp(beta*z)
//' @param weights case weights
//' @param method  ==1 for efron approx
//' @param sort1 sort order of y[,1], within strata
//' @export
// [[Rcpp::export]]

Rcpp::NumericMatrix agscore3(Rcpp::NumericVector tstart,
							 Rcpp::NumericVector tstop,
							 Rcpp::NumericVector event,
							 Rcpp::NumericMatrix covar,
							 Rcpp::IntegerVector strata,
							 Rcpp::NumericVector score,
							 Rcpp::NumericVector weights,
							 Rcpp::IntegerVector sort1,
							 int method)
{
	int i, j, k;

	int person;
	int currentstrata;
	double denom = 0.0;

	double e_denom;
	double risk, dtime;
	double hazard, meanwt;
	double deaths, downwt;
	int dd;

	double cumhaz = 0.0;
	double d2;

	int n = covar.nrow();
	int nvar = covar.ncol();

	Rcpp::NumericMatrix resid(n, nvar);

	int i1 = n - 1;

	Rcpp::NumericVector a(nvar);	
	Rcpp::NumericVector a2(nvar);	
	Rcpp::NumericVector mean(nvar); 
	Rcpp::NumericVector mh1(nvar);	
	Rcpp::NumericVector mh2(nvar);	
	Rcpp::NumericVector mh3(nvar);	
	Rcpp::NumericVector xhaz(nvar); 

	currentstrata = strata[n - 1];

	for (person = n - 1; person >= 0;)
	{
		dtime = tstop[person];
		if (strata[person] != currentstrata)
		{
			// /* first obs of a new strata, finish off prior one */
			for (; i1 >= 0 && sort1[i1] > person; i1--)
			{
				k = sort1[i1];
				for (j = 0; j < nvar; j++)
				{
					resid(k, j) -= score[k] * (cumhaz * covar(k, j) - xhaz[j]);
				}
			}
			// /* rezero */
			cumhaz = 0;
			denom = 0;
			for (j = 0; j < nvar; j++)
			{
				a[j] = 0;
				xhaz[j] = 0;
			}
			currentstrata = strata[person];
		}
		else
		{
			for (; i1 >= 0 && tstart[sort1[i1]] >= dtime; i1--)
			{
				k = sort1[i1]; /* observation to be removed from risk set */
				if (strata[k] != currentstrata)
					break;
				risk = score[k] * weights[k];
				denom -= risk;
				for (j = 0; j < nvar; j++)
				{
					resid(k, j) -= score[k] * (cumhaz * covar(k, j) - xhaz[j]);
					a[j] -= risk * covar(k, j);
				}
			}
		}

		// /* count up over this time point */
		e_denom = 0;
		meanwt = 0;
		deaths = 0;
		// done at definition.
		// for (i = 0; i < nvar; i++)
		// 	a2[i] = 0;

		for (; person >= 0 && tstop[person] == dtime; person--)
		{
			// /*
			// ** this next line is rare: the first obs of the next strata
			// **  has exactly the same time value as the last person
			// **  of the current strata
			// */
			if (strata[person] != currentstrata)
			{
				break;
			}
			for (j = 0; j < nvar; j++)
			{
				resid(person, j) = (covar(person, j) * cumhaz - xhaz[j]) *
								   score[person];
			}
			risk = score[person] * weights[person];
			denom += risk; // /* denominator of xbar(t) and hazard */
			for (i = 0; i < nvar; i++)
			{
				a[i] += risk * covar(person, i); // /* numerator of xbar(t) */
			}

			if (event[person] == 1)
			{
				deaths++;
				e_denom += risk;
				meanwt += weights[person];
				for (i = 0; i < nvar; i++)
				{
					a2[i] = a2[i] + risk * covar(person, i);
				}
			}
		}
		if (deaths > 0)
		{ // /* update all the values */
			if (deaths < 2 || method == 0)
			{
				// /* easier case */
				hazard = meanwt / denom;
				cumhaz += hazard;
				for (i = 0; i < nvar; i++)
				{
					mean[i] = a[i] / denom;
					xhaz[i] += mean[i] * hazard;
					for (j = person + 1; j <= person + deaths; j++)
					{
						resid(j, i) += covar(j, i) - mean[i];
					}
				}
			}
			else
			{
				// /*
				// ** Efron case.  If there are k deaths, we treat it as k
				// **  separate additions to the hazard.  For the second one
				// **  each death has prob (k-1)/k of being present, then (k-2)/k,
				// **  etc.  The idea is that the deaths actually occur in some
				// **  order, we just don't know what that order is.
				// ** Say k=3 and h1, h2, h3 are the jumps in hazard.  The cumhaz
				// **  and xhaz go up as we would expect.
				// ** The deaths get an addition  (x_i - xbar)/3 at each of the
				// **  3 pseudo death times, and also a "look ahead" correction
				// **  since they don't deserve the full increment of hazard.
				// */
				// done at definition now.
				// for (i = 0; i < nvar; i++)
				// {
				// 	mh1[i] = 0;
				// 	mh2[i] = 0;
				// 	mh3[i] = 0;
				// }
				meanwt /= deaths; // /* average weight of a death */
				for (dd = 0; dd < deaths; dd++)
				{
					downwt = dd / deaths;
					d2 = denom - downwt * e_denom;
					hazard = meanwt / d2;
					cumhaz += hazard;
					for (i = 0; i < nvar; i++)
					{
						mean[i] = (a[i] - downwt * a2[i]) / d2;
						xhaz[i] += mean[i] * hazard;
						mh1[i] += hazard * downwt;
						mh2[i] += mean[i] * hazard * downwt;
						mh3[i] += mean[i] / deaths;
					}
				}

				for (j = person + 1; j <= person + deaths; j++)
				{
					for (i = 0; i < nvar; i++)
					{
						resid(j, i) += (covar(j, i) - mh3[i]) +
									   score[j] * (covar(j, i) * mh1[i] - mh2[i]);
					}
				}
			}
		}
	}

	// /*
	// ** finish those in the final stratum
	// */
	for (; i1 >= 0; i1--)
	{
		k = sort1[i1];
		for (j = 0; j < nvar; j++)
			resid(k, j) -= score[k] * (covar(k, j) * cumhaz - xhaz[j]);
	}

	return (resid);
}
