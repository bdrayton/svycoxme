

my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 11
my_nk = 11

my_X = c("X1", "X2", "X3")

ds <- svycoxme::one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta)) %>%
      dplyr::mutate(M2 = rep(c("l", "r"), ceiling(nrow(.)/2))[seq_len(my_k * my_nk)])


debugonce(coxme::coxme)
fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M) + (1|M2), data = ds)

coxme::VarCorr(fit)

coxme_est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

str(fit)

one_dataset

#' the dataset needs
#' IDs - individual, event (for repeated events), and one or more cluster ids.
#' covariates
#' failure or censoring times
#'

#' stuff that might get controlled
#' distributions for covariates, including cluster-level variation
#' number of covariates
#' number of cluster levels and number of clusters
#' number of total observations
#' variance of each set random effects
#' fixed effect coefficents
#' amount of random error
#'


#' The flow:
#' give the function a model formula
#' For each fixed effect, it will seek an expression from a list of expressions also passed in.
#'
#' something similar for random effects.
#'
#' pass in a list of control parameters. this will have elements specifying clustering and sample sizes.
#' As well, fixed effects, error, and random effects variance.
#'
#'
#'
#'

one_dataset <- function(formula, dists, dist_args){

  # decompose formula and check against names in dists.


  vars <- lapply(dists, lazyeval::f_eval, data = dist_args)

  data.frame(vars)

}

ds <- one_dataset(dists = list(X1 = ~rnorm(n, 0, 1),
                               X2 = ~rep(rnorm(k, 0, 1), each = nk),
                               X3 = ~rep(rep(c(1, 0), ceiling(k/2))[seq_len(k)], each = nk),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)]),
                  dist_args = list(k = 10, nk = 10, n = 100))

ds

test_formula <- formula("Y ~ X1 + X2 + X3 + (1|M/M2)")

lme4:::subbars(test_formula)

lme4:::getFixedFormula(test_formula)

lme4:::findbars

coxme:::expand.nested(coxme:::formula1(test_formula)$random[[1]])

coxme:::findIntercept

coxme:::formula1(test_formula)

coxme:::formula2(test_formula)


lme4:::getCovariateFormula(test_formula)


fr.form <- lme4:::subbars(test_formula)
environment(fr.form) <- environment(test_formula)

mf$formula <- fr.form

lme4::lFormula(test_formula)

fr <- eval(mf, parent.frame())

fr <- factorize(fr.form, fr, char.only = TRUE)
attr(fr, "formula") <- formula

n <- nrow(fr)
reTrms <- mkReTrms(lme4:::findbars(lme4:::RHSForm(test_formula)), fr)

lme4::mkReTrms(lme4:::findbars(lme4:::RHSForm(test_formula)), fr)

findbars(test_formula)

lme4:::reOnly(test_formula)


my_fun <- function(k, stuff_using_k){

  lazyeval::f_eval(stuff_using_k, data = list(k = k))

}

my_fun(k = 10, stuff_using_k = ~rnorm(n = k))


f <- ~rnorm(n = k)

lazyeval::f_eval(f)





