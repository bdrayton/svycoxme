

# test new estimation setup.

library(Matrix)

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M)
my_k = 10
my_nk = 4
my_theta = c(M = 1)
my_beta = c(1, -0.7, 0.5)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = my_k, nk = my_nk, n = my_k*my_nk),
                  coefficients = my_beta,
                  random_effect_variance = my_theta
)

# compare use of partial hessian (the random effects part) to the full hessian
# in estimation of theta

my_ests <- est_parameters(my_formula, ds, method = "ppl")

my_ests2 <- est_parameters(my_formula, ds, method = "ppl", control = list(re_only = FALSE))

cbind(my_ests$theta, my_ests2$theta)

attr(my_ests, "theta_ests")

attr(my_ests2, "theta_ests")

# try and get a better numerical approximation for the hessian.



####
my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1 | M2)
my_k = 50
my_nk = 10
my_theta = c(M1 = 2, M2 = 1)
my_beta = c(1, -0.7, 0.5)

ds <- one_dataset(my_formula,
                    dists = list(X1 = ~rnorm(n),
                                 X2 = ~rnorm(k * nk),
                                 X3 = ~rbinom(n, 1, 0.5),
                                 M1 = ~rep(1:k, each = nk),
                                 M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                                 error = ~rexp(n, 10),
                                 stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                    dist_args = list(k = my_k, nk = my_nk, n = my_k * my_nk),
                    coefficients = my_beta,
                    random_effect_variance = my_theta
)

# compare use of partial hessian (the random effects part) to the full hessian
# in estimation of theta

debugonce(est_parameters)
my_ests <- est_parameters(my_formula, ds, method = "ppl")

my_ests2 <- est_parameters(my_formula, ds, method = "ppl", control = list(re_only = FALSE))

cbind(my_ests$theta, my_ests2$theta)

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)


start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

theta_ipl(theta_start, my_formula, make_ppl(parsed_data),
          other_args = list(re_only = TRUE,
                            n_fixed = n_fixed,
                            stat = stat,
                            start_params = start_params))


lp.ppl(params = start_params, formula = my_formula, parsed_data = make_ppl(parsed_data),
   other_args = list(n_fixed = n_fixed,
                     stat = stat,
                     theta = theta_start))

lp_gr.ppl(params = start_params, formula = my_formula, parsed_data = make_ppl(parsed_data),
      other_args = list(n_fixed = n_fixed,
                        stat = stat,
                        theta = theta_start))










my_ests_ppl <- est_parameters(my_formula, ds, method = "ppl", start_params = c(my_ests$beta, my_ests$b), theta_start = my_ests$theta)

my_ests_pplextra <- est_parameters(my_formula, ds, method = "pplextra",
                                   start_params = c(my_ests_ppl$beta, my_ests_ppl$b), theta_start = my_ests_ppl$theta)




coxme_fit <- coxme::coxme(my_formula, data = ds)

coxme::VarCorr(coxme_fit)

coxme::random.effects(coxme_fit)
coxme::fixed.effects(coxme_fit)



### test gradient functions for pplextra

# this stuff is needed for each call to lp and lp_grd, but doesn't change, so
# I calculate it once here, and pass it in.

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

theta_start <- get_start_theta(length(parsed_data$reTrms$flist))

n_fixed <- length(attr(terms(coxme:::formula1(my_formula)$fixed), "order"))

fixed_formula <- lme4:::getFixedFormula(my_formula)

fit0 <- survival::coxph(fixed_formula, data = ds)


start_params <- c(coef(fit0),
                  rnorm(n = length(parsed_data$reTrms$Lind),
                        mean = 0,
                        sd = sqrt(theta_start[parsed_data$reTrms$Lind])))

lp(params = start_params, formula = my_formula, parsed_data = make_pplextra(parsed_data),
   other_args = list(n_fixed = n_fixed,
                     stat = stat,
                     theta = my_theta))



# assumes that the response is of the form Surv(time, stat). Behaviour for other Surv formats is undefined.
stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

# debugonce(lp_gr_beta.pplextra)

lp_gr(params = start_params, formula = my_formula, parsed_data = make_ppl(parsed_data),
      other_args = list(n_fixed = n_fixed,
                        stat = stat,
                        theta = my_theta))

lp_gr(params = start_params, formula = my_formula, parsed_data = make_pplextra(parsed_data),
      other_args = list(n_fixed = n_fixed,
                        stat = stat,
                        theta = my_theta))

# test theta_ipl

theta_ipl(1, formula = my_formula, parsed_data = make_ppl(parsed_data),
          other_args = list(n_fixed = n_fixed,
                            stat = stat,
                            theta = my_theta,
                            start_params = start_params))

theta_ipl(1, formula = my_formula, parsed_data = make_pplextra(parsed_data),
          other_args = list(n_fixed = n_fixed,
                            stat = stat,
                            theta = my_theta,
                            start_params = start_params))

beta_b_est <- optim(par = start_params,
                    fn = lp,
                    gr = lp_gr,
                    formula = my_formula,
                    parsed_data = make_ppl(parsed_data),
                    other_args = list(theta = coxme::VarCorr(coxme_fit)$M,
                                      stat = stat),
                    method = "BFGS",
                    control = list(fnscale = -1, trace = 3))


beta_b_est_pplextra <- optim(par = c(0.5, beta_b_est$par[-1]),
                    fn = lp,
                    gr = lp_gr,
                    formula = my_formula,
                    parsed_data = make_pplextra(parsed_data),
                    other_args = list(theta = coxme::VarCorr(coxme_fit)$M,
                                      stat = stat,
                                      n_fixed = n_fixed),
                    method = "BFGS",
                    control = list(fnscale = -1, trace = 3))

cbind(beta_b_est$par, beta_b_est_pplextra$par)

# does pplextra look right when varying one fixed parameter?

# the true value is 1
ll_pplextra <- lapply(seq(0, 2, by = 0.1), function(B){

  lp.pplextra(c(B, beta_b_est$par[-1]), my_formula, parsed_data,
              list(theta = 1,
                   stat = stat,
                   n_fixed = n_fixed))

})

plot(seq(0, 2, by = 0.1),
     ll_pplextra)



beta_b_est_pplextra <- optim(par = beta_b_est$par,
                             fn = lp,
                             gr = NULL,
                             formula = my_formula,
                             parsed_data = make_pplextra(parsed_data),
                             other_args = list(theta = coxme::VarCorr(coxme_fit)$M,
                                               stat = stat,
                                               n_fixed = n_fixed),
                             method = "BFGS",
                             control = list(fnscale = -1, trace = 3))

cbind(beta_b_est$par, beta_b_est_pplextra$par)

# the ppl extra method seems quite bad. it goes out of bounds too easily, which
# makes it really slow (it's slow anyway). estimates are not improved noticably.
# actually using it wont work, because it takes too long. > 12 hours for one
# estimate_parameters call.











