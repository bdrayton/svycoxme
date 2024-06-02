



n = 5000

X = data.frame(X1 = rnorm(n),
               X2 = rnorm(n),
               M = rep(c(0,1), each = n/2),
               id = seq_len(n))

dset = draw_event_times(formula = Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = X,
                        coefficients = c(X1 = 1, X2 = 0.5),
                        random_effect_variance = c(M = 0),
                        id = id,
                        baseline_hazard = 1,
                        event = "single")

dset = dplyr::left_join(dset, dplyr::select(X, id, M), by = dplyr::join_by(id))


coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r1 = residuals.coxme(coxmefit, data = dset)

# debugonce(residuals2.coxme)
r2 = residuals2.coxme(coxmefit, data = dset)

all.equal(r1, r2)

microbenchmark::microbenchmark(
  residuals.coxme(coxmefit, data = dset),
  residuals2.coxme(coxmefit, data = dset),
  times = 10
)



dset$strat = rep(c(1:2), n/2)

coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M) + strata(strat),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r3 = residuals2.coxme(coxmefit, data = dset)

summary(coxmefit)



n = 100

X = data.frame(X1 = rnorm(n),
               X2 = rnorm(n),
               M = rep(c(0,1), each = n/2),
               id = seq_len(n))

dset = draw_event_times(formula = Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = X,
                        coefficients = c(X1 = 1, X2 = 0.5),
                        random_effect_variance = c(M = 0),
                        id = id,
                        baseline_hazard = 1,
                        event = "single")

dset = dplyr::left_join(dset, dplyr::select(X, id, M), by = dplyr::join_by(id))


coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r1 = residuals.coxme(coxmefit, data = dset, type = "dfbeta")

# debugonce(residuals2.coxme)
r2 = residuals2.coxme(coxmefit, data = dset, type = "dfbeta")

all.equal(r1, r2)


# test timing of new function.
n = 100000

X = data.frame(X1 = rnorm(n),
               X2 = rnorm(n),
               M = rep(c(0,1), each = n/2),
               id = seq_len(n))

dset = draw_event_times(formula = Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = X,
                        coefficients = c(X1 = 1, X2 = 0.5),
                        random_effect_variance = c(M = 0),
                        id = id,
                        baseline_hazard = 1,
                        event = "single")

dset = dplyr::left_join(dset, dplyr::select(X, id, M), by = dplyr::join_by(id))


coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

microbenchmark::microbenchmark(
  residuals.coxme(coxmefit, data = dset),
  times = 10
)







