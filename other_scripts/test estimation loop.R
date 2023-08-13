


my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M)
my_k = 10
my_nk = 4
my_theta = c(M = 2)
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

ests <- optim(par = start_params,
              fn = lp,
              gr = lp_gr,
              formula = my_formula,
              parsed_data = make_ppl(parsed_data),
              other_args = list(n_fixed = n_fixed,
                                stat = stat,
                                theta = theta_start),
              method = "BFGS",
              control = list(fnscale = -1),
              hessian = TRUE)

current_lp <- lp(ests$par, formula = my_formula, parsed_data = make_ppl(parsed_data),
   other_args = list(n_fixed = n_fixed,
                     stat = stat,
                     theta = theta_start))

# debugonce(theta_ipl2)

# this likelihood is not quadratic. it always has a maximum at 0 (in the limit).
# to make this work, need an expression for theta, which I only have for the
# relatively simple shared frailty model. Additionally, using the optim nest,
# I can use a method that keeps theta in bounds.

theta_ipl2(2, formula = my_formula, parsed_data = parsed_data,
           Kbb = ests$hessian[-(1:n_fixed), -(1:n_fixed)], value = -attr(current_lp, "penalty"))

test_thetas <- seq(0.00001, 5, by = 0.1)

theta_ll <- lapply(test_thetas, theta_ipl2, formula = my_formula, parsed_data = parsed_data,
       Kbb = ests$hessian[-(1:n_fixed), -(1:n_fixed)], value = -attr(current_lp, "penalty"))

plot(test_thetas, theta_ll, type = "l")

theta_ll <- lapply(test_thetas, theta_ipl,
       formula = my_formula, parsed_data = make_ppl(parsed_data), other_args = list(n_fixed = n_fixed,
                                                                                            start_params = start_params,
                                                                                            stat = stat,
                                                                                            re_only = TRUE))
plot(test_thetas, theta_ll, type = "l")












