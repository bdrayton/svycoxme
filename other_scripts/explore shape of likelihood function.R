
# i want to better understand the shape of the likelihood function, and why it
# needs to have the local optimum of beta and b to be the right shape.
# I'm hoping I can get some intuition about this by looking at plots with just
# one fixed effect and one random effect variance.
# In addition, I'll use this framework to explore the gradient function.
# Currently it doesn't evaluate to 0 at the maximum, which is weird.


my_formula <- survival::Surv(stat_time, stat)~ X1 + (1 | M1)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               M1 = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1),
                  random_effect_variance = c(M1 = 1)
)

est_parameters(my_formula, data = ds)

test_thetas <- seq(0.001, 3, by = 0.05)

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

test_thetas <- seq(0.01, 3, by = 0.05)

theta_ll <- lapply(test_thetas, theta_ipl2, formula = my_formula, parsed_data = parsed_data,
                   Kbb = ests$hessian[-(1:n_fixed), -(1:n_fixed)], value = -attr(current_lp, "penalty"))

plot(test_thetas, theta_ll, type = "l")

theta_ll <- lapply(test_thetas, theta_ipl,
                   formula = my_formula, parsed_data = make_ppl(parsed_data), other_args = list(n_fixed = n_fixed,
                                                                                                start_params = start_params,
                                                                                                stat = stat,
                                                                                                re_only = FALSE))
theta_ll[[1]]


dets <- sapply(theta_ll, attr, "det(D)")
hess <- sapply(theta_ll, attr, "detKbb")
values <- sapply(theta_ll, attr, "value")

plot(test_thetas, theta_ll, ylim = range(c(dets, hess, values)))
points(test_thetas, dets)
points(test_thetas, hess)
points(test_thetas, values)
points(test_thetas, -0.5 * dets -0.5 * hess)


# what is the penalty doing?

penalties <- sapply(theta_ll, function(ll){

  beta_b <- attr(ll, "beta_b")

  b <- Matrix(beta_b[-seq(n_fixed)], ncol = 1)

  D <- attr(ll, "D")

  t(b) %*% solve(D) %*% b

})

penalties





# modify theta_ipl so that it returns the optimum beta and b


theta_ipl(0.5, formula = my_formula,
          parsed_data = make_ppl(parsed_data),
          other_args = list(n_fixed = n_fixed,
                            start_params = start_params,
                            stat = stat,
                            re_only = TRUE))


# modify theta_ipl so it returns the components
# determinant(D)$modulus
# determinant(Kbb)$modulus
# value









