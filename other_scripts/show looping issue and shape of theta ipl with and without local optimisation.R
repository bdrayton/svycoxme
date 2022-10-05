

# illustrate the looping issue.

my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1| M1:M2)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1, -0.7, 0.5),
                  random_effect_variance = c(M1 = 1, `M1:M2` = 2)
)


ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

parsed_data$reTrms$flist




ests_list <- list()

# debugonce(est_parameters)
ests_list[[1]] <- est_parameters(my_formula, data = ds,
                                 control = control.list(factr = 1e3,
                                                        reltol = 1e-13))

for(i in 1:10) {

  current_params <- ests_list[[i]]

  ests_list[[i + 1]] <- est_parameters(my_formula, data = ds,
                                       start_params = c(ests_list[[i]]$beta, ests_list[[i]]$b),
                                       theta_start = ests_list[[i]]$theta,
                                       control = control.list(factr = 1e3,
                                                              reltol = 1e-13))

}


# extract theta
lapply(ests_list, names)
thetas_df <- sapply(ests_list, "[[", "theta", simplify = "matrix") |>
  data.frame() |>
  mutate(theta = factor(c("theta M1", "theta M2")), .before = everything()) |>
  pivot_longer(cols = X1:X11) |>
  mutate(Iteration = factor(name, levels = paste0("X", 1:11), labels = as.character(1:11)))

thetas_df %>%
  group_by(theta) %>%
  summarise(min(value) - max(value))

#zigzagging?
ggplot(thetas_df, aes(Iteration, value)) + geom_point() +
  facet_grid(rows = vars(theta), scales = "free")


# same with betas
betas_df <- sapply(ests_list, "[[", "beta", simplify = "matrix") |>
  data.frame() |>
  mutate(beta = factor(paste("beta", 1:3)), .before = everything()) |>
  pivot_longer(cols = X1:X11) |>
  mutate(Iteration = factor(name, levels = paste0("X", 1:11), labels = as.character(1:11)))

#zigzagging?
ggplot(betas_df, aes(Iteration, value)) + geom_point() +
  facet_grid(rows = vars(beta), scales = "free")

betas_df %>%
  group_by(beta) %>%
  summarise(min(value) - max(value))


# what's happening with the hessians now?

lapply(ests_list, function(x){
  sqrt(diag(-solve(attr(x, "theta_est")$hessian)))
})










# show what happens with iterative optimisation - theta ipl is not quadratic.

# for fixed beta and b,

ds_sorted <- sortAndIndex(ds, sort_vars = stat_time)

parsed_data <- lme4::lFormula(my_formula, data = ds_sorted)

stat <- Matrix(unclass(parsed_data$fr[,1])[, "status"], ncol = 1)

beta_b_est <- attr(ests_list[[1]], "beta_b_est")

theta_ipl2(theta = ests_list[[1]]$theta, formula = my_formula, parsed_data = parsed_data,
           Kbb = beta_b_est$hessian[-(1:3), -(1:3)], value = beta_b_est$value)

test_thetas <- data.frame(seq(0.001, 5, by = 0.1), ests_list[[1]]$theta[2])

test_theta_ll <- apply(test_thetas, 1, function(thetas){

  theta_ipl2(theta = thetas, formula = my_formula, parsed_data = parsed_data,
             Kbb = beta_b_est$hessian[-(1:3), -(1:3)], value = beta_b_est$value)

})

plot(test_thetas[,1], test_theta_ll)

theta_ipl(ests_list[[1]]$theta, formula = my_formula, parsed_data = make_ppl(parsed_data),
          other_args = list(start_params = beta_b_est$par,
                            re_only = TRUE,
                            n_fixed = 3,
                            stat = stat))

test_theta_ll2 <- apply(test_thetas, 1, function(thetas){

  theta_ipl(thetas, formula = my_formula, parsed_data = make_ppl(parsed_data),
            other_args = list(start_params = beta_b_est$par,
                              re_only = TRUE,
                              n_fixed = 3,
                              stat = stat))

})

plot(test_thetas[,1], test_theta_ll2)
abline(v = ests_list[[1]]$theta[1])


test_thetas <- data.frame(ests_list[[1]]$theta[1],seq(0.1, 7, by = 0.1))

test_theta_ll2 <- apply(test_thetas, 1, function(thetas){

  theta_ipl(thetas, formula = my_formula, parsed_data = make_ppl(parsed_data),
            other_args = list(start_params = beta_b_est$par,
                              re_only = TRUE,
                              n_fixed = 3,
                              stat = stat))

})

plot(test_thetas[-1,2], test_theta_ll2[-1])
abline(v = ests_list[[1]]$theta[2])



# if you can get a se, bootstrap the t stat.
# second order improvements.

# athol anderson has a book too.

# Rao and Wu - stratified weighted version.

# for clustering, need to resample culstering.

# if 5 per stratum, resample 4 instead. fpc stuff.


# need to think carefully about the sampling.
# - population vs sampling uncertainty.

# and then clustering - consider calibration
# - need to calibrate after you make the replicate weights.

# calibration after linearisation


# shudong's thesis - available in the library.

# weights - point estimates might be wrong
#         - variance estimates will be wrong

# consider convergence thresholds (tighten), and step sizes.
# convergence of the numerical derivative.

# It has to do with the interplay between the convergence threshold of the
# likelihood, and the step size for the numerical derivative. For
# smaller convergence threshold, the likelihood will look smooth if the step size
# of the numerical derivative is bigger.















