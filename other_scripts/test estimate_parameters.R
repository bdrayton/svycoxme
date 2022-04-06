
my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 10
my_nk = 10

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

nb <- dplyr::n_distinct(sample_data$M)

start_parameters = c(coef(fit), rep(0, nb))

names(start_parameters) <- c(my_X, paste0("Z", seq_len(nb)))


current_estimates <- estimate_parameters(start_parms = start_parameters, theta = 1, X = my_X, t = t,
                                         cluster = "M", dij = stat, data = sample_data)

max_iter = 10
convergence_threshold = 0.000001

estimate_history <- list(current_estimates)

for (i in 1:max_iter) {

  current_estimates <- estimate_parameters(start_parms = current_estimates$new_parms,
                                           theta = current_estimates$new_theta,
                                           X = my_X,
                                           t = t,
                                           cluster = "M",
                                           dij = stat,
                                           data = sample_data)

  estimate_history[[i+1]] <- current_estimates

  biggest_diff = max(abs(c(estimate_history[[i]]$new_theta - estimate_history[[i+1]]$new_theta,
            estimate_history[[i]]$new_parms - estimate_history[[i+1]]$new_parms)), na.rm = TRUE)

  if(biggest_diff <= convergence_threshold){
    cat("converged in", i, "iterations")
    break
  }

}

tail(estimate_history, 1)

coxme_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

coxme::ranef(coxme_fit)



# get theta estimates

theta_ests <- sapply(estimate_history, "[[", "new_theta")

plot(theta_ests)

# get the fixed effect estimates

parm_ests <- sapply(estimate_history, "[[", "new_parms")

parm_ests %>%
  as.data.frame() %>%
  tibble::rownames_to_column() %>%
  tidyr::pivot_longer(cols = dplyr::starts_with("V")) %>%
  dplyr::filter(grepl("X", rowname)) %>%
  ggplot(aes(name, value)) + geom_point() + facet_grid(rows = vars(rowname), scales = "free")











