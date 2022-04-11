

my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 50
my_nk = 10

max_iter = 100
convergence_threshold = 0.0001

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

nb <- dplyr::n_distinct(sample_data$M)

start_parameters = c(coef(fit), rep(0, nb))

names(start_parameters) <- c(my_X, paste0("Z", seq_len(nb)))

current_estimates <- estimate_parameters(start_parms = start_parameters, theta = 4, X = my_X, t = "t",
                                         cluster = "M", dij = stat, data = sample_data)

estimate_history <- list(current_estimates)

for (i in 1:max_iter) {

  current_estimates <- try(estimate_parameters(start_parms = current_estimates$new_parms,
                                           theta = current_estimates$new_theta,
                                           X = my_X,
                                           t = "t",
                                           cluster = "M",
                                           dij = stat,
                                           data = sample_data))

  if("try-error" %in% class(current_estimates)) {
    cat("fit failed on iteration", i)
    break
  }

  estimate_history[[i+1]] <- current_estimates

  biggest_diff = max(abs(c(estimate_history[[i]]$new_theta - estimate_history[[i+1]]$new_theta,
            estimate_history[[i]]$new_parms - estimate_history[[i+1]]$new_parms)), na.rm = TRUE)

  if(biggest_diff <= convergence_threshold){
    cat("converged in", i, "iterations")
    break
  }

}


# get theta estimates

theta_ests <- sapply(estimate_history, "[[", "new_theta")

plot(theta_ests)


algo_fit <- tail(estimate_history, 1)

coxme_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

data.frame(
  true_vals = c(my_theta, my_beta, attr(sample_data, "random_effects")),
  coxme = c(coxme::VarCorr(coxme_fit)$M,
    coxme::fixef(coxme_fit),
    coxme::ranef(coxme_fit)$M),
  my_ests = unlist(algo_fit))




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
  dplyr::mutate(iteration = as.numeric(gsub("V", "", name))) %>%
  ggplot(aes(iteration, value)) +
    geom_line() +
    facet_grid(rows = vars(rowname),
               scales = "free")











