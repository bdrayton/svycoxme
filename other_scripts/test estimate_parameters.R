

my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 50
my_nk = 10

max_iter = 100
convergence_threshold = 0.00001 # same as coxme

my_X = c("X1", "X2", "X3")



sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

nb <- dplyr::n_distinct(sample_data$M)

start_parameters = c(coef(fit), rep(0, nb))

names(start_parameters) <- c(my_X, paste0("Z", seq_len(nb)))

current_estimates <- estimate_parameters(start_parms = start_parameters, theta = 4, X = my_X, t = "t",
                                         cluster = "M", dij = stat, data = sample_data)

estimate_history <- list(current_estimates)
start_time <- Sys.time()
for (i in 1:max_iter) {

  current_estimates <- try(estimate_parameters(start_parms = current_estimates$new_parms,
                                           theta = current_estimates$new_theta,
                                           X = my_X,
                                           t = "t",
                                           cluster = "M",
                                           dij = stat,
                                           data = sample_data))

  if("try-error" %in% class(current_estimates)) {
    cat("fit failed on iteration", i, "\n")
    break
  }

  estimate_history[[i+1]] <- current_estimates

  biggest_diff = max(abs(c(estimate_history[[i]]$new_theta - estimate_history[[i+1]]$new_theta,
            estimate_history[[i]]$new_parms - estimate_history[[i+1]]$new_parms)), na.rm = TRUE)

  if(biggest_diff <= convergence_threshold){
    cat("converged in", i, "iterations", "\n")
    break
  }

}
end_time <- Sys.time()

end_time - start_time

# get theta estimates

theta_ests <- sapply(estimate_history, "[[", "new_theta")

plot(theta_ests)

coxme_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

coxme::VarCorr(coxme_fit)$M

algo_fit <- tail(estimate_history, 1)


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

r1 <- estimate_parameters_loop(beta = c(1, -0.7, 0.5),
                         theta = 1,
                         k = 50,
                         nk = 10,
                         cluster = "M",
                         max_iter = 100,
                         convergence_threshold = 0.00001,
                         start_theta = 4,
                         X = c("X1", "X2", "X3"))

theta_history <- sapply(r1$estimate_history, "[[", "new_theta")

plot(theta_history)


