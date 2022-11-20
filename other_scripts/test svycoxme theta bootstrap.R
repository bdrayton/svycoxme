# cluster structure in the population.
cluster_str <- data.frame(table(Size = rpois(2.5e5, 2) + 6)) |>
  dplyr::filter(Freq >=10)

cluster_str_list <- split(cluster_str, seq(nrow(cluster_str)))

max_cluster_digits <- max(nchar(as.character(cluster_str$Size)))
max_cluster_freq_digits <- max(nchar(as.character(cluster_str$Freq)))


pop_list <- lapply(cluster_str_list, function(cluster_info){

  k <- cluster_info$Freq
  nk <- as.numeric(as.character(cluster_info$Size))

  k_id <- formatC(k, width = max_cluster_freq_digits, flag = "0")
  nk_id <- formatC(nk, width = max_cluster_digits, flag = "0")

  # the_data <- one_dataset(~X1 + X2 + X3 + stratum + (1 | M1) + (1 | M1:M2),
  the_data <- one_dataset(~X1 + X2 + X3 + stratum + (1 | M1),
                          dists = list(X1 = ~rnorm(n),
                                       X2 = ~rep(rnorm(k), each = nk),
                                       X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                       M1 = ~rep(1:k, each = nk),
                                       M2 = ~rep(rep(c("l","r"), ceiling(nk/2))[seq(nk)], k),
                                       error = ~rexp(n, 10),
                                       stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                                       stratum = ~rep(c(0, 1), c(floor(2/3 * k) * nk, ceiling(1/3 * k) * nk))),
                          dist_args = list(k = k, nk = nk,
                                           n = k * nk),
                          coefficients = c(1, -0.7, 0.5, -0.5),
                          random_effect_variance = list(M1 = 1)
                          # random_effect_variance = list(M1 = 1, `M1:M2` = 1)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M1, sep = "_" ))

})

pop <- Reduce(rbind.data.frame, pop_list)

##### do this 100 times
# extra theta from coxme,
# theta and CI from svycoxme

one_rep <- function(){

  sample_clusters <- dplyr::bind_rows(
    pop |>
      dplyr::filter(stratum == 0) |>
      dplyr::distinct(stratum, id) |>
      dplyr::mutate(prob = 50/dplyr::n()) |>
      dplyr::slice_sample(n = 50),
    pop |>
      dplyr::filter(stratum == 1) |>
      dplyr::distinct(stratum, id) |>
      dplyr::mutate(prob = 50/dplyr::n()) |>
      dplyr::slice_sample(n = 50)
  ) |> dplyr::select(stratum, id, prob)

  sample_data2 <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id")) |>
    # dplyr::mutate(weight = 1/prob, scaled_weight = (1/prob)/(1/mean(prob)),
    #               M1M2 = interaction(M1, M2))
    dplyr::mutate(weight = 1/prob, scaled_weight = (1/prob)/(1/mean(prob)))


  d2 <- survey::svydesign(~id, weights = ~scaled_weight, strata = ~stratum, data = sample_data2)
  d3 <- survey::as.svrepdesign(d2)

  #
  # svycoxme_fit <- svycoxme(survival::Surv(stat_time, stat) ~ X1+X2+X3+stratum+(1|M1)+(1|M1M2),
  #                          design = d3)
  #
  # coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1+X2+X3+stratum+(1|M1)+(1|M1M2),
  #                           data = sample_data2, weights = scaled_weight)
  #

  svycoxme_fit <- svycoxme(survival::Surv(stat_time, stat) ~ X1+X2+X3+stratum+(1|M1),
                           design = d3)

  coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1+X2+X3+stratum+(1|M1),
                            data = sample_data2, weights = scaled_weight)

  svycoxme_theta <- unlist(coxme::VarCorr(svycoxme_fit))

  ci = svycoxme_theta + sqrt(diag(svycoxme_fit$vvar)) %o% qt(c(0.025, 0.975), Inf)

  colnames(ci) <- c("lower", "upper")

  data.frame(coxme_theta = unlist(coxme::VarCorr(coxme_fit)),
             svycoxme_theta = svycoxme_theta, ci)

}

one_rep()

reps <- replicate(100, one_rep(), simplify = FALSE)

library(tidyverse)

reps_df <- Reduce(rbind.data.frame, reps)

reps_df$theta_name = rep(c("M1"), 100)
reps_df$rep = rep(1:100, each = 1)

# hits
reps_df$true_theta = rep(c(1), 100)

reps_df$hits <- with(reps_df, true_theta >= lower & true_theta <= upper)

xtabs(~theta_name + hits, data = reps_df)

reps_df %>%
  filter(theta_name == "M1") %>%
  arrange(desc(svycoxme_theta)) %>%
  mutate(y = row_number()) %>%
  ggplot(aes(x = svycoxme_theta, y = y,xmin = lower, xmax = upper, colour = hits)) +
    geom_point() +
    geom_errorbarh()

# errors
reps_df %>%
  mutate(error1 = coxme_theta - true_theta) %>%
  ggplot(aes(error1, colour = theta_name)) + geom_density()

reps_df %>%
  mutate(error2 = svycoxme_theta - true_theta) %>%
  ggplot(aes(error2, colour = theta_name)) + geom_density()

reps_df %>%
  ggplot(aes(svycoxme_theta, coxme_theta, colour = theta_name)) + geom_point()










