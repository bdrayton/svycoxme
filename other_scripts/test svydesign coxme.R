# test svycoxme.svydesign

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
    dplyr::mutate(weight = 1/prob, scaled_weight = (1/prob)/(1/mean(prob))) |>
    dplyr::ungroup()

  d2 <- survey::svydesign(~id, weights = ~scaled_weight, strata = ~stratum, data = sample_data2)

  # weights are already rescaled.
  svycoxme_fit <- svycoxme(survival::Surv(stat_time, stat) ~ X1+X2+X3+stratum+(1|M1),
                           design = d2, rescale = FALSE)

  coxme_fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1+X2+X3+stratum+(1|M1),
                            data = sample_data2, weights = scaled_weight)



  ci <- rbind(confint(svycoxme_fit), confint(coxme_fit))
  colnames(ci) <- c("lower", "upper")

  point <- c(coef(svycoxme_fit), coef(coxme_fit))

  data.frame(fun  = rep(c("svycoxme", "coxme"), each = 4),
             coef = names(point),
             ci   = ci,
             correct_value = c(1, -0.7, 0.5, -0.5))

}

one_rep()

reps <- replicate(1000, one_rep(), simplify = FALSE)

res <- Reduce(rbind, reps) |> data.frame() |>
  dplyr::mutate(hit = correct_value >= ci.lower & correct_value <= ci.upper)

# coverage at the individual level is a bit worse (0.8 vs 0.9)
# coverage at the cluster level is dismal (0.25-0.3)
res |>
  dplyr::group_by(fun, coef) |>
  dplyr::summarise(mean(hit))



library(ggplot2)

res |>
  tidyr::pivot_longer(cols = c("svycoxme", "coxme")) |>
  ggplot(aes(value, colour = name)) + geom_density() +
  facet_grid(rows = vars(coef))









