
# compare coxph and svycoxph


library(survey)
library(coxme)

# cluster structure in the population.
cluster_str <- data.frame(table(Size = rpois(2.5e5, 2) + 2)) |>
  dplyr::filter(Freq >=10)

cluster_str_list <- split(cluster_str, seq(nrow(cluster_str)))

max_cluster_digits <- max(nchar(as.character(cluster_str$Size)))
max_cluster_freq_digits <- max(nchar(as.character(cluster_str$Freq)))

pop_list <- lapply(cluster_str_list, function(cluster_info){

  k <- cluster_info$Freq
  nk <- as.numeric(as.character(cluster_info$Size))

  k_id <- formatC(k, width = max_cluster_freq_digits)
  nk_id <- formatC(nk, width = max_cluster_digits)

  the_data <- one_dataset(~X1 + X2 + X3 + stratum + (1 | M),
                          dists = list(X1 = ~rnorm(n),
                                       X2 = ~rep(rnorm(k), each = nk),
                                       X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                       M = ~rep(1:k, each = nk),
                                       error = ~rexp(n, 10),
                                       stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                                       stratum = ~rep(c(0, 1), c(floor(2/3 * k) * nk, ceiling(1/3 * k) * nk))),
                          dist_args = list(k = k, nk = nk,
                                           n = k * nk),
                          coefficients = c(1, -0.7, 0.5, -1),
                          random_effect_variance = c(M=0)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))

})

pop <- Reduce(rbind.data.frame, pop_list)

############
# the pop true coefs will be a bit different.

true_coefs = c(1, -0.7, 0.5, -1)



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

  sample_data <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id")) |>
    dplyr::mutate(scaled_weight = (1/prob)/(1/mean(prob)))

  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data)
  d3 <- as.svrepdesign(d2)


  coxph_robust_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, data = sample_data, robust = TRUE, weights = scaled_weight)
  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2)
  svycoxph_rep_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d3)

  data.frame(model = rep(c("coxph_weighted_robust", "svycoxph", "svycoxph_rep"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(coxph_robust_fit_d2),
                          coef(svycoxph_fit_d2),
                          coef(svycoxph_rep_fit_d2)),
             rbind(confint(coxph_robust_fit_d2),
                   confint(svycoxph_fit_d2),
                   confint(svycoxph_rep_fit_d2))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..)

}

one_rep()

reps <- replicate(1000, one_rep(), simplify = FALSE)

reps_df <- Reduce(rbind.data.frame, reps)

prop.table(xtabs(~hit + model + X, data = reps_df), margin = c(2, 3))


