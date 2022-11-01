# demonstrate calibration with coxph
library(survey)
library(ggplot2)


# make a population

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
                          random_effect_variance = c(M = 0)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))

})

pop <- Reduce(rbind.data.frame, pop_list)

pop$fX3 <- factor(pop$X3)

# hist(pop$X1, plot = FALSE)

plot(density(pop$X1))

# cut into 10 bits of roughly equal size
X1_cut <- pop$X1

pop$X1_dec <- with(pop, cut(X1, quantile(X1,  0:10/10), labels = 1:10))

############
# the pop true coefs will be a bit different.
# true_coefs = c(1, -0.7, 0.5, -1)
pop_fit <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, data = pop)
true_coefs <- coef(pop_fit)

library(survey)

pop.totals <- colSums(model.matrix(~X1_dec * X3, data = pop))

one_rep <- function(){
# cluster sampling. Also generates scaled weights for coxph
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
  # d3 <- as.svrepdesign(d2)

  d2cal <- calibrate(d2, ~X1_dec * X3, pop.totals,
                     calfun = "raking", bounds = c(0, Inf), maxit = 1000,
                     aggregate.index = ~id)

  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2)
  svycoxph_fit_d2cal <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2cal)

  data.frame(model = rep(c("svycoxph", "svycoxph_calibrated"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(svycoxph_fit_d2),
                          coef(svycoxph_fit_d2cal)),
             rbind(confint(svycoxph_fit_d2),
                   confint(svycoxph_fit_d2cal))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..) |>
    suppressWarnings()

}

one_rep()
# debugonce(one_rep)

reps <- replicate(1000, try(one_rep()), simplify = FALSE)

reps_df <- Reduce(rbind.data.frame, reps)

prop.table(xtabs(~hit + model + X, data = reps_df), margin = c(2, 3))

# CIs from calibrated weights seem to be overly wide

interval_ratio <- function(df) {

  df$range <-  abs(df$X2.5.. - df$X97.5..)

  df$range[1:4]/df$range[5:8]

}

ratios <- Reduce(rbind, lapply(reps, interval_ratio)) |> data.frame()

ratios |>
  tidyr::pivot_longer(cols = c("X1", "X2", "X3", "X4")) |>
  ggplot(aes(value, colour = name)) + geom_density()

# all approx 1
colMeans(ratios)


# error are about the same
reps_df |>
  ggplot(aes(error, colour = model)) + geom_density()

reps_df |>
  ggplot(aes(error, colour = model)) + geom_density() +
  facet_grid(rows = vars(X), scales = "free")

# small bias, not different by model.
reps_df |>
  dplyr::group_by(X, model) |>
  dplyr::summarise(mean(error))



#### without cluster correction

one_rep <- function(){
  # cluster sampling. Also generates scaled weights for coxph
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
  # d3 <- as.svrepdesign(d2)

  d2cal <- calibrate(d2, ~X1_dec * X3, pop.totals,
                     calfun = "raking", bounds = c(0, Inf), maxit = 1000)

  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2)
  svycoxph_fit_d2cal <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2cal)

  data.frame(model = rep(c("svycoxph", "svycoxph_calibrated"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(svycoxph_fit_d2),
                          coef(svycoxph_fit_d2cal)),
             rbind(confint(svycoxph_fit_d2),
                   confint(svycoxph_fit_d2cal))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..) |>
    suppressWarnings()

}

one_rep()
# debugonce(one_rep)

reps2 <- replicate(1000, try(one_rep()), simplify = FALSE)

reps2_df <- Reduce(rbind.data.frame, reps2)

# worse.
prop.table(xtabs(~hit + model + X, data = reps2_df), margin = c(2, 3))


ratios <- Reduce(rbind, lapply(reps2, interval_ratio)) |> data.frame()

# confidence intervals are narrower for variables that have been
# used in the calibration.
ratios |>
  tidyr::pivot_longer(cols = c("X1", "X2", "X3", "X4")) |>
  ggplot(aes(value, colour = name)) + geom_density()

colMeans(ratios)

# error depends on coef. stratum and x3 have a bias.
reps2_df |>
  ggplot(aes(error, colour = model)) + geom_density()

reps2_df |>
  ggplot(aes(error, colour = model)) + geom_density() +
  facet_grid(rows = vars(X), scales = "free")


##### try dump some data

one_rep <- function(){
  # cluster sampling. Also generates scaled weights for coxph
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

  # 20% clusters mcar
  sample_clusters <- dplyr::slice_sample(sample_clusters, prop = .8)

  sample_data <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id")) |>
    dplyr::mutate(scaled_weight = (1/prob)/(1/mean(prob)))

  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data)
  # d3 <- as.svrepdesign(d2)

  d2cal <- calibrate(d2, ~X1_dec * X3, pop.totals,
                     calfun = "raking", bounds = c(0, Inf), maxit = 1000,
                     aggregate.index = ~id)

  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2)
  svycoxph_fit_d2cal <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2cal)

  data.frame(model = rep(c("svycoxph", "svycoxph_calibrated"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(svycoxph_fit_d2),
                          coef(svycoxph_fit_d2cal)),
             rbind(confint(svycoxph_fit_d2),
                   confint(svycoxph_fit_d2cal))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..) |>
    suppressWarnings()

}

one_rep()
# debugonce(one_rep)

reps <- replicate(1000, try(one_rep()), simplify = FALSE)

reps_df <- Reduce(rbind.data.frame, reps)

prop.table(xtabs(~hit + model + X, data = reps_df), margin = c(2, 3))

# CIs from calibrated weights seem to be overly wide

ratios <- Reduce(rbind, lapply(reps, interval_ratio)) |> data.frame()

ratios |>
  tidyr::pivot_longer(cols = c("X1", "X2", "X3", "X4")) |>
  ggplot(aes(value, colour = name)) + geom_density()

# all approx 1 maybe calibrated is a little wider.
colMeans(ratios)

# error are about the same
reps_df |>
  ggplot(aes(error, colour = model)) + geom_density()

reps_df |>
  ggplot(aes(error, colour = model)) + geom_density() +
  facet_grid(rows = vars(X), scales = "free")

# small bias, not different by model.
reps_df |>
  dplyr::group_by(X, model) |>
  dplyr::summarise(mean(error))

head(reps_df)

reps_df |>
  dplyr::filter(X == "X2") |>
  dplyr::arrange(estimate) |>
  dplyr::group_by(model) |>
  dplyr::mutate(est_rank = dplyr::row_number()) |>
  ggplot(aes( y = est_rank, xmin = X2.5.., xmax = X97.5.., colour = hit)) +
    geom_errorbarh() +
    facet_grid(cols = vars(model)) +
    geom_vline(xintercept = -0.7)


# add in calibrated svyrep design.

one_rep <- function(){
  # cluster sampling. Also generates scaled weights for coxph
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

  d2cal <- calibrate(d2, ~X1_dec * X3, pop.totals,
                     calfun = "raking", bounds = c(0, Inf), maxit = 1000,
                     aggregate.index = ~id)

  d3cal <- calibrate(d3, ~X1_dec * X3, pop.totals,
                     calfun = "raking", bounds = c(0, Inf), maxit = 1000,
                     aggregate.index = ~id)

  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum,
                              design = d2)
  svycoxph_fit_d2cal <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum,
                                 design = d2cal)
  svycoxph_fit_d3 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum,
                              design = d3)
  svycoxph_fit_d3cal <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum,
                                 design = d3cal)


  data.frame(model = rep(c("svycoxph", "svycoxph_calibrated", "svycoxph_rep", "svycoxph_repcal"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(svycoxph_fit_d2),
                          coef(svycoxph_fit_d2cal),
                          coef(svycoxph_fit_d3),
                          coef(svycoxph_fit_d3cal)),
             rbind(confint(svycoxph_fit_d2),
                   confint(svycoxph_fit_d2cal),
                   confint(svycoxph_fit_d3),
                   confint(svycoxph_fit_d3cal))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..) |>
    suppressWarnings()

}

one_rep()
# debugonce(one_rep)

reps4m <- replicate(1000, try(one_rep()), simplify = FALSE)

table(sapply(reps4m, class))

# dump the fails
successes <- reps4m

successes[sapply(reps4m, class) == "try-error"] <- NULL

reps_df <- Reduce(rbind.data.frame, successes)

prop.table(xtabs(~hit + model + X, data = reps_df), margin = c(2, 3))

# CIs from calibrated weights seem to be overly wide

ratios <- Reduce(rbind, lapply(successes, interval_ratio)) |> data.frame()

ratios |>
  tidyr::pivot_longer(cols = c("X1", "X2", "X3", "X4")) |>
  ggplot(aes(value, colour = name)) + geom_density()

# all approx 1 maybe calibrated is a little wider.
colMeans(ratios)

# error are about the same
reps_df |>
  ggplot(aes(error, colour = model)) + geom_density()

reps_df |>
  ggplot(aes(error, colour = model)) + geom_density() +
  facet_grid(rows = vars(X), scales = "free")

# small bias, not different by model.
reps_df |>
  dplyr::group_by(X, model) |>
  dplyr::summarise(mean(error))

head(reps_df)

reps_df |>
  dplyr::filter(X == "X2") |>
  dplyr::arrange(estimate) |>
  dplyr::group_by(model) |>
  dplyr::mutate(est_rank = dplyr::row_number()) |>
  ggplot(aes( y = est_rank, xmin = X2.5.., xmax = X97.5.., colour = hit)) +
  geom_errorbarh() +
  facet_grid(cols = vars(model)) +
  geom_vline(xintercept = -0.7)

















