# cluster structure in the population.
cluster_str <- data.frame(table(Size = rpois(2.5e5, 2) + 2)) |>
  dplyr::filter(Freq >=10)

cluster_str_list <- split(cluster_str, seq(nrow(cluster_str)))

max_cluster_digits <- max(nchar(as.character(cluster_str$Size)))
max_cluster_freq_digits <- max(nchar(as.character(cluster_str$Freq)))

pop_list <- lapply(cluster_str_list, function(cluster_info){

  k <- cluster_info$Freq
  nk <- as.numeric(as.character(cluster_info$Size))

  k_id <- formatC(k, width = max_cluster_freq_digits, flag = "0")
  nk_id <- formatC(nk, width = max_cluster_digits, flag = "0")

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
                          coefficients = c(1, -0.7, 0.5, -0.5),
                          random_effect_variance = list(M = 1)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))

})

pop <- Reduce(rbind.data.frame, pop_list)

#####

debug(survey:::svycoxph.svyrep.design)

one_rep <- function(true_coefs = c(1, -0.7, 0.5)){

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
    dplyr::mutate(weight = 1/prob, scaled_weight = (1/prob)/(1/mean(prob)))

  d2 <- survey::svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data2)
  d3 <- survey::as.svrepdesign(d2)

  mixed_effects_formula <- survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M)
  fixed_effects_formula <- lme4:::getFixedFormula(mixed_effects_formula)

  svycoxph_fit_d2 <- try(survey::svycoxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2))
  svycoxph_fit_d3 <- try(survey::svycoxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, design = d3))
  svycoxme_fit_d3 <- try(svycoxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), design = d3))
  coxme_fit <- try(coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data2))
  coxme_fit_weighted <- try(coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M),
                                         data = sample_data2, weights = weight))
  coxme_fit_scale_weighted <- try(coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M),
                                               data = sample_data2, weights = scaled_weight))

  list(f1 = svycoxph_fit_d2,
       f2 = svycoxph_fit_d3,
       f3 = svycoxme_fit_d3,
       f4 = coxme_fit,
       f5 = coxme_fit_weighted,
       f6 = coxme_fit_scale_weighted)

# # leave processing out of the rep. makes it easier to catch errors and handle them later.
#   data.frame(model = rep(c("svy_coxph", "svrep_coxph", "svrep_coxme"), each = length(true_coefs)),
#              X = names(coef(svycoxph_fit_d2)),
#              true_value = true_coefs,
#              estimate = c(coef(svycoxph_fit_d2), coef(svycoxph_fit_d3), coef(svycoxme_fit_d3)),
#              rbind(confint(svycoxph_fit_d2),
#                    confint(svycoxph_fit_d3),
#                    confint(svycoxme_fit_d3))) |>
#     dplyr::mutate(error = estimate - true_value,
#                   hit = true_value >= X2.5.. & true_value <= X97.5..)

}

one_rep()

my_reps <- replicate(1000, one_rep(), simplify = FALSE)

my_models <- unlist(my_reps, recursive = FALSE)

is_error <- sapply(my_models, function(x) "try-error" %in%  class(x))

table(is_error)

table(unlist(lapply(my_models, class)))

f1_reps <- lapply(my_reps, "[[", "f1")
f2_reps <- lapply(my_reps, "[[", "f2")
f3_reps <- lapply(my_reps, "[[", "f3")
f4_reps <- lapply(my_reps, "[[", "f4")
f5_reps <- lapply(my_reps, "[[", "f5")
f6_reps <- lapply(my_reps, "[[", "f6")



# which models failed?
sapply(f1_reps, function(x) "try-error" %in%  class(x)) |> table()
sapply(f2_reps, function(x) "try-error" %in%  class(x)) |> table()
sapply(f3_reps, function(x) "try-error" %in%  class(x)) |> table()
sapply(f4_reps, function(x) "try-error" %in%  class(x)) |> table()
sapply(f5_reps, function(x) "try-error" %in%  class(x)) |> table()
sapply(f6_reps, function(x) "try-error" %in%  class(x)) |> table()

# drop errors
f3_errors <- sapply(f3_reps, function(x) "try-error" %in%  class(x))
f3_reps[f3_errors] <- NULL

f4_errors <- sapply(f4_reps, function(x) "try-error" %in%  class(x))
f4_reps[f4_errors] <- NULL

f5_errors <- sapply(f5_reps, function(x) "try-error" %in%  class(x))
f5_reps[f5_errors] <- NULL




f1_coef <- Reduce(rbind, lapply(f1_reps, coef)) |> data.frame(model = "svycoxph.svydesign")

f2_coef <- Reduce(rbind, lapply(f2_reps, coef)) |> data.frame(model = "svycoxph.svyrep")

f3_coef <- Reduce(rbind, lapply(f3_reps, coef)) |> data.frame(model = "svycoxme.svyrep")

f4_coef <- Reduce(rbind, lapply(f4_reps, coef)) |> data.frame(model = "coxme")

f5_coef <- Reduce(rbind, lapply(f5_reps, coef)) |> data.frame(model = "coxme.weights")

f6_coef <- Reduce(rbind, lapply(f6_reps, coef)) |> data.frame(model = "coxme.scaled_weights")




plot_data <- rbind.data.frame(f1_coef, f2_coef, f3_coef, f4_coef, f5_coef, f6_coef) |>
  tidyr::pivot_longer(cols = c("X1", "X2", "X3"), names_to = "parameter") |>
  dplyr::filter(dplyr::between(value, -5, 2.5))

library(ggplot2)
plot_data |>
  ggplot(aes(value, colour = model)) + geom_density() + facet_grid(rows = vars(parameter), scales = "free") +
  geom_vline(data = dplyr::filter(plot_data, parameter == "X1"), aes(xintercept = 1)) +
  geom_vline(data = dplyr::filter(plot_data, parameter == "X2"), aes(xintercept = -0.7)) +
  geom_vline(data = dplyr::filter(plot_data, parameter == "X3"), aes(xintercept = 0.5))

plot_data |>
  dplyr::filter(parameter == "X1") |>
  ggplot(aes(value, colour = model)) + geom_density() +
  geom_vline(aes(xintercept = 1))

plot_data |>
  dplyr::filter(parameter == "X2") |>
  ggplot(aes(value, colour = model)) + geom_density() +
  geom_vline(aes(xintercept = -0.7))

plot_data |>
  dplyr::filter(parameter == "X3") |>
  ggplot(aes(value, colour = model)) + geom_density() +
  geom_vline(aes(xintercept = 0.5))


# compare the best fits
plot_data |>
  dplyr::filter(parameter == "X1", model %in% c("svycoxme.svyrep", "coxme.scaled_weights")) |>
  ggplot(aes(value, colour = model)) + geom_density() +
  geom_vline(aes(xintercept = 1))

plot_data |>
  dplyr::filter(parameter == "X2", model %in% c("svycoxme.svyrep", "coxme.scaled_weights")) |>
  ggplot(aes(value, colour = model)) + geom_density() +
  geom_vline(aes(xintercept = -0.7))

plot_data |>
  dplyr::filter(parameter == "X3", model %in% c("svycoxme.svyrep", "coxme.scaled_weights")) |>
  ggplot(aes(value, colour = model)) + geom_density() +
  geom_vline(aes(xintercept = 0.5))





# just look at one svycoxme fit

f3_reps[[1]]$var

my_confit <- function(fit){
  x = coef(fit)

  v = fit$var

  lower = x - 1.96 * sqrt(diag(v))
  upper = x + 1.96 * sqrt(diag(v))

  res <- cbind(lower, upper)

  rownames(res) <- names(x)
  colnames(res) <- c("2.5 %", "97.5 %")

  res

}

is_hit <- function(fit, true_values){

  CI <- confint(fit)

  true_values >= CI[,1] & true_values <= CI[,2]

}

is_hit(f3_reps[[1]], true_values = c(1, -0.7, 0.5))

i = 1

models = c("svycoxph.svydesign", "svycoxph.svyrep", "svycoxme.svyrep", "coxme","coxme.weights", "coxme.scaled_weights")

for(i in seq_along(models)){
hits <- eval(parse(text = glue::glue("lapply(f{i}_reps, is_hit, true_values = c(1, -0.7, 0.5))")))

cat(models[i],": ", mean(Reduce(rbind, hits), na.rm = TRUE), "\n")

}







data.frame(Reduce(rbind,f3_coef)) |>
  tidyr::pivot_longer(cols = c("X1", "X2", "X3"), names_to = "parameter") |>
  ggplot(aes(value)) + geom_density() + facet_grid(rows = vars(parameter))

estimates_list <- lapply(c(f1_reps, f2_reps, f3_reps), function(model) {

  beta <- coef(model)

  data.frame(true_value = c(1, -0.7, 0.5),
        parameter = names(beta),
        estimate = beta,
        confint(model))

  }
)

estimates_df <- Reduce(rbind.data.frame, estimates_list) |>
  dplyr::mutate(error = estimate - true_value,
                hit = true_value >= X2.5.. & true_value <= X97.5..)




estimates_df$model <- rep(c("coxph", "coxph", "coxme"), c(length(f1_reps),
                                                          length(f2_reps),
                                                          length(f3_reps)))

estimates_df$error.method <- rep(c("svy", "svrep", "svrep"), c(length(f1_reps),
                                                          length(f2_reps),
                                                          length(f3_reps)))

head(estimates_df)

xtabs(~hit + parameter + paste(model, error.method, sep = "."), data = estimates_df) |>
  prop.table(margin = c(2,3)) |> round(digits = 3)

library(ggplot2)
estimates_df |> ggplot(aes(estimate, colour = paste(model, error.method, sep = "."))) +
  geom_density() + facet_grid(rows = vars(parameter))


pop |>
  dplyr::group_by(id) |>
  dplyr::summarise(cluster_size = dplyr::n()) |>
  dplyr::count(cluster_size)


# select 2 per cluster so that pr selections vary.

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
  dplyr::group_by(id) |>
  dplyr::mutate(pr_selection = prob * 2/dplyr::n()) |>
  dplyr::slice_sample(n = 2)

d2 <- survey::svydesign(~id, probs = ~pr_selection, strata = ~stratum, data = sample_data2)
d3 <- survey::as.svrepdesign(d2)

svycoxph_fit_d2 <- survey::svycoxph(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2)

svycoxph_fit_d3 <- svycoxme::svycoxph.svyrep.design(survival::Surv(stat_time, stat) ~ X1 + X2 + X3, design = d3)

debugonce(svycoxme.svyrep.design)
svycoxme_fit_d3 <- svycoxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), design = d3)


data.frame(model = rep(c("svy_coxph", "svrep_coxph", "svrep_coxme"), each = length(true_coefs)),
           X = names(coef(svycoxph_fit_d2)),
           true_value = true_coefs,
           estimate = c(coef(svycoxph_fit_d2), coef(svycoxph_fit_d3), coef(svycoxme_fit_d3)),
           rbind(confint(svycoxph_fit_d2),
                 confint(svycoxph_fit_d3),
                 confint(svycoxme_fit_d3))) |>
  dplyr::mutate(error = estimate - true_value,
                hit = true_value >= X2.5.. & true_value <= X97.5..)

