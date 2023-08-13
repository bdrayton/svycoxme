# survey version of coxme

#' formula      passed in via formula. Will need to check how cluster(), strata(), offset() etc are handled.
#' data         passed in via design
#' subset       handled prior to coxme call.
#' weights      passed in via design
#' na.action    pass in via ...
#' init         pass in via ...
#' control = coxme.control passed in via ...
#' ties         pass in via ...
#' varlist      pass in via ...
#' vfixed       pass in vai ...
#' vinit        pass in via ...
#' x            ''
#' y            ''
#' refine.n     ''
#' fixed   depreciated
#' random  depreciated
#' variance depreciated
#' ... passed to coxme.control

?svycoxph.survey.design
?svycoxph.survey.design2
?svycoxph.svyrep.design

# svycoxph.svyrep.design
# formula - needs to handle cluster()?
# design
# subset
# rescale
#

library(survey)
library(coxme)

ds <- one_dataset(~X1 + X2 + X3 + (1 | M),
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M = ~rep(1:k, each = nk),
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 200),
                  coefficients = c(1, 1, 1),
                  random_effect_variance = list(M = 1)
)

coxfit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = ds)

formula(coxfit)
subset
weights
na.action
init



residuals(coxfit)

summary(coxfit)

# generate a population
# 100-1000 clusters of 2-10 observations # consider aligning cluster size/prevalence with interRAI

replicate(10000, mean(round(runif(1000, min = 1, max = 9)))) |> mean()
# i want about 1e6 observations.
# mean cluster size = 4
# need 2.5e5 clusters.

replicate(10000, mean(rpois(1000, 2) + 2)) |> mean()

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
              coefficients = c(1, -0.7, 0.5, -0.5),
              random_effect_variance = list(M = 2)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))

})


pop <- Reduce(rbind.data.frame, pop_list)

# id identifies clusters.
dplyr::n_distinct(pop$id)

distinct_ids <- unique(pop$id)

# sample id, and then use sampled ids to filter data.

chosen_few <- sample(distinct_ids, 150)

sample_data <- pop |> dplyr::filter(id %in% dplyr::all_of(chosen_few))

sample_data$sampling_weight <- dplyr::n_distinct(pop$id) / 150

# try stratified sampling, using X3 to stratify. Take 2/3 from X3 == 0

sample_data2 <- dplyr::bind_rows(
  pop |>
    dplyr::filter(stratum == 1) |>
    dplyr::mutate(sampling_weight = dplyr::n()/50) |>
    dplyr::slice_sample(n = 50),
  pop |>
    dplyr::filter(stratum == 0) |>
    dplyr::mutate(sampling_weight = dplyr::n()/100) |>
    dplyr::slice_sample(n = 100)
)

library(survey)

d1 <- svydesign(~id, weights = ~sampling_weight, data = sample_data)

data.frame(svytable(~stat, design = d1))

table(pop$stat)

one_rep <- function(){

  chosen_few <- sample(distinct_ids, 150)

  sample_data <- pop |> dplyr::filter(id %in% dplyr::all_of(chosen_few))

  sample_data$sampling_weight <- dplyr::n_distinct(pop$id) / 150

  d1 <- svydesign(~id, weights = ~sampling_weight, data = sample_data)

  ci <- confint(svymean(~stat, design = d1))

  0.8 >= ci[1] & 0.8 <= ci[2]

}

one_rep()

in_ci <- replicate(1000, one_rep())

mean(in_ci)

#####

one_rep <- function(){

  sample_data2 <- dplyr::bind_rows(
    pop |>
      dplyr::filter(stratum == 1) |>
      dplyr::mutate(sampling_weight = dplyr::n()/50) |>
      dplyr::slice_sample(n = 50),
    pop |>
      dplyr::filter(stratum == 0) |>
      dplyr::mutate(sampling_weight = dplyr::n()/100) |>
      dplyr::slice_sample(n = 100)
  )

  d2 <- svydesign(~id, weights = ~sampling_weight, strata = ~stratum, data = sample_data2)

  ci <- confint(svymean(~stat, design = d2))

  0.8 >= ci[1] & 0.8 <= ci[2]

}

in_ci <- replicate(1000, one_rep())

mean(in_ci)

# ignore the random effect, and fit the models.

coxfit_d1 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d1)

summary(coxfit_d1)

# fit with coxme
library(coxme)
coxmefit_d1 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data)

summary(coxmefit_d1)

#### same thing with stratified model.

coxph_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3, data = sample_data2)
svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2)
coxme_fit_d2 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data2)

coef(coxph_fit_d2)
coef(svycoxph_fit_d2)
coef(coxme_fit_d2)

confint(coxph_fit_d2)
confint(svycoxph_fit_d2)
confint(coxme_fit_d2)

# draw a stratified sample, fit the 3 models, extract ci hits, coef errors.

true_coefs <- c(1, -0.7, 0.5)

one_rep <- function(){

  sample_data2 <- dplyr::bind_rows(
    pop |>
      dplyr::filter(stratum == 1) |>
      dplyr::mutate(sampling_weight = dplyr::n()/50) |>
      dplyr::slice_sample(n = 50),
    pop |>
      dplyr::filter(stratum == 0) |>
      dplyr::mutate(sampling_weight = dplyr::n()/100) |>
      dplyr::slice_sample(n = 100)
  )

  d2 <- svydesign(~id, weights = ~sampling_weight, strata = ~stratum, data = sample_data2)

  coxph_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3, data = sample_data2)
  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2)
  coxme_fit_d2 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data2)

  data.frame(model = rep(c("coxhp", "svycoxph", "coxme"), each = length(true_coefs)),
             X = names(coef(coxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(coxph_fit_d2), coef(svycoxph_fit_d2), coef(coxme_fit_d2)),
             rbind(confint(coxph_fit_d2),
                   confint(svycoxph_fit_d2),
                   confint(coxme_fit_d2))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..)

}



reps <- replicate(1000, one_rep(), simplify = FALSE)

############

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

  sample_data2 <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id"))

  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data2)

  coxph_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3, data = sample_data2)
  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2)
  coxme_fit_d2 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data2)

  data.frame(model = rep(c("coxph", "svycoxph", "coxme"), each = length(true_coefs)),
             X = names(coef(coxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(coxph_fit_d2), coef(svycoxph_fit_d2), coef(coxme_fit_d2)),
             rbind(confint(coxph_fit_d2),
                   confint(svycoxph_fit_d2),
                   confint(coxme_fit_d2))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..)
}

one_rep()
############

debug(survey:::svycoxph.svyrep.design)

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

  sample_data2 <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id"))

  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data2)
  d3 <- as.svrepdesign(d2)
  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2)
  svycoxph_fit_d3 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d3)

  data.frame(model = rep(c("svycox", "svrepcox"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(svycoxph_fit_d2), coef(svycoxph_fit_d3)),
             rbind(confint(svycoxph_fit_d2),
                   confint(svycoxph_fit_d3))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..)
}

one_rep()

############
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

  sample_data2 <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id"))

  # scale weights
  sample_data2 <- sample_data2 |> dplyr::mutate(scaled_weight = (1/prob)/(1/mean(prob)))

  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data2)
  d3 <- svydesign(~id, weights = ~scaled_weight, strata = ~stratum, data = sample_data2)

  coxph_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3, data = sample_data2)
  coxph_robust_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3, data = sample_data2, robust = TRUE, weights = scaled_weight)
  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d2)
  svycoxph_fit_d3 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3, design = d3)
  coxme_fit_d2 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data2)
  coxme_fit_d3 <- coxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M), data = sample_data2, weights = scaled_weight)

  data.frame(model = rep(c("coxph", "coxph_weighted_robust", "svycoxph", "svycoxph_scaled", "coxme", "coxme_weighted"), each = length(true_coefs)),
             X = names(coef(coxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(coxph_fit_d2), coef(coxph_robust_fit_d2),
                          coef(svycoxph_fit_d2), coef(svycoxph_fit_d3),
                          coef(coxme_fit_d2), coef(coxme_fit_d3)),
             rbind(confint(coxph_fit_d2),
                   confint(coxph_robust_fit_d2),
                   confint(svycoxph_fit_d2),
                   confint(svycoxph_fit_d3),
                   confint(coxme_fit_d2),
                   confint(coxme_fit_d3))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..)

}

one_rep()

reps <- replicate(1000, one_rep(), simplify = FALSE)

reps_df <- Reduce(rbind.data.frame, reps)

prop.table(xtabs(~hit + model + X, data = reps_df), margin = c(2, 3))




# stat ignores the fact that in repeated measures, only the last period of observation can be censored (right?).
# I'm dumping the ranom effects info, but could extract it if needed (it's returned as an attribute by one_dataset)

# sample 50 clusters



one_dataset(~X1 + X2 + X3 + stratum + (1 | M),
            dists = list(X1 = ~rnorm(n),
                         X2 = ~rnorm(k * nk),
                         X3 = ~rbinom(n, 1, 0.5),
                         M = ~rep(1:k, each = nk),
                         error = ~rexp(n, 10),
                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                         stratum = ~rep(c(0, 1), c(floor(2/3 * k) * nk, ceiling(1/3 * k) * nk))),
            dist_args = list(k = 50, nk = 4, n = 200),
            coefficients = c(1, 1, 1, 0.5),
            random_effect_variance = list(M = 1)
)


k = 50
nk = 7

k * nk

floor(k/3)










