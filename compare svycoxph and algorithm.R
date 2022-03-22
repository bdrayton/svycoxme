# compare my algorithm to svycoxph
# for a scenario where random effects are 0, so they should be the same.

library(devtools)

document()
load_all()

theta = 0.5
number_of_clusters = 200
fixed_effects = c(1, -0.7, 0.5)

set.seed(203948)


myPop <- make_population(theta = theta, N = number_of_clusters,
                         fixed_effects = fixed_effects,
                         prCensoring = 0.2)

true_random_effects <- attr(myPop, "b")

n_clusters = 20
n_sample_per_cluster = 50

myNi <- myPop |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::pull(n)

# some clusters are smaller than ni if any of these cluster are chosen, the sample size
# will be less than n * ni. It also means sample sizes vary a bit.
# This behavior is inherited from Wang 2019

mySample <- sample_clusters(myPop, cluster, n_clusters, n_sample_per_cluster, z_columns = glue::glue("Z{1:number_of_clusters}"))

x_columns = paste0("X", 1:3)
z_columns = paste0("Z", sort(unique(mySample$cluster)))

D_theta = diag(n_clusters) * theta

rownames(D_theta) <- z_columns
colnames(D_theta) <- z_columns

cox_start <- coef(survival::coxph(survival::Surv(t, d) ~ X1 + X2 + X3, data = mySample))

# debugonce(estimate_all_parameters)

test_results <- estimate_all_parameters(
  fixed_effects_start = cox_start,
  random_effects_start = rep(0, n_clusters),
  theta_start = 0.2,
  x_columns = x_columns,
  z_columns = z_columns,
  t = t,
  i = cluster,
  wi = wi,
  wji = wji,
  dij = d,
  data = mySample,
  max_iterations = 1000,
  eps = 0.0001)


library(ggplot2)

lapply(test_results, "[[", "fixed_effects") |>
  Reduce(f = rbind) |>
  as.data.frame() |>
  dplyr::mutate(iteration = dplyr::row_number() - 1, .before = everything()) |>
  tidyr::pivot_longer(cols = glue::glue("X{1:3}")) |>
  ggplot(aes(iteration, value, colour = name))  + geom_point() + geom_line()


lapply(test_results, "[[", "theta") |>
  Reduce(f = rbind) |>
  as.data.frame() |>
  dplyr::mutate(iteration = dplyr::row_number() - 1, .before = everything()) |>
  ggplot(aes(iteration, V1))  + geom_point() + geom_line()

lapply(test_results, "[[", "random_effects") |>
  Reduce(f = rbind) |>
  as.data.frame() |>
  dplyr::mutate(iteration = dplyr::row_number() - 1, .before = everything()) |>
  tidyr::pivot_longer(cols = glue::glue("V{1:n_clusters}")) |>
  ggplot(aes(iteration, value, colour = name))  + geom_point() + geom_line()

eps = c()
for (i in seq_along(test_results)[-1]) {
  eps = c(eps, max(abs(unlist(test_results[i-1]) - unlist(test_results[i]))))
}

plot(eps)

library(survey)

des <- svydesign(ids = ~cluster, weights = ~wi, strata = NULL, data = mySample)

cbind(
  true = fixed_effects,
  survey = coef(svycoxph(survival::Surv(t, d) ~ X1 + X2 + X3, design = des)),
  survivial = cox_start,
  svycoxme = tail(test_results,1)[[1]][["fixed_effects"]])







