
# population

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
                          random_effect_variance = list(M = 1)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))

})

pop <- Reduce(rbind.data.frame, pop_list)

# I expect the random effects to have mean = 0 and var = 2.
pop |> dplyr::distinct(id, .keep_all = TRUE) |>
  dplyr::summarise(mean(M_re),
                   var(M_re))

# a random sample from the population

sample <- dplyr::slice_sample(pop, n = 1000)

# calculate a outcome, y, that is linearly dependent on X, plus a specified amount of error

lm_data <- sample |> dplyr::mutate(error = rnorm(dplyr::n(), 0, sqrt(1.5)),
                                   y = X1 * 1 + X2 * -0.7 + X3 * 0.5 + error,
                                   X3 = as.factor(X3))

lm_fit <- lm(y ~ X1 + X2 + X3, data = lm_data)

summary(lm_fit)
plot(lm_fit, which = 1)

lm_res <- residuals(lm_fit)

mean(lm_res)
sum(lm_res)
plot(lm_res)
hist(lm_res)
var(lm_res)

plot(lm_data$error, lm_res)
abline(0, 1, col = "green")

res_fit <- lm(lm_data$error ~ lm_res)
summary(res_fit)
abline(coef(res_fit), col = "red")

# differences between error and residual

error_res_diff <- lm_data$error - lm_res

mean(error_res_diff)
hist(error_res_diff)

# intercept seems to be hard to fit correctly. Should be 0, but only get close with lots of data (n = ~10000)

# what do residuals behave with lme?
# added M_re into the linear term
lme_data <- sample |> dplyr::mutate(error = rnorm(dplyr::n(), 0, sqrt(1.5)),
                                   y = X1 * 1 + X2 * -0.7 + X3 * 0.5 + M_re + error,
                                   X3 = as.factor(X3))
library(lme4)
lme_fit <- lmer(y ~ X1 + X2 + X3 + (1 | M), data = lme_data)

summary(lme_fit)
plot(lme_fit, which = 1)

lme_res <- residuals(lme_fit)

mean(lme_res)
sum(lme_res)
plot(lme_res)
hist(lme_res)
var(lme_res)

plot(lme_data$error, lme_res)
abline(0, 1, col = "green")

res_fit <- lm(lme_data$error ~ lme_res)
summary(res_fit)
abline(coef(res_fit), col = "red")

# differences between error and residual

error_res_diff <- lm_data$error - lm_res

mean(error_res_diff)
hist(error_res_diff)













