# evaluate the impact of changing numerical gradient step size on
# hessian stability.

# see results in:
# C:\Users\Bradley\OneDrive - The University of Auckland\PhD\outputs\scratchpads and notes\theta estimation\optimisation struggles.docx


my_formula <- survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1| M1:M2)

ds <- one_dataset(my_formula,
                  dists = list(X1 = ~rnorm(n),
                               X2 = ~rnorm(k * nk),
                               X3 = ~rbinom(n, 1, 0.5),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  coefficients = c(1, -0.7, 0.5),
                  random_effect_variance = c(M1 = 1, `M1:M2` = 1.5)
)

ds$`M1:M2` <- with(ds, interaction(M1, M2))

coxfit <- coxme::coxme(survival::Surv(stat_time, stat)~X1 + X2 + X3 + (1 | M1) + (1| `M1:M2`), data = ds)

coxme::VarCorr(coxfit)

ests_list <- list()

# debugonce(est_parameters)
ests_list[[1]] <- est_parameters(my_formula, data = ds,
                                 control = control.list(factr = 1e3,
                                                        reltol = 1e-13,
                                                        ndeps = c(1e-2, 1e-2)))

for(i in 1:10) {

  ests_list[[i + 1]] <- est_parameters(my_formula, data = ds,
                                       start_params = c(ests_list[[i]]$beta, ests_list[[i]]$b),
                                       theta_start = ests_list[[i]]$theta,
                                       control = control.list(factr = 1e3,
                                                              reltol = 1e-13,
                                                              ndeps = c(1e-2, 1e-2)))

}

# saveRDS(object = ests_list, file = "C:/Users/Bradley/OneDrive - The University of Auckland/PhD/outputs/scratchpads and notes/theta estimation/ests_ndeps=1e-6.rds")


# extract theta
library(tidyverse)

# lapply(ests_list, names)
thetas_df <- sapply(ests_list, "[[", "theta", simplify = "matrix") |>
  data.frame() |>
  mutate(theta = factor(c("theta M1", "theta M2")), .before = everything()) |>
  pivot_longer(cols = X1:X11) |>
  mutate(Iteration = factor(name, levels = paste0("X", 1:11), labels = as.character(1:11)))

thetas_df %>%
  group_by(theta) %>%
  summarise(min(value) - max(value))

#zigzagging?
ggplot(thetas_df, aes(Iteration, value)) + geom_point() +
  facet_grid(rows = vars(theta), scales = "free")


# same with betas
betas_df <- sapply(ests_list, "[[", "beta", simplify = "matrix") |>
  data.frame() |>
  mutate(beta = factor(paste("beta", 1:3)), .before = everything()) |>
  pivot_longer(cols = X1:X11) |>
  mutate(Iteration = factor(name, levels = paste0("X", 1:11), labels = as.character(1:11)))

#zigzagging?
ggplot(betas_df, aes(Iteration, value)) + geom_point() +
  facet_grid(rows = vars(beta), scales = "free")

betas_df %>%
  group_by(beta) %>%
  summarise(min(value) - max(value))


# what's happening with the hessians now?

theta_var_ests <- lapply(ests_list, function(x){
  sqrt(diag(-solve(attr(x, "theta_est")$hessian)))
})

theta_var_df <- Reduce(rbind, theta_var_ests) |>
  as.data.frame() |>
  pivot_longer(cols = everything())

theta_var_df |> ggplot(aes(value)) + geom_dotplot() + facet_wrap(vars(name), scales = "free")

theta_var_df |>
  group_by(name) |>
  summarise(diff(range(value)))





