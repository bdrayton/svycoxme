
old_seed <- set.seed(9867)
n_samp = 100
k = 2000
nk = 10
n = nk * k
Z1 = rep(rep(1:5, c(700, 600, 400, 200, 100)), each = nk)
Z2 = rnorm(k)
Z3 = runif(k)
mu_X1 = 0.5 * (Z1 + 1)
mu_X2 = 0.5 * (Z2 + Z3)
mu_X3 = Z3
X1 = rnorm(n, mean = mu_X1, sd = 1)
X2 = rep(rnorm(k, mean = mu_X2, sd = 1), each = nk)
X3 = rep(rbinom(k, 1, mu_X3), each = nk)
X = cbind(X1, X2, X3)
beta = matrix(c(1, -0.7, 0.5))
theta = 0.5
b = rep(rnorm(k, sd = sqrt(theta)), each = nk)
baseline_hazard = 1.5
hazard = baseline_hazard * exp(X %*% beta + b)
event_time = - (log(runif(n)) / hazard)
censoring_time = runif(n, 0, 0.9)
stat = event_time <= censoring_time
stat_time = ifelse(stat, event_time, censoring_time)
group_id = rep(1:k, each = nk)
obs_id = seq(n)
pop = data.frame(X1, X2, X3, Z1, Z2 = rep(Z2, each = nk),
                 Z3 = rep(Z3, each = nk), stat_time, stat,
                 group_id, obs_id)

pop$sampled = pop$group_id %in% sample(unique(pop$group_id), size = n_samp)

usethis::use_data(pop, overwrite = TRUE, compress = "xz")

samp_srcs = pop[which(pop$sampled), , drop = FALSE]
samp_srcs$fpc = length(unique(group_id))
samp_srcs$weight = samp_srcs$fpc/n_samp

usethis::use_data(samp_srcs, overwrite = TRUE, compress = "xz")

set.seed(old_seed)
