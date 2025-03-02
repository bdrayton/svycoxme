
library(survival)


n = 5000

X = data.frame(X1 = rnorm(n),
               X2 = rnorm(n),
               M = rep(c(0,1), each = n/2),
               id = seq_len(n))

dset = draw_event_times(formula = Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = X,
                        coefficients = c(X1 = 1, X2 = 0.5),
                        random_effect_variance = c(M = 0),
                        id = id,
                        baseline_hazard = 1,
                        event = "single")

dset = dplyr::left_join(dset, dplyr::select(X, id, M), by = dplyr::join_by(id))


coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r1 = residuals.coxme(coxmefit, data = dset)

# debugonce(residuals2.coxme)
r2 = residuals2.coxme(coxmefit, data = dset)

all.equal(r1, r2)

microbenchmark::microbenchmark(
  residuals.coxme(coxmefit, data = dset),
  residuals2.coxme(coxmefit, data = dset),
  times = 10
)


# try with different start times.

library(survival)

k = 200 # clusters
nk = 5 # cluster size

N = nk * k # total obs

# clusters sampled
n = 100

baserate = 0.5
beta = c(0.5, -0.25)
theta = 0
b = rnorm(k, sd = theta)

Z1 = rnorm(N) # obs level
Z2 = rep(rnorm(k), each = nk)  # cluster level
group_id = rep(1:k, each = nk)

Z1_mean = tapply(Z1, group_id, FUN = mean)
Z2_mean = tapply(Z2, group_id, FUN = mean)

mu_X1 = 0.5 * Z1

mu_X2 = 0.5 * (Z1_mean + Z2_mean)

X1 = rnorm(N, mean = mu_X1, sd = 0.25)
X2 = rep(rnorm(k, mean = mu_X2, sd = 0.5), each = nk)

rate = baserate * exp(cbind(X1, X2)%*% matrix(beta) + rep(b, each = nk))

event_time = rexp(N, rate)

dset = data.frame(event_time, stat = 1, X1, X2, Z1, Z2, group_id)

coxmefit = coxme::coxme(Surv(event_time, stat) ~ X1 + X2 + (1 | group_id),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r1 = residuals.coxme(coxmefit, data = dset)

# debugonce(residuals2.coxme)
r2 = residuals2.coxme(coxmefit, data = dset)

all.equal(r1, r2)

dset$start_time = 0

coxmefit = coxme::coxme(Surv(start_time, event_time, stat) ~ X1 + X2 + (1 | group_id),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r1 = residuals.coxme(coxmefit, data = dset)

# debugonce(residuals2.coxme)
r2 = residuals2.coxme(coxmefit, data = dset)

all.equal(r1, r2)

# counting time.
dset <- dset %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(t_stop = cumsum(event_time)) %>%
  dplyr::mutate(t_start = dplyr::lag(t_stop, default = 0), .before = t_stop) %>%
  dplyr::ungroup()

coxmefit = coxme::coxme(Surv(t_start, t_stop, stat) ~ X1 + X2 + (1 | group_id),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

r1 = residuals.coxme(coxmefit, data = dset)

# debugonce(residuals2.coxme)
r2 = residuals2.coxme(coxmefit, data = dset)

all.equal(r1, r2)




