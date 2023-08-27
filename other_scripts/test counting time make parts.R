# test make_parts on counting process time survival models

# should work with time varying covariates and counting process time.

# just start with counting process time.

library(survival)

# Create a simple data set for a time-dependent model
test2 <- list(
  # start=c(1,2,5,2,1,7,3,4,8,8),
    start=c(0,0,0,0,0,0,0,0,0,0),
              stop=c(2,3,6,7,8,9,9,9,14,17),
              event=c(1,1,1,1,1,1,1,0,0,0),
              x=c(1,0,0,1,0,1,1,1,0,0),
              x2 = rnorm(10)) %>%
  as.data.frame()

fit <- coxph(Surv(start, stop, event) ~ x + x2, test2, x = TRUE)

with(test2, Surv(stop, event)) |> str()

library(ggplot2)

test2$y <- rev(seq(nrow(test2)))

ggplot(test2, aes(x = start, xend = stop, y = y, yend = y)) + geom_segment()

# reordering
test2 <- test2[order(test2$stop, test2$start), ]
test2$y <- rev(seq(nrow(test2)))
ggplot(test2, aes(x = start, xend = stop, y = y, yend = y)) + geom_segment()

n = nrow(test2)

start_test <- test2$start <= rep(test2$start, each = n)
stop_test <- test2$stop >= rep(test2$stop, each = n)

in_risk_set_matrix <- matrix(start_test & stop_test, nrow = n)

exp_risk_score = predict(fit, type = "risk")

X <- fit$x

fast_risk_sets(matrix(exp_risk_score))

n = nrow(response)

colSums(in_risk_set_matrix * exp_risk_score)

# this is S1_hat, a n * p matrix
exp_risk_score_X <- exp_risk_score * X

# also an n * p matrix
# at_risk_X <- fast_risk_sets(exp_risk_score_X)

apply(exp_risk_score_X, 2, function(X_j){

  colSums(in_risk_set_matrix * X_j)

})

XtX_i <- matrix(apply(X, 1, tcrossprod), ncol = nrow(X))

exp_risk_score_XtX <- t(XtX_i) * exp_risk_score

# I think I can give this to fast_risk_sets.
# at_risk_XtX <- fast_risk_sets(exp_risk_score_XtX)

apply(exp_risk_score_XtX, 2, function(X_j){

  colSums(in_risk_set_matrix * X_j)

})


#### compare parts when using counting or right time.

# debugonce(make_parts.coxph)
parts1 <- make_parts(fit, data = test2, weights = rep(1, n))

test2 <- list(
  # start=c(1,2,5,2,1,7,3,4,8,8),
  start=c(0,0,0,0,0,0,0,0,0,0),
  stop=c(2,3,6,7,8,9,9,9,14,17),
  event=c(1,1,1,1,1,1,1,0,0,0),
  x=c(1,0,0,1,0,1,1,1,0,0),
  x2 = rnorm(10)) %>%
  as.data.frame()

fit <- coxph(Surv(stop, event) ~ x + x2, test2, x = TRUE)

parts2 <- make_parts(fit, data = test2, weights = rep(1, n))

all.equal(parts1, parts2)


# compare residuals with therneau residuals with each time type.

# right time (effectively)
test2 <- list(
  start=c(0,0,0,0,0,0,0,0,0,0),
  # start=c(1,2,5,2,1,7,3,4,8,8),
  stop=c(2,3,6,7,8,9,9,9,14,17),
  event=c(1,1,1,1,1,1,1,0,0,0),
  x=c(1,0,0,1,0,1,1,1,0,0),
  x2 = rnorm(10)) %>%
  as.data.frame()

fit <- coxph(Surv(start, stop, event) ~ x + x2, test2, x = TRUE, ties = 'breslow')

parts3 <- make_parts(fit, data = test2, weights = rep(1, n))

score_mine <- calc_ui(parts3)
score_therneau <- resid(fit, type = "score")

all(score_mine - score_therneau < .Machine$double.neg.eps*10)

### when using counting time:

test2 <- list(
  # start=c(0,0,0,0,0,0,0,0,0,0),
  start=c(1,2,5,2,1,7,3,4,8,8),
  stop=c(2,3,6,7,8,9,9,9,14,17),
  event=c(1,1,1,1,1,1,1,0,0,0),
  x=c(1,0,0,1,0,1,1,1,0,0),
  x2 = rnorm(10)) %>%
  as.data.frame()

fit <- coxph(Surv(start, stop, event) ~ x + x2, test2, x = TRUE, ties = 'breslow')

parts3 <- make_parts(fit, data = test2, weights = rep(1, n))

score_mine <- calc_ui(parts3)
score_therneau <- resid(fit, type = "score")

all(score_mine - score_therneau < .Machine$double.neg.eps*10)

score_mine - score_therneau










