

# compare residuals using agscore from coxph and mine.

# get a model fit.

library(survival)

k = 200 # clusters
nk = 5 # cluster size

N = nk * k # total obs

# clusters sampled
n = 100

baserate = 0.5
beta = c(0.5, -0.25)
theta = 0.5
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
dset$start_time = 0

dset <- dset %>%
  dplyr::group_by(group_id) %>%
  dplyr::mutate(t_stop = cumsum(event_time)) %>%
  dplyr::mutate(t_start = dplyr::lag(t_stop, default = 0), .before = t_stop) %>%
  dplyr::ungroup()

#
# coxfit = coxph(Surv(t_start, t_stop, stat) ~ X1 + X2, data = dset,
#                x = TRUE, y = TRUE, ties = "breslow")
#
# coxmefit = coxme(Surv(t_start, t_stop, stat) ~ X1 + X2 + (1 | group_id), data = dset,
#                  x = TRUE, y = TRUE, ties = "breslow")
#

coxfit = coxph(Surv(start_time, event_time, stat) ~ X1 + X2, data = dset,
               x = TRUE, y = TRUE, ties = "breslow")

coxmefit = coxme(Surv(event_time, stat) ~ X1 + X2 + (1 | group_id), data = dset,
                 x = TRUE, y = TRUE, ties = "breslow")



coef(coxfit)
coef(coxmefit)

coxfit$residuals


type = "score"
otype = "score"

  n <- length(coxfit$residuals)
  rr <- coxfit$residuals
  y <- coxfit$y
  x <- coxfit[["x"]]
  vv <- drop(coxfit$naive.var)
  if (is.null(vv))
    vv <- drop(coxfit$var)
  weights <- coxfit$weights
  if (is.null(weights))
    weights <- rep(1, n)
  strat <- coxfit$strata
  method <- coxfit$method

  Terms <- coxfit$terms
  strats <- attr(Terms, "specials")$strata
  ny <- ncol(y)
  status <- y[, ny, drop = TRUE]
  nvar <- ncol(x)
  # set up strata, order.
  if (is.null(strat)) {
      ord <- order(y[, ny - 1], -status)
      newstrat <- integer(n)
      istrat <- integer(n)
  } else {
    istrat <- as.integer(strat)
    ord <- order(istrat, y[, ny - 1], -status)
    newstrat <- c(diff(as.numeric(istrat[ord])) !=
                    0, 1)
  }
  newstrat[n] <- 1
  x <- x[ord, ]
  y <- y[ord, ]
  score <- exp(coxfit$linear.predictors)[ord]
  # score <- exp(coxmefit$linear.predictor)[ord]
  istrat <- istrat[ord]
  if (ny == 3) {
    if (is.null(strat))
      sort1 <- order(y[, 1])
    else sort1 <- order(istrat, y[, 1])
  }

  storage.mode(y) <- storage.mode(x) <- "double"
  storage.mode(newstrat) <- "integer"
  storage.mode(score) <- storage.mode(weights) <- "double"
  resid <- .Call(survival:::Cagscore3, y, x, istrat, score, weights[ord],
                 as.integer(method == "efron"), sort1 - 1L)

  resid <- .Call(survival:::Ccoxscore2, y, x, istrat, score, weights[ord],
                 as.integer(method == "efron"))

  if (nvar > 1) {
    rr <- matrix(0, n, nvar)
    rr[ord, ] <- resid
    dimnames(rr) <- list(names(coxfit$residuals), names(coxfit$coefficients))
  }

  my_tstart = y[, 1, drop = TRUE]
  my_tend = y[, 2, drop = TRUE]
  my_status = y[, 3, drop = TRUE]

 resid2 = svycoxme::agscore3(my_tstart, my_tend, my_status, covar = x, strata = istrat,
                     score = score, weights = weights[ord], sort1 = sort1 - 1L,
                     method = as.integer(method == "efron"))

 all.equal(resid, resid2)

 resid3 = residuals.coxme(coxmefit, data = dset)

 identical(resid, resid3)


 # fit a mixed effects model. modify to give results of the random effects model.
 coxmefit = coxme::coxme(Surv(start_time, event_time, status) ~ X1 + X2 + (1 | group_id),
                         data = dset, x = TRUE, y = TRUE, ties = "breslow")

 coxmefit$coefficients <- coef(coxfit)

 coxmefit$frail <- list(M = c(`0` = 0.0, `1` = 0.0))

 coxme::random.effects(coxmefit)

 coef(coxmefit)
 coef(coxfit)

 parts <- make_parts.coxme(coxmefit, dset)

 lapply(parts, is.matrix)

 rr <- with(parts,{
   C_calc_ui(time_start = time_start,
             time_stop = time_stop,
             stat = stat,
             weights = weights,
             exp_risk_score = exp_risk_score,
             S0 = S0,
             S1_X = S1_X,
             X = X,
             weighted = TRUE)
 })

 rr



########

method_1 = function() {

     type = "score"
     otype = "score"

     n <- length(coxfit$residuals)
     rr <- coxfit$residuals
     y <- coxfit$y
     x <- coxfit[["x"]]
     vv <- drop(coxfit$naive.var)
     if (is.null(vv))
       vv <- drop(coxfit$var)
     weights <- coxfit$weights
     if (is.null(weights))
       weights <- rep(1, n)
     strat <- coxfit$strata
     method <- coxfit$method

     Terms <- coxfit$terms
     strats <- attr(Terms, "specials")$strata
     ny <- ncol(y)
     status <- y[, ny, drop = TRUE]
     nvar <- ncol(x)
     # set up strata, order.
     if (is.null(strat)) {
       ord <- order(y[, ny - 1], -status)
       newstrat <- integer(n)
       istrat <- integer(n)
     } else {
       istrat <- as.integer(strat)
       ord <- order(istrat, y[, ny - 1], -status)
       newstrat <- c(diff(as.numeric(istrat[ord])) !=
                       0, 1)
     }
     newstrat[n] <- 1
     x <- x[ord, ]
     y <- y[ord, ]
     score <- exp(coxfit$linear.predictors)[ord]
     istrat <- istrat[ord]
     if (ny == 3) {
       if (is.null(strat))
         sort1 <- order(y[, 1])
       else sort1 <- order(istrat, y[, 1])
     }

     storage.mode(y) <- storage.mode(x) <- "double"
     storage.mode(newstrat) <- "integer"
     storage.mode(score) <- storage.mode(weights) <- "double"

     resid = svycoxme::agscore3(my_tstart, my_tend, my_status, covar = x, strata = istrat,
                                 score = score, weights = weights[ord], sort1 = sort1 - 1L,
                                 method = as.integer(method == "efron"))

     resid

}

method_2 = function() {

     parts <- make_parts.coxme(coxmefit, dset)

     lapply(parts, is.matrix)

     rr <- with(parts,{
       C_calc_ui(time_start = time_start,
                 time_stop = time_stop,
                 stat = stat,
                 weights = weights,
                 exp_risk_score = exp_risk_score,
                 S0 = S0,
                 S1_X = S1_X,
                 X = X,
                 weighted = TRUE)
     })

     rr

 }


# timing test.

n = 10000

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

  coxfit = coxph(Surv(start_time, stop_time, status) ~ X1 + X2, data = dset, x = TRUE, y = TRUE, ties = "breslow")

  # fit a mixed effects model. modify to give results of the random effects model.
  coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                          data = dset, x = TRUE, y = TRUE, ties = "breslow")

  coxmefit$coefficients <- coef(coxfit)

  coxmefit$frail <- list(M = c(`0` = 0.0, `1` = 0.0))


all.equal(method_1(), method_2() )




microbenchmark::microbenchmark(method_1(), method_2(),
                               setup =  expression({
  dset = draw_event_times(formula = Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                          data = X,
                          coefficients = c(X1 = 1, X2 = 0.5),
                          random_effect_variance = c(M = 0),
                          id = id,
                          baseline_hazard = 1,
                          event = "single")

  dset = dplyr::left_join(dset, dplyr::select(X, id, M), by = dplyr::join_by(id))

  coxfit = coxph(Surv(start_time, stop_time, status) ~ X1 + X2, data = dset, x = TRUE, y = TRUE, ties = "breslow")

  # fit a mixed effects model. modify to give results of the random effects model.
  coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M),
                          data = dset, x = TRUE, y = TRUE, ties = "breslow")

  coxmefit$coefficients <- coef(coxfit)

  coxmefit$frail <- list(M = c(`0` = 0.0, `1` = 0.0))
}))



n = 1000

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

dset$strat = rep(c(1:2), n/2)


# can coxme have strata
coxmefit = coxme::coxme(Surv(start_time, stop_time, status) ~ X1 + X2 + (1 | M) + strata(strat),
                        data = dset, x = TRUE, y = TRUE, ties = "breslow")

coxmefit$strata


