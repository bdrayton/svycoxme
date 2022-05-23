#'
#' Purpose: Run modified version of coxme.fit to learn more about the estimation
#' of theta and beta and b.
#'
#' Inside coxme.fit, optim is called to maximise theta - the so called outer loop.
#' the function that is maximised finds optimum beta and b internally.
#'
#' I have modified coxme.fit so that it returns
#'
#' coxme.fit is run inside coxme.
#'
#' first I just want to get back what the loklig function returns as theta varies.
#'
#'
#'
#'

devtools::install("C:/Users/Bradley/OneDrive - The University of Auckland/PhD/code scratch pads/coxme",
                  build_vignettes = FALSE)

my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 10
my_nk = 10

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

# debugonce(coxme::coxme)

real_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = sample_data)

real_fit

fit <- coxme:::coxme_mod(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = sample_data)

mfit <- do.call('optim', c(list(par= fit$theta, fn=fit$fn, gr=NULL),
                           fit$control, fit$logpar))

theta <- mfit$par

exp(theta) == coxme::VarCorr(real_fit)


logliks <- sapply(seq(0.01, 1, by = 0.001), function(x){
  fit$fn(x,
         varlist = fit$logpar$varlist,
         vparm = fit$logpar$vparm,
         ntheta = fit$logpar$ntheta,
         ncoef = fit$logpar$ncoef,
         kfun = fit$logpar$kfun,
         init = fit$logpar$init,
         fit0 = fit$logpar$fit0,
         iter = fit$logpar$iter,
         timedep = fit$logpar$timedep)
})


plot(seq(0.01, 1, by = 0.001), logliks,)

seq(0.01, 1, by = 0.001)[min(logliks) == logliks] |> exp()


library(survival)

library(coxme)

fit <- coxme(Surv(t, stat) ~ X1 + X2 + X3 + (1 | M), data = sortAndIndex(sample_data, t))


d1 <- add_Z(sample_data, "M")

d2 <- sortAndIndex(d1, t)

d3 <- calcLinearPredictor(d2, my_X, attr(d1, "Z_names"), parms = c(fixef(fit), ranef(fit)$M))

d4 <- d3 |>
  dplyr::arrange(desc(index)) |>
  dplyr::mutate(cumsum_A = cumsum(A)) |>
  dplyr::arrange(index) |>
  dplyr::mutate(d_Hazard = stat/cumsum_A,
                cumsum_Hazard = cumsum(d_Hazard))

d4$cumsum_Hazard

cbind(d3$A, predict(fit, type = "risk"))





