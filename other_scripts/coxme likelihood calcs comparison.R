#'  This document aims to illuminate where my method diverges from coxme.
#'
#'  We know that my likelihood function, gradient function and hessian maxtix function
#'  return the same values as coxme when evaluated at the coxme estimates of b and beta.
#'
#'  Theta is a different story. My theta estimate is different (worse).
#'
#'  here I want to examine some of the intermediate steps is coxme estimation method.
#'
#'  Briefly, coxme fits a series of models. First, coxph is used to get starting values
#'  for the fixed effects. and random effects start at 0.
#'
#'  next, these starting values are passed to optim, along with a likilood function for theta.
#'  This theta function finds the mles for b and beta and then returns a likelihood for the
#'  theta likelihood function, evaluated as these maximised b and beta.
#'
#'  This theta likelihood seems to be where our approaches diverge, so I will be looking at
#'  the code and data around that point for clues into what is happening.
#'


my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 10
my_nk = 10

ds <- one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + X3 + (1|M), data = ds)

coxme_est_theta <- coxme::VarCorr(fit)$M

coxme_est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

parm_names <- paste0(rep(c("X", "Z"), c(3, my_k)), c(1:3, seq_len(my_k)))

names(coxme_est_parms) <- parm_names

my_loglik <- lp(parms = coxme_est_parms,
                X = c("X1", "X2", "X3"),
                cluster = "M",
                t = t, dij = stat,
                theta = coxme_est_theta,
                data = ds)

# debugonce(coxme::coxme)
debugonce(coxme:::coxme.fit)

# this is the key part:

# .Call(Ccoxfit6b, as.integer(c(iter, iter)), as.double(init), ikmat@blocks, ikmat@rmat)
# ilik <- fit$loglik[2] - 0.5 * (sum(log(diag(gkmat))) +
#                                  fit$hdet)

#' it's in the "logfun" function for theta, which is passed to optim. So the questions are:
#' what is fit$loglik[2] , and the other parts? do I have functions that return the same numbers?
#'
#' what is fit0? why the -(1 + ilik - fit0) step?
#'
#' it also gets called for the potential start values for theta.
#'
#' fit$loglik[2] - this I think is the loglikelihood evaluated at maximums of beta and b
#'
#' sum(log(diag(gkmat))) if the trace of the logged penalty matrix.
#'
#' fit$hdet log determinant of the sparse portion of H
#'
#' fit$loglik
#'
#'
#'


# fit0 = -275.091572412296


# $theta
# -7.82404601086
# 4e-04
#
# $beta
# [1]  0.0027984241 -0.0006796661  0.0027857531  0.0010365816 -0.0029880198  0.0011195241  0.0005011004 -0.0025541616 -0.0030972577  0.0010777219  0.5786779052 -0.2329656985
# [13] -0.0298623774
#
# $loglik
# [1] -276.6478 -275.0354
#
# $hdet
# [1] 78.26614
#
# $iter
# [1] 5


my_parms =  c(0.5786779052, -0.2329656985,
              -0.0298623774, 0.0027984241, -0.0006796661,  0.0027857531,  0.0010365816,
              -0.0029880198,  0.0011195241,  0.0005011004, -0.0025541616,
              -0.0030972577, 0.0010777219)

lp(parms = my_parms,
   X = c("X1", "X2", "X3"),
   cluster = "M",
   t = t, dij = stat,
   theta = exp(-7.82404601086),
   data = ds)


################

# theta  -0.446287102628

# $beta
# [1]  1.14170937 -0.41150386  1.38062518  0.10912330 -1.17788400  0.77794480 -0.07114515 -0.91425244 -1.27330540  0.43868820  0.70600728 -0.36854523  0.04664341
#
# $loglik
# [1] -276.6478 -251.7455
#
# $hdet
# [1] 19.70085
#
# $iter
# [1] 5

# with these betas and thetas, can i recover the loglik?

my_parms =  c(0.70600728, -0.36854523,  0.04664341, 1.14170937, -0.41150386,
              1.38062518,  0.10912330, -1.17788400,  0.77794480, -0.07114515,
              -0.91425244, -1.27330540,  0.43868820  )

lp(parms = my_parms,
   X = c("X1", "X2", "X3"),
   cluster = "M",
   t = t, dij = stat,
   theta = exp(-0.446287102628),
   data = ds)

# it's not the same. perhaps there is some scaling?








