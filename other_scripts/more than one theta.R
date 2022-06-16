

# to extend beyond the shared frailty model, code must allow length(theta) > 1,
# and also various random effects structures. I'll be seeking inspiration from lme4 for this

# some scenarios
# nested random effects
# non-nested random effects
#

# what happens when you give the formula parser in lme4 a formula with a Surv left hand side?

library(survival)
library(lme4)
library(lme4pureR)

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, subset=Days>=2))

# step by step

parsed_formula <- lFormula(Reaction ~ Days + (Days|Subject), sleepstudy, subset = Days >= 2)

# contains:
#  model.frame
#  model.frame for X variables (intercept and Days)
#  list with
#    transposed Zt. blocks have one row for each Z term. one column per observation. one block per subject.
#    theta, c(1, 0, 1)
#    Lind - c(1, 2, 3) repeated 18 times
#    Gp 0 36 dont know what this is.
#    lower - presumably lower bounds.
#    Lambdat 36 * 36 sparse Matrix. block of 2 * 2. they're upper triangle blocks. the lower triangle is sparse, diag is 1, 1 and the upper off-diagonal is 0.
#    flist - list with one term per random effect. factor with the random effect groups. this has an attribute called "assign" with value 1.
#    cnms - list with the names of X variables for each random effect. maybe for slopes?
#    Ztlist - list of matrices. the first (an only in this case) is called "Days | Subject and is 36 * 144. Same structure as Zt. Zt would be an rbind of Ztlist i think.
#  REML - logical flag
#  formula - the call formula
#  wmsgs - count of warning messages.


# Look at another dataset
ds <- one_dataset(~ X1 + (1 | M1) + (1| M2),
                  dists = list(X1 = ~rnorm(n),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 10, nk = 4, n = 40),
                  coefficients = c(X1 = 1),
                  random_effect_variance = list(M1 = 1, M2 = 2)
)

lFormula(stat_time ~ X1 + (X1 | M1), ds)
lFormula(stat_time ~ X1 + (1 | M1), ds)
lFormula(stat_time ~ X1 + (1 | M1) + (1 | M2), ds)
lFormula(stat_time ~ X1 + (1 | M1) + (1 | M1:M2), ds)

lFormula(stat_time ~ X1 + (X1 || M1), ds)




## ----phiToTheta----------------------------------------------------------
phiToTheta <- function(phi) {
  theta5 <- -(phi[2]*phi[3])/phi[4]
  c(phi[1:4], theta5, phi[5])
}

## ----compute deviance function modular-----------------------------------
# lf <- lFormula(stat_time ~ X1 + (1 | M1), ds)
lf <- lFormula(Reaction ~ Days + (1 | Subject), data = subset(sleepstudy, Days >= 2))
devf <- do.call(mkLmerDevfun, lf)

## ----wrapper modular-----------------------------------------------------
devfWrap <- function(phi) devf(phiToTheta(phi))

## ----opt modular---------------------------------------------------------
opt <- minqa::bobyqa(par = lf$reTrms$theta[-5],
              fn = devfWrap,
              lower = lf$reTrms$lower[-5])

devfWrap(lf$reTrms$theta)


## optimisation in R
## minqa::bobyqa gives the same value as REMLcrit.
## but what does par mean? It's theta: random-effects parameter estimates: these are parameterized as the relative Cholesky factors of each random effect term

lf <- lFormula(Reaction ~ Days + (1 | Subject), data = subset(sleepstudy, Days >= 2))

devfun <- pls(lf$X, lf$fr$Reaction, lf$reTrms$Zt, lf$reTrms$Lambdat, thfun = function(theta) theta[lf$reTrms$Lind])

opt <- minqa::bobyqa(par = lf$reTrms$theta,
              fn = devfun,
              lower = lf$reTrms$lower)

fit <- lmer(Reaction ~ Days + (1 | Subject), data = subset(sleepstudy, Days >= 2))
REMLcrit(fit)

getME(fit, name = "theta")

## lets look at other random effect structures

test_form <- stat_time ~ X1 + (X1 | M1)

ds <- one_dataset(test_form,
                  dists = list(X1 = ~rnorm(n),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 10, nk = 4, n = 40),
                  coefficients = c(X1 = 1),
                  random_effect_variance = list(M1 = 1, M2 = 2)
)

lf <- lFormula(test_form, data = ds)
debugonce(pls)
devfun <- pls(lf$X, lf$fr[,1, drop = TRUE], lf$reTrms$Zt, lf$reTrms$Lambdat, thfun = function(theta) theta[lf$reTrms$Lind])

devenv <- environment(devfun)
devenv_parent <- parent.env(devenv)

parent.env(devenv_parent)

child.env

get("Whalf", devenv_parent)

lapply(ls(envir = parent.env(devenv)),
       function(object_name) {
         get(object_name, parent.env(devenv)) |> class()
       })

ls(envir = devenv)

get("u", devenv)






