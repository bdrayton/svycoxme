

# to extend beyond the shared frailty model, code must allow length(theta) > 1,
# and also various random effects structures. I'll be seeking inspiration from lme4 for this

# some scenarios
# nested random effects
# non-nested random effects
#


# what happens when you give the formula parser in lme4 a formula with a Surv left hand side?

library(survival)

# modify functions as needed.
devtools::install("C:\\Users\\Bradley\\Documents\\PhD_local\\lme4pureR")
library(lme4pureR)

library(lme4)
library(minqa)

# using the sleep study data, examine what the functions do

n <- nrow(sleepstudy)
x <- sleepstudy$Days
z <- rnorm(n)
X <- cbind(1, x)
ZZ <- cbind(1, z)
grp <- as.factor(sleepstudy$Subject)

RE <- mkRanefStructures(list(grp), list(ZZ))

Z <- t(RE$Zt)

head(Z)

fixed_effects <- rnorm(ncol(X))
random_effects <- rnorm(ncol(Z))

# y <- as.numeric(X%*%fixed_effects + Z%*%random_effects + rnorm(n))
y <- sleepstudy$Reaction
m <- lmer.fit(y,X,ZZ,grp)
m$par
Lambdat <- RE$Lambdat
Lambdat

summary(lmer(y ~ x + (1|grp)))




(fm1 <- lmer(Reaction ~ Days + (Days | Subject), sleepstudy))

(fm1 <- lmer(Reaction ~ Days + (Days|Subject), sleepstudy, subset=Days>=2))
## independent model
(fm2 <- lmer(Reaction ~ Days + (1|Subject) + (0+Days|Subject), sleepstudy, subset=Days>=2))

# fit a coxme model to these data
# see help(sleepstudy) in lme4 for information about the data. Models here are taken from the examples.

sleep_study <- sleepstudy

# everyone "fails"
sleep_study$stat <- 1

sleep_study2 <- subset(sleep_study, Days >= 2)

coxme::coxme(Surv(Reaction, stat) ~ Days + (1|Subject), data = sleep_study2)

coxme::coxme(Surv(Reaction, stat) ~ Days + (Days|Subject), data = sleep_study2)

coxme::coxme(Surv(Reaction, stat) ~ Days + (1|Subject) + (0+Days|Subject), data = sleep_study2)

# can I pass Surv to lmer? I expect it to fail, but I want to know where it fails.

lmer(survival::Surv(Reaction, stat) ~ Days + (1|Subject), data = sleep_study2)


parsed_formula <- lmer_formula_parsing(Surv(Reaction, stat) ~ Days + (1|Subject), data = sleep_study2)
str(parsed_formula)

parsed_formula <- lmer_formula_parsing(Surv(Reaction, stat) ~ Days + (Days|Subject), data = sleep_study2)
str(parsed_formula)

parsed_formula <- lmer_formula_parsing(Surv(Reaction, stat) ~ Days + (1|Subject) + (0+Days|Subject), data = sleep_study2)
str(parsed_formula)

# lmer control, and what my version will be is to be determined.

# The next part, where the deviance function is defined, is where lmer() fails.
# It fails because there are some operations on the Surv object that return a matrix, when a scalar is expected.

debugonce(lmer)
lmer(Reaction ~ Days + (1|Subject), data = sleep_study2)










