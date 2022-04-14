



my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 50
my_nk = 10

max_iter = 100
convergence_threshold = 0.0001

my_X = c("X1", "X2", "X3")



sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

library(coxme)

debugonce(coxme)

coxme_fit <- coxme::coxme(survival::Surv(t, stat) ~ X1 + X2 + strata(X3) + (1 | M), data = sample_data)





formula <- survival::Surv(t, stat) ~ X1 + X2 + strata(X3) + (1 | M)

head(model.frame(formula, data = sample_data))

m <- model.frame(coxme:::subbar(formula), data = sample_data)

flist <- coxme:::formula1(formula)

special <- c("strata", "cluster")
Terms <- terms(flist$fixed, special)
attr(Terms, "intercept") <- 1
strats <- attr(Terms, "specials")$strata
cluster <- attr(Terms, "specials")$cluster
if (length(cluster)) {
  stop("A cluster() statement is invalid in coxme")
}

if (length(strats)) {
  temp <- untangle.specials(Terms, "strata", 1)
  if (length(temp$vars) == 1)
    strata.keep <- m[[temp$vars]]
  else strata.keep <- strata(m[, temp$vars], shortlabel = T)
  strats <- as.numeric(strata.keep)
  X <- model.matrix(Terms[-temp$terms], m)[, -1, drop = F]
}


coxme:::formula2(flist$random[[1]])





