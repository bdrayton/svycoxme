
data(cancer, package = "survival")

invisible(lung)

fit2 <- coxme::coxme(survival::Surv(time, status) ~ ph.ecog + age + (1|inst), lung)

resid(fit2, data = lung, type = "dfbeta")

