
# more tests for score residuals.
library(survival)
cgd0 <- cgd0

newcgd <- tmerge(data1=cgd0[, 1:13], data2=cgd0, id=id, tstop=futime)
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime1))

newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime2))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime3))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime4))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime5))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime6))
newcgd <- tmerge(newcgd, cgd0, id=id, infect = event(etime7))
newcgd <- tmerge(newcgd, newcgd, id, enum=cumtdc(tstart))

newcgd_sorted <- newcgd[order(newcgd$tstop, newcgd$tstart), ]

newcgd_sorted$y = nrow(newcgd):1

newcgd$y = nrow(newcgd):1

library(ggplot2)
ggplot(newcgd_sorted, aes(x = tstart, xend = tstop, y = y, yend = y, shape = factor(infect))) + geom_segment() +
  geom_point(aes(x = tstop, y = y))

fit <- coxph(Surv(tstart, tstop, infect) ~ treat + inherit + steroids,
             data = newcgd_sorted, cluster = id, ties = 'breslow')

# debugonce(make_parts.coxph)
my_parts <- make_parts(fit, data = newcgd_sorted, weights = rep(1, nrow(newcgd)))

# debugonce(calc_ui.coxph_parts)
my_score <- calc_ui(my_parts)

score <- resid(fit, type = "score")

all.equal(as.matrix(my_score), as.matrix(score))

my_score |> colSums()
score |> colSums()

##### coxme, svycoxme

#### counting time with coxme

my_des <- svydesign(~id, weights = ~1, data = newcgd_sorted)
my_des_jackknife <- as.svrepdesign(my_des, type = "JK1")

svycoxme_fit_jackknife <- svycoxme(Surv(tstart, tstop, infect) ~ treat + inherit + steroids + (1|id),
                                   des = my_des_jackknife)

fit <- coxme::coxme(Surv(tstart, tstop, infect) ~ treat + inherit + steroids + (1|id),
                    data = newcgd_sorted)

# debugonce(make_parts.coxme)
parts3 <- make_parts(fit, data = newcgd_sorted, weights = rep(1, nrow(newcgd_sorted)))
debugonce(calc_ui)
score_mine <- calc_ui(parts3)

debugonce(get_information)
vv <- get_information(fit)

uivv <- score_mine %*% vv

sandwich_var <- survey::svyrecvar(uivv, my_des$cluster, my_des$strata,
                  my_des$fpc, postStrata = my_des$postStrata)

diag(sandwich_var)[1:3]

vcov(svycoxme_fit_jackknife) |> diag()

myb <- unlist(coxme::ranef(fit), use.names = FALSE)

names(myb) <- 1:4

debugonce(ScoreAll)
ScoreAll(beta0  = coxme::fixef(fit),
         b0 = myb,
         mydata_sample1 = test2 ,
         theta0 = unlist(coxme::VarCorr(fit)))






