
# more tests for score residuals.

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
ggplot(newcgd, aes(x = tstart, xend = tstop, y = y, yend = y, shape = factor(infect))) + geom_segment() +
  geom_point(aes(x = tstop, y = y))

fit <- coxph(Surv(tstart, tstop, infect) ~ treat + inherit + steroids,
             data = newcgd_sorted, cluster = id, ties = 'breslow')

debugonce(make_parts.coxph)
my_parts <- make_parts(fit, data = newcgd_sorted, weights = rep(1, nrow(newcgd)))

debugonce(calc_ui.coxph_parts)
my_score <- calc_ui(my_parts)

score <- resid(fit, type = "score")

all.equal(as.matrix(my_score), as.matrix(score))

my_score |> colSums()
score |> colSums()


