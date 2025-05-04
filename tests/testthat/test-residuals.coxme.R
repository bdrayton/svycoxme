# test-residuals.coxme

# Compare the output from
# residuals.coxme and residuals.coxph
# when the coxme.object has been doctored
# so that the random effects are 0 and
# the coefficients are equal to those from coxph.
# The residuals should be equal, down to negligible decimal differences.

# test residuals against coxph
# if the linear predictors and variance are the same, the residuals should be the
# same

test_that("residuals are the same: right-censored time",{

  des <- svydesign(ids = ~group_id, weights = ~1, data = samp_srcs, fpc = ~fpc)

  fit_svycoxph <- svycoxph(Surv(stat_time, stat) ~ X1,
                           design = des)

  fit_svycoxme <- svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id),
                           design = des)

  # cook the fit
  fit_svycoxme$linear.predictor <- fit_svycoxph$linear.predictors

  fit_svycoxme$var <- vcov(fit_svycoxph)

  dfbeta_svycoxph <- resid(fit_svycoxph, type = 'dfbeta', weighted = TRUE) |> as.matrix()
  dfbeta_svycoxme <- resid(fit_svycoxme, data = des$variables, type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxph) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxph, dfbeta_svycoxme)

})


# check with tied event times
test_that("residuals are the same: right-censored time; ties",{

  dset2 <- samp_srcs

  dset2$stat_time <- round(dset2$stat_time, 3)

  des <- svydesign(ids = ~group_id, weights = ~1, data = dset2, fpc = ~fpc)

  fit_svycoxph <- svycoxph(Surv(stat_time, stat) ~ X1,
                           design = des)

  fit_svycoxme <- svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id),
                           design = des)

  # cook the fit
  fit_svycoxme$linear.predictor <- fit_svycoxph$linear.predictors

  fit_svycoxme$var <- vcov(fit_svycoxph)

  dfbeta_svycoxph <- resid(fit_svycoxph, type = 'dfbeta', weighted = TRUE) |> as.matrix()
  dfbeta_svycoxme <- resid(fit_svycoxme, data = des$variables, type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxph) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxph, dfbeta_svycoxme)

})

# check with counting process time.
test_that("residuals are the same: counting-process time",{

  dset2 = samp_srcs

  dset2$start_time = 0

  des <- svydesign(ids = ~group_id, weights = ~1, data = dset2, fpc = ~fpc)

  fit_svycoxph <- svycoxph(Surv(start_time, stat_time, stat) ~ X1,
                           design = des)

  fit_svycoxme <- svycoxme(Surv(start_time, stat_time, stat) ~ X1 + (1 | group_id),
                           design = des)

  # cook the fit
  fit_svycoxme$linear.predictor <- fit_svycoxph$linear.predictors

  fit_svycoxme$var <- vcov(fit_svycoxph)

  dfbeta_svycoxph <- resid(fit_svycoxph, type = 'dfbeta', weighted = TRUE) |> as.matrix()
  dfbeta_svycoxme <- resid(fit_svycoxme, data = des$variables, type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxph) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxph, dfbeta_svycoxme)

})


# check with tied event times
test_that("residuals are the same: counting process time; ties",{

  dset2 <- samp_srcs

  dset2$stat_time <- round(dset2$stat_time, 3) + 0.0001
  dset2$start_time <- 0

  des <- svydesign(ids = ~group_id, weights = ~1, data = dset2, fpc = ~fpc)

  fit_svycoxph <- svycoxph(Surv(stat_time, stat) ~ X1,
                           design = des)

  fit_svycoxme <- svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id),
                           design = des)

  # cook the fit
  fit_svycoxme$linear.predictor <- fit_svycoxph$linear.predictors

  fit_svycoxme$var <- vcov(fit_svycoxph)

  dfbeta_svycoxph <- resid(fit_svycoxph, type = 'dfbeta', weighted = TRUE) |> as.matrix()
  dfbeta_svycoxme <- resid(fit_svycoxme, data = des$variables, type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxph) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxph, dfbeta_svycoxme)

})


# compare right-censored and counting process time
test_that("residuals are the same: right-censored vs counting-process time",{

  dset2 = samp_srcs

  dset2$start_time = 0

  des <- svydesign(ids = ~group_id, weights = ~1, data = dset2, fpc = ~fpc)

  fit_svycoxme_rc <- svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id),
                              design = des)

  fit_svycoxme_cp <- svycoxme(Surv(start_time, stat_time, stat) ~ X1 + (1 | group_id),
                           design = des)

  dfbeta_svycoxme_rc <- resid(fit_svycoxme_rc, data = des$variables,
                              type = 'dfbeta', weighted = TRUE)

  dfbeta_svycoxme_cp <- resid(fit_svycoxme_cp, data = des$variables,
                              type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxme_rc) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme_cp) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxme_rc, dfbeta_svycoxme_cp)

})


# check with tied event times
test_that("residuals are the same: right-censored vs counting-process time; ties",{

  dset2 <- samp_srcs

  dset2$stat_time <- round(dset2$stat_time, 3) + 0.0001
  dset2$start_time = 0

  des <- svydesign(ids = ~group_id, weights = ~1, data = dset2, fpc = ~fpc)

  fit_svycoxme_rc <- svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id),
                              design = des)

  fit_svycoxme_cp <- svycoxme(Surv(start_time, stat_time, stat) ~ X1 + (1 | group_id),
                              design = des)

  dfbeta_svycoxme_rc <- resid(fit_svycoxme_rc, data = des$variables,
                              type = 'dfbeta', weighted = TRUE)

  dfbeta_svycoxme_cp <- resid(fit_svycoxme_cp, data = des$variables,
                              type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxme_rc) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme_cp) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxme_rc, dfbeta_svycoxme_cp)

})


test_that("residuals are the same: right-censored time; missing data",{

  dmiss = samp_srcs

  i_miss = sample(seq_along(dmiss$stat_time), floor(0.1 * nrow(dmiss)))

  dmiss$stat_time[i_miss] <- NA

  des <- svydesign(ids = ~group_id, weights = ~1, data = dmiss, fpc = ~fpc)

  fit_svycoxph <- svycoxph(Surv(stat_time, stat) ~ X1,
                           design = des)

  fit_svycoxme <- svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id),
                           design = des)

  # cook the fit
  fit_svycoxme$linear.predictor <- fit_svycoxph$linear.predictors

  fit_svycoxme$var <- vcov(fit_svycoxph)

  dfbeta_svycoxph <- resid(fit_svycoxph, type = 'dfbeta', weighted = TRUE) |> as.matrix()
  dfbeta_svycoxme <- resid(fit_svycoxme, data = des$variables, type = 'dfbeta', weighted = TRUE)

  dimnames(dfbeta_svycoxph) <- list(NULL, NULL)
  dimnames(dfbeta_svycoxme) <- list(NULL, NULL)

  expect_equal(dfbeta_svycoxph, dfbeta_svycoxme)

})














