
# test for inadvertent changes to the output from print, coef, and summary.

test_that('compare summary(fit1) to saved value', {
  des <- svydesign(ids = ~group_id, weights = ~weight, data = samp_srcs)

  fit1 = svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id), design = des)

  expect_equal(summary(fit1), summary_fit1)
})


test_that("compare coef(fit) to saved value", {
  des <- svydesign(ids = ~group_id, weights = ~weight, data = samp_srcs)

  fit1 = svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id), design = des)

  expect_equal(coef(fit), coef_fit1)
})


test_that("compare print(fit) to saved value", {
  des <- svydesign(ids = ~group_id, weights = ~weight, data = samp_srcs)

  fit1 = svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id), design = des)

  expect_equal(print(fit), print_fit1)
})

# test the use of futures

# test residuals against coxph




