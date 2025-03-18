# generate objects for tests
# save for use in tests.

# When you run this script, which you'll need to do if
# you change the output of print, coef, or summary (on purpose)
# or samp_scrs changes, you'll need to check these objects to
# make sure they look right. The tests referencing these objects
# are there to catch accidental changes to the output.

data("samp_srcs")

des <- svydesign(ids = ~group_id, weights = ~weight, data = samp_srcs)

fit1 = svycoxme(Surv(stat_time, stat) ~ X1 + (1 | group_id), design = des)

coef_fit1 <- coef(fit1)
summary_fit1 <- summary(fit1)
print_fit1 <- print(fit1)

usethis::use_data(coef_fit1
                  ,summary_fit1
                  ,print_fit1
                  ,internal = TRUE)





