# the simplest example of calibration using dfbetas.

# sample needs one more variable than population

# sample variable needs to be predicted using other sample variables.

# # try and use imputed variables instead of the hybrid.

# from Thomas's book

library(survey)

data(api)
clus1<-svydesign(id=~dnum, weights=~pw,
                  data=apiclus1, fpc=~fpc)
# fit a model for api 2000 using variables available
# at phase 2. These are x variables.

m0 <- svyglm(api00 ~ ell + mobility + emer, clus1)

# idk if I need complete data
apipop_complete <- apipop[complete.cases(apipop[, c("ell", "mobility", "emer")]), ]
(pop_totals <- model.matrix(~api99+ell+mobility+emer, data = apipop) |> colSums())

# calibrate the weights using pop totals
var_cal <- calibrate(clus1, formula = ~api99+ell+mobility+emer,
                     pop=pop_totals, bounds = c(0.1, 10))

# fit a model for api 1999 using the calibrated design
m1 <- svyglm(api00 ~ ell + mobility + emer, design = var_cal)

# fit the same model to the population
popmodel <- glm(api99 ~ ell + mobility + emer, data = apipop,
             na.action = na.exclude)

# extract the dfbetas from this pop-level fit.
inffun <- dfbeta(popmodel)

# add the influence functions to the design
index <- match(apiclus1$snum, apipop$snum)
clus1if <- update(clus1,
                  ifint =      inffun[index, 1],
                  ifell =      inffun[index, 2],
                  ifmobility = inffun[index, 3],
                  ifemer =     inffun[index, 4])

# calibrate using the influence functions
if_cal <- calibrate(clus1if,
                    formula = ~ifint + ifell + ifmobility + ifemer,
                    pop = c(pop_totals[1], ifint = 0, ifell = 0, ifmobility = 0,
                            ifemer = 0))

# fit a model for api2000 using the IF calibrated design
m2 <- svyglm(api00 ~ ell + mobility + emer, design = if_cal)

# Uncalibrated design
(m0_summary <- coef(summary(m0)))

# Calibrated with population variables
(m1_summary <- coef(summary(m1)))

# Calibrated with population influence functions
(m2_summary <- coef(summary(m2)))

# compare standard errors

std_err <- Reduce(cbind, lapply(list(m0_summary, m1_summary, m2_summary), function(df) df[, "Std. Error"]))

# factor of reduction
round(matrix(std_err[ , 1], nr = 4, nc = 3) / std_err , 1)

##############
# use this calibration approach with cox models and the wilms tumour data.

# also consider
# data(nwts2ph, package = "addhazard")
data(nwtsco, package = "addhazard")
nwts <- nwtsco

# set up strata
nwts <- nwts |>
  dplyr::mutate(instit = factor(instit, c(0, 1), labels = c("FH", "UH")),
         stage4 = stage == 4,
         stage = factor(stage, levels = 1:4,
                        labels = c("I-II", "I-II",
                                   "III-IV", "III-IV")),
         age_bin = factor(0 + age<1,
                          labels = c("<1Yr", ">=1Yrs")),
         strata =interaction(instit, stage, age_bin, relaps))

## collapse smallest 13 strata together.
s_counts <- nwts |>
  dplyr::count(strata) |>
  dplyr::arrange(desc(n)) |>
  dplyr::mutate(rank = dplyr::row_number()) |>
  dplyr::mutate(strat4 = dplyr::if_else(rank <= 3, rank, 4L)) |>
  dplyr::group_by(strat4) |>
  dplyr::mutate(weight = dplyr::case_when(strat4 == 1 ~ n/120,
                            strat4 == 2 ~ n/160,
                            strat4 == 3 ~ n/120,
                            TRUE ~ 1))

nwts_join <- dplyr::left_join(nwts, dplyr::select(s_counts, strata, strat4, weight)) |>
  dplyr::mutate(id = dplyr::row_number())

# nwts_pop <- dplyr::select(nwts_join, instit, age_bin, stage, relaps,
#                    trel, id, histol, weight, strat4, tumdiam)

nwts_pop <- nwts_join

# draw the complex sample.
nwts_samp <- dplyr::bind_rows(
  dplyr::filter(nwts_pop, strat4 == 1) |>
    dplyr::slice_sample(n = 120),
  dplyr::filter(nwts_pop, strat4 == 2) |>
    dplyr::slice_sample(n = 160),
  dplyr::filter(nwts_pop, strat4 == 3) |>
    dplyr::slice_sample(n = 120),
  dplyr::filter(nwts_pop, strat4 == 4))

nwts_pop$in.subsample <- nwts_pop$id %in% nwts_samp$id

# fit imptation model
impmodel <- glm(histol~instit + age_bin + stage4*study, data = nwts_pop, subset = in.subsample, family = binomial )

# predict histology for everyone
nwts_pop$imphist <- predict(impmodel, newdata = nwts_pop, type = "response")

# replace predicted histology with actual histology where available (i.e. sampled children)
nwts_pop$imphist[nwts_pop$in.subsample] <- nwts$histol[nwts_pop$in.subsample]

# fit a model to the population data supplemented with the imputed histology
ifmodel <- coxph(Surv(trel, relaps) ~ imphist * age + stage * tumdiam, data = nwts_pop)

# append influence functions to the pop data
inffun <- resid(ifmodel, "dfbeta")
colnames(inffun) <- paste0("if", seq(ncol(inffun)))

nwts_pop_if <- cbind(nwts_pop, inffun)
# design with pop data
if_design_tp <- twophase(id = list(~1, ~1), subset = ~in.subsample,
                      strata = list(NULL, ~strat4), data = nwts_pop_if)

if_design <- svydesign(~1, strata = ~strat4, weights = ~weight,
                         data = nwts_pop_if, subset = ~in.subsample)

pop_totals <- colSums(model.matrix(~ if1 + if2 + if3 + if4 + if5 + strat4,
                                   data = nwts_pop_if))

if_cal <- calibrate(if_design, phase = 2, calfun = "raking",
                    formula = ~ if1 + if2 + if3 + if4 + if5 + strat4,
                    population = pop_totals)

if_cal_tp <- calibrate(if_design_tp, phase = 2, calfun = "raking",
                    formula = ~ if1 + if2 + if3 + if4 + if5 + strat4)

# design with sample data. same as if_design
nwts_design <- svydesign(~1, strata = ~strat4, weights = ~weight,
                         data = nwts_samp)

# fit a model using the uncalibrated design
m0 <- svycoxph(Surv(trel, relaps) ~  histol * age + stage * tumdiam, design = if_design)

cal_form <- ~ imphist * age + stage * tumdiam
pop_totals <- colSums(model.matrix(cal_form, data = nwts_pop))

# calibrate the weights using pop totals
var_cal <- calibrate(if_design, formula = cal_form,
                     pop=pop_totals) # dropped , bounds = c(0.1, 10)

# fit a model using the variable calibrated design
m1 <- svycoxph(Surv(trel, relaps) ~ histol * age + stage * tumdiam, design = var_cal)
# fit a model using the influence function calibrated design with imputation
m1a <- svycoxph(Surv(trel, relaps) ~ histol * age + stage * tumdiam, design = if_cal)
# same, but with the twophase design
m1a_tp <- svycoxph(Surv(trel, relaps) ~ histol * age + stage * tumdiam, design = if_cal_tp)

# regular cox regression on the population data
m2 <- coxph(Surv(trel, relaps) ~ histol * age + stage * tumdiam, data = nwts_pop_if)

# Uncalibrated design
(m0_summary <- coef(summary(m0)))

# Calibrated with population variables # really biased?
(m1_summary <- coef(summary(m1)))

# Calibrated with population influence functions using imputed data
(m1a_summary <- coef(summary(m1a)))

m2_summary <- coef(summary(m2))

### influence function calibration without imputation (instit swapped for histol)
popmodel <- coxph(Surv(trel, relaps) ~ instit * age + stage * tumdiam,
                  data = nwts_pop, na.action = na.exclude)

# extract the dfbetas from this pop-level fit.
inffun <- resid(popmodel, "dfbeta")

# add the influence functions to the design
index <- match(nwts_samp$id, nwts_pop$id)
nwts_design_if <- update(nwts_design,
                  if1 = inffun[index, 1],
                  if2 = inffun[index, 2],
                  if3 = inffun[index, 3],
                  if4 = inffun[index, 4],
                  if5 = inffun[index, 5])

# calibrate using the influence functions
if_cal <- calibrate(nwts_design_if,
                    formula = ~if1 + if2 + if3 + if4 + if5,
                    pop = c(pop_totals[1], if1 = 0, if2 = 0, if3 = 0, if4 = 0, if5 = 0))

# fit a model using the IF calibrated design
m3 <- svycoxph(Surv(trel, relaps) ~ histol * age + stage * tumdiam, design = if_cal)

(m3_summary <- coef(summary(m3)))
m2_summary
m1a_summary

# the imputed CI is narrower than the population fit. FPC correction? the CI for coxph is wrong without it.
confint(m3)
confint(m2)
confint(m1a)
confint(m1a_tp)



fpc_props <- colSums(model.matrix(~factor(strat4)-1, data = nwts_pop))/nrow(nwts_pop)

fpc_df <- data.frame(strat4 = 1:4, prop = 0.99)

nwts_pop_fpc <- dplyr::left_join(nwts_pop, fpc_df, by = "strat4")

des <- svydesign(~1, probs = 1, strata = ~strat4, data = nwts_pop_fpc,
                 fpc = ~prop)

m4 <- svycoxph(Surv(trel, relaps) ~ histol * age + stage * tumdiam, design = des)

confint(m4)
confint(m3)
confint(m2)
confint(m1a)

# any different using twophase?




