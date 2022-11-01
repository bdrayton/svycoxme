
# compare svycoxph and svycoxph calibrated.
library(survey)
library(dplyr)


# replicate example in Breslow

# see for code:
# http://staff.pubhealth.ku.dk/~pd/survey-2009/misc/calibrate-copenhagen-2x2.pdf
#
# https://search.r-project.org/CRAN/refmans/addhazard/html/nwtsco.html
# another National Wilms Tumor Study

#
# data(nwtco, package="survival")
# nrow(nwtco)
# names(nwtco)
#
# # Sample all relapses, all patients with unfavorable histology
# # by local pathologist, 10% of remainder
#
# sample_from <- nwtco$seqno[!(nwtco$rel == 1 | nwtco$instit == 1)]
# keep <- c(nwtco$seqno[nwtco$rel == 1 | nwtco$instit == 1],
#           sample(sample_from, size = 0.1 * length(sample_from)))
#
# nwtco$incc2 <- nwtco$seqno %in% keep
#
#
# dccs2<-twophase(id=list(~seqno,~seqno),
#                 strata=list(NULL,~interaction(rel,instit)),
#                 data=nwtco, subset=~incc2)
#
# gccs8<-calibrate(dccs2, phase=2,
#                  formula=~interaction(rel,stage,instit))
#
# ###
# set up sampling.

# this is the phase 1 dataset. assumed to be drawn from an unobsered super-population.
data(nwtsco, package = "addhazard")
nwts <- nwtsco

# set up strata

nwts <- nwts %>%
  mutate(instit = factor(instit, c(0, 1), labels = c("FH", "UH")),
         stage = factor(stage, levels = 1:4,
                        labels = c("I-II", "I-II",
                                   "III-IV", "III-IV")),
         age_bin = factor(0 + age<1,
                          labels = c("<1Yr", ">=1Yrs")),
         strata =interaction(instit, stage, age_bin, relaps))

## collapse smallest 13 strata together.
s_counts <- nwts %>%
  count(strata) %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number()) %>%
  mutate(strat4 = if_else(rank <= 3, rank, 4L)) %>%
  group_by(strat4) %>%
  mutate(weight = case_when(strat4 == 1 ~ n/120,
                            strat4 == 2 ~ n/160,
                            strat4 == 3 ~ n/120,
                            TRUE ~ 1))

nwts_join <- left_join(nwts, select(s_counts, strata, strat4, weight)) %>%
  mutate(id = row_number())

nwts_pop <- select(nwts_join, instit, age_bin, stage, relaps,
                   trel, id, histol, weight, strat4, tumdiam)

# draw the complex sample.
nwts_samp <- bind_rows(
  filter(nwts_pop, strat4 == 1) %>%
    dplyr::slice_sample(n = 120),
  filter(nwts_pop, strat4 == 2) %>%
    dplyr::slice_sample(n = 160),
  filter(nwts_pop, strat4 == 3) %>%
    dplyr::slice_sample(n = 120),
  filter(nwts_pop, strat4 == 4))

nwts_design <- svydesign(~1, strata = ~strat4, weights = ~weight,
                         data = nwts_samp)

impmodel2 <- glm(histol~instit*age_bin*stage*relaps,
                 data=nwts_design$variables,
                 family=binomial,
                 weights = weights(nwts_design))

library(MuMIn)

options(na.action = "na.fail")
dredged_models <- MuMIn::dredge(impmodel2)
summary(MuMIn::get.models(dredged_models, 1)[[1]])
options(na.action = "na.omit")


# best model according to dredge
impmodel <- svyglm(formula = histol ~ age_bin + instit + relaps +
                     stage + age_bin:stage + instit:stage,
                   family = quasibinomial,
                   design = nwts_design)

# impute histol for everyone in the population
nwts_pop$imphist <- predict(impmodel, newdata=nwts_pop, type="response")

# replace imputed values with known values where available
nwts_pop$imphist[nwts_samp$id] <- nwts_samp$histol

# fit the phase one model with imputed data
ifmodel <- coxph(Surv(trel, relaps) ~ imphist,
                 data = nwts_pop)

inffun <- matrix(resid(ifmodel, type = "dfbeta"), ncol = 1)

colnames(inffun) <- c("if1")

nwts_if <- cbind(nwts_pop, inffun)

names(nwts_if)

# pairs(nwts_if[, c('imphist', 'age', 'stage', 'tumdiam', 'if1', 'if2', 'if3', 'if4', 'if5', 'if6')])

# with(nwts_if, plot(if1, imphist))

nwts_if <- mutate(nwts_if, in.sample = id %in% nwts_samp$id)

if_design <- twophase(id = list(~1, ~1), subset = ~in.sample,
                      strata = list(NULL, ~strat4),
                      data = nwts_if)


if_cal <- calibrate(if_design, phase=2, calfun="raking",
                    ~if1+strat4)


m1 <- svycoxph(Surv(trel, relaps)~histol*age_bin+stage*tumdiam,
               design=if_design)

m1a <- svycoxph(Surv(trel, relaps)~histol*age_bin+stage*tumdiam,
                design=if_cal)

round(coef(summary(m1)),2)
round(coef(summary(m1a)),2)




# sample simdat

a_sample <- bind_rows(
  filter(simdat, strat4 == 1) %>%
    dplyr::slice_sample(n = 120),
  filter(simdat, strat4 == 2) %>%
    dplyr::slice_sample(n = 160),
  filter(simdat, strat4 == 3) %>%
    dplyr::slice_sample(n = 120),
  filter(simdat, strat4 == 4))




nwts$seqno <- seq_len(nrow(nwts))
sample_from <- nwts$seqno[!(nwts$relaps == 1 | nwts$instit == 1)]
keep <- c(nwtco$seqno[nwts$rel == 1 | nwts$instit == 1],
          sample(sample_from, size = 0.1 * length(sample_from)))

# sample size for strata that don't have pr_selection = 1
n_sample <- ceiling(0.1 * length(sample_from))

# derive weights
nwts <- nwts %>%
  mutate(pr_1 = relaps == 1 | instit ==1) %>%
  group_by(pr_1) %>%
  mutate(n = n()) %>%
  ungroup() %>%
  mutate(pr_selection = case_when(pr_1 ~ 1,
                                  TRUE ~ n_sample/n))

nwts$in.subsample <- nwts$seqno %in% keep

nwts_design <- svydesign(~1, probs = ~pr_selection, strata = ~pr_1, data = nwts)

impmodel <- glm(histol~instit,
                data=nwts, subset=in.subsample, family=binomial)
nwts$imphist <- predict(impmodel, newdata=nwts, type="response")
nwts$imphist[nwts$in.subsample] <- nwts$histol[nwts$in.subsample]
ifmodel <- coxph(Surv(trel,relaps)~imphist*age+I(stage>2)*tumdiam,
                 data=nwts)

summary(ifmodel)

inffun <- resid(ifmodel, "dfbeta")
colnames(inffun) <- paste("if",1:6,sep="")

head(inffun)

nwts_if <- cbind(nwts, inffun)

# pairs(nwts_if[, c('imphist', 'age', 'stage', 'tumdiam', 'if1', 'if2', 'if3', 'if4', 'if5', 'if6')])

# with(nwts_if, plot(if1, imphist))


nwts_design2 <- twophase(id = list(~seqno, ~seqno), subset = ~in.subsample,
                  strata = list(NULL, ~interaction(instit, relaps)),
                  data = nwts)

if_design <- twophase(id = list(~1, ~1), subset = ~in.subsample,
                      strata = list(NULL, ~interaction(instit, relaps)),
                      data = nwts_if)

if_cal <- calibrate(if_design, phase=2, calfun="raking",
                    ~if1+if2+if3+if4+if5+if6+relaps*instit)
m1 <- svycoxph(Surv(trel, relaps)~histol*age+I(stage>2)*tumdiam,
               design=nwts_design)
m1a <- svycoxph(Surv(trel, relaps)~histol*age+I(stage>2)*tumdiam,
                design=nwts_design2)
m2 <- svycoxph(Surv(trel, relaps)~histol*age+I(stage>2)*tumdiam,
               design=if_cal)
m3 <- coxph(Surv(trel, relaps)~imphist*age+I(stage>2)*tumdiam,
            data=nwts)
m4 <- coxph(Surv(trel, relaps)~histol*age+I(stage>2)*tumdiam,
            data=nwts)

round(cbind(coef(m1), coef(m1a), coef(m2), coef(m3), coef(m4)), 3)

round(cbind(
  coef(summary(m1))[,"se(coef)"],
  coef(summary(m1a))[,"se(coef)"],
  coef(summary(m2))[,"se(coef)"],
  coef(summary(m3))[,"se(coef)"],
  coef(summary(m4))[,"se(coef)"]), 3)

# look at the influence functions

with(filter(nwts_if, imphist ==1), plot(age, if2))
with(filter(nwts_if, imphist ==0), plot(age, if2))




data(api)
clus1<-svydesign(id=~dnum, weights=~pw, data=apiclus1, fpc=~fpc)

m0 <- svyglm(api00~ell+mobility+emer, clus1)
var_cal <- calibrate(clus1, formula=~api99+ell+mobility+emer,
                       pop=c(6194,3914069, 141685, 106054, 70366),
                       bounds=c(0.1,10))
m1<-svyglm(api00~ell+mobility+emer, design=var_cal)

popmodel <- glm(api99~ell+mobility+emer, data=apipop,
                    na.action=na.exclude)

inffun <- dfbeta(popmodel)

head(inffun)

index <- match(apiclus1$snum, apipop$snum)
clus1if <- update(clus1, ifint = inffun[index,1],
                    ifell=inffun[index,2], ifmobility=inffun[index,3],
                    ifemer=inffun[index,4])
if_cal <- calibrate(clus1if,
                      formula=~ifint+ifell+ifmobility+ifemer,
                      pop=c(6194,0,0,0,0))

m2<-svyglm(api00~ell+mobility+emer, design=if_cal)
coef(summary(m0))
coef(summary(m1))
coef(summary(m2))


nrow(nwtco)
names(nwtco)

table(nwtco$histol, exclude = NULL)
table(nwtco$instit, exclude = NULL)

# make strata from event-free survival, age, stage, institutional histology

nwtco$instit <- factor(nwtco$instit, levels = c(1, 2),
                       labels = c("Favourable", "Unfavourable"))

table(nwtco$instit)

nwtco$stage <- factor(nwtco$stage, levels = 1:4,
                      labels = c("I-II", "I-II", "III-IV", "III-IV"))

table(nwtco$stage)

nwtco$age_bin <- factor(0 + nwtco$age<12,
                        labels = c("<1Yr",
                                   ">=1Yrs"))

# rel is the relapse indicator.
table(nwtco$rel)

nwtco$strata <- with(nwtco, interaction(instit, stage, age_bin, rel))

c(table(nwtco$strata))

## collapse smallest 13 strata together.

library(dplyr)

s_counts <- nwtco %>%
  count(strata) %>%
  arrange(desc(n)) %>%
  mutate(rank = row_number()) %>%
  mutate(strat4 = if_else(rank <= 3, rank, 4L),
         weight = case_when(strat4 == 1 ~ 1783/120,
                            strat4 == 2 ~ 971/160,
                            strat4 == 3 ~ 419/120,
                            TRUE ~ 1))

simdat <- left_join(nwtco, select(s_counts, strata, strat4, weight))

# sample simdat

a_sample <- bind_rows(
  filter(simdat, strat4 == 1) %>%
    dplyr::slice_sample(n = 120),
  filter(simdat, strat4 == 2) %>%
    dplyr::slice_sample(n = 160),
  filter(simdat, strat4 == 3) %>%
    dplyr::slice_sample(n = 120),
  filter(simdat, strat4 == 4))

# pretend that histol isn't know, and predict it using sample data

des <- svydesign(ids = ~1, strata = ~strat4, weights = ~weight,
          data = a_sample)


# 1 = favourable histology

m1 <- svyglm(I(histol == 2) ~ instit, family = quasibinomial(),
             design = des)
m2 <- svyglm(I(histol == 2) ~ instit + as.factor(rel),
             family = quasibinomial(), design = des)
m3 <- svyglm(I(histol == 2) ~ instit + rel + as.factor(strat4),
             family = quasibinomial(), design = des)

AIC(m1, m2, m3)

# use prediction equations to impute values for everyone

preds <- (predict(m2, newdata = simdat[, c("instit", "rel")], type = "response") > 0.5) + 1

# sensitivity and specificity pretty close to numbers in the paper.
xtabs(~ preds + simdat$histol)

simdat2 <- simdat

simdat2$histol = preds

# fit model to full data

imp_fit <- coxph(Surv(edrel, rel) ~ age_bin + stage + histol + instit, data = simdat2)

inflnc <- resid(imp_fit, type = "dfbeta")

plot(density(inflnc[,4]))



# cluster structure in the population.
cluster_str <- data.frame(table(Size = rpois(2.5e5, 2) + 2)) |>
  dplyr::filter(Freq >=10)

cluster_str_list <- split(cluster_str, seq(nrow(cluster_str)))

max_cluster_digits <- max(nchar(as.character(cluster_str$Size)))
max_cluster_freq_digits <- max(nchar(as.character(cluster_str$Freq)))

pop_list <- lapply(cluster_str_list, function(cluster_info){

  k <- cluster_info$Freq
  nk <- as.numeric(as.character(cluster_info$Size))

  k_id <- formatC(k, width = max_cluster_freq_digits)
  nk_id <- formatC(nk, width = max_cluster_digits)

  the_data <- one_dataset(~X1 + stratum + (1 | M),
                          dists = list(X1 = ~rnorm(n),
                                       M = ~rep(1:k, each = nk),
                                       error = ~rexp(n, 10),
                                       stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                                       stratum = ~rep(c(0, 1), c(floor(2/3 * k) * nk, ceiling(1/3 * k) * nk))),
                          dist_args = list(k = k, nk = nk,
                                           n = k * nk),
                          coefficients = c(1, 0),
                          random_effect_variance = c(M=0)
  )

  dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))

})

pop <- Reduce(rbind.data.frame, pop_list)

# keep the first obs for each cluster

pop_1 <- pop |>
  dplyr::group_by(id)|>
  dplyr::filter(dplyr::row_number() == 1) |>
  dplyr::ungroup()

# add a variable that is related to the model variable.

pop_1$X4 = rnorm(nrow(pop_1)) + 0.5 * pop_1$X1

# add cut up version of X1 for calibration
pop_1 <- pop_1 |> dplyr::mutate(X1_cat = cut(X1, quantile(X1, 0:10/10), labels = 1:10))

# population level coefs

pop_fit <- coxph(Surv(stat_time, stat) ~ X1 + stratum, data = pop_1)

# beta tilde from breslow et al
# not normally observed
pop_coef <- coef(pop_fit)

# beta hat, from a sample, n, has variance from beta tilde (from finite to superpop),
# plus the variance of beta hat - beta tilde which is the uncertainty due to incomplete
# data for the cohort "phase two" component of variance.

# weighting should reduce phase 2 variance. Auxiliary data must be correlated with data in the model.

# fit separate models for predicting each partly missing variable. best for just one or two variables.

# let's say X4 is available for everyone.

with(pop_1[sample.int(n = nrow(pop_1), size = 20000), ] , plot(X4, X1))

fit <- lm(X1 ~ X4, data = pop_1)
abline(coef(fit), col = "red")

#### draw a complex sample
sample_clusters <- dplyr::bind_rows(
  pop_1 |>
    dplyr::filter(stratum == 0) |>
    dplyr::distinct(stratum, id) |>
    dplyr::mutate(prob = 200/dplyr::n()) |>
    dplyr::slice_sample(n = 200),
  pop_1 |>
    dplyr::filter(stratum == 1) |>
    dplyr::distinct(stratum, id) |>
    dplyr::mutate(prob = 200/dplyr::n()) |>
    dplyr::slice_sample(n = 200)
) |> dplyr::select(stratum, id, prob)

sample_data <- dplyr::left_join(sample_clusters, pop_1, by = c("stratum", "id")) |>
  dplyr::mutate(scaled_weight = (1/prob)/(1/mean(prob)))


#### 1. Develop weighted regression modls from phase 2 data for prediction of partly missing information.
des <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data)

pmod <- svyglm(X1 ~ X4, family = gaussian(), design = des)

summary(pmod)

# predict X1 for the entire cohort

pop_1$X1_pred <- predict(pmod, newdata = pop_1[, "X4", drop = FALSE])

# relationship in the sample
svyplot(X1~X4, design = des)
abline(coef(svyglm(X1 ~ X4, family = gaussian(), design = des)))

# predicted vs true value in the population.
with(pop_1[sample.int(n = nrow(pop_1), size = 10000), ] , plot(X1, X1_pred))
summary(fit <- lm(X1_pred ~ X1, data = pop_1))
abline(coef(fit), col = "red")
abline(h = 0)
# fit a model using the imputed values

# cut X1_pred into deciles

pop_1 <- pop_1 |> dplyr::mutate(X1_pred_cat = cut(X1_pred, quantile(X1_pred, 0:10/10), labels = 1:10))

pop_fit_imp <- coxph(Surv(stat_time, stat) ~ X1_pred_cat, data = pop_1)

summary(pop_fit_imp)

# extract the df_betas

influence <- resid(pop_fit_imp, type = "dfbeta")

head(influence)

# use df_betas to calibrate weights.

pop_totals <- colSums(influence + 1)

cal_names(~ X1_cat, des)

des_cal <- calibrate(des, ~X1_cat, pop_totals, calfun = "raking")






############
# the pop true coefs will be a bit different.
# true_coefs = c(1, -0.7, 0.5, -1)


table(pop$X1_cat)

# fit a svycoxph model. extract dfbetas. are they correlated with the vars?












one_rep <- function(){

  sample_clusters <- dplyr::bind_rows(
    pop |>
      dplyr::filter(stratum == 0) |>
      dplyr::distinct(stratum, id) |>
      dplyr::mutate(prob = 50/dplyr::n()) |>
      dplyr::slice_sample(n = 50),
    pop |>
      dplyr::filter(stratum == 1) |>
      dplyr::distinct(stratum, id) |>
      dplyr::mutate(prob = 50/dplyr::n()) |>
      dplyr::slice_sample(n = 50)
  ) |> dplyr::select(stratum, id, prob)

  sample_data <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id")) |>
    dplyr::mutate(scaled_weight = (1/prob)/(1/mean(prob)))

  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data)
  d3 <- as.svrepdesign(d2)


  coxph_robust_fit_d2 <- coxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, data = sample_data, robust = TRUE, weights = scaled_weight)
  svycoxph_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d2)
  svycoxph_rep_fit_d2 <- svycoxph(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum, design = d3)

  data.frame(model = rep(c("coxph_weighted_robust", "svycoxph", "svycoxph_rep"), each = length(true_coefs)),
             X = names(coef(svycoxph_fit_d2)),
             true_value = true_coefs,
             estimate = c(coef(coxph_robust_fit_d2),
                          coef(svycoxph_fit_d2),
                          coef(svycoxph_rep_fit_d2)),
             rbind(confint(coxph_robust_fit_d2),
                   confint(svycoxph_fit_d2),
                   confint(svycoxph_rep_fit_d2))) |>
    dplyr::mutate(error = estimate - true_value,
                  hit = true_value >= X2.5.. & true_value <= X97.5..)

}

one_rep()

reps <- replicate(1000, one_rep(), simplify = FALSE)

reps_df <- Reduce(rbind.data.frame, reps)

prop.table(xtabs(~hit + model + X, data = reps_df), margin = c(2, 3))


# replicate example in Breslow
data(nwtco, package="survival")


data(airquality)

## ignoring missingness, using model-based standard error
summary(lm(log(Ozone)~Temp+Wind, data=airquality))

## Without covariates to predict missingness we get
## same point estimates, but different (sandwich) standard errors
daq<-estWeights(airquality, formula=~1,subset=~I(!is.na(Ozone)))
summary(svyglm(log(Ozone)~Temp+Wind,design=daq))

## Reweighting based on weather, month
d2aq<-estWeights(airquality, formula=~Temp+Wind+Month,
                 subset=~I(!is.na(Ozone)))
summary(svyglm(log(Ozone)~Temp+Wind,design=d2aq))

