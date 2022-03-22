

library(devtools)


document()
load_all()


theta = 0.5
number_of_clusters = 200
fixed_effects = c(1, -0.7, 0.5)

set.seed(203948)
myPop <- make_population(theta = theta, N = 200, fixed_effects = fixed_effects, prCensoring = 0.2)

true_random_effects <- attr(myPop, "b")

n_clusters = 20

myNi <- myPop |>
  dplyr::group_by(cluster) |>
  dplyr::summarise(n = dplyr::n()) |>
  dplyr::pull(n)

# some clusters are smaller than ni if any of these cluster are chosen, the sample size
# will be less than n * ni. It also means sample sizes vary a bit.
# This behavior is inherited from Wang 2019

mySample <- sample_clusters(myPop, cluster, n_clusters, 100, z_columns = glue::glue("Z{1:number_of_clusters}"))

names(mySample)

# The next step is to fit the Wang likelihood. For this I need weights.
# cluster level weights depend on the number of clusters sampled, and the total
# number of clusters.

x_columns = paste0("X", 1:3)
z_columns = paste0("Z", sort(unique(mySample$cluster)))

D_theta = diag(n_clusters) * theta

rownames(D_theta) <- z_columns
colnames(D_theta) <- z_columns

cox_start <- coef(survival::coxph(survival::Surv(t, d) ~ X1 + X2 + X3, data = mySample))

test_results <- estimate_all_parameters(
  fixed_effects_start = cox_start,
  random_effects_start = rep(0, n_clusters),
  theta_start = 0.2,
  x_columns = x_columns,
  z_columns = z_columns,
  t = t,
  i = cluster,
  wi = wi,
  wji = wji,
  dij = d,
  data = mySample,
  max_iterations = 1000,
  eps = .Machine$double.eps)


library(ggplot2)

lapply(test_results, "[[", "fixed_effects") |>
  Reduce(f = rbind) |>
  as.data.frame() |>
  dplyr::mutate(iteration = dplyr::row_number() - 1, .before = everything()) |>
  tidyr::pivot_longer(cols = glue::glue("X{1:3}")) |>
  ggplot(aes(iteration, value, colour = name))  + geom_point() + geom_line()


lapply(test_results, "[[", "theta") |>
  Reduce(f = rbind) |>
  as.data.frame() |>
  dplyr::mutate(iteration = dplyr::row_number() - 1, .before = everything()) |>
  ggplot(aes(iteration, V1))  + geom_point() + geom_line()

lapply(test_results, "[[", "random_effects") |>
  Reduce(f = rbind) |>
  as.data.frame() |>
  dplyr::mutate(iteration = dplyr::row_number() - 1, .before = everything()) |>
  tidyr::pivot_longer(cols = glue::glue("V{1:n_clusters}")) |>
  ggplot(aes(iteration, value, colour = name))  + geom_point() + geom_line()






test_results[[length(test_results)]][["random_effects"]] |> var() |> sqrt()



debugonce(calcLinearPredictor)

minuslp(c(fixed_effects, rep(0, n_clusters)), X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta, data = mySample)

# This gives the maximised estimates of b and beta, for a given theta. This needs to be iterated with
# a step to give the maximised theta for a given b and beta.
res <- optim(c(fixed_effects, rep(0, n_clusters)), fn = lp, gr = lp_grd,
             X = x_columns, Z = z_columns,
             t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
             data = mySample, hessian = TRUE, control = list(fnscale = -1),
             method = "BFGS")

H <- res$hessian

random_effects_index <- (length(fixed_effects)+1):nrow(H)

H22 <- H[random_effects_index, random_effects_index]

bhat <- res$par[random_effects_index]

# need cluster level weights in order corresponding to the order of random effects
# currently i'm assuming the clusters are ordered correctly, but a sensible
# sorting step would be prudent.

# this will only be needed one per run, as wi doesn't change, so should be passed in
# as a parmeter, not calculated on the fly.
wi <- mySample |>
  dplyr::distinct(cluster, wi) |>
  dplyr::pull(wi)

-0.5 * diag(wi) %*% solve(D_theta)

debugonce(calcPenalty)
calcPenalty(mySample, D = D_theta, Z = z_columns, i = cluster, wi = wi)


theta_hat = (sum(wi * bhat * bhat) + sum(wi^2 * diag(solve(H22)))) / sum(wi)

sortedIndexedData <- sortAndIndex(mySample, sort_vars = t)

addedLP <- calcLinearPredictor(sortedIndexedData, X = x_columns, Z = z_columns, parms = c(fixed_effects, rep(0, 50)))

names(addedLP)

calcRiskSets(data = addedLP, vars = z_columns, varCol = "Zi", wi = wi, wji = wji, A = A, index = index)

varColName <- rlang::sym(varCol)

addedLP |>
  tidyr::pivot_longer(cols = dplyr::all_of(z_columns), names_to = "Zi") |>
  dplyr::mutate(
    wi_wji_A = wi * wji * A,
    "Zi_wi_wji_A" := wi_wji_A * value,
    "Zi" := forcats::as_factor(Zi)
  ) |>
  dplyr::arrange(dplyr::desc(index)) |>
  dplyr::group_by(Zi) |>
  dplyr::mutate(
    cumsum_wi_wji_A = cumsum(wi_wji_A),
    "cumsum_Zi_wi_wji_A" := cumsum(Zi_wi_wji_A),
    cumsum_wi_wji_A_squared = cumsum_wi_wji_A^2
  ) |>
  dplyr::arrange(index)

dlp_b(c(fixed_effects, rep(0, 50)), X = x_columns, Z = z_columns,
      t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta, data = mySample)

dlp_beta(c(fixed_effects, rep(0, 50)), X = x_columns, Z = z_columns,
      t = t, wi = wi, wji = wji, dij = d, data = mySample)

lp_grd(c(fixed_effects, rep(0, 50)), X = x_columns, Z = z_columns,
       t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta, data = mySample)

#################################
# the call with hessian = TRUE is really slow, so I will have to implement
# optimisation that uses the derivatives and second derivatives.
# These will need careful testing.


system.time(
res2 <- optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
              X = x_columns, Z = z_columns,
              t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
              data = mySample, control = list(fnscale = -1)))

system.time(
  res2 <- optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
                X = x_columns, Z = z_columns,
                t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
                data = mySample, control = list(fnscale = -1),
                method = "BFGS"))

microbenchmark::microbenchmark(
  optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
        X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample, control = list(fnscale = -1), hessian = TRUE),
  optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
        X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample, control = list(fnscale = -1), hessian = TRUE,
        method = "BFGS"), times = 10)


res <- optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
             X = x_columns, Z = z_columns,
             t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
             data = mySample, control = list(fnscale = -1), hessian = TRUE,
             method = "BFGS")

myHess <- res$hessian

dim(myHess)

myHess22 <- myHess[4:53, 4:53]






debugonce(dlp_b_b)
# tests for the second derivatives.

resbb <- dlp_b_b(c(fixed_effects, rep(0, n_clusters)), X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample)

plot_data_1 <- tidyr::pivot_wider(resbb, names_from = Zm, values_from = ll) |>
  tibble::column_to_rownames("Zi") |> as.matrix()


library(plotly)
fig <- plot_ly(x = 1:n_clusters,
               y = 1:n_clusters,
               z = plot_data_1) |> add_surface()

fig

fig <- plot_ly(x = 1:50,
               y = 1:50,
               z = H22) |> add_surface()
fig


fig <- plot_ly(x = 1:53,
               y = 1:53,
               z = myHess) |> add_surface()

fig

# this is the penalty in the log likelihood
-0.5 * matrix(bhat, nr = 1) %*% diag(wi) %*% solve(D_theta) %*% matrix(bhat)

# first derivative (wrt b)
-0.5 * diag(wi) %*% solve(D_theta) %*% matrix(bhat)

# second derivative (wrt b)
-0.5 * diag(wi) %*% solve(D_theta)




microbenchmark::microbenchmark(
  optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
        X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample, control = list(fnscale = -1)),
  optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
        X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample, control = list(fnscale = -1),
        method = "BFGS"), times = 10)


stats::nlminb(c(fixed_effects, rep(0, 50)),
      objective = minuslp,
      gradient = minus_lp_grd,
      X = x_columns, Z = z_columns,
      t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
      data = mySample)


# do i get the same results when using numeric and analytic gradients?

r1 <- optim(c(fixed_effects, rep(0, 50)), fn = lp,
        X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample, control = list(fnscale = -1),
        method = "BFGS")

r2 <- optim(c(fixed_effects, rep(0, 50)), fn = lp, gr = lp_grd,
        X = x_columns, Z = z_columns,
        t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
        data = mySample, control = list(fnscale = -1),
        method = "BFGS")

true_values <- c(fixed_effects, true_random_effects[as.numeric(gsub("Z", "", z_columns))])

comp1 <- cbind(nograd = r1$par, grad = r2$par, diff = r1$par - r2$par,
      eps1 = abs(r1$par - true_values),
      eps2 = abs(r2$par - true_values)) |> tibble::as_tibble()

plot(density(comp1$eps1))
lines(density(comp1$eps2))

plot(true_values, type = "p")
points(comp1$nograd, col = "red")
points(comp1$grad, col = "blue")

# With the gradient the approximation is better, but the difference is small.
# Much faster though.



# more trying with second derivatives.

parameter_values <- c(fixed_effects, rep(0, n_clusters))

d1 <- mySample
d2 <- sortAndIndex(d1, sort_vars = t)

d3 <- calcLinearPredictor(d2, x_columns, z_columns, parameter_values)

debugonce(calcRiskSets)
d4 <- calcRiskSets(d3, z_columns, "Zi")


dplyr::inner_join(
d4 |> dplyr::select(index, Zi, cumsum_Zi_wi_wji_A) |> tidyr::pivot_wider(names_from = Zi, values_from = cumsum_Zi_wi_wji_A),
d4 |> dplyr::select(index, Zi, value) |> tidyr::pivot_wider(names_from = Zi, values_from = value),
by = 'index') |> View()

# by hand linear predictors
lphand <- d2[, x_columns] |> as.matrix() %*% fixed_effects + d1[, z_columns] |> as.matrix() %*% rep(0, n_clusters)

ZCrossProducts <- calcCrossProducts(d2, z_columns, z_columns, "Zi", "Zm")


debugonce(dlp_b_b)

resbb <- dlp_b_b(c(fixed_effects, rep(0, n_clusters)), X = x_columns, Z = z_columns,
                 t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
                 data = mySample)

myZ <- d3[,z_columns] |> as.matrix()

w_i <- d3$wi

wji <- d3$wji

A <- d3$A

(wji * myZ * A)[rev(seq_along(A)), ]

d3 |>
  dplyr::select(index, cluster, wji, dplyr::starts_with("Z"), A) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("Z"), ~ wji * .x * A,
                              .names = "wji_A_{col}")) |>
  dplyr::arrange(dplyr::desc(index)) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("wji_A_"), ~ cumsum(.x),
                              .names = "cumsum_{col}")) |>
  dplyr::arrange(index) |>
  dplyr::group_by(cluster)


################# Try and get this hessian correct.

theta = 0.5
population_n_clusters = 200
fixed_effects = c(1, -0.7, 0.5)

set.seed(203)
myPop <- make_population(theta = theta, N = population_n_clusters,
                         fixed_effects = fixed_effects, prCensoring = 0.2)

true_random_effects <- attr(myPop, "b")

sample_n_clusters = 20
sample_n_per_cluster = 100


oneRep <- function(...){
#draw a sample

mySample <- sample_clusters(myPop, cluster, sample_n_clusters, sample_n_per_cluster,
                            z_columns = glue::glue("Z{1:sample_n_clusters}"))


#set up info for estimation

x_columns = paste0("X", 1:3)
z_columns = paste0("Z", sort(unique(mySample$cluster)))

D_theta = diag(sample_n_clusters) * theta

rownames(D_theta) <- z_columns
colnames(D_theta) <- z_columns

cox_start <- coef(survival::coxph(survival::Surv(t, d) ~ X1 + X2 + X3, data = mySample))
start_values <- c(cox_start, rep(0, length(z_columns)))

# calculate hessian for fixed effects numerically
r2 <- optim(start_values, fn = lp, gr = lp_grd,
            X = x_columns, Z = z_columns,
            t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
            data = mySample, control = list(fnscale = -1), hessian = TRUE,
            method = "BFGS")

hes_num <- r2$hessian[1:3, 1:3]


## calculate hessian for fixed effects analytically

sortedIndexedData <- sortAndIndex(mySample, sort_vars = t )

addedLP <- calcLinearPredictor(sortedIndexedData, X = x_columns, Z = z_columns, parms = start_values)

addedCumsums <- calcRiskSets(addedLP, x_columns, "Xr")

mySample2 <- mySample |>
  dplyr::mutate(w_kl = wi * wji)

XX <- apply(mySample2[, x_columns], 1, function(X){
  matrix(X) %*% t(matrix(X))
}, simplify = FALSE)


XXcumsum_parts <- lapply(seq_along(XX), function(kl){

  d <- addedCumsums[kl*3,]

  d$wi * d$wji * XX[[kl]] * d$A

})

XXcumsum <- lapply(seq_along(XXcumsum_parts), function(i){
  Reduce(`+`, XXcumsum_parts[i:length(XXcumsum_parts)])
})

hessian_parts <- lapply(seq_along(XXcumsum), function(i){

  d <- addedCumsums[i*3 - 2:0, ]

  d$wi * d$wji * d$d * (d$cumsum_Xr_wi_wji_A^2/d$cumsum_wi_wji_A_squared - XXcumsum[[i]]/d[1, ]$cumsum_wi_wji_A)

})

hes_ana <- Reduce(`+`, hessian_parts)

list(hes_num = hes_num, hes_ana = hes_ana)

}


all_reps <- lapply(1:1000, oneRep)

saveRDS(all_reps, "all_reps.rds")

all_hes_num <- lapply(all_reps, "[[", 1)

long_hes_num <- Reduce(rbind, all_hes_num) |>
  tibble::as_tibble() |>
  dplyr::mutate(Method = "Numerical",
                Xr = rep(c("X1", "X2", "X3"), length(all_hes_num))) |>
  tidyr::pivot_longer(cols = c(X1, X2, X3), names_to = "Xs")

all_hes_ana <- lapply(all_reps, "[[", 2)

long_hes_ana <- Reduce(rbind, all_hes_ana) |>
  tibble::as_tibble() |>
  dplyr::rename(X1 = V1, X2 = V2, X3 = V3) |>
  dplyr::mutate(Method = "Analytical",
                Xr = rep(c("X1", "X2", "X3"), length(all_hes_ana))) |>
  tidyr::pivot_longer(cols = c(X1, X2, X3), names_to = "Xs")


library(ggplot2)

dplyr::bind_rows(long_hes_num, long_hes_ana) |>
  ggplot(aes(x = value, group = Method, colour = Method)) +
  geom_density() +
  facet_grid(vars(Xr), vars(Xs), scales = "free")

long_hes_num |>
  ggplot(aes(x = value, group = Method, colour = Method)) +
  geom_density() +
  facet_wrap(facets = vars(Xr, Xs), scales = "free")

long_hes_ana |>
  ggplot(aes(x = value, group = Method, colour = Method)) +
  geom_density() +
  facet_wrap(facets = vars(Xr, Xs), scales = "free")





###########

try2 <- addedLP |>
  dplyr::select(cluster, index, wi, wji, all_of(x_columns), A) |>
  dplyr::mutate(dplyr::across(all_of(x_columns), ~ wi * wji * .x * A, .names = "wi_wji_A_{col}")) |>
  dplyr::arrange(dplyr::desc(index)) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("wi_wji_A_X"), ~cumsum(.x), .names = "{col}_cumsum")) |>
  dplyr::arrange(index) |>
  dplyr::select(cluster, index, wi_wji_A_X1_cumsum, wi_wji_A_X2_cumsum, wi_wji_A_X3_cumsum) |>
  tidyr::pivot_longer(cols = c(wi_wji_A_X1_cumsum, wi_wji_A_X2_cumsum, wi_wji_A_X3_cumsum)) |>
  dplyr::mutate(Xr = gsub("wi_wji_A_(X[1-3])_cumsum", "\\1", name)) |>
  dplyr::select(-name)

try1 <- addedCumsums |>
  dplyr::select(cluster, index , Xr, cumsum_Xr_wi_wji_A) |>
  dplyr::mutate(Xr = as.character(Xr))

dplyr::left_join(try1, try2, by = c("index", "Xr")) |>
  dplyr::mutate(equal = value == cumsum_Xr_wi_wji_A,
                diff = value - cumsum_Xr_wi_wji_A) |>
  dplyr::filter(!equal) |>
  dplyr::summarise(max(diff))

addedLP |>
  dplyr::select(cluster, index, wi, wji, d, all_of(x_columns), A) |>
  dplyr::mutate(dplyr::across(all_of(x_columns), ~ wi * wji * .x * A, .names = "wi_wji_A_{col}")) |>
  dplyr::arrange(dplyr::desc(index)) |>
  dplyr::mutate(dplyr::across(dplyr::starts_with("wi_wji_A_X"), ~cumsum(.x), .names = "{col}_cumsum"),
                wi_wji_A_cumsum = cumsum(wi * wji * A)) |>
  dplyr::arrange(index) |>
  dplyr::mutate(lX1 = wi * wji * d * (X1 - wi_wji_A_X1_cumsum/wi_wji_A_cumsum),
                lX2 = wi * wji * d * (X2 - wi_wji_A_X2_cumsum/wi_wji_A_cumsum),
                lX3 = wi * wji * d * (X3 - wi_wji_A_X3_cumsum/wi_wji_A_cumsum)) |>
  dplyr::summarise(X1 = sum(lX1),
                   X2 = sum(lX2),
                   X3 = sum(lX3))
