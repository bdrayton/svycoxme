
library(tidyverse)

my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 200

one_rep <- function(){

  sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

  my_parms <- c(my_beta, b <- attr(sample_data, "random_effects"))

  D = my_theta * diag(length(b))

  Z_matrix <- model.matrix( ~ as.factor(M) - 1, data = sample_data)

  colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

  data_with_Z <- dplyr::bind_cols(sample_data, data.frame(Z_matrix))

  d2 <- sortAndIndex(data_with_Z, t)

  d3 <- calcLinearPredictor(data = d2, X = c("X1"), Z = colnames(Z_matrix), parms = my_parms[c(-2, -3)])

  d4 <- calcRiskSets(d3, vars = "X1", varCol = "Xr")

  d5 <- calcCrossProducts(d2, "X1", "X1", "Xr", "Xs")


  d6 <- left_join(d4, d5, by = "index")

  d7 <- d6 %>%
    arrange(desc(index)) %>%
    mutate(XrXs_A = XrXs * A,
           cumsum_XrXs_A = cumsum(XrXs_A)) %>%
    arrange(index)

  by_hand_hess_X1 <- d7 %>%
    mutate(
      ll_parts = stat * ((cumsum_Xr_A^2)/(cumsum_A)^2 - cumsum_XrXs_A/cumsum_A)
    ) %>%
    summarise(ll = sum(ll_parts)) %>%
    pull(ll)


  fit <- survival::coxph(survival::Surv(t, stat) ~ X1, data = sample_data)

  start_parameters = c(coef(fit), rep(0, length(b)))

  names(start_parameters) <- c("X1", paste0("Z", seq_len(length(b))))

  fit_optim <- optim(par = start_parameters,
                     fn = lp,
                     gr = lp_grd,
                     X = c("X1"),
                     t = t,
                     cluster = "M",
                     dij = stat,
                     D = D,
                     data = sample_data,
                     method = "BFGS",
                     control = list(fnscale = -1),
                     hessian = TRUE)

  data.frame(analytical = by_hand_hess_X1,
             numerical = fit_optim$hessian["X1", "X1"])

}

one_rep()

tests <- replicate(100, one_rep(), simplify = "matrix")

tests_df <- matrix(unlist(tests),nc = 2, byrow = TRUE) %>% as.data.frame()

colnames(tests_df)  <- c("analytical", "numerical")


tests_df %>%
  mutate(diff = analytical - numerical) %>%
  ggplot(aes(x = diff)) + geom_density()

tests_df %>%
  mutate(diff = analytical - numerical) %>%
  summarise(mean(diff),
            mean(analytical),
            mean(numerical))


###################################


my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 10



sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

my_parms <- c(my_beta, b <- attr(sample_data, "random_effects"))

D = my_theta * diag(length(b))

Z_matrix <- model.matrix( ~ as.factor(M) - 1, data = sample_data)

Z_names <- paste0("Z", seq(ncol(Z_matrix)))


colnames(Z_matrix) <- paste0("Z", seq(ncol(Z_matrix)))

data_with_Z <- dplyr::bind_cols(sample_data, data.frame(Z_matrix))

d2 <- sortAndIndex(data_with_Z, t)

d3 <- calcLinearPredictor(data = d2, X = c("X1", "X2", "X3"), Z = colnames(Z_matrix), parms = my_parms)

d4 <- calcRiskSets(d3, vars = c("X1", "X2", "X3"), varCol = "Xr")

d5 <- calcCrossProducts(d2, c("X1", "X2", "X3"), c("X1", "X2", "X3"), "Xr", "Xs")

d6 <- dplyr::left_join(d4, d5, by = c("index", "Xr"))

d7 <- d6 %>%
  dplyr::arrange(desc(index)) %>%
  dplyr::group_by(Xr, Xs) %>%
  dplyr::mutate(XrXs_A = XrXs * A,
                cumsum_XrXs_A = cumsum(XrXs_A)) %>%
  dplyr::arrange(index) %>%
  dplyr::group_by(Xr, Xs)

by_hand_hess_11 <- d7 %>%
  dplyr::mutate(
    ll_parts = stat * ((cumsum_Xr_A^2)/(cumsum_A)^2 - cumsum_XrXs_A/cumsum_A)
  ) %>%
  dplyr::summarise(ll = sum(ll_parts), .groups = "drop")

# dif wrt b

d4 <- calcRiskSets(d3, vars = Z_names, varCol = "Zr")

d5 <- calcCrossProducts(d2, Z_names, Z_names, "Zr", "Zs")

d6 <- dplyr::left_join(d4, d5, by = c("index", "Zr"))

d7 <- d6 %>%
  dplyr::arrange(desc(index)) %>%
  dplyr::group_by(Zr, Zs) %>%
  dplyr::mutate(ZrZs_A = ZrZs * A,
                cumsum_ZrZs_A = cumsum(ZrZs_A)) %>%
  dplyr::arrange(index) %>%
  dplyr::group_by(Zr, Zs)

by_hand_hess_22 <- d7 %>%
  dplyr::mutate(
    ll_parts = stat * ((cumsum_Zr_A^2)/(cumsum_A)^2 - cumsum_ZrZs_A/cumsum_A)
  ) %>%
  dplyr::summarise(ll = sum(ll_parts), .groups = "drop")

penalty = solve(D)

unpenalised <- by_hand_hess_22 %>%
  tidyr::pivot_wider(names_from = Zs, values_from = ll) %>%
  tibble::column_to_rownames("Zr") %>%
  as.matrix()

penalised = unpenalised - penalty

fit <- survival::coxph(survival::Surv(t, stat) ~ X1 + X2 + X3, data = sample_data)

start_parameters = c(coef(fit), rep(0, length(b)))

names(start_parameters) <- c("X1", "X2", "X3", paste0("Z", seq_len(length(b))))

fit_optim <- optim(par = start_parameters,
                   fn = lp,
                   gr = lp_grd,
                   X = c("X1", "X2", "X3"),
                   t = t,
                   cluster = "M",
                   dij = stat,
                   D = D,
                   data = sample_data,
                   method = "BFGS",
                   control = list(fnscale = -1),
                   hessian = TRUE)

diag(fit_optim$hessian)[c(-1, -2, -3)] |> sum()

sum(diag(penalised))

diag(fit_optim$hessian)[c(-1, -2, -3)]

diag(penalised)

### d beta d b

d4x <- calcRiskSets(d3, vars = c("X1", "X2", "X3"), varCol = "Xr")

d4z <- calcRiskSets(d3, vars = Z_names, varCol = "Zr")

d4 <- dplyr::left_join(
  d4x %>% dplyr::select(index, stat, A, cumsum_A, Xr, cumsum_Xr_A),
  d4z %>%  dplyr::select(index, Zr, cumsum_Zr_A), by = "index")

d5 <- calcCrossProducts(d2, c("X1", "X2", "X3"), Z_names, "Xr", "Zr")

d6 <- dplyr::left_join(d4, d5, by = c("index", "Xr", "Zr"))

d7 <- d6 %>%
  dplyr::arrange(desc(index)) %>%
  dplyr::group_by(Xr, Zr) %>%
  dplyr::mutate(XrZr_A = XrZr * A,
                cumsum_XrZr_A = cumsum(XrZr_A)) %>%
  dplyr::arrange(index) %>%
  dplyr::group_by(Xr, Zr)

by_hand_hess_21 <- d7 %>%
  dplyr::mutate(
    ll_parts = stat * ((cumsum_Xr_A * cumsum_Zr_A)/(cumsum_A)^2 - cumsum_XrZr_A/cumsum_A)) %>%
  dplyr::summarise(ll = sum(ll_parts), .groups = "drop")

by_hand_hess_21 %>%
  tidyr::pivot_wider(names_from = Zr, values_from = ll)

by_hand_hess_21 %>%
  tidyr::pivot_wider(names_from = Xr, values_from = ll)








diag(fit_optim$hessian)[c(-1, -2, -3)] |> sum()

sum(diag(penalised))

diag(fit_optim$hessian)[c(-1, -2, -3)]

diag(penalised)






