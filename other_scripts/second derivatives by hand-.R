

my_beta = c(1, -0.7, 0.5)
my_theta = 0.2
my_k = 10
my_nk = 20





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

tests <- replicate(100, one_rep(), simplify = "matrix")

tests_df <- matrix(unlist(tests),nc = 2, byrow = TRUE) %>% as.data.frame()

colnames(tests_df)  <- c("analytical", "numerical")


tests_df %>%
  mutate(diff = analytical - numerical) %>%
  ggplot(aes(x = diff)) + geom_density()

tests_df %>%
  mutate(diff = analytical - numerical) %>%
  summarise(mean(diff))

chol(diag(diag(-fit_optim$hessian)))

fit_coxme$hmat











