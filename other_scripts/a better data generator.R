
my_beta = c(1, -0.7, 0.5)
my_theta = 1
my_k = 11
my_nk = 11

my_X = c("X1", "X2", "X3")

ds <- svycoxme::one_dataset(list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta)) %>%
      dplyr::mutate(M2 = rep(c("l", "r"), ceiling(nrow(.)/2))[seq_len(my_k * my_nk)])


debugonce(coxme::coxme)
fit <- coxme::coxme(survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1|M) + (1|M2), data = ds)

coxme::VarCorr(fit)

coxme_est_parms <- c(coxme::fixef(fit), coxme::ranef(fit)$M)

str(fit)

one_dataset

#' the dataset needs
#' IDs - individual, event (for repeated events), and one or more cluster ids.
#' covariates
#' failure or censoring times
#'

#' stuff that might get controlled
#' distributions for covariates, including cluster-level variation
#' number of covariates
#' number of cluster levels and number of clusters
#' number of total observations
#' variance of each set random effects
#' fixed effect coefficents
#' amount of random error
#'


#' The flow:
#' give the function a model formula
#' For each fixed effect, it will seek an expression from a list of expressions also passed in.
#'
#' something similar for random effects.
#'
#' pass in a list of control parameters. this will have elements specifying clustering and sample sizes.
#' As well, fixed effects, error, and random effects variance.
#'
#'
#'
#'

get_term.labels <- function(formula_terms){

  attr(formula_terms, "term.labels")

}

get_group_label <- function(random_term){

  term_list <- coxme:::formula2(random_term)

  as.character(term_list$group)

}

check_var_labels <- function(var_labels, dist_labels){

  test <- var_labels %in% dist_labels

  not_in_dist_labels <- paste0(var_labels[!test], collapse = ", ")

  if(!all(test))
    stop(paste0("All variables in formula must have be present in dists.\n  These variables are in your formula but not in dists: ", not_in_dist_labels))

}


one_dataset <- function(formula, dists, dist_args, coefficients, random_effects_variance, seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  flist <- coxme:::formula1(formula)

  fixed_terms <- stats::terms(flist$fixed)

  fixed_term_labels <- get_term.labels(fixed_terms)

  random_term_labels <- sapply(flist$random, get_group_label)
  names(random_term_labels) <- random_term_labels

  term_labels <- c(fixed_term_labels, random_term_labels)

  check_var_labels(term_labels, names(dists))

  vars <- lapply(dists[term_labels], lazyeval::f_eval, data = dist_args)

  vars_df <- data.frame(vars)

  lp_random_effects <- sapply(random_term_labels, function(re){

    re_data <- vars_df[, re, drop = TRUE]

    re_factor <- factor(re_data, levels = unique(re_data))

    group_counts <- table(re_factor)

    b <- rnorm(length(group_counts), mean = 0, sd = sqrt(random_effects_variance[[re]]))

    b_rep <- dplyr::left_join(data.frame(re_name = re_factor),
                              data.frame(re_name = names(group_counts), re = b), by = "re_name")

    b_rep$re

  }, simplify = "matrix")

  risk_score <- .rowSums(cbind(as.matrix(vars_df[, fixed_term_labels]) %*% coefficients, lp_random_effects),
                         m = nrow(lp_random_effects), n = ncol(lp_random_effects) + 1)

  error <- lazyeval::f_eval(dists$error, data = dist_args)

  vars_df$stat_time <- exp(-risk_score) * error

  vars_df$stat <- lazyeval::f_eval(dists$stat, data = dist_args)

  vars_df

}

# debugonce(one_dataset)

ds <- one_dataset(formula = survival::Surv(stat_time, stat) ~ X1 + X2 + X3 + (1 | M1),
                  dists = list(X1 = ~rnorm(n, 0, 1),
                               X2 = ~rep(rnorm(k, 0, 1), each = nk),
                               X3 = ~rep(rep(c(1, 0), ceiling(k/2))[seq_len(k)], each = nk),
                               M1 = ~rep(1:k, each = nk),
                               M2 = ~rep(c("l", "r"), ceiling(n/2))[seq_len(n)],
                               error = ~rexp(n, 10),
                               stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                  dist_args = list(k = 50, nk = 4, n = 50*4),
                  random_effects_variance = list(M1 = 1, M2 = 2),
                  coefficients = c(1, -0.7, 0.5))






