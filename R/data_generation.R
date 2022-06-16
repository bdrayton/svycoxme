
#' extract term.labels attribute
#'
#' \link{\code{terms}} returns the formula it is given with a decomposition of attributes attached.
#' This function is a wrapper for `attr(., "term.labels")`
#'
#'
#'


get_group_label <- function(random_term){

  term_list <- coxme:::formula2(random_term)

  as.character(term_list$group)

}


#' Check formula against dists names
#'
#' Checks that the variables in the right hand side of the formula all have matching
#' expressions is the dists list.
#'
#' @param var_labels Names of the variables in the right hand side of the formula
#' @param dist_labels Names of the expressions in dists
#'
#' @export


check_var_labels <- function(var_labels, dist_labels){

  test <- var_labels %in% dist_labels

  not_in_dist_labels <- paste0(var_labels[!test], collapse = ", ")

  if(!all(test))
    stop(paste0("All variables in formula must have be present in dists.\n  These variables are in your formula but not in dists: ", not_in_dist_labels))

}



#'  create a dataset for testing cox model stuff.
#'
#' @param formula Model equation for the model you plan to fit to the data. Only the right hand side is used at the moment.
#' @param dists expressions for generating the variables in the formula
#' @param dist_args any values that are needed for the expressions in dists
#' @param coefficients vector of fixed effect coefficients. If there are more coefficients than fixed effects, only those up to the number of fixed effects will be used and the rest ignored.
#' @param random_effect_variance list of thetas
#' @param seed passed to set.seed. Use if you need the random processes to be reproducible.
#'
#'
#'
#' @export

one_dataset <- function(formula, dists, dist_args, coefficients = c(), random_effect_variance = list(), seed = NULL){

  if(!is.null(seed)) set.seed(seed)

  flist <- coxme:::formula1(formula)

  fixed_terms <- stats::terms(flist$fixed)

  fixed_term_labels <- attr(fixed_terms, "term.labels")

  if(length(fixed_term_labels) > length(coefficients))
    stop("Each fixed effect in the formula must have a corresponding coefficient")

  all_dists_are_formula <- all(sapply(dists, class) == "formula")

  if(!all_dists_are_formula)
    stop("All terms in dists must be left hand formulas e.g. ~rnorm(n)")

  random_term_labels <- sapply(flist$random, get_group_label, simplify = "vector")
  # naming the random terms means the results from sapply will inherit the names.
  names(random_term_labels) <- random_term_labels

  if(length(random_term_labels) > length(random_effect_variance))
    stop("Each random effect in the formula must have a corresponding random_effect_variance")

  # when there are no random terms, random_term_labels is list(). c(vector(), list()) results in a list, which I don't want. unlisting converts list() this to NULL, which makes term_labels a vector.
  term_labels <- c(fixed_term_labels, unlist(random_term_labels))

  check_var_labels(term_labels, names(dists))

  vars <- lapply(dists[term_labels], lazyeval::f_eval, data = dist_args)

  vars_df <- list2DF(vars)

  lp_random_effects <- sapply(random_term_labels, function(re){

    re_data <- vars_df[, re, drop = TRUE]

    re_factor <- factor(re_data, levels = unique(re_data))

    group_counts <- table(re_factor)

    b <- rnorm(length(group_counts), mean = 0, sd = sqrt(random_effect_variance[[re]]))

    b_rep <- dplyr::left_join(data.frame(re_name = re_factor),
                              data.frame(re_name = names(group_counts), re = b), by = "re_name")

    b_rep$re

  }, simplify = "matrix")

  # transpose non-zero length coefficients
  if(length(coefficients))
    coefficients <- t(coefficients)

  risk_score_parts_list <- list(fixed = tcrossprod(as.matrix(vars_df[, fixed_term_labels]), coefficients),
                                random = lp_random_effects)

  risk_score_parts_list <- risk_score_parts_list[lapply(risk_score_parts_list, length)>0]

  risk_score_parts_matrix <- Reduce(cbind, risk_score_parts_list)

  error <- lazyeval::f_eval(dists$error, data = dist_args)

  if(!is.null(risk_score_parts_matrix)) {
    risk_score <- .rowSums(risk_score_parts_matrix,
                           m = nrow(risk_score_parts_matrix), n = ncol(risk_score_parts_matrix))
  } else {
    risk_score <- rep(0, length(error))
    vars_df <- data.frame(stat_time = numeric(length(error)))
  }

  vars_df$stat_time <- exp(-risk_score) * error

  vars_df$stat <- lazyeval::f_eval(dists$stat, data = dist_args)

  vars_df

}
