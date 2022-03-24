#'
#'
#' #'
#' #' test if a value is in bounds.
#' #'
#' #'
#' #' @param b a numeric vector
#' #' @param sqrtTheta a numeric scalar
#' #'
#'
#'
#' check_limits <- function(b, sqrtTheta) {
#'   dplyr::between(b, -1.5 * sqrtTheta, 1.5 * sqrtTheta)
#' }
#'
#' #'
#' #' calculate cluster size for a set of clusters
#' #'
#' #' @param b vector of random numbers from N(0, \theta)
#' #' @param theta variance of the normal distribution that b is drawn from
#' #' @param base constant used is the calculation of cluster size
#' #'
#' #' Calculates cluster sizes according to the method specified in Wang 2019.
#' #'
#' #' @export
#' #'
#' #'
#'
#' make_cluster_size <- function(b, theta, base = 75, ... = ...) {
#'   sqrtTheta <- sqrt(theta)
#'
#'   inLimits <- check_limits(b, sqrtTheta)
#'
#'   Ni <- dplyr::if_else(
#'     inLimits,
#'     base * exp(b),
#'     base * exp(sign(b) * 1.5 * sqrtTheta)
#'   ) |> round()
#' }
#'
#'
#' #' Make a finite population
#' #'
#' #' @param theta variance of the random effects
#' #' @param N number of clusters in the population
#' #' @param ... arguments passed to other functions
#' #'
#' #' @details cluster sizes are determined via make_cluster_size. Everything is set up
#' #' as in Wang 2019.
#' #'
#' #' @returns Returns a tibble with the columns
#' #' \itemize{
#' #' \item cluster factor with numeric levels identifying cluster membership
#' #' \item t failure or censoring time
#' #' \item d indicator of failure. 1 indicates failure, 0 censoring
#' #' \item X1:X3 fixed effect covariates. X1 ~ N(0, 1) and varies across individuals and
#' #' clusters. X2 ~ N(0, 1) varies across clusters, but constant within clusters.
#' #' \item X3 ~ binomial(0.5, N) Binary 'treatment allocation' at the cluster level.
#' #' \item Z1:ZN columns indicating cluster membership.
#' #' }
#' #'
#' #' an additional attribute has been attached to the tibble: 'b', which holds the
#' #' cluster-level random effects/frailties
#' #'
#' #'
#'
#'
#'
#' make_population <- function(theta, N, fixed_effects, random_effects = NULL,  prCensoring, ...) {
#'   sqrt_theta <- sqrt(theta)
#'
#'   if(is.null(random_effects)) {
#'
#'       random_effects <- rnorm(N, mean = 0, sd = sqrt_theta)
#'
#'   }
#'
#'   names(random_effects) <- as.character(1:N)
#'
#'   cluster_size_multiplier <- rnorm(N, mean = 0, sd = sqrt_theta)
#'
#'   Ni <- make_cluster_size(b = cluster_size_multiplier, theta = theta, ... = ...)
#'
#'   Z <- model.matrix(~ Z - 1, data = data.frame(Z = factor(rep(1:N, Ni))))
#'
#'   X <- matrix(c(
#'     rnorm(sum(Ni), 0, 1),
#'     rep(rnorm(N, 0, 1), Ni),
#'     rep(rbinom(N, 1, 0.5), Ni)
#'   ),
#'   ncol = 3, dimnames = list(
#'     rownames = NULL,
#'     colnames = c("X1", "X2", "X3")
#'   )
#'   )
#'
#'   eij <- rexp(sum(Ni), 0.1)
#'
#'   tij <- exp(-X %*% fixed_effects - Z %*% random_effects) * eij
#'
#'   dij <- sample(c(0, 1), sum(Ni), replace = TRUE, prob = c(prCensoring, 1 - prCensoring))
#'
#'   myData <- tibble::tibble(
#'     cluster = rep(1:N, Ni),
#'     t = tij[, 1],
#'     d = dij,
#'     tibble::as_tibble(X),
#'     tibble::as_tibble(Z)
#'   )
#'
#'   attr(myData, "b") <- random_effects
#'
#'   myData
#'
#' }
#'
#' #' draw a sample of clusters and all members is those clusters
#' #'
#' #' @param population tibble of the complete population to sample from
#' #' @param cluster variable in population identifying cluster membership
#' #' @param n number of clusters to sample
#' #' @param ni number of observations to sample from within each cluster
#' #' @param z_columns all the z columns in the population tibble
#' #' @param prop named vector of selection probabilities for clusters. Defaults to NULL i.e. simple
#' #' random sampling
#' #'
#' #' @importFrom rlang :=
#'
#'
#'
#' sample_clusters <- function(population, cluster, n, ni, z_columns, prob = NULL) {
#'
#'   cluster_label <- rlang::as_label(rlang::enquo(cluster))
#'
#'   distinct_clusters <- dplyr::distinct(population, {{ cluster }}) |>
#'     dplyr::pull({{ cluster }})
#'
#'   if(is.null(prob)) {
#'     sampled_clusters <- sample(distinct_clusters, n, prob = prob)
#'     cluster_weights <- tibble::tibble({{ cluster }} := distinct_clusters,
#'                                       wi = n/length(distinct_clusters))
#'   } else {
#'     sampled_clusters <- sample(names(prob), n, prob = prob)
#'     cluster_weights <- tibble::tibble({{ cluster }} := names(prob),
#'                                       wi = 1/prob1)
#'   }
#'
#'   sample <- population |>
#'     dplyr::filter({{ cluster }} %in% sampled_clusters) |>
#'     dplyr::group_by({{ cluster }}) |>
#'     dplyr::slice_sample(n = ni) |>
#'     dplyr::mutate(wji = n/dplyr::n()) |>
#'     dplyr::ungroup()
#'
#'   # we don't want Z columns for unsampled clusters. Z columns could be suffixed
#'   # 1:nClusters, but I'm leaving as-is so that the column name aligns with the
#'   # cluster number.
#'
#'   z_to_drop <- z_columns[!(z_columns %in% glue::glue("Z{sampled_clusters}"))]
#'
#'   sample <- sample |>
#'     dplyr::select(-all_of(z_to_drop))
#'
#'   weighted_sample <- dplyr::left_join(sample, cluster_weights, by = cluster_label)
#'
#'   weighted_sample
#'
#' }
#'
#'
#'
#'
#'
#'
#'
#'
#'
#'
