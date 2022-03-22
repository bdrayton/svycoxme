
#' estimate_all_parameters
#'
#' estimates \eqn{\beta}, b and \eqn{\theta} using the method described in Wang
#' 2019 and Ripatti et al 2001
#'
#' @export


estimate_all_parameters <- function(
fixed_effects_start,
random_effects_start,
theta_start,
x_columns,
z_columns,
t,
i,
wi,
wji,
dij,
data,
max_iterations = 1000,
eps = .Machine$double.eps) {

  on.exit(return(all_estimates))

  random_effects_index <- seq_along(random_effects_start) + length(fixed_effects_start)

  cluster_weights <- data |>
    dplyr::distinct({{ i }}, {{ wi }}) |>
    dplyr::pull({{ wi }})

  current_estimates = list(fixed_effects = fixed_effects_start,
                           random_effects = random_effects_start,
                           theta = theta_start)

  current_eps = 1

  D_theta = diag(length(z_columns)) * theta_start

  rownames(D_theta) <- z_columns
  colnames(D_theta) <- z_columns

  all_estimates <- list(current_estimates)

  for (iteration in 1:max_iterations) {

    if (current_eps <= eps) break

    params <- c(current_estimates[["fixed_effects"]], current_estimates[["random_effects"]])

    res <- optim(params, fn = lp, gr = lp_grd,
                 X = x_columns, Z = z_columns,
                 t = t, i = cluster, wi = wi, wji = wji, dij = d, D = D_theta,
                 data = mySample, hessian = TRUE, control = list(fnscale = -1),
                 method = "BFGS")

    H22 <- res$hessian[random_effects_index, random_effects_index]

    bhat <- res$par[random_effects_index]

    theta_hat <- (sum(cluster_weights * bhat * bhat) +
                   sum(cluster_weights^2 * diag(solve(H22)))) / sum(cluster_weights)


    # This is a hack to force theta to be positive.
    theta_hat <- abs(theta_hat)

    diag(D_theta) <- theta_hat

    former_estimates <- current_estimates

    current_estimates <- list(fixed_effects = res$par[-random_effects_index],
                              random_effects = bhat,
                              theta = theta_hat)

    all_estimates[[iteration + 1]] <- current_estimates

    current_eps <- max(abs(unlist(former_estimates) - unlist(current_estimates)))

  }

}






