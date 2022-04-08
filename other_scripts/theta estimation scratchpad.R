

my_beta = c(1, -0.7, 0.5)
my_theta = 2
my_k = 10
my_nk = 10

my_X = c("X1", "X2", "X3")

sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

b <- attr(sample_data, "random_effects")

my_params <- c(my_beta[seq_along(my_X)], b)

# my_hessian <- ppl_hessian(parms = my_params, X = my_X, t = t,
#                           cluster = "M", dij = stat, data = sample_data,
#                           theta = my_theta)

my_K_ppl <- bb(parms = my_params, X = my_X, t = t,
               cluster = "M", dij = stat, data = sample_data,
               theta = my_theta, return_matrix = TRUE)

est_theta(b, my_K_ppl)

thetas_to_try <- seq(from = 0.3, to = 10, by = 0.01)

ll <- sapply(thetas_to_try, est_theta3, b = b, K_ppl = my_K_ppl)

plot(thetas_to_try, ll, type = "l")

optim(par = 10, fn = est_theta2, b = b, K_ppl = my_K_ppl,
      method = "Brent", lower = 0.1, upper = 50)

optim(par = 10, fn = est_theta3, b = b, K_ppl = my_K_ppl,
      method = "Brent", lower = 0.1, upper = 50)

est_theta2(1, b = b, K_ppl = my_K_ppl)

est_theta3(1, b = b, K_ppl = my_K_ppl)




inv_D = solve(D)

inv_D*inv_D

1/theta

# when D I*theta, then,

optim(par = 10, fn = est_theta2, b = b, K_ppl = my_K_ppl,
      method = "Brent", lower = 0.1, upper = 50)


est_theta2(0.2, b = b, K_ppl = my_K_ppl)

svycoxme::est_theta2(0.2, b = b, K_ppl = my_K_ppl)

#
# est_theta2 <- function(par, b, K_ppl){
#
#   theta = par
#
#   I <- diag(length(b))
#
#   D <- theta * I
#
#   inv_D <- solve(D)
#
#   c(-0.5 * (sum(diag(inv_D)) + sum(diag(solve(K_ppl) %*% inv_D %*% inv_D  )) - t(b) %*% inv_D %*% inv_D %*% b))
#
# }
#

t(b) %*% inv_D %*% inv_D %*% b


#############simulation

my_beta = c(1, -0.7, 0.5)
my_theta = 2
my_k = 50
my_nk = 10

my_X = c("X1", "X2", "X3")



res <- lapply(1:100, function(i){

  sample_data <- one_dataset(control = list(k = my_k, nk = my_nk, beta = my_beta, theta = my_theta))

  b <- attr(sample_data, "random_effects")

  my_params <- c(my_beta[seq_along(my_X)], b)

  my_K_ppl <- bb(parms = my_params, X = my_X, t = t,
                 cluster = "M", dij = stat, data = sample_data,
                 theta = my_theta, return_matrix = TRUE)

  # when D I*theta, then,

  c(
  optim(par = 0.2, fn = est_theta2, b = b,
        K_ppl = my_K_ppl,
        method = "Brent", lower = 0, upper = 100)$par,
  (b %*% b + sum(diag(solve(my_K_ppl))))/length(b))

})

df <- Reduce(rbind, res) |> as.data.frame(row.names = "")

names(df) = c("optim", "fun")

library(ggplot2)

# values
df %>%
  tidyr::pivot_longer(cols = c("optim", "fun")) %>%
  ggplot(aes(value, colour = name)) + geom_boxplot()

# residuals
df %>%
  tidyr::pivot_longer(cols = c("optim", "fun")) %>%
  ggplot(aes(value - my_theta, colour = name)) + geom_boxplot()




sum(diag(-solve(my_K_ppl) %*% inv_D %*% inv_D))

-sum(diag(inv_D %*% inv_D) * diag(solve(my_K_ppl)))


est_theta(b, my_K_ppl)














