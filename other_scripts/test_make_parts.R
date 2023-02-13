

# update, added weights.

devtools::load_all("C:/Users/Bradley/Documents/PhD_local/svycoxme")

cluster_str <- data.frame(table(Size = rpois(2500, 2) + 6)) |>
  dplyr::filter(Freq >=10)

cluster_str_list <- split(cluster_str, seq(nrow(cluster_str)))

max_cluster_digits <- max(nchar(as.character(cluster_str$Size)))
max_cluster_freq_digits <- max(nchar(as.character(cluster_str$Freq)))

set.seed(949742)

pop_list <- lapply(cluster_str_list, function(cluster_info){
  
  k <- cluster_info$Freq
  nk <- as.numeric(as.character(cluster_info$Size))
  
  k_id <- formatC(k, width = max_cluster_freq_digits, flag = "0")
  nk_id <- formatC(nk, width = max_cluster_digits, flag = "0")
  
  the_data <- one_dataset(~X1 + X2 + (1 | M1),
                          dists = list(X1 = ~rnorm(n),
                                       X2 = ~rnorm(n),
                                       M1 = ~rep(1:k, each = nk),
                                       error = ~rexp(n, 10),
                                       stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n)),
                          dist_args = list(k = k, nk = nk,
                                           n = k * nk),
                          coefficients = c(X1 = 1, X2 = 1),
                          random_effect_variance = list(M1 = 0)
                          
  )
  
  dplyr::mutate(the_data, id = paste(nk_id,k_id, M1, sep = "_" ))
  
})

pop <- Reduce(rbind.data.frame, pop_list)
#
# # weights
# n_clusters <- dplyr::n_distinct(pop$id)
#
# n_clusters_in_samp <- 50
#
# # sample data
# samp_cluster_ids <- unique(pop$id)[sample.int(n_clusters, n_clusters_in_samp)]
#
# my_samp <- pop[pop$id %in% samp_cluster_ids, ]
#
# my_samp$weight <- (n_clusters_in_samp/n_clusters)^-1

# simple random sample
my_samp <- pop[sample.int(nrow(pop), 500),]

# sort samp by time.
my_samp <- my_samp[order(my_samp$stat_time),]


my_samp$weights = nrow(pop)/500

N_hat <- sum(my_samp$weights)

coxfit_weighted <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2, data = my_samp, weights = weights)
coxfit_unweighted <- survival::coxph(survival::Surv(stat_time, stat) ~ X1 + X2, data = my_samp)

# The same?
cbind(
  coef(coxfit_weighted),
  coef(coxfit_unweighted))

parts_weighted <- make_parts(coxfit_weighted, data = my_samp)
parts_unweighted <- make_parts(coxfit_unweighted, data = my_samp)

parts_unweighted$weights
parts_weighted$weights

parts <- parts_weighted


# calculate things by hand and compare.

names(parts)

all(parts$stat == my_samp$stat)
all(parts$time == my_samp$time)
all(parts$weights == my_samp$weights)

# 
parts$S0

beta <- coef(coxfit_weighted)

exp_Xbeta <- exp(as.matrix(my_samp[c("X1", "X2")]) %*% beta)

w_exp_Xbeta <- my_samp$weights * exp_Xbeta

# S0 is the cumulative sum of this rev(this vector)

S0 <- rev(cumsum(rev(w_exp_Xbeta)))

# seems ok
max(abs(parts$S0 - S0))
plot(S0, parts$S0)

# s1 has two columns, one for each X. 

s1_1 <- rev(cumsum(rev(my_samp$X1 * w_exp_Xbeta)))
s1_2 <- rev(cumsum(rev(my_samp$X2 * w_exp_Xbeta)))

# seems fine
max(abs(parts$S1[,1] - s1_1))
plot(s1_1, parts$S1[,1])

max(abs(parts$S1[,2] - s1_2))
plot(s1_2, parts$S1[,2])

# calc diff increases with size of term.
plot(parts$S1[,1] - s1_1, parts$S1[,1])
lines(lowess(parts$S1[,1] - s1_1, parts$S1[,1]))

# weight was wrong, so lets see if it's fixed now.

















































