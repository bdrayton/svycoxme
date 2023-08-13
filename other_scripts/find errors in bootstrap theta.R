
#' the theta bootstrap is failing a lot. 
#' 
#' here I run many iterations and catch everything so I can diagnose what is going wrong.
#' 
#' for theta = 2, I expect ~20% errors.
#' 

nreps = 100 # reps per theta
my_theta = 2

true_coefs = c(1, -0.7, 0.5, -1)

library(survival)
devtools::load_all(path = "C:/Users/bdra011/Documents/PhD_local/svycoxme")

make_pop <- function(theta, clusters = 2.5e5, lambda = 2, cluster_minimum = 2){
  
  cluster_str <- data.frame(table(Size = rpois(clusters, lambda) + cluster_minimum)) %>%
    dplyr::filter(Freq >=10)
  
  cluster_str_list <- split(cluster_str, seq(nrow(cluster_str)))
  
  max_cluster_digits <- max(nchar(as.character(cluster_str$Size)))
  max_cluster_freq_digits <- max(nchar(as.character(cluster_str$Freq)))
  
  pop_list <- lapply(cluster_str_list, function(cluster_info){
    
    k <- cluster_info$Freq
    nk <- as.numeric(as.character(cluster_info$Size))
    
    k_id <- formatC(k, width = max_cluster_freq_digits)
    nk_id <- formatC(nk, width = max_cluster_digits)
    
    the_data <- one_dataset(~X1 + X2 + X3 + stratum + (1 | M),
                            dists = list(X1 = ~rnorm(n),
                                         X2 = ~rep(rnorm(k), each = nk),
                                         X3 = ~rep(rbinom(k, 1, 0.5), each = nk),
                                         M = ~rep(1:k, each = nk),
                                         error = ~rexp(n, 10),
                                         stat = ~sample(rep(c(0, 1), round(n * c(0.2, 0.8))), n),
                                         stratum = ~rep(c(0, 1), c(floor(2/3 * k) * nk, ceiling(1/3 * k) * nk))),
                            dist_args = list(k = k, nk = nk,
                                             n = k * nk),
                            coefficients = true_coefs,
                            random_effect_variance = c(M=theta)
    )
    
    dplyr::mutate(the_data, id = paste(nk_id,k_id, M, sep = "_" ))
    
  })
  
  pop <- Reduce(rbind.data.frame, pop_list)
  
  attr(pop, "cluster_str") <- cluster_str
  
  pop
  
}


one_rep <- function(rep = 0){
  
  sample_clusters <- dplyr::bind_rows(
    pop %>%
      dplyr::filter(stratum == 0) %>%
      dplyr::distinct(stratum, id) %>%
      dplyr::mutate(prob = 50/dplyr::n()) %>%
      dplyr::slice_sample(n = 50),
    pop %>%
      dplyr::filter(stratum == 1) %>%
      dplyr::distinct(stratum, id) %>%
      dplyr::mutate(prob = 50/dplyr::n()) %>%
      dplyr::slice_sample(n = 50)
  ) %>% dplyr::select(stratum, id, prob)
  
  sample_data <- dplyr::left_join(sample_clusters, pop, by = c("stratum", "id")) %>%
    dplyr::mutate(scaled_weight = (1/prob)/(1/mean(prob)))
  
  d2 <- svydesign(~id, probs = ~prob, strata = ~stratum, data = sample_data)
  d3 <- as.svrepdesign(d2)
  
  svycoxme_fit <- svycoxme(Surv(stat_time, stat) ~ X1 + X2 + X3 + stratum + (1 | M), design = d3)
  
  svycoxme_fit
  
}

# set cluster
cl <- parallel::makeCluster(parallel::detectCores()-1)

# Set a different seed on each member of the cluster (just in case)
parallel::clusterSetRNGStream(cl)

# load svycoxme in all the nodes.
parallel::clusterEvalQ(cl, {
  devtools::load_all(path = "C:/Users/bdra011/Documents/PhD_local/svycoxme")
  library(survival)
})

parallel::clusterExport(cl, c("one_rep"))

pop <- make_pop(theta = my_theta, clusters = 1e4, lambda = 100, cluster_minimum = 0)

parallel::clusterExport(cl, c("pop"))

fits <- parallel::parLapply(cl, seq(nreps), function(i) try(one_rep(i)))

# dump errors
not_error <- sapply(fits, function(one_fit) !("try-error" %in% class(one_fit)))

# fits <- fits[not_error]

parallel::stopCluster(cl)









