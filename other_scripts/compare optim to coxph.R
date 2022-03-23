


library(ggplot2)
library(dplyr)

res_df <- readr::read_rds("C:/Users/bdra011/OneDrive - The University of Auckland/PhD/outputs/simulations/data/results_optim_vs_coxph.rds")

res_df %>%
  group_by(method) %>%
  summarise(mean(X1),
            mean(X2),
            mean(X3))

res_df %>%
  tidyr::pivot_longer(cols = c("X1", "X2", "X3")) %>%
  ggplot(aes(name, value, colour = method)) + geom_point(position = position_dodge(width = 0.5))






