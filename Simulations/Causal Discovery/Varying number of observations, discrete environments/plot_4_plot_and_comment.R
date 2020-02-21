# Plot 4 varying n_obs plots and comments

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

res_1_2pi <- readRDS(file = "Simulations/Causal Discovery/Plot4/res_varying_n_obs_1_2pi.RDS")
res_1_2pi_ell1 <- readRDS(file = "Simulations/Causal Discovery/Plot4/res_varying_n_obs_1_2pi_ell1.RDS")
res_1_10pi <- readRDS(file = "Simulations/Causal Discovery/Plot4/res_varying_n_obs_1_10pi.RDS")
res_10_2pi <- readRDS(file = "Simulations/Causal Discovery/Plot4/res_varying_n_obs_10_2pi.RDS")
res_2pi_reg <- readRDS(file = "Simulations/Causal Discovery/Plot4/res_varying_n_obs_reg.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_1_2pi)/50, colSums(res_1_2pi_ell1)/50, 
                                 colSums(res_1_10pi)/50, colSums(res_10_2pi)/50, colSums(res_2pi_reg)/50), 
                  "Number_of_observations" = rep(c(20, 50, 100, 150, 200, 300, 400), 5), 
                  "SD" = c(apply(res_1_2pi, 2, sd), apply(res_1_2pi_ell1, 2, sd), 
                           apply(res_1_10pi, 2, sd), apply(res_10_2pi, 2, sd), apply(res_2pi_reg, 2, sd)),
                  "Exp" = c(rep("2pi", 7), rep("2pi_l1", 7), rep("10pi", 7), rep("c10_2pi", 7), 
                            rep("2pi_reg", 7)))
ggplot(data = tmp, aes(x = Number_of_observations, y = Accuracy, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Causal Discovery Accuracy with discrete environments") + 
  xlab("Number of observations") + 
  ylab("Accuracy of Causal Discovery") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.2)
