###### Plots and comments

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

# Discrete environments

res_200 <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_200.RDS")
res_50 <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50.RDS")
res_50_ell <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50_ell1.RDS")
res_50_10pi <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50_10pi.RDS")
res_50_reg <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50_reg.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_200)/50, colSums(res_50)/50, 
                                 colSums(res_50_ell)/50, colSums(res_50_10pi)/50, colSums(res_50_reg)/50), 
                  "Influence_of_H" = rep(seq(0, 20, length.out = 11), 5), 
                  "SD" = c(apply(res_200, 2, sd), apply(res_50, 2, sd), 
                           apply(res_50_ell, 2, sd), apply(res_50_10pi, 2, sd), apply(res_50_reg, 2, sd)),
                  "Exp" = c(rep("200_obs", 11), rep("50_obs", 11), rep("50_obs_l1", 11), rep("50_obs_l1_10pi", 11), 
                            rep("50_obs_reg", 11)))
ggplot(data = tmp, aes(x = Influence_of_H, y = Accuracy, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Causal Discovery Accuracy with discrete environments") + 
  xlab("Influence of H") + 
  ylab("Accuracy of Causal Discovery") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.2)


# Continuous environments

res_200 <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_200-kopi.RDS")
res_50 <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_50-kopi.RDS")
res_50_ell <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_50_ell1-kopi.RDS")
res_50_10pi <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_50_10pi-kopi.RDS")
res_50_reg <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_reg-kopi.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_200)/50, colSums(res_50)/50, 
                                 colSums(res_50_ell)/50, colSums(res_50_10pi)/50, colSums(res_50_reg)/50), 
                  "Influence_of_H" = rep(seq(0, 20, length.out = 11), 5), 
                  "SD" = c(apply(res_200, 2, sd), apply(res_50, 2, sd), 
                           apply(res_50_ell, 2, sd), apply(res_50_10pi, 2, sd), apply(res_50_reg, 2, sd)),
                  "Exp" = c(rep("200_obs", 11), rep("50_obs", 11), rep("50_obs_l1", 11), rep("50_obs_l1_10pi", 11), 
                            rep("50_obs_reg", 11)))
ggplot(data = tmp, aes(x = Influence_of_H, y = Accuracy, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Causal Discovery Accuracy with continuous environments") + 
  xlab("Influence of H") + 
  ylab("Accuracy of Causal Discovery") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.2)

