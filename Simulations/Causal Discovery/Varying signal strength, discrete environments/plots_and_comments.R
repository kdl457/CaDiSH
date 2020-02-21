##### Plotting simulation results of causal disovery
setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

# Discrete environments

res_fixed_200 <- readRDS(file = "Simulations/Causal Discovery/res_matrix_fixed_200.RDS")
res_fixed_50 <- readRDS(file = "Simulations/Causal Discovery/res_matrix_fixed_50.RDS")
res_fixed_50_ell_1 <- readRDS(file = "Simulations/Causal Discovery/res_matrix_fixed_50_ell_1.RDS")
res_fixed_50_10pi <- readRDS(file = "Simulations/Causal Discovery/res_matrix_fixed_50_10pi.RDS")
res_fixed_reg <- readRDS(file = "Simulations/Causal Discovery/res_matrix_fixed_reg.RDS")
plot(seq(0.1, 10, length.out = 19), colSums(res_fixed_200)/50, col = "red", ylim = c(0, 1))
points(seq(0.1, 10, length.out = 19), colSums(res_fixed_50)/50, col = "blue")
points(seq(0.1, 10, length.out = 19), colSums(res_fixed_50_ell_1)/50, col = "green")
points(seq(0.1, 10, length.out = 19), colSums(res_fixed_50_10pi)/50, col = "purple")
# Makes sense since assumption of H varying slower than environments is violated.

tmp_fixed <- data.frame("Accuracy" = c(colSums(res_fixed_200)/50, colSums(res_fixed_50)/50, 
                            colSums(res_fixed_50_ell_1)/50, colSums(res_fixed_50_10pi)/50, colSums(res_fixed_reg)/50),
                  "Alpha_val" = rep(seq(0.1, 10, length.out = 19), 5), 
                  "SD" = c(apply(res_fixed_200, 2, sd), apply(res_fixed_50, 2, sd), 
                           apply(res_fixed_50_ell_1, 2, sd), apply(res_fixed_50_10pi, 2, sd),
                           apply(res_fixed_reg, 2, sd)),
                  "Exp" = c(rep("200_obs", 19), rep("50_obs", 19), rep("50_obs_l1", 19), rep("50_obs_l1_10pi", 19),
                            rep("reg", 19)))

ggplot(data = tmp_fixed, aes(x = Alpha_val, y = Accuracy, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Causal Discovery Accuracy with fixed coefficients in each environment") + 
  xlab("Regression coefficient of the hidden variable in Y") + 
  ylab("Accuracy of Causal Discovery") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=0.2)




