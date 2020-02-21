# Plot 3: plots and comments

# Continuous environments

res_200 <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_200.RDS")
res_50 <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_50.RDS")
res_50_ell <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_50_ell.RDS")
res_50_reg <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_50_reg.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_200)/50, colSums(res_50)/50, 
                                 colSums(res_50_ell)/50, colSums(res_50_reg)/50), 
                  "Frequency_of_H" = rep(seq(0, 60, length.out = 11), 4), 
                  "SD" = c(apply(res_200, 2, sd), apply(res_50, 2, sd), 
                           apply(res_50_ell, 2, sd), apply(res_50_reg, 2, sd)),
                  "Exp" = c(rep("200_obs", 11), rep("50_obs", 11), rep("50_obs_l1", 11), rep("50_obs_reg", 11)))
ggplot(data = tmp, aes(x = Frequency_of_H, y = Accuracy, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Causal Discovery Accuracy with continuous environments") + 
  xlab("Frequency of H") + 
  ylab("Accuracy of Causal Discovery") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.2)

