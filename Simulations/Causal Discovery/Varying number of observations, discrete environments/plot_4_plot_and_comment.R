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
                  "Method" = c(rep("l_2", 7), rep("l_1", 7), rep("l_2", 14), rep("LM", 7)),
                  "Setting" = c(rep("2pi", 14), rep("999pi", 7), rep("Cy_10", 7), rep("2pi", 7)))

g <- ggplot(data = tmp, aes(x = factor(tmp$Number_of_observations), y = Accuracy, color = Method, shape = Setting)) + 
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Number of observations"))) + 
  ylab("Prob. of recov. the causal parent") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])), expression(paste("l"[2])), "LM"), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) + 
  scale_shape_discrete(labels = c(expression(paste(omega[H] == 2* pi)), 
                                  expression(paste(omega[h] == 10 * pi)), expression(c[YH2] == 10)), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))
g$labels$colour <- "Method"
g$labels$shape <- "Setting"


ggsave(filename = "Graphics/causal_n_obs_discrete.pdf", plot = g, width = 5.8, height = 3)

