###### Plots and comments

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

# Discrete environments

res_200 <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_200.RDS")
res_50 <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50.RDS")
res_50_ell <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50_ell1.RDS")
res_50_10pi <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50_10pi.RDS")
res_50_reg <- readRDS(file = "Simulations/Causal Discovery/Plot2/res_varying_influence_50_reg.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_50)/50, 
                                 colSums(res_50_ell)/50, colSums(res_50_10pi)/50, colSums(res_50_reg)/50), 
                  "Influence_of_H" = rep(seq(0, 20, length.out = 11), 4), 
                  "SD" = c(apply(res_50, 2, sd), 
                           apply(res_50_ell, 2, sd), apply(res_50_10pi, 2, sd), apply(res_50_reg, 2, sd)),
                  "Method" = c(rep("l_2", 11), rep("l_1", 11), rep("l_2", 11), rep("LM", 11)), 
                  "Setting" = c(rep("50_obs", 22), rep("50_obs_10_pi", 11) , rep("50_obs", 11)))

g <- ggplot(data = tmp, aes(x = factor(tmp$Influence_of_H), y = Accuracy, color = Method, shape = Setting)) + 
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Influence of H", ", " ,c[YH2]))) + 
  ylab("Prob. of recov. the causal parent") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])),expression(paste("l"[2])), "LM"), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_shape_discrete(labels = c(expression(paste(2* pi)), expression(paste(10 * pi))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))
g$labels$colour <- "Method"
g$labels$shape <- expression(paste(omega[H]))

ggsave(filename = "Graphics/causal_influence_H_disc.pdf", plot = g, width = 5.8, height = 3)

# Continuous environments

res_200 <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_200-kopi.RDS")
res_50 <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_50-kopi.RDS")
res_50_ell <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_50_ell1-kopi.RDS")
res_50_10pi <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_50_10pi-kopi.RDS")
res_50_reg <- readRDS(file = "Simulations/Causal Discovery/Plot2b/res_varying_cont_reg-kopi.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_50)/50, 
                                 colSums(res_50_ell)/50, colSums(res_50_10pi)/50, colSums(res_50_reg)/50), 
                  "Influence_of_H" = rep(seq(0, 20, length.out = 11), 4), 
                  "SD" = c(apply(res_50, 2, sd), 
                           apply(res_50_ell, 2, sd), apply(res_50_10pi, 2, sd), apply(res_50_reg, 2, sd)),
                  "Method" = c(rep("l_2", 11), rep("l_1", 11), rep("l_2", 11), rep("LM", 11)), 
                  "Setting" = c(rep("50_obs", 22), rep("50_obs_10_pi", 11) , rep("50_obs", 11)))

g <- ggplot(data = tmp, aes(x = factor(tmp$Influence_of_H), y = Accuracy, color = Method, shape = Setting)) + 
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Influence of H", ", " ,c[YH2]))) + 
  ylab("Prob. of recov. the causal parent") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])),expression(paste("l"[2])), "LM"), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_shape_discrete(labels = c(expression(paste(2* pi)), expression(paste(10 * pi))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))
g$labels$colour <- "Method"
g$labels$shape <- expression(paste(omega[H]))



ggsave(filename = "Graphics/causal_influence_H_cont.pdf", plot = g, width = 5.8, height = 3)

