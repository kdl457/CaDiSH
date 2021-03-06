#### Plots for random graph with random parents - using AUC on models of size 1

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

res_ell_2 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1.RDS")
res_ell_1 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_ell_1.RDS")
res_reg <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_reg.RDS")

tmp <- data.frame("Average_AUC" = c(colSums(res_ell_2)/20, colSums(res_ell_1)/20, 
                                    colSums(res_reg)/20), 
                  "Frequency_of_H" = rep(seq(from = 0, to = 40, length.out = 11), 3), 
                  "SD" = c(apply(res_ell_2, 2, sd), apply(res_ell_1, 2, sd), apply(res_reg, 2, sd)),
                  "Exp" = c(rep("ell_2", 11), rep("ell_1", 11), rep("lin_reg", 11)))

g <- ggplot(data = tmp, aes(x = factor(tmp$Frequency_of_H), y = Average_AUC, color = Exp)) + 
  geom_errorbar(aes(ymin = Average_AUC - SD, ymax=Average_AUC + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Frequency of H", ", " ,omega[H]))) + 
  ylab("Average AUC") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), legend.position = c(0.795,0.88), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])),expression(paste("l"[2])),"LM" ), 
                       guide = guide_legend(direction = "horizontal", title.position = "top", title.hjust = 0.5)) 
g$labels$colour <- "Method"

ggsave(filename = "Graphics/AUC_PA_Y_random_models_1.pdf", plot = g, width = 5.8, height = 3)


# With environments on all but Y

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

res_ell_2 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_env_all.RDS")
res_ell_1 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_ell_1_env_all.RDS")
res_reg <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_reg_env_all.RDS")

tmp <- data.frame("Average_AUC" = c(colSums(res_ell_2)/20, colSums(res_ell_1)/20, 
                                    colSums(res_reg)/20), 
                  "Frequency_of_H" = rep(seq(from = 0, to = 40, length.out = 11), 3), 
                  "SD" = c(apply(res_ell_2, 2, sd), apply(res_ell_1, 2, sd), apply(res_reg, 2, sd)),
                  "Exp" = c(rep("ell_2", 11), rep("ell_1", 11), rep("lin_reg", 11)))

ggplot(data = tmp, aes(x = Frequency_of_H, y = Average_AUC, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Average AUC with randomized graph and random parents of Y, models of size 1") + 
  xlab("Frequency of H") + 
  ylab("Average AUC") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Average_AUC - SD, ymax=Average_AUC + SD), width=.2)

# With small noise coefficients

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

res_ell_2 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_small_noise.RDS")
res_ell_1 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_ell_1_small_noise.RDS")
res_reg <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_reg_small_noise.RDS")

tmp <- data.frame("Average_AUC" = c(colSums(res_ell_2)/20, colSums(res_ell_1)/20, 
                                    colSums(res_reg)/20), 
                  "Frequency_of_H" = rep(seq(from = 0, to = 40, length.out = 11), 3), 
                  "SD" = c(apply(res_ell_2, 2, sd), apply(res_ell_1, 2, sd), apply(res_reg, 2, sd)),
                  "Exp" = c(rep("ell_2", 11), rep("ell_1", 11), rep("lin_reg", 11)))

ggplot(data = tmp, aes(x = Frequency_of_H, y = Average_AUC, fill = Exp, color = Exp)) + 
  geom_point(size = 3)  + 
  ggtitle("Average AUC with randomized graph and random parents of Y, models of size 1") + 
  xlab("Frequency of H") + 
  ylab("Average AUC") + 
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  geom_errorbar(aes(ymin = Average_AUC - SD, ymax=Average_AUC + SD), width=.2)


