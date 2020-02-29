##### Plotting simulation results of causal disovery
setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

# Discrete environments

res_fixed_200 <- readRDS(file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_200.RDS")
res_fixed_50 <- readRDS(file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_50.RDS")
res_fixed_50_ell_1 <- readRDS(file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_50_ell_1.RDS")
res_fixed_50_10pi <- readRDS(file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_50_10pi.RDS")
res_fixed_reg <- readRDS(file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_reg.RDS")

tmp_fixed <- data.frame("Accuracy" = c(colSums(res_fixed_200)/50, colSums(res_fixed_50)/50, 
                            colSums(res_fixed_50_ell_1)/50, colSums(res_fixed_50_10pi)/50, colSums(res_fixed_reg)/50),
                  "Alpha_val" = rep(seq(0.1, 10, length.out = 19), 5), 
                  "SD" = c(apply(res_fixed_200, 2, sd), apply(res_fixed_50, 2, sd), 
                           apply(res_fixed_50_ell_1, 2, sd), apply(res_fixed_50_10pi, 2, sd),
                           apply(res_fixed_reg, 2, sd)),
                  "Method" = c(rep("l_2", 38), rep("l_1", 19), rep("l_2", 19), rep("LM", 19)), 
                  "Setting" = c(rep("999_obs", 19), rep("50_obs", 38), rep("50_obs_10_pi", 19) , rep("50_obs", 19)))

g <- ggplot(data = tmp_fixed, aes(x = tmp_fixed$Alpha_val, y = Accuracy, color = Method, shape = Setting)) + 
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Influence of H", ", " ,c[YH2]))) + 
  ylab("Prob. of recov. the causal parent") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])),expression(paste("l"[2])), "LM"), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) +
  scale_shape_discrete(labels = c(expression(omega[H] == paste(2* pi)), expression(omega[H] == paste(10 * pi)), 
                                  "200 obs."), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))
g$labels$colour <- "Method"
g$labels$shape <- "Setting"

ggsave(filename = "Graphics/causal_disc_fixed.pdf", plot = g, width = 5.8, height = 3)


