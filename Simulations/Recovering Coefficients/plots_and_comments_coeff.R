##### Plotting simulation results of coefficient recovery
setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

# H(t) = sin(2 pi t)

res_full <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_50_new.RDS")
res_sin_1 <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_1_new.RDS")
res_sin_2 <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_2_new.RDS")
res_sin_3 <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_3_new.RDS")


tmp <- data.frame("Recovery_Rate" = c(colSums(res_full)/50, colSums(res_sin_1)/50, 
                                 colSums(res_sin_2)/50, colSums(res_sin_3)/50), 
                  "H_Y_val" = rep(seq(0.1, 20, length.out = 11), 4), 
                  "SD" = c(apply(res_full, 2, sd), apply(res_sin_1, 2, sd), 
                           apply(res_sin_2, 2, sd), apply(res_sin_3, 2, sd)),
                  "Exp" = c(rep("50_obs_full", 11), rep("50_obs_1", 11), rep("50_obs_2", 11), rep("50_obs_3", 11)))

g <- ggplot(data = tmp, aes(x = tmp$H_Y_val, y = Recovery_Rate, color = Exp)) + 
  geom_errorbar(aes(ymin = Recovery_Rate - SD, ymax=Recovery_Rate + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Influence of H", ", " ,c[YH2]))) + 
  ylab("Prob. of CaDiSH beating tvlm") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("X"^1)),expression(paste("X"^2)),
                                  expression(paste("X"^3)), expression(paste("Average"))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) 
g$labels$colour <- "Predictor"

ggsave(filename = "Graphics/coeff_recov_sin.pdf", plot = g, width = 5.8, height = 3)

# H(t) = log(t + 0.1)

res_full_log <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_50_log_new.RDS")
res_log_1 <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_1_log_new.RDS")
res_log_2 <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_2_log_new.RDS")
res_log_3 <- readRDS(file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_3_log_new.RDS")

tmp_log <- data.frame("Recovery_Rate" = c(colSums(res_full_log)/50, colSums(res_log_1)/50, 
                                      colSums(res_log_2)/50, colSums(res_log_3)/50), 
                  "H_Y_val" = rep(seq(0.1, 20, length.out = 11), 4), 
                  "SD" = c(apply(res_full_log, 2, sd), apply(res_log_1, 2, sd), 
                           apply(res_log_2, 2, sd), apply(res_log_3, 2, sd)),
                  "Exp" = c(rep("50_obs_full_log", 11), rep("50_obs_1", 11), rep("50_obs_2", 11), rep("50_obs_3", 11)))

g <- ggplot(data = tmp_log, aes(x = tmp_log$H_Y_val, y = Recovery_Rate, color = Exp)) + 
  geom_errorbar(aes(ymin = Recovery_Rate - SD, ymax=Recovery_Rate + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Influence of H", ", " ,c[YH2]))) + 
  ylab("Prob. of CaDiSH beating tvlm") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("X"^1)),expression(paste("X"^2)),
                                  expression(paste("X"^3)), expression(paste("Average"))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) 
g$labels$colour <- "Predictor"

ggsave(filename = "Graphics/coeff_recov_log.pdf", plot = g, width = 5.8, height = 3)
