# Plot 3: plots and comments

# Continuous environments

res_200 <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_200.RDS")
res_50 <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_50.RDS")
res_50_ell <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_50_ell.RDS")
res_50_reg <- readRDS(file = "Simulations/Causal Discovery/Plot3/res_varying_freq_50_reg.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_50)/50, 
                                 colSums(res_50_ell)/50, colSums(res_50_reg)/50, colSums(res_200)/50), 
                  "Frequency_of_H" = rep(seq(0, 60, length.out = 11), 4), 
                  "SD" = c(apply(res_50, 2, sd), apply(res_50_ell, 2, sd), 
                           apply(res_50_reg, 2, sd), apply(res_200, 2, sd)),
                  "Method" = c(rep("l_2", 11), rep("l_1", 11), rep("LM", 11), rep("l_2", 11)),
                  "Setting" = c(rep("50_obs", 33), rep("999_obs", 11)))


g <- ggplot(data = tmp, aes(x = factor(tmp$Frequency_of_H), y = Accuracy, color = Method, shape = Setting)) + 
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Frequency of H", ", " ,omega[H]))) + 
  ylab("Prob. of recov. the causal parent") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])), expression(paste("l"[2])),"LM" ), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5)) + 
  scale_shape_discrete(labels = c("50", "200"), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))
g$labels$colour <- "Method"
g$labels$shape <- "Obs."

ggsave(filename = "Graphics/causal_frequency_H.pdf", plot = g, width = 5.8, height = 3)
