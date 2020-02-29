#### Plots for random graph with 1 parent

setwd("~/Dropbox/KU/Thesis/R")
library(ggplot2)

res_ell_2 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_matrix_1_parent.RDS")
res_ell_1 <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_matrix_1_parent_ell_1.RDS")
res_reg <- readRDS(file = "Simulations/Causal Discovery/Randomized Graph/res_matrix_1_reg.RDS")

tmp <- data.frame("Accuracy" = c(colSums(res_ell_2)/50, colSums(res_ell_1)/50, 
                                 colSums(res_reg)/50), 
                  "Frequency_of_H" = rep(seq(from = 0, to = 40, length.out = 11), 3), 
                  "SD" = c(apply(res_ell_2, 2, sd), apply(res_ell_1, 2, sd), apply(res_reg, 2, sd)),
                  "Exp" = c(rep("ell_2", 11), rep("ell_1", 11), rep("lin_reg", 11)))

g <- ggplot(data = tmp, aes(x = factor(tmp$Frequency_of_H), y = Accuracy, color = Exp)) + 
  geom_errorbar(aes(ymin = Accuracy - SD, ymax=Accuracy + SD), width=.5, position = "dodge") + 
  geom_point(size = 2, position = position_dodge(width = 0.5))  + 
  xlab(expression(paste("Frequency of H", ", " ,omega[H]))) + 
  ylab("Prob. of recov. the causal parent") + 
  theme_minimal() + 
  theme(axis.title = element_text(size = 12), legend.position = c(0.795,0.82), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  scale_color_discrete(labels = c(expression(paste("l"[1])),expression(paste("l"[2])),"LM" ), 
                       guide = guide_legend(direction = "horizontal", title.position = "top", title.hjust = 0.5)) 
g$labels$colour <- "Method"

ggsave(filename = "Graphics/random_graph_PA_Y_1.pdf", plot = g, width = 5.8, height = 3)

