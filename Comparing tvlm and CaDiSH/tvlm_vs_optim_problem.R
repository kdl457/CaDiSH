#### Comparing optim problem to tvlm - an example.

setwd("~/Dropbox/KU/Thesis/R")
source(file = "Functions/source_file.R")
library(tvReg)
library(reshape2)
library(ggplot2)

set.seed(1)
n_sim <- 200
data <- sim_data(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_2pi, 
                 hidden_y = 10, y_var = TRUE)
data <- data.frame(data)
data$Y <- (data$Y - mean(data$Y)) / sd(data$Y)
data$X1 <- (data$X1 - mean(data$X1)) / sd(data$X1)
data$X2 <- (data$X2 - mean(data$X2)) / sd(data$X2)
data$X3 <- (data$X3 - mean(data$X3)) / sd(data$X3)

tv_model_1 <- tvLM(formula = Y~X1, data = data)
tv_model_2 <- tvLM(formula = Y~X2, data = data)
tv_model_3 <- tvLM(formula = Y~X3, data = data)

#### Choosing lambda

# lambda_LOOCV <- choose_lambda(data = data, lambda_val = seq(20, 60, length.out = 20), n_pred = 3)
# lambda_LOOCV$LOOCV_val 40

res_X1 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
      data = data[,c(1,2)], pen_fun = beta_pen_2_sq, lambda = 40,
      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
      control = list(maxit = 500))

res_X2 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                data = data[,c(1,3)], pen_fun = beta_pen_2_sq, lambda = 40,
                gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                control = list(maxit = 500))

res_X3 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                data = data[,c(1,4)], pen_fun = beta_pen_2_sq, lambda = 40,
                gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                control = list(maxit = 500))

X1_coeff <- rep(-1, n_sim)
X1_coeff[((n_sim/4) + 1):(3*n_sim/4)] <- 1

X3_coeff <- data$h
X3_coeff[((n_sim/4) + 1):(n_sim/2)] <- -X3_coeff[((n_sim/4) + 1):(n_sim/2)]
X3_coeff[((n_sim/2) + 1):((3*n_sim)/4)] <- 3*X3_coeff[((n_sim/2) + 1):((3*n_sim)/4)]

plot_data_frame <- data.frame("t" = data$t, "X1_tvlm_estimate" =  tv_model_1$tvcoef[,2], 
                              "X2_tvlm_estimate" =  tv_model_2$tvcoef[,2],
                              "X3_tvlm_estimate" =  tv_model_3$tvcoef[,2], 
                              "X1_obj_estimate" = res_X1$par[(n_sim + 1): (2*n_sim)],
                              "X2_obj_estimate" = res_X2$par[(n_sim + 1): (2*n_sim)],
                              "X3_obj_estimate" = res_X3$par[(n_sim + 1): (2*n_sim)],
                              "X1_true_coefficient" = X1_coeff, 
                              "X2_true_coefficient" = data$h, 
                              "X3_true_coefficient" = X3_coeff)

sum((plot_data_frame$X1_tvlm_estimate - plot_data_frame$X1_true_coefficient)^2) # 277.87
sum((plot_data_frame$X1_obj_estimate - plot_data_frame$X1_true_coefficient)^2) # 13.79

temp_hat_1 <- paste("SSR(hat(beta^{1})(t)) ==", 278)
temp_tilde_1 <- paste("SSR(tilde(beta^{1})(t)) ==", 13.8)

plot_data_frame_1 <- melt(plot_data_frame[,c(1,2,5,8)], id = "t", variable.name = "Method")
p_1 <- ggplot(data = plot_data_frame_1, aes(x = t, y = value, color = Method)) + geom_line() + 
  labs(y = expression(paste("Estimates of", " ",  beta^{1}, (t))), x = "Time",  fill = "Estimates") + 
  scale_fill_discrete(labels = c("tvlm estimate", "Objective estimate", "True coefficient")) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  annotate("text", x = 0.9, y = 0.85, label = temp_hat_1, cex = 3, parse = TRUE) + 
  annotate("text", x = 0.9, y = 0.60, label = temp_tilde_1, cex = 3, parse = TRUE) + 
  scale_color_discrete(labels = c(expression(paste(hat(beta^{1}), (t))),
                                  expression(paste(tilde(beta^{1}), (t))),
                                  expression(paste(beta^{1}, (t)))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave(filename = "Graphics/tvlm_obj_est_X1.pdf", plot = p_1, width = 5.8, height = 3)

sum((plot_data_frame$X2_tvlm_estimate - plot_data_frame$X2_true_coefficient)^2) # 0.97
sum((plot_data_frame$X2_obj_estimate - plot_data_frame$X2_true_coefficient)^2) # 0.32

temp_hat_2 <- paste("SSR(hat(beta^{2})(t)) ==", 0.97)
temp_tilde_2 <- paste("SSR(tilde(beta^{2})(t)) ==", 0.32)

plot_data_frame_2 <- melt(plot_data_frame[,c(1,3,6,9)], id = "t", variable.name = "Method")

p_2 <- ggplot(data = plot_data_frame_2, aes(x = t, y = value, color = Method)) + geom_line() + 
  labs(y = expression(paste("Estimates of", " ",  beta^{2}, (t))), x = "Time",  fill = "Estimates") + 
  scale_fill_discrete(labels = c("tvlm estimate", "Objective estimate", "True coefficient")) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  annotate("text", x = 0.9, y = 0.85, label = temp_hat_2, cex = 3, parse = TRUE) + 
  annotate("text", x = 0.9, y = 0.60, label = temp_tilde_2, cex = 3, parse = TRUE) + 
  scale_color_discrete(labels = c(expression(paste(hat(beta^{2}), (t))),
                                  expression(paste(tilde(beta^{2}), (t))),
                                  expression(paste(beta^{2}, (t)))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave(filename = "Graphics/tvlm_obj_est_X2.pdf", plot = p_2, width = 5.8, height = 3)


sum((plot_data_frame$X3_tvlm_estimate - plot_data_frame$X3_true_coefficient)^2) # 107.7
sum((plot_data_frame$X3_obj_estimate - plot_data_frame$X3_true_coefficient)^2) # 99.7

temp_hat_3 <- paste("SSR(hat(beta^{3})(t)) ==", 108)
temp_tilde_3 <- paste("SSR(tilde(beta^{3})(t)) ==", 100)

plot_data_frame_3 <- melt(plot_data_frame[,c(1,4,7,10)], id = "t", variable.name = "Method")
p_3 <- ggplot(data = plot_data_frame_3, aes(x = t, y = value, color = Method)) + geom_line() + 
  labs(y = expression(paste("Estimates of", " ",  beta^{3}, (t))), x = "Time",  fill = "Estimates") + 
  scale_fill_discrete(labels = c("tvlm estimate", "Objective estimate", "True coefficient")) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 12), 
        legend.background = element_rect(linetype = 1, size = 0.5, colour = 1)) + 
  annotate("text", x = 0.9, y = 0.85, label = temp_hat_3, cex = 3, parse = TRUE) + 
  annotate("text", x = 0.9, y = 0.40, label = temp_tilde_3, cex = 3, parse = TRUE) + 
  scale_color_discrete(labels = c(expression(paste(hat(beta^{3}), (t))),
                                  expression(paste(tilde(beta^{3}), (t))),
                                  expression(paste(beta^{3}, (t)))), 
                       guide = guide_legend(title.position = "top", title.hjust = 0.5))

ggsave(filename = "Graphics/tvlm_obj_est_X3.pdf", plot = p_3, width = 5.8, height = 3)

######### For X1:
#### Showing how SSR of estimates and true coefficients changes with choice of bandwidth

h <- seq(from = 0.03, to = 0.1, length.out = 100)

SSR_h <- c()
for (i in 1:length(h)) {
  SSR_h[i] <- sum((tvLM(formula = Y~X1, data = data, bw = h[i])$tvcoef[,2] - plot_data_frame$X1_true_coefficient)^2)
}

#### Showing SSR of estimates and true coefficients changes with choice of lambda

lambda <- seq(from = 0.01, to = 80, length.out = 400)

SSR_lambda <- c()
for (i in 1:length(lambda)) {
  res <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                  data = data[,c(1,2)], pen_fun = beta_pen_2_sq, lambda = lambda[i],
                  gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                  control = list(maxit = 500))$par
  res <- res[(n_sim + 1): (2*n_sim)]
  SSR_lambda[i] <- sum((res - plot_data_frame$X1_true_coefficient)^2)
}

SSR_data_1 <- data.frame("SSR_lambda" = SSR_lambda, "SSR_h" = SSR_h, "Lambda" = lambda, "h" = h)
library(ggplot2)
library(gridExtra)

p_SSR_h_1 <- ggplot(data = SSR_data_1, aes(x = h, y = SSR_h)) + geom_point(color = "red", size = 0.5) + 
  labs(y = expression(paste("SSR"[h])), x = "Bandwidth h") + ylim(c(10,40)) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 12)) + 
  geom_hline(yintercept =  min(SSR_h), color = "black", size = 0.7, linetype = "dashed")

p_SSR_lambda_1 <- ggplot(data = SSR_data_1, aes(x = lambda, y = SSR_lambda)) + geom_point(color = "blue", size = 0.5) + 
  labs(y = expression(paste("SSR"[lambda])), x = expression(paste("Roughness penalty", " ", lambda))) + 
  theme_minimal() + ylim(c(10,40)) + 
  theme(axis.title = element_text(size = 12)) + 
  geom_hline(yintercept =  min(SSR_h), color = "black", size = 0.7, linetype = "dashed")

g_1 <- grid.arrange(p_SSR_h_1, p_SSR_lambda_1, nrow = 1)

ggsave(filename = "Graphics/SSR_h_lambda_X1.pdf", plot = g_1, width = 5.8, height = 3)


#### For X2:

SSR_h_2 <- c()
for (i in 1:length(h)) {
  SSR_h_2[i] <- sum((tvLM(formula = Y~X2, data = data, bw = h[i])$tvcoef[,2] - plot_data_frame$X2_true_coefficient)^2)
}

SSR_lambda_2 <- c()
for (i in 1:length(lambda)) {
  res <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
               data = data[,c(1,3)], pen_fun = beta_pen_2_sq, lambda = lambda[i],
               gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
               control = list(maxit = 500))$par
  res <- res[(n_sim + 1): (2*n_sim)]
  SSR_lambda_2[i] <- sum((res - plot_data_frame$X2_true_coefficient)^2)
}

SSR_data_2 <- data.frame("SSR_lambda" = SSR_lambda_2, "SSR_h" = SSR_h_2, "Lambda" = lambda, "h" = h)

p_SSR_h_2 <- ggplot(data = SSR_data_2, aes(x = h, y = SSR_h)) + geom_point(color = "red", size = 0.5) + 
  labs(y = expression(paste("SSR"[h])), x = "Bandwidth h") + ylim(c(0.3,1.5)) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 12)) + 
  geom_hline(yintercept =  min(SSR_h_2), color = "black", size = 0.7, linetype = "dashed")

p_SSR_lambda_2 <- ggplot(data = SSR_data_2, aes(x = lambda, y = SSR_lambda)) + geom_point(color = "blue", size = 0.5) + 
  labs(y = expression(paste("SSR"[lambda])), x = expression(paste("Roughness penalty", " ", lambda))) + 
  theme_minimal() + ylim(c(0.3,1.5)) + 
  theme(axis.title = element_text(size = 12)) + 
  geom_hline(yintercept =  min(SSR_h_2), color = "black", size = 0.7, linetype = "dashed")

g_2 <- grid.arrange(p_SSR_h_2, p_SSR_lambda_2, nrow = 1)

ggsave(filename = "Graphics/SSR_h_lambda_X2.pdf", plot = g_2, width = 5.8, height = 3)
#### For X3:

SSR_h_3 <- c()
for (i in 1:length(h)) {
  SSR_h_3[i] <- sum((tvLM(formula = Y~X3, data = data, bw = h[i])$tvcoef[,2] - plot_data_frame$X3_true_coefficient)^2)
}

SSR_lambda_3 <- c()
for (i in 1:length(lambda)) {
  res <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
               data = data[,c(1,4)], pen_fun = beta_pen_2_sq, lambda = lambda[i],
               gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
               control = list(maxit = 500))$par
  res <- res[(n_sim + 1): (2*n_sim)]
  SSR_lambda_3[i] <- sum((res - plot_data_frame$X3_true_coefficient)^2)
}

SSR_data_3 <- data.frame("SSR_lambda" = SSR_lambda_3, "SSR_h" = SSR_h_3, "Lambda" = lambda, "h" = h)

p_SSR_h_3 <- ggplot(data = SSR_data_3, aes(x = h, y = SSR_h)) + geom_point(color = "red", size = 0.5) + 
  labs(y = expression(paste("SSR"[h])), x = "Bandwidth h") + ylim(c(95,130)) + 
  theme_minimal() +
  theme(axis.title = element_text(size = 12)) + 
  geom_hline(yintercept =  min(SSR_h_3), color = "black", size = 0.7, linetype = "dashed")

p_SSR_lambda_3 <- ggplot(data = SSR_data_3, aes(x = lambda, y = SSR_lambda)) + geom_point(color = "blue", size = 0.5) + 
  labs(y = expression(paste("SSR"[lambda])), x = expression(paste("Roughness penalty", " ", lambda))) + 
  theme_minimal() + ylim(c(95,130)) + 
  theme(axis.title = element_text(size = 12)) + 
  geom_hline(yintercept =  min(SSR_h_3), color = "black", size = 0.7, linetype = "dashed")

g_3 <- grid.arrange(p_SSR_h_3, p_SSR_lambda_3, nrow = 1)

ggsave(filename = "Graphics/SSR_h_lambda_X3.pdf", plot = g_3, width = 5.8, height = 3)
