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

plot_data_frame_1 <- melt(plot_data_frame[,c(1,2,5,8)], id = "t", variable.name = "Method")
p_1 <- ggplot(data = plot_data_frame_1, aes(x = t, y = value, color = Method)) + geom_line() + 
  ggtitle("Comparison of tvlm and our CaDiSH method on X1") + xlab("Time") + 
  ylab("Regression coefficient for X1") + labs(fill = "Estimates") + 
  scale_fill_discrete(labels = c("tvlm estimate", "Objective estimate", "True coefficient")) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 30), axis.title = element_text(size = 22), 
        legend.text = element_text(size = 14), legend.title = element_text(size = 14)) +
  annotate("text", x = 0.9, y = 0.85, label = "SSR score for tvlm = 277.9", cex = 4.5) + 
  annotate("text", x = 0.9, y = 0.75, label = "SSR score for our method = 13.8", cex = 4.5)

ggsave(filename = "Graphics/tvlm_obj_est_X1.pdf", plot = p_1, width = 7, height = 4, dpi = 72)

sum((plot_data_frame$X2_tvlm_estimate - plot_data_frame$X2_true_coefficient)^2) # 0.97
sum((plot_data_frame$X2_obj_estimate - plot_data_frame$X2_true_coefficient)^2) # 0.32

plot_data_frame_2 <- melt(plot_data_frame[,c(1,3,6,9)], id = "t", variable.name = "Method")
ggplot(data = plot_data_frame_2, aes(x = t, y = value, color = Method)) + geom_line() + 
  ggtitle("Comparison of tvlm and our CaDiSH method on X2") + xlab("Time") +
  ylab("Regression coefficient for X2") + 
  scale_fill_discrete(labels = c("tvlm estimate", "Objective estimate", "True coefficient")) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  annotate("text", x = 0.9, y = 0.85, label = "SSR score for tvlm = 0.97") + 
  annotate("text", x = 0.9, y = 0.75, label = "SSR score for our method = 0.32")

sum((plot_data_frame$X3_tvlm_estimate - plot_data_frame$X3_true_coefficient)^2) # 107.7
sum((plot_data_frame$X3_obj_estimate - plot_data_frame$X3_true_coefficient)^2) # 99.7

plot_data_frame_3 <- melt(plot_data_frame[,c(1,4,7,10)], id = "t", variable.name = "Method")
ggplot(data = plot_data_frame_3, aes(x = t, y = value, color = Method)) + geom_line() + 
  ggtitle("Comparison of tvlm and our CaDiSH method on X3") + xlab("time") + 
  ylab("Regression coefficient for X1") + 
  scale_fill_discrete(labels = c("tvlm estimate", "Objective estimate", "True coefficient")) + 
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, size = 22), axis.title = element_text(size = 14)) +
  annotate("text", x = 0.9, y = 0.85, label = "SSR score for tvlm = 107.7") + 
  annotate("text", x = 0.9, y = 0.70, label = "SSR score for our method = 99.7")

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

par(mfrow = c(1, 2))
plot(h, SSR_h, pch = 19, col = "red", xlab = "bandwidth h", ylab = "SSR value of tvlm estimate", 
     cex.lab=1.5, cex.axis=1.4, cex.main=1.5, cex.sub=1.5)
abline(h = min(SSR_h), col = "red", lwd = 3)
plot(lambda, SSR_lambda, pch = 19, col = "blue", xlab = "lambda", 
     ylab = "SSR value of objective function estimate", 
     cex.lab=1.5, cex.axis=1.4, cex.main=1.5, cex.sub=1.5)
abline(h = min(SSR_h), col = "red", lwd = 3)
mtext("Comparing tvlm to our novel method using X1 as predictor", side = 3, line = -2, outer = TRUE, cex = 2)


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

par(mfrow = c(1, 2))
plot(h, SSR_h_2, pch = 19, col = "red", xlab = "bandwidth h", ylab = "SSR value of tvlm estimate", 
     cex.lab=1.5, cex.axis=1.4, cex.main=1.5, cex.sub=1.5)
abline(h = min(SSR_h_2), col = "red", lwd = 3)
plot(lambda, SSR_lambda_2, pch = 19, col = "blue", xlab = "lambda", 
     ylab = "SSR value of objective function estimate", ylim = c(0, min(SSR_h_2) + 0.05), 
     cex.lab=1.5, cex.axis=1.4, cex.main=1.5, cex.sub=1.5)
abline(h = min(SSR_h_2), col = "red", lwd = 3)
mtext("Comparing tvlm to our novel method using X2 as predictor", side = 3, line = -2, outer = TRUE, cex = 2)

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

par(mfrow = c(1, 2))
plot(h, SSR_h_3, pch = 19, col = "red", xlab = "bandwidth h", ylab = "SSR value of tvlm estimate", 
     cex.lab=1.5, cex.axis=1.4, cex.main=1.5, cex.sub=1.5)
abline(h = min(SSR_h_3), col = "red", lwd = 3)
plot(lambda, SSR_lambda_3, ylim = c(97.7, 107), pch = 19, col = "blue", xlab = "lambda", 
     ylab = "SSR value of objective function estimate", 
     cex.lab=1.5, cex.axis=1.4, cex.main=1.5, cex.sub=1.5)
abline(h = min(SSR_h_3), col = "red", lwd = 3)
mtext("Comparing tvlm to our novel method using X3 as predictor", side = 3, line = -2, outer = TRUE, cex = 2)


