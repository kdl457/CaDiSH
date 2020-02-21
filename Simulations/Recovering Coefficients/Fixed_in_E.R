###### Simulations for recovering time-varying coefficients
setwd("~/Dropbox/KU/Thesis/R")
source(file = "Functions/source_file.R")
library(tvReg)

set.seed(1)
SSR_res <- c()
SSR_res_1 <- c()
SSR_res_2 <- c()
SSR_res_3 <- c()
SSR_h_1_vec <- c()
SSR_h_2_vec <- c()
SSR_h_3_vec <- c()
SSR_lambda_1_vec <- c()
SSR_lambda_2_vec <- c()
SSR_lambda_3_vec <- c()
SSR_h_list <- list()
SSR_lambda_list <- list()
SSR_h_res <- list()
SSR_lambda_res <- list()
coeff_accuracy <- c()
coeff_acc_1 <- c()
coeff_acc_2 <- c()
coeff_acc_3 <- c()
B <- 20 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 1.2 # Fixed lambda. Found by LOOCV with h_Y = 5.
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
coeff_matrix_fixed_50 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))
coeff_matrix_fixed_1 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))
coeff_matrix_fixed_2 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))
coeff_matrix_fixed_3 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
    for (i in 1:B) {
      c_1y <- coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)]
      c_23 <- coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)]
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_fixed_env_2(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_2pi, 
                                 alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                                 alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                                 c_23 = c_23, h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                                 n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                                 hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                                 alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                                 c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                                 c_1y = c_1y, 
                                 n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)], c_y2 = 5)
      data <- data.frame(data)
      data$Y <- (data$Y - mean(data$Y)) / sd(data$Y)
      data$X1 <- (data$X1 - mean(data$X1)) / sd(data$X1)
      data$X2 <- (data$X2 - mean(data$X2)) / sd(data$X2)
      data$X3 <- (data$X3 - mean(data$X3)) / sd(data$X3)
      res_X1 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,2)], pen_fun = beta_pen_2_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      res_X2 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,3)], pen_fun = beta_pen_2_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      res_X3 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,4)], pen_fun = beta_pen_2_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      
      # Recovering time-varying coefficients - versus existing tvLM method
      # Fitting tvLM models
      X1_tvlm_estimate <- tvLM(formula = Y~X1, data = data)$tvcoef[,2]
      X2_tvlm_estimate <- tvLM(formula = Y~X2, data = data)$tvcoef[,2]
      X3_tvlm_estimate <- tvLM(formula = Y~X3, data = data)$tvcoef[,2]
      
      # Calculating true coefficients
      X1_coeff <- rep(1, n_sim)
      X1_coeff[1:(n_sim/4)] <- 1 / (c_1y[1])
      X1_coeff[(3 * n_sim/4) : n_sim] <- 1 / (c_1y[1])
      X1_coeff[((n_sim/4) + 1):(n_sim/2)] <- 1 / (c_1y[2])
      X1_coeff[((n_sim/2) + 1):(3 * n_sim/4)] <- 1 / (c_1y[3])
      
      X2_coeff <- 5 + h_y_vec[h_y] * data$h
      
      X3_coeff <- 5 + h_y_vec[h_y] * data$h
      X3_coeff[1:(n_sim/4)] <- c_23[1] * X3_coeff[1:(n_sim/4)]
      X3_coeff[(3 * n_sim/4) : n_sim] <- c_23[1] * X3_coeff[(3 * n_sim/4) : n_sim]
      X3_coeff[((n_sim/4) + 1):(n_sim/2)] <- c_23[2] * X3_coeff[((n_sim/4) + 1):(n_sim/2)]
      X3_coeff[((n_sim/2) + 1):(3 * n_sim/4)] <- c_23[3] * X3_coeff[((n_sim/2) + 1):(3 * n_sim/4)]
      
      # Calcuating SSR values and determining if tvlm or our method has the lowest for each predictor
      
      SSR_h_1 <- sum((X1_tvlm_estimate - X1_coeff)^2)
      SSR_h_2 <- sum((X2_tvlm_estimate - X2_coeff)^2)
      SSR_h_3 <- sum((X3_tvlm_estimate - X3_coeff)^2)
      SSR_lambda_1 <- sum((res_X1$par[(n_sim + 1): (2*n_sim)] - X1_coeff)^2)
      SSR_lambda_2 <- sum((res_X2$par[(n_sim + 1): (2*n_sim)] - X2_coeff)^2)
      SSR_lambda_3 <- sum((res_X3$par[(n_sim + 1): (2*n_sim)] - X3_coeff)^2)
      
      SSR_res[3 * i - 2] <- which.min(c(SSR_h_1, SSR_lambda_1))
      SSR_res_1[i] <- SSR_res[3 * i - 2]
      SSR_res[3 * i - 1] <- which.min(c(SSR_h_2, SSR_lambda_2))
      SSR_res_2[i] <- SSR_res[3 * i - 1]
      SSR_res[3 * i] <- which.min(c(SSR_h_3, SSR_lambda_3))
      SSR_res_3[i] <- SSR_res[3 * i]
      
      SSR_h_1_vec[i] <- SSR_h_1
      SSR_h_2_vec[i] <- SSR_h_2
      SSR_h_3_vec[i] <- SSR_h_3
      SSR_lambda_1_vec[i] <- SSR_lambda_1
      SSR_lambda_2_vec[i] <- SSR_lambda_2
      SSR_lambda_3_vec[i] <- SSR_lambda_3
      
    }
    SSR_h_list[[j]] <- matrix(c(SSR_h_1_vec, SSR_h_2_vec, SSR_h_3_vec), nrow = B)
    SSR_lambda_list[[j]] <- matrix(c(SSR_lambda_1_vec, SSR_lambda_2_vec, SSR_lambda_3_vec), nrow = B)
    coeff_accuracy[j] <- sum(SSR_res == 2)/(3*B)
    coeff_acc_1[j] <- sum(SSR_res_1 == 2)/(B)
    coeff_acc_2[j] <- sum(SSR_res_2 == 2)/(B)
    coeff_acc_3[j] <- sum(SSR_res_3 == 2)/(B)
  }
  SSR_h_res[[h_y]] <- SSR_h_list
  SSR_lambda_res[[h_y]] <- SSR_lambda_list
  coeff_matrix_fixed_50[,h_y] <- coeff_accuracy
  coeff_matrix_fixed_1[,h_y] <- coeff_acc_1
  coeff_matrix_fixed_2[,h_y] <- coeff_acc_2
  coeff_matrix_fixed_3[,h_y] <- coeff_acc_3
}

saveRDS(object = coeff_matrix_fixed_50, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_50_new.RDS")
saveRDS(object = coeff_matrix_fixed_1, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_1_new.RDS")
saveRDS(object = coeff_matrix_fixed_2, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_2_new.RDS")
saveRDS(object = coeff_matrix_fixed_3, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_3_new.RDS")
saveRDS(object = SSR_h_res, file = "Simulations/Recovering Coefficients/SSR_h_res_sin.RDS")
saveRDS(object = SSR_lambda_res, file = "Simulations/Recovering Coefficients/SSR_lambda_res_sin.RDS")

##### Repeat with H(t) = log(t + 0.1)

set.seed(1)
SSR_res <- c()
SSR_res_1 <- c()
SSR_res_2 <- c()
SSR_res_3 <- c()
SSR_h_1_vec <- c()
SSR_h_2_vec <- c()
SSR_h_3_vec <- c()
SSR_lambda_1_vec <- c()
SSR_lambda_2_vec <- c()
SSR_lambda_3_vec <- c()
SSR_h_list <- list()
SSR_lambda_list <- list()
SSR_h_res_log <- list()
SSR_lambda_res_log <- list()
coeff_accuracy <- c()
coeff_acc_1 <- c()
coeff_acc_2 <- c()
coeff_acc_3 <- c()
B <- 20 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 1.2 # Fixed lambda. Found by LOOCV with h_Y = 5.
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
coeff_matrix_fixed_50_log <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))
coeff_matrix_fixed_1_log <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))
coeff_matrix_fixed_2_log <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))
coeff_matrix_fixed_3_log <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
    for (i in 1:B) {
      c_1y <- coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)]
      c_23 <- coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)]
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_fixed_env_2(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_log, 
                                   alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                                   alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                                   c_23 = c_23, h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                                   n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                                   hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                                   alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                                   c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                                   c_1y = c_1y, 
                                   n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)], c_y2 = 5)
      data <- data.frame(data)
      data$Y <- (data$Y - mean(data$Y)) / sd(data$Y)
      data$X1 <- (data$X1 - mean(data$X1)) / sd(data$X1)
      data$X2 <- (data$X2 - mean(data$X2)) / sd(data$X2)
      data$X3 <- (data$X3 - mean(data$X3)) / sd(data$X3)
      res_X1 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,2)], pen_fun = beta_pen_2_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      res_X2 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,3)], pen_fun = beta_pen_2_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      res_X3 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,4)], pen_fun = beta_pen_2_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      
      # Recovering time-varying coefficients - versus existing tvLM method
      # Fitting tvLM models
      X1_tvlm_estimate <- tvLM(formula = Y~X1, data = data)$tvcoef[,2]
      X2_tvlm_estimate <- tvLM(formula = Y~X2, data = data)$tvcoef[,2]
      X3_tvlm_estimate <- tvLM(formula = Y~X3, data = data)$tvcoef[,2]
      
      # Calculating true coefficients
      X1_coeff <- rep(1, n_sim)
      X1_coeff[1:(n_sim/4)] <- 1 / (c_1y[1])
      X1_coeff[(3 * n_sim/4) : n_sim] <- 1 / (c_1y[1])
      X1_coeff[((n_sim/4) + 1):(n_sim/2)] <- 1 / (c_1y[2])
      X1_coeff[((n_sim/2) + 1):(3 * n_sim/4)] <- 1 / (c_1y[3])
      
      X2_coeff <- 5 + h_y_vec[h_y] * data$h
      
      X3_coeff <- 5 + h_y_vec[h_y] * data$h
      X3_coeff[1:(n_sim/4)] <- c_23[1] * X3_coeff[1:(n_sim/4)]
      X3_coeff[(3 * n_sim/4) : n_sim] <- c_23[1] * X3_coeff[(3 * n_sim/4) : n_sim]
      X3_coeff[((n_sim/4) + 1):(n_sim/2)] <- c_23[2] * X3_coeff[((n_sim/4) + 1):(n_sim/2)]
      X3_coeff[((n_sim/2) + 1):(3 * n_sim/4)] <- c_23[3] * X3_coeff[((n_sim/2) + 1):(3 * n_sim/4)]
      
      # Calcuating SSR values and determining if tvlm or our method has the lowest for each predictor
      
      SSR_h_1 <- sum((X1_tvlm_estimate - X1_coeff)^2)
      SSR_h_2 <- sum((X2_tvlm_estimate - X2_coeff)^2)
      SSR_h_3 <- sum((X3_tvlm_estimate - X3_coeff)^2)
      SSR_lambda_1 <- sum((res_X1$par[(n_sim + 1): (2*n_sim)] - X1_coeff)^2)
      SSR_lambda_2 <- sum((res_X2$par[(n_sim + 1): (2*n_sim)] - X2_coeff)^2)
      SSR_lambda_3 <- sum((res_X3$par[(n_sim + 1): (2*n_sim)] - X3_coeff)^2)
      
      SSR_res[3 * i - 2] <- which.min(c(SSR_h_1, SSR_lambda_1))
      SSR_res_1[i] <- SSR_res[3 * i - 2]
      SSR_res[3 * i - 1] <- which.min(c(SSR_h_2, SSR_lambda_2))
      SSR_res_2[i] <- SSR_res[3 * i - 1]
      SSR_res[3 * i] <- which.min(c(SSR_h_3, SSR_lambda_3))
      SSR_res_3[i] <- SSR_res[3 * i]
      
      SSR_h_1_vec[i] <- SSR_h_1
      SSR_h_2_vec[i] <- SSR_h_2
      SSR_h_3_vec[i] <- SSR_h_3
      SSR_lambda_1_vec[i] <- SSR_lambda_1
      SSR_lambda_2_vec[i] <- SSR_lambda_2
      SSR_lambda_3_vec[i] <- SSR_lambda_3
      
    }
    SSR_h_list[[j]] <- matrix(c(SSR_h_1_vec, SSR_h_2_vec, SSR_h_3_vec), nrow = B)
    SSR_lambda_list[[j]] <- matrix(c(SSR_lambda_1_vec, SSR_lambda_2_vec, SSR_lambda_3_vec), nrow = B)
    coeff_accuracy[j] <- sum(SSR_res == 2)/(3*B)
    coeff_acc_1[j] <- sum(SSR_res_1 == 2)/(B)
    coeff_acc_2[j] <- sum(SSR_res_2 == 2)/(B)
    coeff_acc_3[j] <- sum(SSR_res_3 == 2)/(B)
  }
  SSR_h_res_log[[h_y]] <- SSR_h_list
  SSR_lambda_res_log[[h_y]] <- SSR_lambda_list
  coeff_matrix_fixed_50_log[,h_y] <- coeff_accuracy
  coeff_matrix_fixed_1_log[,h_y] <- coeff_acc_1
  coeff_matrix_fixed_2_log[,h_y] <- coeff_acc_2
  coeff_matrix_fixed_3_log[,h_y] <- coeff_acc_3
}

saveRDS(object = coeff_matrix_fixed_50_log, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_50_log_new.RDS")
saveRDS(object = coeff_matrix_fixed_1_log, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_1_log_new.RDS")
saveRDS(object = coeff_matrix_fixed_2_log, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_2_log_new.RDS")
saveRDS(object = coeff_matrix_fixed_3_log, file = "Simulations/Recovering Coefficients/coeff_matrix_fixed_3_log_new.RDS")
saveRDS(object = SSR_h_res_log, file = "Simulations/Recovering Coefficients/SSR_h_res_log.RDS")
saveRDS(object = SSR_lambda_res_log, file = "Simulations/Recovering Coefficients/SSR_lambda_res_log.RDS")
