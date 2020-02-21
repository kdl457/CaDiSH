##### Simulation from SCMs with fixed coefficients in each environment

setwd("~/Dropbox/KU/Thesis/R")
source(file = "Functions/source_file.R")
library(tvReg)

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 40 # Fixed lambda. LOOCV could be used but takes way too long. Maybe do it for the first.
n_sim <- 200 # Number of observations from each SCM
h_y_vec <- seq(0.1, 10, length.out = 19) # Values of h_y - large h_y = large effect of H on beta_2 in Y
res_matrix_fixed_200 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
for (i in 1:B) {
  cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
  data <- sim_data_fixed_env(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_2pi, 
                             alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                             alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                             c_23 = coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)], h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                             n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                             hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                             alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                             c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                             c_1y = coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)], 
                             n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)])
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
  
  # Causal discovery
  parent_res[i] <- which.min(c(res_X1$value, res_X2$value, res_X3$value))
}
CD_accuracy[j] <- sum(parent_res == 2)/B
}
res_matrix_fixed_200[,h_y] <- CD_accuracy
}

res_matrix_fixed_200

saveRDS(object =  res_matrix_fixed_200, file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_200.RDS")

###### Compare with n_sim = 50

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 40 # Fixed lambda. LOOCV could be used but takes way too long. Maybe do it for the first.
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0.1, 10, length.out = 19) # Values of h_y - large h_y = large effect of H on beta_2 in Y
res_matrix_fixed_50 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_fixed_env(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_2pi, 
                                 alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                                 alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                                 c_23 = coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)], h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                                 n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                                 hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                                 alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                                 c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                                 c_1y = coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)], 
                                 n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)])
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
      
      # Causal discovery
      parent_res[i] <- which.min(c(res_X1$value, res_X2$value, res_X3$value))
    }
    CD_accuracy[j] <- sum(parent_res == 2)/B
  }
  res_matrix_fixed_50[,h_y] <- CD_accuracy
}

res_matrix_fixed_50

saveRDS(object =  res_matrix_fixed_50, file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_50.RDS")

####### Compare beta_pen_sq

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 0.04 # Fixed lambda. LOOCV could be used but takes way too long. Maybe do it for the first.
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0.1, 10, length.out = 19) # Values of h_y - large h_y = large effect of H on beta_2 in Y
res_matrix_fixed_50_1 <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_fixed_env(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_2pi, 
                                 alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                                 alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                                 c_23 = coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)], h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                                 n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                                 hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                                 alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                                 c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                                 c_1y = coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)], 
                                 n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)])
      data <- data.frame(data)
      data$Y <- (data$Y - mean(data$Y)) / sd(data$Y)
      data$X1 <- (data$X1 - mean(data$X1)) / sd(data$X1)
      data$X2 <- (data$X2 - mean(data$X2)) / sd(data$X2)
      data$X3 <- (data$X3 - mean(data$X3)) / sd(data$X3)
      res_X1 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,2)], pen_fun = beta_pen_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      res_X2 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,3)], pen_fun = beta_pen_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      res_X3 <- optim(par = rep(0, 2*n_sim), fn = obj_fun_1,
                      data = data[,c(1,4)], pen_fun = beta_pen_sq, lambda = lambda,
                      gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                      control = list(maxit = 500))
      
      # Causal discovery
      parent_res[i] <- which.min(c(res_X1$value, res_X2$value, res_X3$value))
    }
    CD_accuracy[j] <- sum(parent_res == 2)/B
  }
  res_matrix_fixed_50_1[,h_y] <- CD_accuracy
}

res_matrix_fixed_50_1

saveRDS(object =  res_matrix_fixed_50_1, file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_50_ell_1.RDS")


###### n_sim = 50, f_hidden = sin(10 * pi * t)

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 40 # Fixed lambda. LOOCV could be used but takes way too long. Maybe do it for the first.
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0.1, 10, length.out = 19) # Values of h_y - large h_y = large effect of H on beta_2 in Y
res_matrix_fixed_50_10pi <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_fixed_env(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_10pi, 
                                 alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                                 alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                                 c_23 = coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)], h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                                 n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                                 hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                                 alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                                 c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                                 c_1y = coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)], 
                                 n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)])
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
      
      # Causal discovery
      parent_res[i] <- which.min(c(res_X1$value, res_X2$value, res_X3$value))
    }
    CD_accuracy[j] <- sum(parent_res == 2)/B
  }
  res_matrix_fixed_50_10pi[,h_y] <- CD_accuracy
}

res_matrix_fixed_50_10pi

saveRDS(object =  res_matrix_fixed_50_10pi, file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_50_10pi.RDS")

###### Compare with linear regression

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 40 # Fixed lambda. LOOCV could be used but takes way too long. Maybe do it for the first.
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0.1, 10, length.out = 19) # Values of h_y - large h_y = large effect of H on beta_2 in Y
res_matrix_fixed_reg <- matrix(NA, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  coeff_var_mat <- matrix(runif(n = 17 * B * n_CI, min = - 5, max = 5), nrow = 17 * B)
  h_var_mat <- matrix(runif(n = 4 * B * n_CI, min = - 10, max = 10), nrow = 4 * B)   
  noise_var_mat <- matrix(runif(n = 8 * B * n_CI, min = 0, max = 5), nrow = 8 * B)  
  for (j in 1:n_CI) {
    coeff_var <- coeff_var_mat[,j]
    h_var <- h_var_mat[,j]
    noise_var <- noise_var_mat[,j]
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_fixed_env(sim_length = n_sim, breaks = c(0, 1/4, 1/2, 3/4, 1), h_func = f_hidden_sin_2pi, 
                                 alpha_3 = coeff_var[17*(i-1) + 1], h_3 = h_var[4*(i-1) + 1], n_3 = noise_var[8*(i-1) + 1], 
                                 alpha_2 = coeff_var[(17*(i-1) + 2):(17*(i-1) + 4)], 
                                 c_23 = coeff_var[(17*(i-1) + 5):(17*(i-1) + 7)], h_2 = h_var[(4*(i-1) + 2):(4*(i-1) + 4)], 
                                 n_2 = noise_var[(8*(i-1) + 2):(8*(i-1) + 4)], alpha_y = coeff_var[17*(i-1) + 8], 
                                 hidden_y = h_y_vec[h_y], n_y = noise_var[8*(i-1) + 5], 
                                 alpha_1 = coeff_var[(17*(i-1) + 9):(17*(i-1) + 11)], 
                                 c_13 = coeff_var[(17*(i-1) + 12):(17*(i-1) + 14)], 
                                 c_1y = coeff_var[(17*(i-1) + 15):(17*(i-1) + 17)], 
                                 n_1 = noise_var[(8*(i-1) + 6):(8*(i-1) + 8)])
      data <- data.frame(data)
      
      # Causal discovery
      parent_res[i] <- which.min(c(summary(lm(data$Y ~ data$X1))$coefficients[,4][2],
                                   summary(lm(data$Y ~ data$X2))$coefficients[,4][2],
                                   summary(lm(data$Y ~ data$X3))$coefficients[,4][2]))
    }
    CD_accuracy[j] <- sum(parent_res == 2)/B
  }
  res_matrix_fixed_reg[,h_y] <- CD_accuracy
}

res_matrix_fixed_reg

saveRDS(object =  res_matrix_fixed_reg, file = "Simulations/Causal Discovery/Plot1/res_matrix_fixed_reg.RDS")