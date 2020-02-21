##### Plot 2b: Varying the influence of H with continuous environments

setwd("~/Dropbox/KU/Thesis/R")
source(file = "Functions/source_file.R")

##########
########## n_sim = 50, beta_pen_2_sq, w_0 = 2*pi
##########

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 3.4 # Fixed lambda. This is found using LOOCV on 1 dataset
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
w_0 <- 2*pi # Period.
res_varying_cont_50 <- matrix(0, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  c_y2_var <- runif(B * n_CI, 0, - 0.5 + 10 + 10 - 0.5)
  for (k in 1:length(c_y2_var)) {
    if( c_y2_var[k] < (9.5) ){
      c_y2_var[k] <- -10 + c_y2_var[k]
    }else{
      c_y2_var[k] <- 0.5 + c_y2_var[k] - 9.5
    } 
  }
  c_y2_mat <- matrix(c_y2_var, nrow = B)
  
  for (j in 1:n_CI) {
    phases <- runif(n = 7 * B, -1, 1)
    coeff_var <- runif(n = 2 * B, min = - 5, max = 5)
    h_var <- runif(n = B, min = - 10, max = 10)
    noise_var <- runif(n = 4 * B, min = 0, max = 5)
    c_y2 <- c_y2_mat[,j]
    
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_varying_2(sim_length = n_sim, alpha = 5, w_0 = w_0, alpha_3 = coeff_var[2*(i-1) + 1], 
                               h_3 = h_var[i], n_3 = noise_var[4*(i-1) + 1], n_2 = noise_var[4*(i-1) + 2], 
                               alpha_y = coeff_var[2*(i-1) + 2], hidden_y = h_y_vec[h_y], 
                               n_y = noise_var[4*(i-1) + 3], n_1 = noise_var[4*(i-1) + 4], 
                               phases = phases[(7*(i-1) + 1) : (7*i)], c_y2 = c_y2[i])
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
  res_varying_cont_50[,h_y] <- CD_accuracy
}

res_varying_cont_50
saveRDS(object =  res_varying_cont_50, file = "Simulations/Causal Discovery/res_varying_cont_50.RDS")


##########
########## n_sim = 50, beta_pen_sq, w_0 = 2*pi
##########

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 0.1 # Fixed lambda. This is found using LOOCV on 1 dataset
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
w_0 <- 2*pi # Period.
res_varying_cont_50_ell1 <- matrix(0, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  c_y2_var <- runif(B * n_CI, 0, - 0.5 + 10 + 10 - 0.5)
  for (k in 1:length(c_y2_var)) {
    if( c_y2_var[k] < (9.5) ){
      c_y2_var[k] <- -10 + c_y2_var[k]
    }else{
      c_y2_var[k] <- 0.5 + c_y2_var[k] - 9.5
    } 
  }
  c_y2_mat <- matrix(c_y2_var, nrow = B)
  
  for (j in 1:n_CI) {
    phases <- runif(n = 7 * B, -1, 1)
    coeff_var <- runif(n = 2 * B, min = - 5, max = 5)
    h_var <- runif(n = B, min = - 10, max = 10)
    noise_var <- runif(n = 4 * B, min = 0, max = 5)
    c_y2 <- c_y2_mat[,j]
    
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_varying_2(sim_length = n_sim, alpha = 5, w_0 = w_0, alpha_3 = coeff_var[2*(i-1) + 1], 
                                 h_3 = h_var[i], n_3 = noise_var[4*(i-1) + 1], n_2 = noise_var[4*(i-1) + 2], 
                                 alpha_y = coeff_var[2*(i-1) + 2], hidden_y = h_y_vec[h_y], 
                                 n_y = noise_var[4*(i-1) + 3], n_1 = noise_var[4*(i-1) + 4], 
                                 phases = phases[(7*(i-1) + 1) : (7*i)], c_y2 = c_y2[i])
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
  res_varying_cont_50_ell1[,h_y] <- CD_accuracy
}

res_varying_cont_50_ell1
saveRDS(object =  res_varying_cont_50_ell1, file = "Simulations/Causal Discovery/res_varying_cont_50_ell1.RDS")

##########
########## n_sim = 50, beta_pen_2_sq, w_0 = 10*pi
##########

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 3.4 # Fixed lambda. This is found using LOOCV on 1 dataset
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
w_0 <- 10*pi # Period.
res_varying_cont_50_10pi <- matrix(0, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  c_y2_var <- runif(B * n_CI, 0, - 0.5 + 10 + 10 - 0.5)
  for (k in 1:length(c_y2_var)) {
    if( c_y2_var[k] < (9.5) ){
      c_y2_var[k] <- -10 + c_y2_var[k]
    }else{
      c_y2_var[k] <- 0.5 + c_y2_var[k] - 9.5
    } 
  }
  c_y2_mat <- matrix(c_y2_var, nrow = B)
  
  for (j in 1:n_CI) {
    phases <- runif(n = 7 * B, -1, 1)
    coeff_var <- runif(n = 2 * B, min = - 5, max = 5)
    h_var <- runif(n = B, min = - 10, max = 10)
    noise_var <- runif(n = 4 * B, min = 0, max = 5)
    c_y2 <- c_y2_mat[,j]
    
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_varying_2(sim_length = n_sim, alpha = 5, w_0 = w_0, alpha_3 = coeff_var[2*(i-1) + 1], 
                                 h_3 = h_var[i], n_3 = noise_var[4*(i-1) + 1], n_2 = noise_var[4*(i-1) + 2], 
                                 alpha_y = coeff_var[2*(i-1) + 2], hidden_y = h_y_vec[h_y], 
                                 n_y = noise_var[4*(i-1) + 3], n_1 = noise_var[4*(i-1) + 4], 
                                 phases = phases[(7*(i-1) + 1) : (7*i)], c_y2 = c_y2[i])
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
  res_varying_cont_50_10pi[,h_y] <- CD_accuracy
}

res_varying_cont_50_10pi
saveRDS(object =  res_varying_cont_50_10pi, file = "Simulations/Causal Discovery/res_varying_cont_50_10pi.RDS")

##########
########## n_sim = 200, beta_pen_2_sq, w_0 = 2*pi
##########

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 40 # Fixed lambda. This is found using LOOCV on 1 dataset
n_sim <- 200 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
w_0 <- 2*pi # Period.
res_varying_cont_200 <- matrix(0, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  c_y2_var <- runif(B * n_CI, 0, - 0.5 + 10 + 10 - 0.5)
  for (k in 1:length(c_y2_var)) {
    if( c_y2_var[k] < (9.5) ){
      c_y2_var[k] <- -10 + c_y2_var[k]
    }else{
      c_y2_var[k] <- 0.5 + c_y2_var[k] - 9.5
    } 
  }
  c_y2_mat <- matrix(c_y2_var, nrow = B)
  
  for (j in 1:n_CI) {
    phases <- runif(n = 7 * B, -1, 1)
    coeff_var <- runif(n = 2 * B, min = - 5, max = 5)
    h_var <- runif(n = B, min = - 10, max = 10)
    noise_var <- runif(n = 4 * B, min = 0, max = 5)
    c_y2 <- c_y2_mat[,j]
    
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_varying_2(sim_length = n_sim, alpha = 5, w_0 = w_0, alpha_3 = coeff_var[2*(i-1) + 1], 
                                 h_3 = h_var[i], n_3 = noise_var[4*(i-1) + 1], n_2 = noise_var[4*(i-1) + 2], 
                                 alpha_y = coeff_var[2*(i-1) + 2], hidden_y = h_y_vec[h_y], 
                                 n_y = noise_var[4*(i-1) + 3], n_1 = noise_var[4*(i-1) + 4], 
                                 phases = phases[(7*(i-1) + 1) : (7*i)], c_y2 = c_y2[i])
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
  res_varying_cont_200[,h_y] <- CD_accuracy
}

res_varying_cont_200
saveRDS(object =  res_varying_cont_200, file = "Simulations/Causal Discovery/res_varying_cont_200.RDS")

##########
########## n_sim = 50, w_0 = 2*pi, regression
##########

set.seed(1)
parent_res <- c()
CD_accuracy <- c()
B <- 50 # Number of simulations for each specific alpha
n_CI <- 50 # Number of repetitions of simulations for each specific alpha
lambda <- 3.4 # Fixed lambda. This is found using LOOCV on 1 dataset
n_sim <- 50 # Number of observations from each SCM
h_y_vec <- seq(0, 20, length.out = 11) # Values of h_y - large h_y = large effect of H on beta_2 in Y
w_0 <- 2*pi # Period.
res_varying_cont_reg <- matrix(0, nrow = n_CI, ncol = length(h_y_vec))

for (h_y in 1:length(h_y_vec)) {
  c_y2_var <- runif(B * n_CI, 0, - 0.5 + 10 + 10 - 0.5)
  for (k in 1:length(c_y2_var)) {
    if( c_y2_var[k] < (9.5) ){
      c_y2_var[k] <- -10 + c_y2_var[k]
    }else{
      c_y2_var[k] <- 0.5 + c_y2_var[k] - 9.5
    } 
  }
  c_y2_mat <- matrix(c_y2_var, nrow = B)
  
  for (j in 1:n_CI) {
    phases <- runif(n = 7 * B, -1, 1)
    coeff_var <- runif(n = 2 * B, min = - 5, max = 5)
    h_var <- runif(n = B, min = - 10, max = 10)
    noise_var <- runif(n = 4 * B, min = 0, max = 5)
    c_y2 <- c_y2_mat[,j]
    
    for (i in 1:B) {
      cat(paste("h_y", h_y_vec[h_y], "Repetition", j, "Number",i))
      data <- sim_data_varying_2(sim_length = n_sim, alpha = 5, w_0 = w_0, alpha_3 = coeff_var[2*(i-1) + 1], 
                                 h_3 = h_var[i], n_3 = noise_var[4*(i-1) + 1], n_2 = noise_var[4*(i-1) + 2], 
                                 alpha_y = coeff_var[2*(i-1) + 2], hidden_y = h_y_vec[h_y], 
                                 n_y = noise_var[4*(i-1) + 3], n_1 = noise_var[4*(i-1) + 4], 
                                 phases = phases[(7*(i-1) + 1) : (7*i)], c_y2 = c_y2[i])
      data <- data.frame(data)
      
      # Causal discovery
      parent_res[i] <- which.min(c(summary(lm(data$Y ~ data$X1))$coefficients[,4][2],
                                   summary(lm(data$Y ~ data$X2))$coefficients[,4][2],
                                   summary(lm(data$Y ~ data$X3))$coefficients[,4][2]))
    }
    CD_accuracy[j] <- sum(parent_res == 2)/B
  }
  res_varying_cont_reg[,h_y] <- CD_accuracy
}

res_varying_cont_reg
saveRDS(object =  res_varying_cont_reg, file = "Simulations/Causal Discovery/res_varying_cont_reg.RDS")