#### Script for simulating graph and data with Y and 6 covariates. Fixing one parent of Y
setwd("~/Dropbox/KU/Thesis/R")
source(file = "Functions/source_file.R")

# Number of parents of Y = 1
# Varying number of children of H
# Varying number of children of E
set.seed(1)
# Idea: Fix frequency omega_E = 8 pi, and vary omega_H from 0 to 30
omega_E <- 8 * pi
n_obs <- 50
t <- seq(from = 0, to = 1, length.out = n_obs)

lambda <- 5.5 # Found by LOOCV on one dataset
omega_H <- seq(from = 0, to = 40, length.out = 11)
B <- 50 # Number of simulations for each value of omega
n_CI <- 50 # Number of repetitions used for confidence intervals

parent_res <- c() # Is the correct parent recovered?
CD_accuracy <- c() # Probability of recovering correct parent
parent_res_ell_1 <- c()
CD_accuracy_ell_1 <- c()
res_matrix_1_parent <- matrix(NA, nrow = n_CI, ncol = length(omega_H))
res_matrix_1_parent_ell_1 <- matrix(NA, nrow = n_CI, ncol = length(omega_H))

for (om_H in 1:length(omega_H)) {
  Hidden_var <- sin(omega_H[om_H] * t)
  for (h in 1:n_CI) {
    cat(paste("Omega_H", omega_H[om_H], "Number", h))
    for (b in 1:B) {
Y_rank <- as.numeric(sample(c(2:7), size = 1)) # Causal ranking of Y. 7 variables in total. Y must have one parent
CH_h <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Which variable are affected by H.
CH_E <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Which variables are affected by E.

parent_list <- list()
for (j in 2:7) { # Create parent sets
  parent_list[[j]] <- as.numeric(sample(c(1: (j - 1)), sample(c(0: (j - 1)), size = 1)))
}

parent_list[[Y_rank]] <- as.numeric(sample(c(1: (Y_rank - 1)), 1)) # Fixing PA_Y to size 1

# Creating list of variables
X_list <- list("X_1" = c(), "X_2" = c(), "X_3" = c(), "X_4" = c(), "X_5" = c(), "X_6" = c(), "X_7" = c())

# Simulating without effects of parents
for (j in 1:7) {
  phase <- runif(1, -1, 1)
  c_1 <- runif(1, -5, 5)
  c_N <- runif(1, 0.2, 5)
  c_H <- runif(1, -5, 5)
  c_E <- runif(1, -5, 5)
  noise <- rnorm(n = n_obs)
  X_list[[j]] <- (j %in% CH_h) * c_H * Hidden_var + (j %in% CH_E)  *  c_E * sin(omega_E * t + phase) + 
    c_1 + c_N * noise   
}

# Adding parent effects for all

for (j in c(1:7)) {
  if (j %in% c(1:7)[c(1:7) != Y_rank]) {
    k <- length(parent_list[[j]])
    if (k > 0) {
      parent_coeff <- runif(k, 0, - 0.5 + 10 + 10 - 0.5)
      c_E <- runif(k, -5, 5)
      for (i in 1:k) {
        if( parent_coeff[i] < (9.5) ){
          parent_coeff[i] <- -10 + parent_coeff[i]
        }else{
          parent_coeff[i] <- 0.5 + parent_coeff[i] - 9.5
        } 
      }
      for (i in 1:k) {
        X_list[[j]] <- X_list[[j]] + parent_coeff[i] * X_list[[parent_list[[j]][i]]] + 
          (j %in% CH_E)  *  c_E[k] * sin(omega_E * t + phase) * X_list[[parent_list[[j]][i]]]
      }  
    }
  } else {
    PA_Y <- length(parent_list[[Y_rank]])
    parent_coeff <- runif(2 * PA_Y, 0, - 0.5 + 10 + 10 - 0.5)
    for (i in 1:(2*PA_Y)) {
      if( parent_coeff[i] < (9.5) ){
        parent_coeff[i] <- -10 + parent_coeff[i]
      }else{
        parent_coeff[i] <- 0.5 + parent_coeff[i] - 9.5
      } 
    }
    for (i in 1:PA_Y) {
      X_list[[Y_rank]] <- X_list[[Y_rank]] + parent_coeff[2*i - 1] * Hidden_var * 
        X_list[[parent_list[[Y_rank]][i]]] + 
        parent_coeff[2 * i] * X_list[[parent_list[[Y_rank]][i]]]
    }
  }
}

# Normalize data

for (i in 1:7) {
  X_list[[i]] <- (X_list[[i]] - mean(X_list[[i]])) / sd(X_list[[i]])
}

# Create results using CaDiSH

result_vec <- c()

for (i in c(1:7)[c(1:7) != Y_rank]) {
  result <- optim(par = rep(0, 2*n_obs), fn = obj_fun_1,
                  data = data.frame("Y" = X_list[[Y_rank]], "X" = X_list[[i]]), pen_fun = beta_pen_2_sq, 
                  lambda = lambda, gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                  control = list(maxit = 500))$value
  result_vec[i] <- result
}

result_vec_ell_1 <- c()

for (i in c(1:7)[c(1:7) != Y_rank]) {
  result <- optim(par = rep(0, 2*n_obs), fn = obj_fun_1,
                  data = data.frame("Y" = X_list[[Y_rank]], "X" = X_list[[i]]), pen_fun = beta_pen_sq, 
                  lambda = lambda, gr = grad_pen_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                  control = list(maxit = 500))$value
  result_vec_ell_1[i] <- result
}

ranking <- order(result_vec) # NA is for Y. Ranks last
correct <- ranking %in% parent_list[[Y_rank]]

ranking_ell_1 <- order(result_vec_ell_1) # NA is for Y. Ranks last
correct_ell_1 <- ranking_ell_1 %in% parent_list[[Y_rank]]

parent_res[b] <- correct[1]
parent_res_ell_1[b] <- correct_ell_1[1]
    }
    CD_accuracy[h] <- sum(parent_res)/B
    CD_accuracy_ell_1[h] <- sum(parent_res_ell_1)/B
  }
  res_matrix_1_parent[,om_H] <- CD_accuracy
  res_matrix_1_parent_ell_1[,om_H] <- CD_accuracy_ell_1
}
res_matrix_1_parent
res_matrix_1_parent_ell_1
saveRDS(object = res_matrix_1_parent, file = "Simulations/Causal Discovery/Randomized Graph/res_matrix_1_parent.RDS")
saveRDS(object = res_matrix_1_parent_ell_1, file = "Simulations/Causal Discovery/Randomized Graph/res_matrix_1_parent_ell_1.RDS")


####### Recovering parent using linear regression

set.seed(1)
# Idea: Fix frequency omega_E = 8 pi, and vary omega_H from 0 to 30
omega_E <- 8 * pi
n_obs <- 50
t <- seq(from = 0, to = 1, length.out = n_obs)

lambda <- 5.5 # Found by LOOCV on one dataset
omega_H <- seq(from = 0, to = 40, length.out = 11)
B <- 50 # Number of simulations for each value of omega
n_CI <- 50 # Number of repetitions used for confidence intervals

parent_res <- c() # Is the correct parent recovered?
CD_accuracy <- c() # Probability of recovering correct parent
parent_res_ell_1 <- c()
CD_accuracy_ell_1 <- c()
res_matrix_1_reg <- matrix(NA, nrow = n_CI, ncol = length(omega_H))

for (om_H in 1:length(omega_H)) {
  Hidden_var <- sin(omega_H[om_H] * t)
  for (h in 1:n_CI) {
    cat(paste("Omega_H", omega_H[om_H], "Number", h))
    for (b in 1:B) {
      Y_rank <- as.numeric(sample(c(2:7), size = 1)) # Causal ranking of Y. 7 variables in total. Y must have one parent
      CH_h <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Which variable are affected by H.
      CH_E <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Which variables are affected by E.
      
      parent_list <- list()
      for (j in 2:7) { # Create parent sets
        parent_list[[j]] <- as.numeric(sample(c(1: (j - 1)), sample(c(0: (j - 1)), size = 1)))
      }
      
      parent_list[[Y_rank]] <- as.numeric(sample(c(1: (Y_rank - 1)), 1)) # Fixing PA_Y to size 1
      
      # Creating list of variables
      X_list <- list("X_1" = c(), "X_2" = c(), "X_3" = c(), "X_4" = c(), "X_5" = c(), "X_6" = c(), "X_7" = c())
      
      # Simulating without effects of parents
      for (j in 1:7) {
        phase <- runif(1, -1, 1)
        c_1 <- runif(1, -5, 5)
        c_N <- runif(1, 0.2, 5)
        c_H <- runif(1, -5, 5)
        c_E <- runif(1, -5, 5)
        noise <- rnorm(n = n_obs)
        X_list[[j]] <- (j %in% CH_h) * c_H * Hidden_var + (j %in% CH_E)  *  c_E * sin(omega_E * t + phase) + 
          c_1 + c_N * noise   
      }
      
      # Adding parent effects for all but Y
      
      for (j in c(1:7)[c(1:7) != Y_rank]) {
        k <- length(parent_list[[j]])
        if (k > 0) {
          parent_coeff <- runif(k, 0, - 0.5 + 10 + 10 - 0.5)
          c_E <- runif(1, -5, 5)
          for (i in 1:k) {
            if( parent_coeff[i] < (4.5) ){
              parent_coeff[i] <- -10 + parent_coeff[i]
            }else{
              parent_coeff[i] <- 0.5 + parent_coeff[i] - 9.5
            } 
          }
          for (i in 1:k) {
            X_list[[j]] <- X_list[[j]] + parent_coeff[i] * X_list[[parent_list[[j]][i]]] + 
              (j %in% CH_E)  *  c_E * sin(omega_E * t + phase) * X_list[[parent_list[[j]][i]]]
          }  
        }
      }
      
      # Adding parent effects for Y - depends on H
      
      PA_Y <- length(parent_list[[Y_rank]])
      parent_coeff <- runif(2 * PA_Y, 0, - 0.5 + 10 + 10 - 0.5)
      for (i in 1:PA_Y) {
        if( parent_coeff[i] < (4.5) ){
          parent_coeff[i] <- -10 + parent_coeff[i]
        }else{
          parent_coeff[i] <- 0.5 + parent_coeff[i] - 9.5
        } 
      }
      for (i in 1:PA_Y) {
        X_list[[Y_rank]] <- X_list[[Y_rank]] + parent_coeff[2*i - 1] * Hidden_var * X_list[[parent_list[[Y_rank]][i]]] + 
          parent_coeff[2 * i] * X_list[[parent_list[[Y_rank]][i]]]
      }
      
      # Normalize data
      
      for (i in 1:7) {
        X_list[[i]] <- (X_list[[i]] - mean(X_list[[i]])) / sd(X_list[[i]])
      }
      
      # Create results using linear regression
      
      
      result_vec <- c()
      for (i in c(1:7)[c(1:7) != Y_rank]) {
        result <- summary(lm(X_list[[Y_rank]] ~ X_list[[i]]))$coefficients[,4][2]
        result_vec[i] <- result
      }
      
      ranking <- order(result_vec) # NA is for Y. Ranks last
      correct <- ranking %in% parent_list[[Y_rank]]
      
      parent_res[b] <- correct[1]
    }
    CD_accuracy[h] <- sum(parent_res)/B
  }
  res_matrix_1_reg[,om_H] <- CD_accuracy
}
res_matrix_1_reg
saveRDS(object = res_matrix_1_reg, file = "Simulations/Causal Discovery/Randomized Graph/res_matrix_1_reg.RDS")

