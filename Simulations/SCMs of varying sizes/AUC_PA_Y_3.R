#### Using AUC
#### Script for simulating graph and data with Y and 6 covariates. Fixing 3 parents of Y

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
B <- 10 # Number of simulations for each value of omega
n_CI <- 10 # Number of repetitions used for confidence intervals

AUC_res <- c() # Is the correct parent recovered?
AUC_average <- c() # Probability of recovering correct parent
AUC_res_ell_1 <- c()
AUC_average_ell_1 <- c()
AUC_res_reg <- c()
AUC_average_reg <- c()
res_AUC_3_parents <- matrix(NA, nrow = n_CI, ncol = length(omega_H))
res_AUC_3_parents_ell_1 <- matrix(NA, nrow = n_CI, ncol = length(omega_H))
res_AUC_3_parents_reg <- matrix(NA, nrow = n_CI, ncol = length(omega_H))

for (om_H in 1:length(omega_H)) {
  Hidden_var <- sin(omega_H[om_H] * t)
  for (h in 1:n_CI) {
    cat(paste("Omega_H", omega_H[om_H], "Number", h))
    for (b in 1:B) {
      Y_rank <- as.numeric(sample(c(4:7), size = 1)) # Causal ranking of Y. 7 variables in total. Y must have 3 parents
      CH_h <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Which variable are affected by H.
      CH_E <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Which variables are affected by E.
      
      parent_list <- list()
      for (j in 2:7) { # Create parent sets
        parent_list[[j]] <- as.numeric(sample(c(1: (j - 1)), sample(c(0: (j - 1)), size = 1)))
      }
      
      parent_list[[Y_rank]] <- as.numeric(sample(c(1: (Y_rank - 1)), 3)) # Fixing PA_Y to size 3
      
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
      
      
      ######
      ###### Using ell_2
      ######
      
      models_1 <- combn(c(1:7)[c(1:7) != Y_rank], m = 1) # All possible models of size 1
      scores_1 <- c()  

      for (i in 1:ncol(models_1)) {
        scores_1[i] <- optim(par = rep(0, 2*n_obs), fn = obj_fun_1,
                        data = data.frame("Y" = X_list[[Y_rank]], "X" = X_list[[models_1[i]]]), pen_fun = beta_pen_2_sq, 
                        lambda = lambda, gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                        control = list(maxit = 500))$value
      }
      
      models_2 <- combn(c(1:7)[c(1:7) != Y_rank], m = 2) # All possible models of size 2
      scores_2 <- c()
      
      for (i in 1:ncol(models_2)) {
        scores_2[i] <- optim(par = rep(0, 3*n_obs), fn = obj_fun_2,
                             data = data.frame("Y" = X_list[[Y_rank]], 
                                               "X_1" = X_list[[models_2[,i][1]]],
                                               "X_2" = X_list[[models_2[,i][1]]],
                                               "X_3" = X_list[[models_2[,i][1]]]), 
                             pen_fun = beta_pen_2_sq, lambda = lambda, gr = grad_pen_2_sq_2, lower = -10, 
                             upper = 10, method = "L-BFGS-B", control = list(maxit = 500))$value
      }
      
      models_3 <- combn(c(1:7)[c(1:7) != Y_rank], m = 3) # All possible models of size 3
      scores_3 <- c()
      
      for (i in 1:ncol(models_3)) {
        scores_3[i] <- optim(par = rep(0, 4*n_obs), fn = obj_fun_3,
                               data = data.frame("Y" = X_list[[Y_rank]], 
                                                 "X_1" = X_list[[models_3[,i][1]]],
                                                 "X_2" = X_list[[models_3[,i][1]]],
                                                 "X_3" = X_list[[models_3[,i][1]]]), 
                               pen_fun = beta_pen_2_sq, lambda = lambda, gr = grad_pen_2_sq_3, lower = -10, 
                               upper = 10, method = "L-BFGS-B", control = list(maxit = 500))$value
      }
      
      ranking_scores <- order(c(scores_1, scores_2, scores_3)) # Ranking of values of objective function

      ranking_models <- list() # Creating list of 15 best models
      for (i in 1:15) {
        loc <- ranking_scores[i]
        if (loc < 7) {
          ranking_models[[i]] <- models_1[,loc]
        } else {
          if (loc > 21) {
            ranking_models[[i]] <- models_3[,loc - 21]
          } else {
            ranking_models[[i]] <- models_2[,loc - 6]
          }
        }
      }
      
      prop_vars <- c() # How often does each predictor occur in the top 15 models
      for (i in c(1:7)[c(1:7) != Y_rank]) {
        prop <- 0
        for (j in 1:15) {
          prop <- prop + (i %in% ranking_models[[j]])
        }
        prop_vars[i] <- prop / 15
      }
      
      ranking_vars <- order(prop_vars, decreasing = TRUE) # Rank predictors according to prop_vars. Y is last
      
      # Calculate AUC scores
    
      AUC <- 0 
      F_count <- 0
      for (i in 1:6) {
        AUC <- (ranking_vars[i] %in% parent_list[[Y_rank]]) * (3 - F_count) + AUC
        F_count <- F_count + (1 - (ranking_vars[i] %in% parent_list[[Y_rank]]))
      }
      
      AUC_res[b] <- AUC / (3^2)
      
      ######
      ###### Using ell_1
      ######
      
      models_1 <- combn(c(1:7)[c(1:7) != Y_rank], m = 1) # All possible models of size 1
      scores_1 <- c()  
      
      for (i in 1:ncol(models_1)) {
        scores_1[i] <- optim(par = rep(0, 2*n_obs), fn = obj_fun_1,
                             data = data.frame("Y" = X_list[[Y_rank]], "X" = X_list[[models_1[i]]]), pen_fun = beta_pen_sq, 
                             lambda = lambda, gr = grad_pen_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                             control = list(maxit = 500))$value
      }
      
      models_2 <- combn(c(1:7)[c(1:7) != Y_rank], m = 2) # All possible models of size 2
      scores_2 <- c()
      
      for (i in 1:ncol(models_2)) {
        scores_2[i] <- optim(par = rep(0, 3*n_obs), fn = obj_fun_2,
                             data = data.frame("Y" = X_list[[Y_rank]], 
                                               "X_1" = X_list[[models_2[,i][1]]],
                                               "X_2" = X_list[[models_2[,i][1]]],
                                               "X_3" = X_list[[models_2[,i][1]]]), 
                             pen_fun = beta_pen_sq, lambda = lambda, gr = grad_pen_sq_2, lower = -10, 
                             upper = 10, method = "L-BFGS-B", control = list(maxit = 500))$value
      }
      
      models_3 <- combn(c(1:7)[c(1:7) != Y_rank], m = 3) # All possible models of size 3
      scores_3 <- c()
      
      for (i in 1:ncol(models_3)) {
        scores_3[i] <- optim(par = rep(0, 4*n_obs), fn = obj_fun_3,
                             data = data.frame("Y" = X_list[[Y_rank]], 
                                               "X_1" = X_list[[models_3[,i][1]]],
                                               "X_2" = X_list[[models_3[,i][1]]],
                                               "X_3" = X_list[[models_3[,i][1]]]), 
                             pen_fun = beta_pen_sq, lambda = lambda, gr = grad_pen_sq_3, lower = -10, 
                             upper = 10, method = "L-BFGS-B", control = list(maxit = 500))$value
      }
      
      ranking_scores <- order(c(scores_1, scores_2, scores_3)) # Ranking of values of objective function
      
      ranking_models <- list() # Creating list of 15 best models
      for (i in 1:15) {
        loc <- ranking_scores[i]
        if (loc < 7) {
          ranking_models[[i]] <- models_1[,loc]
        } else {
          if (loc > 21) {
            ranking_models[[i]] <- models_3[,loc - 21]
          } else {
            ranking_models[[i]] <- models_2[,loc - 6]
          }
        }
      }
      
      prop_vars <- c() # How often does each predictor occur in the top 15 models
      for (i in c(1:7)[c(1:7) != Y_rank]) {
        prop <- 0
        for (j in 1:15) {
          prop <- prop + (i %in% ranking_models[[j]])
        }
        prop_vars[i] <- prop / 15
      }
      
      ranking_vars <- order(prop_vars, decreasing = TRUE) # Rank predictors according to prop_vars. Y is last
      
      # Calculate AUC scores
      
      AUC <- 0 
      F_count <- 0
      for (i in 1:6) {
        AUC <- (ranking_vars[i] %in% parent_list[[Y_rank]]) * (3 - F_count) + AUC
        F_count <- F_count + (1 - (ranking_vars[i] %in% parent_list[[Y_rank]]))
      }
      
      AUC_res_ell_1[b] <- AUC / (3^2)
      
      ######
      ###### Using linear regression
      ######
      
      models_1 <- combn(c(1:7)[c(1:7) != Y_rank], m = 1) # All possible models of size 1
      scores_1 <- c()  
      
      for (i in 1:ncol(models_1)) {
        scores_1[i] <- summary(lm(X_list[[Y_rank]] ~ X_list[[models_1[,i][1]]]))$fstatistic[1]
      }
      
      models_2 <- combn(c(1:7)[c(1:7) != Y_rank], m = 2) # All possible models of size 2
      scores_2 <- c()
      
      for (i in 1:ncol(models_2)) {
        scores_2[i] <- summary(lm(X_list[[Y_rank]] ~ X_list[[models_2[,i][1]]] + 
                                    X_list[[models_2[,i][2]]]))$fstatistic[1]
      }
      
      models_3 <- combn(c(1:7)[c(1:7) != Y_rank], m = 3) # All possible models of size 3
      scores_3 <- c()
      
      for (i in 1:ncol(models_3)) {
        scores_3[i] <- summary(lm(X_list[[Y_rank]] ~ X_list[[models_3[,i][1]]] + X_list[[models_3[,i][2]]] + 
                                    X_list[[models_3[,i][3]]]))$fstatistic[1]
      }
      
      ranking_scores <- order(c(scores_1, scores_2, scores_3), decreasing = TRUE) # Ranking of values of f_stat
      
      ranking_models <- list() # Creating list of 15 best models
      for (i in 1:15) {
        loc <- ranking_scores[i]
        if (loc < 7) {
          ranking_models[[i]] <- models_1[,loc]
        } else {
          if (loc > 21) {
            ranking_models[[i]] <- models_3[,loc - 21]
          } else {
            ranking_models[[i]] <- models_2[,loc - 6]
          }
        }
      }
      
      prop_vars <- c() # How often does each predictor occur in the top 15 models
      for (i in c(1:7)[c(1:7) != Y_rank]) {
        prop <- 0
        for (j in 1:15) {
          prop <- prop + (i %in% ranking_models[[j]])
        }
        prop_vars[i] <- prop / 15
      }
      
      ranking_vars <- order(prop_vars, decreasing = TRUE) # Rank predictors according to prop_vars. Y is last
      
      # Calculate AUC scores
      
      AUC <- 0 
      F_count <- 0
      for (i in 1:6) {
        AUC <- (ranking_vars[i] %in% parent_list[[Y_rank]]) * (3 - F_count) + AUC
        F_count <- F_count + (1 - (ranking_vars[i] %in% parent_list[[Y_rank]]))
      }
      
      AUC_res_reg[b] <- AUC / (3^2)
      
  
      
    }
    AUC_average[h] <- sum(AUC_res) / B
    AUC_average_ell_1[h] <- sum(AUC_res_ell_1) / B
    AUC_average_reg[h] <- sum(AUC_res_reg) / B
  }
  res_AUC_3_parents[,om_H] <- AUC_average
  res_AUC_3_parents_ell_1[,om_H] <- AUC_average_ell_1
  res_AUC_3_parents_reg[,om_H] <- AUC_average_reg
}
res_AUC_3_parents
res_AUC_3_parents_ell_1
res_AUC_3_parents_reg

saveRDS(object = res_AUC_3_parents, file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_3_parents.RDS")
saveRDS(object = res_AUC_3_parents_ell_1, file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_3_parents_ell_1.RDS")
saveRDS(object = res_AUC_3_parents_reg, file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_3_parents_reg.RDS")
