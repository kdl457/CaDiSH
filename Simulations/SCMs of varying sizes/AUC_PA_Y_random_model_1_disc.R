#### Using AUC but only using models of size 1.
#### Script for simulating graph and data with Y and 6 covariates. Varying number of parents of Y
#### DISCRETE ENVIRONMENTS AFFECTING ALL BUT Y.

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

lambda <- 0.1 # Found by LOOCV on one dataset
omega_H <- seq(from = 0, to = 40, length.out = 11)
B <- 50 # Number of simulations for each value of omega
n_CI <- 50 # Number of repetitions used for confidence intervals

AUC_res <- c() # Is the correct parent recovered?
AUC_average <- c() # Probability of recovering correct parent
AUC_res_ell_1 <- c()
AUC_average_ell_1 <- c()
AUC_res_reg <- c()
AUC_average_reg <- c()
res_AUC_ran_parents_1_disc <- matrix(NA, nrow = n_CI, ncol = length(omega_H))
res_AUC_ran_parents_1_ell_1_disc <- matrix(NA, nrow = n_CI, ncol = length(omega_H))
res_AUC_ran_parents_1_reg_disc <- matrix(NA, nrow = n_CI, ncol = length(omega_H))

for (om_H in 1:length(omega_H)) {
  Hidden_var <- sin(omega_H[om_H] * t)
  for (h in 1:n_CI) {
    cat(paste("Omega_H", omega_H[om_H], "Number", h))
    for (b in 1:B) {
      Y_rank <- as.numeric(sample(c(4:7), size = 1)) # Causal ranking of Y. 7 variables in total. Y must have 3 parents
      CH_h <- as.numeric(sample(c(1:7)[c(1:7) != Y_rank], size = sample(1:6, size = 1))) # Variable are affected by H.
      CH_E <- as.numeric(c(1:7)[c(1:7) != Y_rank]) # All are affected by E.
      
      parent_list <- list()
      for (j in 2:7) { # Create parent sets
        parent_list[[j]] <- as.numeric(sample(c(1: (j - 1)), sample(c(0: (j - 1)), size = 1)))
      }
      
      if (length(parent_list[[Y_rank]]) == 0) { # Making sure that Y has at least 1 parent
        parent_list[[Y_rank]] <- as.numeric(sample(c(1: (Y_rank - 1)), sample(c(1: (Y_rank - 1)), size = 1)))
      }
      
      # Creating list of variables
      X_list <- list("X_1" = c(rep(NA, n_obs)), "X_2" = c(rep(NA, n_obs)), "X_3" = c(rep(NA, n_obs)), 
                     "X_4" = c(rep(NA, n_obs)), "X_5" = c(rep(NA, n_obs)), "X_6" = c(rep(NA, n_obs)), 
                     "X_7" = c(rep(NA, n_obs)))
      
      # Simulating without effects of parents
      
      for (j in 1:7) {
        if (j %in% c(1:7)[c(1:7) != Y_rank]) {
        c_1 <- runif(3, -5, 5)
        c_N <- runif(3, 0.2, 5)
        c_H <- runif(3, -5, 5)
        noise <- rnorm(n = n_obs)
        X_list[[j]][which(t<0.25 | t >= 0.75)] <- (j %in% CH_h) * c_H[1] * Hidden_var[which(t<0.25 | t >= 0.75)] + 
          c_1[1] + c_N[1] * noise[which(t<0.25 | t >= 0.75)]
        X_list[[j]][which(t>=0.25 | t < 0.5)] <- (j %in% CH_h) * c_H[1] * Hidden_var[which(t>=0.25 | t < 0.5)] + 
          c_1[1] + c_N[1] * noise[which(t>=0.25 | t < 0.5)]
        X_list[[j]][which(t>=0.5 | t < 0.75)] <- (j %in% CH_h) * c_H[1] * Hidden_var[which(t>=0.5 | t < 0.75)] + 
          c_1[1] + c_N[1] * noise[which(t>=0.5 | t < 0.75)]
        } else {
          c_1 <- runif(1, -5, 5)
          c_N <- runif(1, 0.2, 5)
          noise <- rnorm(n = n_obs)
          X_list[[j]] <- c_1 + c_N * noise   
        }
      }
      
      # Adding parent effects for all. Only Y is not intervened on
      
      for (j in c(1:7)) {
        if (j %in% c(1:7)[c(1:7) != Y_rank]) {
          k <- length(parent_list[[j]])
          if (k > 0) {
            parent_coeff <- runif(3*k, 0, - 0.5 + 10 + 10 - 0.5)
            for (i in 1:(3*k)) {
              if( parent_coeff[i] < (9.5) ){
                parent_coeff[i] <- -10 + parent_coeff[i]
              }else{
                parent_coeff[i] <- 0.5 + parent_coeff[i] - 9.5
              } 
            }
            for (i in 1:k) {
              X_list[[j]][which(t<0.25 | t >= 0.75)] <- X_list[[j]][which(t<0.25 | t >= 0.75)] + 
                parent_coeff[3*i - 2] * X_list[[parent_list[[j]][i]]][which(t<0.25 | t >= 0.75)]
              X_list[[j]][which(t>=0.25 | t < 0.5)] <- X_list[[j]][which(t>=0.25 | t < 0.5)] + 
                parent_coeff[3*i - 1] * X_list[[parent_list[[j]][i]]][which(t>=0.25 | t < 0.5)]
              X_list[[j]][which(t>=0.5 | t < 0.75)] <- X_list[[j]][which(t>=0.5 | t < 0.75)] + 
                parent_coeff[3*i] * X_list[[parent_list[[j]][i]]][which(t>=0.5 | t < 0.75)]
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
                             data = data.frame("Y" = X_list[[Y_rank]], "X" = X_list[[models_1[i]]]), 
                             pen_fun = beta_pen_2_sq, lambda = lambda, gr = grad_pen_2_sq, lower = -10, upper = 10, 
                             method = "L-BFGS-B", control = list(maxit = 500))$value
      }
      
      ranking_scores <- order(scores_1) # Ranking of values of objective function
      
      ranking_models <- list() # Rank models according to scores
      for (i in 1:6) {
        ranking_models[[i]] <- models_1[,ranking_scores[i]]
      }
      
      # Calculate AUC scores
      
      AUC <- 0 
      F_count <- 0
      for (i in 1:6) {
        AUC <- (ranking_models[[i]] %in% parent_list[[Y_rank]]) * (6 - length(parent_list[[Y_rank]]) - F_count) + AUC
        F_count <- F_count + (1 - (ranking_models[[i]] %in% parent_list[[Y_rank]]))
      }
      
      if (length(parent_list[[Y_rank]]) != 6) {
        AUC_res[b] <- AUC / ((6 * length(parent_list[[Y_rank]])) - (length(parent_list[[Y_rank]])^2))
      } else {
        AUC_res[b] <- AUC
      }
      
      ######
      ###### Using ell_1
      ######
      
      models_1 <- combn(c(1:7)[c(1:7) != Y_rank], m = 1) # All possible models of size 1
      scores_1 <- c()  
      
      for (i in 1:ncol(models_1)) {
        scores_1[i] <- optim(par = rep(0, 2*n_obs), fn = obj_fun_1,
                             data = data.frame("Y" = X_list[[Y_rank]], "X" = X_list[[models_1[i]]]), 
                             pen_fun = beta_pen_sq, lambda = lambda, gr = grad_pen_sq, lower = -10, upper = 10, 
                             method = "L-BFGS-B", control = list(maxit = 500))$value
      }
      
      ranking_scores <- order(scores_1) # Ranking of values of objective function
      
      ranking_models <- list() # Rank models according to scores
      for (i in 1:6) {
        ranking_models[[i]] <- models_1[,ranking_scores[i]]
      }
      
      # Calculate AUC scores
      
      AUC <- 0 
      F_count <- 0
      for (i in 1:6) {
        AUC <- (ranking_models[[i]] %in% parent_list[[Y_rank]]) * (6 - length(parent_list[[Y_rank]]) - F_count) + AUC
        F_count <- F_count + (1 - (ranking_models[[i]] %in% parent_list[[Y_rank]]))
      }
      
      if (length(parent_list[[Y_rank]]) != 6) {
        AUC_res_ell_1[b] <- AUC / ((6 * length(parent_list[[Y_rank]])) - (length(parent_list[[Y_rank]])^2))
      } else {
        AUC_res_ell_1[b] <- AUC
      }
      
      ######
      ###### Using linear regression
      ######
      X_list_pred <- X_list
      X_list_pred[Y_rank] <- NULL
      
      models_1 <- combn(c(1:7)[c(1:7) != Y_rank], m = 1) # All possible models of size 1
      scores_1 <- c()  
      
      for (i in 1:ncol(models_1)) {
        scores_1[i] <- summary(lm(X_list[[Y_rank]] ~ X_list[[models_1[,i][1]]]))$fstatistic[1]
      }
      
      ranking_scores <- order(scores_1, decreasing = TRUE) # Ranking of values of f_stat
      
      ranking_models <- list() # Creating list of 15 best models
      for (i in 1:6) {
        ranking_models[[i]] <- models_1[,ranking_scores[i]]
      }
      
      # Calculate AUC scores
      
      AUC <- 0 
      F_count <- 0
      for (i in 1:6) {
        AUC <- (ranking_models[[i]] %in% parent_list[[Y_rank]]) * (6 - length(parent_list[[Y_rank]]) - F_count) + AUC
        F_count <- F_count + (1 - (ranking_models[[i]] %in% parent_list[[Y_rank]]))
      }
      
      if (length(parent_list[[Y_rank]]) != 6) {
        AUC_res_reg[b] <- AUC / ((6 * length(parent_list[[Y_rank]])) - (length(parent_list[[Y_rank]])^2))
      } else {
        AUC_res_reg[b] <- AUC
      }      
      
    }
    AUC_average[h] <- sum(AUC_res) / B
    AUC_average_ell_1[h] <- sum(AUC_res_ell_1) / B
    AUC_average_reg[h] <- sum(AUC_res_reg) / B
  }
  res_AUC_ran_parents_1_disc[,om_H] <- AUC_average
  res_AUC_ran_parents_1_ell_1_disc[,om_H] <- AUC_average_ell_1
  res_AUC_ran_parents_1_reg_disc[,om_H] <- AUC_average_reg
}
colSums(res_AUC_ran_parents_1_disc)/50
colSums(res_AUC_ran_parents_1_ell_1_disc)/50
colSums(res_AUC_ran_parents_1_reg_disc)/50

saveRDS(object = res_AUC_ran_parents_1_disc, 
        file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_disc.RDS")
saveRDS(object = res_AUC_ran_parents_1_ell_1_disc, 
        file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_ell_1_disc.RDS")
saveRDS(object = res_AUC_ran_parents_1_reg_disc, 
        file = "Simulations/Causal Discovery/Randomized Graph/res_AUC_ran_parents_1_reg_disc.RDS")
