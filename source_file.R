####### Beta penalty functions

library(stocks)

beta_pen <- function(x) {
  return(abs(diffs(x, 1)))
}

beta_pen_sq <- function(x) {
  return(sum(diffs(x, 1)^2))
}

beta_pen_2 <- function(x) {
  return(sum(abs(diffs(x, 1)[-1] - diffs(x, 1)[-(length(x) - 1)])))
}

beta_pen_2_sq <- function(x) {
  return(sum((diffs(x, 1)[-1] - diffs(x, 1)[-(length(x) - 1)])^2))
}

beta_pen_harm <- function(x) {
  return(sum(((diffs(x, 1)[-c(1,2)] - diffs(x, 1)[-c(1,(length(x) - 1))]) - 
                (diffs(x, 1)[-c(1,(length(x) - 1))] - diffs(x, 1)[-c((length(x) - 2),(length(x) - 1))]))^2))
}

##### Objective function and gradients for two variables Y and X

obj_fun_1 <- function(beta, data, pen_fun, lambda, ...) {
  n_sim <- length(data[,1])
  alpha_pen <- lambda * pen_fun(beta[1:n_sim], ...)
  beta_pen <- lambda * pen_fun(beta[(n_sim + 1): (2*n_sim)], ...)
  SSR <- sum((data[,1] - data[,2] * beta[(n_sim + 1): (2*n_sim)] - beta[1:n_sim])^2)
  return((SSR + beta_pen + alpha_pen) / nrow(data))
}

grad_pen_sq <- function(beta, data, pen_fun, lambda){
  Y <- data[,1]
  X <- data[,2]
  n <- length(beta) / 2
  X <- matrix(c(rep(1, length(X)), X), nrow = n)
  beta_mat <- matrix(beta, nrow = n)
  grad_mat <- matrix(rep(NA, length(beta)), nrow = n)
  for (i in 1:nrow(grad_mat)) {
    grad_mat[i,1] <- - 2 * Y[i] + 2 * beta_mat[i,1] + 2 * beta_mat[i,2] * X[i,2]
    grad_mat[i,2] <- - 2 * Y[i] * X[i,2] + 2 * beta_mat[i,2] * (X[i,2]^2) + 
      2 * X[i,2] * beta_mat[i,1] * X[i,1]
  }
  grad_mat[1,1] <- grad_mat[1,1] + 2 * lambda * (beta_mat[1,1] - beta_mat[2,1])
  grad_mat[1,2] <- grad_mat[1,2] + 2 * lambda * (beta_mat[1,2] - beta_mat[2,2])
  grad_mat[n,1] <- grad_mat[n,1] + 2 * lambda * (beta_mat[n,1] - beta_mat[n-1,1])
  grad_mat[n,2] <- grad_mat[n,2] + 2 * lambda * (beta_mat[n,2] - beta_mat[n-1,2])
  for (i in 2:(n-1)) {
    grad_mat[i,1] <- grad_mat[i,1] + 2 * lambda * (2 * beta_mat[i,1] - beta_mat[i-1,1] - beta_mat[i+1,1])
    grad_mat[i,2] <- grad_mat[i,2] + 2 * lambda * (2 * beta_mat[i,2] - beta_mat[i-1,2] - beta_mat[i+1,2])
  }
  c(grad_mat[,1], grad_mat[,2])
}

grad_pen_2_sq <- function(beta, data, pen_fun, lambda){
  Y <- data[,1]
  X <- data[,2]
  n <- length(beta) / 2
  X <- matrix(c(rep(1, length(X)), X), nrow = n)
  beta_mat <- matrix(beta, nrow = n)
  grad_mat <- matrix(rep(NA, length(beta)), nrow = n)
  for (i in 1:nrow(grad_mat)) {
    grad_mat[i,1] <- - 2 * Y[i] + 2 * beta_mat[i,1] + 2 * beta_mat[i,2] * X[i,2]
    grad_mat[i,2] <- - 2 * Y[i] * X[i,2] + 2 * beta_mat[i,2] * (X[i,2]^2) + 
      2 * X[i,2] * beta_mat[i,1] * X[i,1]
  }
  grad_mat[1,1] <- grad_mat[1,1] + 2 * lambda * (beta_mat[1,1] - 2 * beta_mat[2,1] + beta_mat[3,1])
  grad_mat[1,2] <- grad_mat[1,2] + 2 * lambda * (beta_mat[1,2] - 2 * beta_mat[2,2] + beta_mat[3,2])
  grad_mat[2,1] <- grad_mat[2,1] + 
    2 * lambda * (- 2 * beta_mat[1,1] + 5 * beta_mat[2,1] - 4 * beta_mat[3,1] + beta_mat[4,1])
  grad_mat[2,2] <- grad_mat[2,2] + 
    2 * lambda * (- 2 * beta_mat[1,2] + 5 * beta_mat[2,2] - 4 * beta_mat[3,2] + beta_mat[4,2])
  grad_mat[n-1,1] <- grad_mat[n-1,1] + 
    2 * lambda * (beta_mat[n-3,1] - 4 * beta_mat[n-2,1] + 5 * beta_mat[n-1,1] - 2 * beta_mat[n,1])
  grad_mat[n-1,2] <- grad_mat[n-1,2] + 
    2 * lambda * (beta_mat[n-3,2] - 4 * beta_mat[n-2,2] + 5 * beta_mat[n-1,2] - 2 * beta_mat[n,2])
  grad_mat[n,1] <- grad_mat[n,1] + 2 * lambda * (beta_mat[n-2,1] - 2 * beta_mat[n-1,1] + beta_mat[n,1])
  grad_mat[n,2] <- grad_mat[n,2] + 2 * lambda * (beta_mat[n-2,2] - 2 * beta_mat[n-1,2] + beta_mat[n,2])
  for (i in 3:(n-2)) {
    grad_mat[i,1] <- grad_mat[i,1] + 
      2 * lambda * (beta_mat[i-2,1] - 4 * beta_mat[i-1,1] + 6 * beta_mat[i,1] - 4 * beta_mat[(i+1),1] + beta_mat[i+2,1])
    grad_mat[i,2] <- grad_mat[i,2] + 
      2 * lambda * (beta_mat[i-2,2] - 4 * beta_mat[i-1,2] + 6 * beta_mat[i,2] - 4 * beta_mat[(i+1),2] + beta_mat[i+2,2])
  }
  c(grad_mat[,1], grad_mat[,2])
}


########### Optimization function

grad_optim <- function(data, par = NA, pen_fun, grad_fun, lambda, lower = -10, upper = 10, maxit = 500){
  n_sim <- length(data[,1])
  cols <- names(data)
  cols <- cols[! cols %in% c("Y", "h", "t")]
  n <- length(cols)
  res_list <- list()
  for (i in 1:n) {
    if (is.na(par)[1] == TRUE) {
      model <- lm(as.formula(paste0("data$Y ~ data$X", i)))
      param <- c(rep(as.numeric(model$coefficients[1]), n_sim), rep(as.numeric(model$coefficients[2]), n_sim))
      tmp <- optim(par = param, fn = obj_fun_1, data = data[,which(names(data) %in% c("Y", paste0("X",i)))], 
                   pen_fun = pen_fun, lambda = lambda, gr = grad_fun, lower = min(param,lower), 
                   upper = max(param,upper), method = "L-BFGS-B", control = list(maxit = maxit))
    } else {
      tmp <- optim(par = c(par[1: n_sim], par[((i * n_sim) + 1): ((i + 1) * n_sim)]), fn = obj_fun_1, 
                   data = data[,which(names(data) %in% c("Y", paste0("X",i)))], pen_fun = pen_fun, lambda = lambda, 
                   gr = grad_fun, lower = min(param,lower), upper = max(param,upper), method = "L-BFGS-B", 
                   control = list(maxit = maxit))
    }
    tmp$SSR <- sum((data[,1] - data[,which(names(data) %in% c("Y", paste0("X",i)))][,2] * 
                      tmp$par[(n_sim + 1): (2*n_sim)] - tmp$par[1:n_sim])^2)
    res_list[[i]] <- tmp
  }
  return(res_list)
}


########### Simulating data

########### Hidden variable assignments

f_hidden <- function(t){
  2 * t - 1
}

f_hidden_10 <- function(t){
  10 * t - 1
}

f_hidden_poly <- function(t){
  (2/3) * t^2 + (4/3) * t - 1 + t^3
}

f_hidden_sin <- function(t){
  sin(20 * t)
}

f_hidden_sin_5 <- function(t){
  5 * sin(20 * t)
}

f_hidden_sin_2pi <- function(t){
  sin(2 * pi * t)
}

f_hidden_sin_10pi <- function(t){
  sin(10 * pi * t)
}

f_hidden_log <- function(t){
  log(t + 0.1)
}

############ Function for simulating data

# In this simulation, X2 is an invariant prediction set under different environments. Y has the
# same structural assignment throughout time which can be affected by H if the hidden_Y variable 
# is non-zero. The estimate of the regression coefficient of X2 on Y should be constant if not.

# sim_length = number of timepoints where the variables are simulated
# breaks = the upper and lower bounds of t along with 4 breaking points, where envir. changes
# h_func = the structural assignment of the hidden variable.
# hidden_y = the linear coefficient of the effect of the hidden variable on y.
# y_var an indicator if the coef. of x^2 in the structural assignment of y should vary with H.

sim_data <- function(sim_length, breaks, h_func, hidden_y, y_var = TRUE){
  t_points <- seq(from = breaks[1], to = breaks[length(breaks)], length.out = sim_length)
  hidden_t <- h_func(t_points)
  x3 <- 1 + hidden_t + 0.5 * rnorm(n = sim_length, mean = 0, sd = 1)
  x2 <- c(rep(NA, sim_length))
  for (i in 1:sim_length) {
    if (t_points[i] < breaks[2]) {
      x2[i] <- 2 + x3[i] + hidden_t[i] + rnorm(n = 1, mean = 0, sd = 1) 
    } else {
      if (t_points[i] < breaks[3]) {
        x2[i] <- 1 - x3[i] + 2 * hidden_t[i] + rnorm(n = 1, mean = 0, sd = 1)
      } else {
        if (t_points[i] < breaks[4]) {
          x2[i] <- 3 * x3[i] - 0.5 * hidden_t[i] + rnorm(n = 1, mean = 0, sd = 1)
        } else {
          x2[i] <- 2 + x3[i] + hidden_t[i] + rnorm(n = 1, mean = 0, sd = 1) 
        }
      }
    }
  }
  if (y_var == FALSE) {
    y <- 1 + x2 + hidden_y * hidden_t + 0.5 * rnorm(n = sim_length, mean = 0, sd = 1)  
  } else {
    y <- 1 + hidden_y * hidden_t * x2 + 0.5 * rnorm(n = sim_length, mean = 0, sd = 1)  
  }
  x1 <- c(rep(NA, sim_length))
  for (i in 1:sim_length) {
    if (t_points[i] < breaks[2]) {
      x1[i] <- 1 + x3[i] - y[i] + 0.7 * rnorm(n = 1, mean = 0, sd = 1) 
    } else {
      if (t_points[i] < breaks[3]) {
        x1[i] <- 1 + 2 * x3[i] + y[i] + 0.5 * rnorm(n = 1, mean = 0, sd = 1)
      } else {
        if (t_points[i] < breaks[4]) {
          x1[i] <- - 1 - x3[i] + y[i] + rnorm(n = 1, mean = 0, sd = 1)
        } else {
          x1[i] <- 1 + x3[i] - y[i] + 0.7 * rnorm(n = 1, mean = 0, sd = 1)
        }
      }
    }
  }
  out.list <- list(Y = y, X1 = x1, X2 = x2, X3 = x3, h = hidden_t, t = t_points)
  return(out.list)
}

# Simulating where each coefficient is constant under each environment

# numbers: alpha_3, h_3, n_3, alpha_y, hidden_y, n_y
# vectors of length three: alpha_2, c_23, h_2, n_2, alpha_1, c_13, c_1y, n_1

sim_data_fixed_env <- function(sim_length, breaks, h_func, 
                               alpha_3, h_3, n_3,
                               alpha_2, c_23, h_2, n_2,
                               alpha_y, hidden_y, n_y,
                               alpha_1, c_13, c_1y, n_1){
  t_points <- seq(from = breaks[1], to = breaks[length(breaks)], length.out = sim_length)
  hidden_t <- h_func(t_points)
  x3 <- alpha_3 + h_3 * hidden_t + n_3 * rnorm(n = sim_length, mean = 0, sd = 1)
  x2 <- c(rep(NA, sim_length))
  for (i in 1:sim_length) {
    if (t_points[i] < breaks[2]) {
      x2[i] <- alpha_2[1] + c_23[1] * x3[i] + h_2[1] * hidden_t[i] + n_2[1] * rnorm(n = 1, mean = 0, sd = 1) 
    } else {
      if (t_points[i] < breaks[3]) {
        x2[i] <- alpha_2[2] + c_23[2] * x3[i] + h_2[2] * hidden_t[i] + n_2[2] * rnorm(n = 1, mean = 0, sd = 1)
      } else {
        if (t_points[i] < breaks[4]) {
          x2[i] <- alpha_2[3] + c_23[3] * x3[i] + h_2[3] * hidden_t[i] + n_2[3] * rnorm(n = 1, mean = 0, sd = 1)
        } else {
          x2[i] <- alpha_2[1] + c_23[1] * x3[i] + h_2[1] * hidden_t[i] + n_2[3] * rnorm(n = 1, mean = 0, sd = 1) 
        }
      }
    }
  }
  y <- alpha_y + hidden_y * hidden_t * x2 + n_y * rnorm(n = sim_length, mean = 0, sd = 1)  
  x1 <- c(rep(NA, sim_length))
  for (i in 1:sim_length) {
    if (t_points[i] < breaks[2]) {
      x1[i] <- alpha_1[1] + c_13[1] * x3[i] + c_1y[1] * y[i] + n_1[1] * rnorm(n = 1, mean = 0, sd = 1) 
    } else {
      if (t_points[i] < breaks[3]) {
        x1[i] <- alpha_1[2] + c_13[2] * x3[i] + c_1y[2] * y[i] + n_1[2] * rnorm(n = 1, mean = 0, sd = 1)
      } else {
        if (t_points[i] < breaks[4]) {
          x1[i] <- alpha_1[3] + c_13[3] * x3[i] + c_1y[3] * y[i] + n_1[3] * rnorm(n = 1, mean = 0, sd = 1)
        } else {
          x1[i] <- alpha_1[1] + c_13[1] * x3[i] + c_1y[1] * y[i] + n_1[1] * rnorm(n = 1, mean = 0, sd = 1)
        }
      }
    }
  }
  out.list <- list(Y = y, X1 = x1, X2 = x2, X3 = x3, h = hidden_t, t = t_points)
  return(out.list)
}

# Simulating where each coefficient is constant under each environment
# with seperated influence of H and signal strength

# numbers: alpha_3, h_3, n_3, c_y2, alpha_y, hidden_y, n_y
# vectors of length three: alpha_2, c_23, h_2, n_2, alpha_1, c_13, c_1y, n_1

sim_data_fixed_env_2 <- function(sim_length, breaks, h_func, 
                               alpha_3, h_3, n_3,
                               alpha_2, c_23, h_2, n_2,
                               alpha_y, c_y2, hidden_y, n_y,
                               alpha_1, c_13, c_1y, n_1){
  t_points <- seq(from = breaks[1], to = breaks[length(breaks)], length.out = sim_length)
  hidden_t <- h_func(t_points)
  x3 <- alpha_3 + h_3 * hidden_t + n_3 * rnorm(n = sim_length, mean = 0, sd = 1)
  x2 <- c(rep(NA, sim_length))
  for (i in 1:sim_length) {
    if (t_points[i] < breaks[2]) {
      x2[i] <- alpha_2[1] + c_23[1] * x3[i] + h_2[1] * hidden_t[i] + n_2[1] * rnorm(n = 1, mean = 0, sd = 1) 
    } else {
      if (t_points[i] < breaks[3]) {
        x2[i] <- alpha_2[2] + c_23[2] * x3[i] + h_2[2] * hidden_t[i] + n_2[2] * rnorm(n = 1, mean = 0, sd = 1)
      } else {
        if (t_points[i] < breaks[4]) {
          x2[i] <- alpha_2[3] + c_23[3] * x3[i] + h_2[3] * hidden_t[i] + n_2[3] * rnorm(n = 1, mean = 0, sd = 1)
        } else {
          x2[i] <- alpha_2[1] + c_23[1] * x3[i] + h_2[1] * hidden_t[i] + n_2[3] * rnorm(n = 1, mean = 0, sd = 1) 
        }
      }
    }
  }
  y <- alpha_y + c_y2*x2 + hidden_y * hidden_t * x2 + n_y * rnorm(n = sim_length, mean = 0, sd = 1)  
  x1 <- c(rep(NA, sim_length))
  for (i in 1:sim_length) {
    if (t_points[i] < breaks[2]) {
      x1[i] <- alpha_1[1] + c_13[1] * x3[i] + c_1y[1] * y[i] + n_1[1] * rnorm(n = 1, mean = 0, sd = 1) 
    } else {
      if (t_points[i] < breaks[3]) {
        x1[i] <- alpha_1[2] + c_13[2] * x3[i] + c_1y[2] * y[i] + n_1[2] * rnorm(n = 1, mean = 0, sd = 1)
      } else {
        if (t_points[i] < breaks[4]) {
          x1[i] <- alpha_1[3] + c_13[3] * x3[i] + c_1y[3] * y[i] + n_1[3] * rnorm(n = 1, mean = 0, sd = 1)
        } else {
          x1[i] <- alpha_1[1] + c_13[1] * x3[i] + c_1y[1] * y[i] + n_1[1] * rnorm(n = 1, mean = 0, sd = 1)
        }
      }
    }
  }
  out.list <- list(Y = y, X1 = x1, X2 = x2, X3 = x3, h = hidden_t, t = t_points)
  return(out.list)
}


##### Simulating where coefficients are varying with time - environments change with high frequency, and not in jumps

sim_data_varying <- function(sim_length, alpha, w_0, alpha_3, h_3, n_3, n_2,
                             alpha_y, hidden_y, n_y, n_1, phases){
  t_points <- seq(from = 0, to = 1, length.out = sim_length)
  hidden_t <- sin(w_0 * t_points + phases[1])
  E_21 <- sin(alpha * w_0 * t_points + phases[2])
  E_22 <- sin(alpha * w_0 * t_points + phases[3])
  E_23 <- sin(alpha * w_0 * t_points + phases[4])
  E_11 <- sin(alpha * w_0 * t_points + phases[5])
  E_12 <- sin(alpha * w_0 * t_points + phases[6])
  E_13 <- sin(alpha * w_0 * t_points + phases[7])
  x3 <- alpha_3 + h_3 * hidden_t + n_3 * rnorm(n = sim_length, mean = 0, sd = 1)
  x2 <- E_21 + E_22 * x3 + E_23 * hidden_t + n_1 * rnorm(n = 1, mean = 0, sd = 1)
  y <- alpha_y + hidden_y * hidden_t * x2 + n_y * rnorm(n = sim_length, mean = 0, sd = 1)  
  x1 <- E_11 + E_12 * x3 + E_13 * y + n_1 * rnorm(n = 1, mean = 0, sd = 1)
  out.list <- list(Y = y, X1 = x1, X2 = x2, X3 = x3, h = hidden_t, t = t_points)
  return(out.list)
}

#### With seperated influence of H and signal strength

sim_data_varying_2 <- function(sim_length, alpha, w_0, alpha_3, h_3, n_3, n_2,
                             alpha_y, c_y2, hidden_y, n_y, n_1, phases){
  t_points <- seq(from = 0, to = 1, length.out = sim_length)
  hidden_t <- sin(w_0 * t_points + phases[1])
  E_21 <- sin(alpha * w_0 * t_points + phases[2])
  E_22 <- sin(alpha * w_0 * t_points + phases[3])
  E_23 <- sin(alpha * w_0 * t_points + phases[4])
  E_11 <- sin(alpha * w_0 * t_points + phases[5])
  E_12 <- sin(alpha * w_0 * t_points + phases[6])
  E_13 <- sin(alpha * w_0 * t_points + phases[7])
  x3 <- alpha_3 + h_3 * hidden_t + n_3 * rnorm(n = sim_length, mean = 0, sd = 1)
  x2 <- E_21 + E_22 * x3 + E_23 * hidden_t + n_1 * rnorm(n = 1, mean = 0, sd = 1)
  y <- alpha_y + c_y2 * x2 + hidden_y * hidden_t * x2 + n_y * rnorm(n = sim_length, mean = 0, sd = 1)  
  x1 <- E_11 + E_12 * x3 + E_13 * y + n_1 * rnorm(n = 1, mean = 0, sd = 1)
  out.list <- list(Y = y, X1 = x1, X2 = x2, X3 = x3, h = hidden_t, t = t_points)
  return(out.list)
}

#### With seperated influence of H and signal strength and with different frequencies for H and E

sim_data_varying_3 <- function(sim_length, w_H, w_E, alpha_3, h_3, n_3, n_2,
                               alpha_y, c_y2, hidden_y, n_y, n_1, phases){
  t_points <- seq(from = 0, to = 1, length.out = sim_length)
  hidden_t <- sin(w_H * t_points + phases[1])
  E_21 <- sin(w_E * t_points + phases[2])
  E_22 <- sin(w_E * t_points + phases[3])
  E_23 <- sin(w_E * t_points + phases[4])
  E_11 <- sin(w_E * t_points + phases[5])
  E_12 <- sin(w_E * t_points + phases[6])
  E_13 <- sin(w_E * t_points + phases[7])
  x3 <- alpha_3 + h_3 * hidden_t + n_3 * rnorm(n = sim_length, mean = 0, sd = 1)
  x2 <- E_21 + E_22 * x3 + E_23 * hidden_t + n_1 * rnorm(n = 1, mean = 0, sd = 1)
  y <- alpha_y + c_y2 * x2 + hidden_y * hidden_t * x2 + n_y * rnorm(n = sim_length, mean = 0, sd = 1)  
  x1 <- E_11 + E_12 * x3 + E_13 * y + n_1 * rnorm(n = 1, mean = 0, sd = 1)
  out.list <- list(Y = y, X1 = x1, X2 = x2, X3 = x3, h = hidden_t, t = t_points)
  return(out.list)
}


######## Chossing roughness penalty lambda by LOOCV

LOOCV_lambda <- function(data, lambda){
  MSE <- c(rep("NA", length(lambda)))
  for (j in seq_along(lambda)) {
    sq_error <- 0
    fit_i <- 0
    for (i in (2:(nrow(data) - 1))) {
      training_data <- data[-i,]
      beta_hat <- optim(par = rep(0, 2*nrow(data) - 2), fn = obj_fun_1,
                        data = training_data, pen_fun = beta_pen_2_sq, lambda = lambda[j],
                        gr = grad_pen_2_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                        control = list(maxit = 500))$par
      alpha_tilde <- (beta_hat[i] + beta_hat[i - 1]) / 2
      beta_tilde <- (beta_hat[nrow(data) + i - 1] + beta_hat[nrow(data) + i - 2]) / 2
      fit_i <- (data[i,1] - beta_tilde * data[i,2] - alpha_tilde)^2 + fit_i
    }
    MSE[j] <- fit_i / (nrow(data) - 2)
  }
  return(as.numeric(MSE))
}

choose_lambda <- function(data, lambda_val, n_pred){
  LOOCV_val <- 0
  for (i in 1:n_pred) {
    LOOCV_val <- LOOCV_val + LOOCV_lambda(data[,c(1,i + 1)], lambda_val)
  }
  return(list("Opt. lambda" = lambda_val[which.min(LOOCV_val)], "lambda_val" = lambda_val, "LOOCV_val" = LOOCV_val))
}


##### The same but with penalty beta_pen_1

LOOCV_lambda_1 <- function(data, lambda){
  MSE <- c(rep("NA", length(lambda)))
  for (j in seq_along(lambda)) {
    sq_error <- 0
    fit_i <- 0
    for (i in (2:(nrow(data) - 1))) {
      training_data <- data[-i,]
      beta_hat <- optim(par = rep(0, 2*nrow(data) - 2), fn = obj_fun_1,
                        data = training_data, pen_fun = beta_pen_sq, lambda = lambda[j],
                        gr = grad_pen_sq, lower = -10, upper = 10, method = "L-BFGS-B", 
                        control = list(maxit = 500))$par
      alpha_tilde <- (beta_hat[i] + beta_hat[i - 1]) / 2
      beta_tilde <- (beta_hat[nrow(data) + i - 1] + beta_hat[nrow(data) + i - 2]) / 2
      fit_i <- (data[i,1] - beta_tilde * data[i,2] - alpha_tilde)^2 + fit_i
    }
    MSE[j] <- fit_i / (nrow(data) - 2)
  }
  return(as.numeric(MSE))
}

choose_lambda_1 <- function(data, lambda_val, n_pred){
  LOOCV_val <- 0
  for (i in 1:n_pred) {
    LOOCV_val <- LOOCV_val + LOOCV_lambda_1(data[,c(1,i + 1)], lambda_val)
  }
  return(list("Opt. lambda" = lambda_val[which.min(LOOCV_val)], "lambda_val" = lambda_val, "LOOCV_val" = LOOCV_val))
}


###### Chossing bandwidth h by LOOCV

LOOCV_h <- function(data, h){
  MSE <- c(rep("NA", length(h)))
  for (j in seq_along(h)) {
    sq_error <- 0
    fit_i <- 0
    for (i in (2:(nrow(data) - 1))) {
      training_data <- data[-i,]
      model <- tvLM(formula = data[,1] ~ data[,2], bw = h[j])
      beta_hat <- model$tvcoef
      alpha_tilde <- (beta_hat[i, 1] + beta_hat[i - 1, 1]) / 2
      beta_tilde <- (beta_hat[i, 2] + beta_hat[i - 1, 2]) / 2
      fit_i <- (data[i,1] - beta_tilde * data[i,2] - alpha_tilde)^2 + fit_i
    }
    MSE[j] <- fit_i / (nrow(data) - 2)
  }
  return(as.numeric(MSE))
}

choose_h <- function(data, h_val, n_pred){
  LOOCV_val <- 0
  for (i in 1:n_pred) {
    LOOCV_val <- LOOCV_val + LOOCV_h(data[,c(1,i + 1)], h_val)
  }
  return(list("Opt. h" = h_val[which.min(LOOCV_val)], "h_val" = h_val, "LOOCV_val" = LOOCV_val))
}
