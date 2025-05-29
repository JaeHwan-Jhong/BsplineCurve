rm(list = ls())

library(dplyr)
library(ggplot2)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(cowplot)
library(rworldmap)
library(sphereplot)

## data time setting
start_date <- as.POSIXct("2018-08-09 00:00:00")
end_date <- as.POSIXct("2019-08-09 23:59:59")

## pre-preocessing
df = read.csv("./Cranes.csv")
data = df[, c("individual.local.identifier", "timestamp", "location.long","location.lat")]
data = data[order(data$timestamp), ]

## time data extraction
time = data["timestamp"]
time = as.matrix(time)

## POSIXct transform
timestamps <- as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S")
data <- data %>%
  filter(timestamps >= start_date & timestamps <= end_date)

## time pre-processing 
time = data["timestamp"]
time = as.matrix(time)
raw_data <- time
timestamps <- as.POSIXct(raw_data, format = "%Y-%m-%d %H:%M:%S")
start_time <- timestamps[1]
time_diffs <- as.numeric(difftime(timestamps, start_time, units = "hours"))

df_time <- data.frame(
  OriginalTime = timestamps,
  TimeElapsed = time_diffs
)

t = df_time$TimeElapsed
t = as.matrix(t)

# make t range 0~365
t = t / max(t) * 365
data$time = t

# define y
y1 = data["location.long"]
y2 = data["location.lat"]

y = NULL
y = cbind(y1, y2)
y = as.matrix(y)

n = length(t)


#################
##smooth-spline##
#################
library(stats)
set.seed(42)

# set repeats
num_repeats <- 50

mse_values <- numeric(num_repeats)
mae_values <- numeric(num_repeats)
mxdv_values <- numeric(num_repeats)

for (i in 1:num_repeats) {
  train_indices <- sample(length(t), floor(length(t) * 0.8))  
  test_indices <- setdiff(seq_along(t), train_indices)
  
  train_t <- t[train_indices]
  test_t <- t[test_indices]
  
  train_y <- y[train_indices, ]
  test_y <- y[test_indices, ]
  
  ss_lon <- smooth.spline(train_t, train_y[,1])
  ss_lat <- smooth.spline(train_t, train_y[,2])
  
  # prediction
  pred_lon <- predict(ss_lon, x = test_t)
  pred_lat <- predict(ss_lat, x = test_t)
  
  pred_combined <- cbind(pred_lon$y, pred_lat$y)
  
  point_errors <- sqrt(rowSums((test_y - pred_combined)^2))
  
  mse_values[i] <- mean(point_errors^2)
  mae_values[i] <- mean(point_errors)
  mxdv_values[i] <- max(point_errors)
}

# results
mse_mean <- mean(mse_values)
mse_sd <- sd(mse_values) / sqrt(num_repeats)

mae_mean <- mean(mae_values)
mae_sd <- sd(mae_values) / sqrt(num_repeats)

mxdv_mean <- mean(mxdv_values)
mxdv_sd <- sd(mxdv_values) / sqrt(num_repeats)

# print results
cat("MSE mean:", mse_mean, "Standard Error:", mse_sd, "\n")
cat("MAE mean:", mae_mean, "Standard Error:", mae_sd, "\n")
cat("MXDV mean:", mxdv_mean, "Standard Error:", mxdv_sd, "\n")


################
### k-smooth ###
################
set.seed(42)  

num_repeats <- 50

mse_values <- numeric(num_repeats)
mae_values <- numeric(num_repeats)
mxdv_values <- numeric(num_repeats)

library(stats)

for (i in 1:num_repeats) {
  train_indices <- sample(length(t), floor(length(t) * 0.8))  
  test_indices <- setdiff(seq_along(t), train_indices)
  
  train_t <- t[train_indices]
  test_t <- t[test_indices]
  
  train_y <- y[train_indices, ]
  test_y <- y[test_indices, ]
  
  # k-smooth
  ks_lon <- ksmooth(train_t, train_y[,1], kernel = "normal", bandwidth = 3, x.points = test_t)
  ks_lat <- ksmooth(train_t, train_y[,2], kernel = "normal", bandwidth = 3, x.points = test_t)

  pred_ksmooth <- cbind(ks_lon$y, ks_lat$y)
  
  point_errors_ksmooth <- sqrt(rowSums((test_y - pred_ksmooth)^2))
  
  mse_values[i] <- mean(point_errors_ksmooth^2)
  mae_values[i] <- mean(point_errors_ksmooth)
  mxdv_values[i] <- max(point_errors_ksmooth)
}

## results
mse_mean <- mean(mse_values)
mse_sd <- sd(mse_values) / sqrt(num_repeats)

mae_mean <- mean(mae_values)
mae_sd <- sd(mae_values) / sqrt(num_repeats)

mxdv_mean <- mean(mxdv_values)
mxdv_sd <- sd(mxdv_values) / sqrt(num_repeats)

# print results
cat("Kernel Smoothing - MSE mean:", mse_mean, "| Standard Error:", mse_sd, "\n")
cat("Kernel Smoothing - MAE mean:", mae_mean, "| Standard Error:", mae_sd, "\n")
cat("Kernel Smoothing - MXDV mean:", mxdv_mean, "| Standard Error:", mxdv_sd, "\n")


###################
## curve fitting ##
###################
set.seed(42)  
order = 3
dimension = 200

num_repeats <- 50

mse_values <- numeric(num_repeats)
mae_values <- numeric(num_repeats)
mxdv_values <- numeric(num_repeats)

n = length(t)
fixed_indices <- c(1, n)
remaining_indices <- 2:(n-1)

## compute test errors to individual (lon + lat)
for (i in 1:num_repeats) {

  n_train_target <- floor(0.8 * n)
  n_train_remaining <- n_train_target - length(fixed_indices)
  
  additional_train <- sort(sample(remaining_indices, size = n_train_remaining))
  train_indices <- sort(c(fixed_indices, additional_train))
  test_indices <- setdiff(1:n, train_indices)
  
  train_y <- y[train_indices, ]
  train_t <- t[train_indices]
  
  test_y <- y[test_indices, ]
  test_t <- t[test_indices]
  
  knots <- knots_quantile(train_t, dimension, order)
  B <- bsplines(train_t, knots, order)
  D <- bspline_jump(knots, order)
  
  Error_fit1 <- bspline.curve.admm_lambdas(train_y[,1], D, B, 
                                           lambdas = NULL, 
                                           lam_max = 10, 
                                           lam_min = 1e-5, 
                                           n_lambda = 200, 
                                           max_iter = 500, 
                                           epsilon = 1e-8, 
                                           eta_c = 1)
  
  Error_fit2 <- bspline.curve.admm_lambdas(train_y[,2], D, B, 
                                           lambdas = NULL, 
                                           lam_max = 10, 
                                           lam_min = 1e-5, 
                                           n_lambda = 200, 
                                           max_iter = 500, 
                                           epsilon = 1e-8, 
                                           eta_c = 1)
  
  B_test <- bsplines(test_t, knots, order)
  
  predict1 <- B_test %*% Error_fit1[[which.min(Error_fit1$bic)]]$xi
  predict2 <- B_test %*% Error_fit2[[which.min(Error_fit2$bic)]]$xi
  predict <- cbind(predict1, predict2)
  
  point_errors_lon <- sqrt(rowSums((test_y[,1] - predict1)^2))
  point_errors_lat <- sqrt(rowSums((test_y[,2] - predict2)^2))
  point_errors <- cbind(point_errors_lat, point_errors_lon)
  
  mse_values[i] <- mean(point_errors^2)
  mae_values[i] <- mean(point_errors)
  mxdv_values[i] <- max(point_errors)
  
}

# results
mse_mean <- mean(mse_values)
mse_sd <- sd(mse_values) / sqrt(num_repeats)

mae_mean <- mean(mae_values)
mae_sd <- sd(mae_values) / sqrt(num_repeats)

mxdv_mean <- mean(mxdv_values)
mxdv_sd <- sd(mxdv_values) / sqrt(num_repeats)

# print results
cat("B-spline ADMM (lon+lat individual) - MSE mean:", mse_mean, "| Standard Error:", mse_sd)
cat("B-spline ADMM (lon+lat individual) - MAE mean:", mae_mean, "| Standard Error:", mae_sd)
cat("B-spline ADMM (lon+lat individual) - MXDV mean:", mxdv_mean, "| Standard Error:", mxdv_sd)


## compute all test errors
mse_all_values <- numeric(num_repeats)
mae_all_values <- numeric(num_repeats)
mxdv_all_values <- numeric(num_repeats)

for (i in 1:num_repeats) {
  n_train_target <- floor(0.8 * n)
  n_train_remaining <- n_train_target - length(fixed_indices)
  
  additional_train <- sort(sample(remaining_indices, size = n_train_remaining))
  train_indices <- sort(c(fixed_indices, additional_train))
  test_indices <- setdiff(1:n, train_indices)
  
  train_y <- y[train_indices, ]
  train_t <- t[train_indices]
  
  test_y <- y[test_indices, ]
  test_t <- t[test_indices]
  
  knots <- knots_quantile(train_t, dimension, order)
  B <- bsplines(train_t, knots, order)
  D <- bspline_jump(knots, order)
  
  Error_fit <- bspline.curve.admm_lambdas(train_y, D, B, 
                                          lambdas = NULL, 
                                          lam_max = 10, 
                                          lam_min = 1e-5, 
                                          n_lambda = 200, 
                                          max_iter = 500, 
                                          epsilon = 1e-8, 
                                          eta_c = 1)
  
  B_test <- bsplines(test_t, knots, order)
  all_predict <- B_test %*% Error_fit[[which.min(Error_fit$bic)]]$xi
  
  all_point_errors <- sqrt(rowSums((test_y - all_predict)^2))
  mse_all_values[i] <- mean(all_point_errors^2)
  mae_all_values[i] <- mean(all_point_errors)
  mxdv_all_values[i] <- max(all_point_errors)
}

mse_all_mean <- mean(mse_all_values)
mse_all_sd <- sd(mse_all_values) / sqrt(num_repeats)

mae_all_mean <- mean(mae_all_values)
mae_all_sd <- sd(mae_all_values) / sqrt(num_repeats)

mxdv_all_mean <- mean(mxdv_all_values)
mxdv_all_sd <- sd(mxdv_all_values) / sqrt(num_repeats)

cat("B-spline ADMM (all train_y) - MSE mean:", mse_all_mean, "| Standard Error:", mse_all_sd)
cat("B-spline ADMM (all train_y) - MAE mean:", mae_all_mean, "| Standard Error:", mae_all_sd)
cat("B-spline ADMM (all train_y) - MXDV mean:", mxdv_all_mean, "| Standard Error:", mxdv_all_sd)


#################
### gam-model ###
#################
library(mgcv)

# set iterations
n_iterations <- 50
set.seed(42)

mse_values <- numeric(n_iterations)
mae_values <- numeric(n_iterations)
mxdv_values <- numeric(n_iterations)
n = length(t)


for(iter in 1:n_iterations) {
  
  train_idx <- sample(1:n, size = round(0.2 * n), replace = FALSE)
  test_idx <- setdiff(1:n, train_idx)
  
  t_train <- t[train_idx]
  y_train <- y[train_idx, ]
  t_test <- t[test_idx]
  y_test <- y[test_idx, ]
  
  train_data <- data.frame(
    t = t_train,
    lon = y_train[,1],
    lat = y_train[,2]
  )
  test_data <- data.frame(t = t_test)
  
  k_candidates <- c(5, 7, 9, 11, 13, 15)
  n_train <- nrow(train_data)
  
  bic_values_lon <- numeric(length(k_candidates))
  bic_values_lat <- numeric(length(k_candidates))
  models_lon <- list()
  models_lat <- list()
  
  for(i in seq_along(k_candidates)) {
    k_val <- k_candidates[i]
    
    model_lon <- gam(lon ~ s(t, k = k_val), data = train_data, method = "REML")
    model_lat <- gam(lat ~ s(t, k = k_val), data = train_data, method = "REML")
    
    models_lon[[i]] <- model_lon
    models_lat[[i]] <- model_lat
    
    bic_values_lon[i] <- AIC(model_lon, k = log(n_train))
    bic_values_lat[i] <- AIC(model_lat, k = log(n_train))
  }
  
  # select best k
  best_index_lon <- which.min(bic_values_lon)
  best_index_lat <- which.min(bic_values_lat)
  
  best_model_lon <- models_lon[[best_index_lon]]
  best_model_lat <- models_lat[[best_index_lat]]
  
  pred_gam_lon <- predict(best_model_lon, newdata = test_data)
  pred_gam_lat <- predict(best_model_lat, newdata = test_data)
  pred_gam <- cbind(pred_gam_lon, pred_gam_lat)
  
  gam_point_errors <- sqrt(rowSums((y_test - pred_gam)^2))
  
  mse_values[iter] <- mean(gam_point_errors^2)
  mae_values[iter] <- mean(gam_point_errors)
  mxdv_values[iter] <- max(gam_point_errors)
}

avg_mse <- mean(mse_values)
sd_mse <- sd(mse_values) / sqrt(num_repeats)

avg_mae <- mean(mae_values)
sd_mae <- sd(mae_values) / sqrt(num_repeats)

avg_mxdv <- mean(mxdv_values)
sd_mxdv <- sd(mxdv_values) / sqrt(num_repeats)

# print results
cat("MSE mean:", avg_mse, "Standard Error:", sd_mse, "\n")
cat("MAE mean:", avg_mae, "Standard Error:", sd_mae, "\n")
cat("MXDV mean:", avg_mxdv, "Standard Error:", sd_mxdv, "\n")


