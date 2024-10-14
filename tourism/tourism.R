library(tidyverse)
library(tsibble)
library(hts)
library(fabletools)
# remotes::install_github("ShanikaLW/ihts", auth_token = "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")
library(ihts)
library(quadprog)
library(foreach)
library(parallel)
library(doParallel)

#################### Functions ################################

### Function to generate ETS forecasts ###
get_forecasts_ets <- function(data, initial_train_size = 96, h = NULL){
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  forecasts <- list()
  residuals <- list()
  actuals <- list()
  
  for(i in 1:(n - initial_train_size)){
    train_end <- initial_train_size + i - 1
    
    train_data <- data[1:train_end, ]
    forecast_list <- matrix(NA, nrow = h, ncol = p)
    residuals_list <- matrix(NA, nrow = train_end, ncol = p)
    actual_data_list <- matrix(NA, nrow = h, ncol = p)
    test_set <- data[-(1:train_end), , drop = FALSE]
    idx <- min(h, nrow(test_set))
    actual_data_list[1:idx, ] <- test_set[1:idx, ]
    
    forecast_results <- foreach(j = 1:p, .packages = c("forecast"),
                                .combine = 'c') %dopar% {
      series <- ts(train_data[, j], frequency = 12)
      fit <- ets(series)
      fcast <- forecast(fit, h = h)
      
      list(forecasts = fcast$mean, 
           residuals = train_data[, j] - fitted(fit))
    }
    
    for (j in 1:p) {
      forecast_list[, j] <- forecast_results[[2 * j - 1]]
      residuals_list[, j] <- forecast_results[[2 * j]]
    }
    
    forecasts[[i]] <- forecast_list
    residuals[[i]] <- residuals_list
    actuals[[i]] <- actual_data_list
  }
  
  stopCluster(cl)
  
  return(list(
    forecasts = forecasts, 
    residuals = residuals, 
    actuals = actuals
  ))
}



### Function to generate ARIMA forecasts ###
get_forecasts_arima <- function(data, initial_train_size = 96, h = NULL){
  
  n <- dim(data)[1]
  p <- dim(data)[2]
  
  num_cores <- detectCores() - 1
  cl <- makeCluster(num_cores)
  registerDoParallel(cl)
  
  forecasts <- list()
  residuals <- list()
  actuals <- list()
  
  for(i in 1:(n - initial_train_size)){
    train_end <- initial_train_size + i - 1
    
    train_data <- data[1:train_end, ]
    forecast_list <- matrix(NA, nrow = h, ncol = p)
    residuals_list <- matrix(NA, nrow = train_end, ncol = p)
    actual_data_list <- matrix(NA, nrow = h, ncol = p)
    test_set <- data[-(1:train_end), , drop = FALSE]
    idx <- min(h, nrow(test_set))
    actual_data_list[1:idx, ] <- test_set[1:idx, ]
    
    forecast_results <- foreach(j = 1:p, .packages = c("forecast"),
                                .combine = 'c') %dopar% {
      series <- ts(train_data[, j], frequency = 12)
      fit <- auto.arima(series)
      fcast <- forecast(fit, h = h)
      
      list(forecasts = fcast$mean, 
           residuals = fit[["residuals"]])
    }
    
    for (j in 1:p) {
      forecast_list[, j] <- forecast_results[[2 * j - 1]]
      residuals_list[, j] <- forecast_results[[2 * j]]
    }
    
    forecasts[[i]] <- forecast_list
    residuals[[i]] <- residuals_list
    actuals[[i]] <- actual_data_list
  }
  
  stopCluster(cl)
  
  return(list(
    forecasts = forecasts, 
    residuals = residuals, 
    actuals = actuals
  ))
}


### Function for unconstrained reconciliation ###
reconcile_f <- function(basef_list, nodes, groups, resid_list, smat) {
  
  ols_results <- list()
  wlss_results <- list()
  wlsv_results <- list()
  shrink_results <- list()
  
  # Iterate over each set of forecasts and residuals
  for (i in seq_along(basef_list)) {
    basef <- basef_list[[i]]
    resid <- resid_list[[i]]
    
    # Perform the reconciliation for each method
    ols_results[[i]] <- combinef(fcasts = basef, nodes = NULL, groups = groups,
                                 weights = NULL, nonnegative = FALSE,
                                 algorithms = "lu", keep = "all",
                                 parallel = TRUE)
    
    wlss_weight <- as.vector(smat %*% matrix(1, nrow = dim(smat)[2]))
    wlss_results[[i]] <- combinef(fcasts = basef, nodes = NULL, groups = groups,
                                  weights = 1 / wlss_weight, nonnegative = FALSE,
                                  algorithms = "lu", keep = "all",
                                  parallel = TRUE)
    
    wlsv_weight <- diag(crossprod(resid) / nrow(resid))
    wlsv_results[[i]] <- combinef(fcasts = basef, nodes = NULL, groups = groups,
                                  weights = 1 / wlsv_weight, nonnegative = FALSE,
                                  algorithms = "lu", keep = "all", 
                                  parallel = TRUE)
    
    shrink_results[[i]] <- MinT(fcasts = basef, nodes = NULL, groups = groups,
                                nonnegative = FALSE, residual = resid,
                                covariance = "shr", algorithms = "lu",
                                keep = "all", parallel = TRUE)
  }
  reconciled_results <- list(
    basef = basef_list, 
    ols_U = ols_results,
    wlss_U = wlss_results,
    wlsv_U = wlsv_results,
    shrink_U = shrink_results
  )
  return(reconciled_results)
}


### Zhang et al. method (I) ###

qp_i <- function(basef, immu_set, smat, weights){
  
  quad_sol <- matrix(0, dim(basef)[1], dim(basef)[2])
  
  Dmat <- t(smat) %*% solve(weights) %*% smat
  Amat <- as.matrix(smat[immu_set, ])
  for (i in 1:dim(basef)[1]){
    dvec <- as.numeric(basef[i,]) %*% solve(weights) %*% smat
    bvec <- as.vector(basef[i,immu_set])
    if (length(immu_set)==1){
      quad_sol[i,] <- smat %*% quadprog::solve.QP(Dmat=t(Dmat), dvec=dvec, 
                                                  Amat=Amat, bvec=bvec, 
                                                  meq = length(immu_set))$solution
    }
    else {
      quad_sol[i,] <- smat %*% quadprog::solve.QP(Dmat=t(Dmat), dvec=dvec, 
                                                  Amat=t(Amat), bvec=bvec, 
                                                  meq = length(immu_set))$solution
    }
  }
  return(quad_sol)
}


reconcile_f_qpi <- function(basef, immu_set, smat, resid) {
  
  # Function to process a single basef and resid pair
  process_pair <- function(basef, resid) {
    ols_weight <- diag(nrow(smat))
    ols_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = ols_weight)
    
    wlss_weight <- diag(as.vector(smat %*% rep(1, ncol(smat))))
    wlss_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = wlss_weight)
    
    wlsv_weight <- diag(diag(crossprod(resid) / nrow(resid)))
    wlsv_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = wlsv_weight)
    
    tar <- hts:::lowerD(resid)
    shrink <- hts:::shrink.estim(resid, tar)
    shrink_weight <- shrink[[1]]
    shrink_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = shrink_weight)
    
    list(
      basef = basef,
      ols_qpi = ols_qpi,
      wlss_qpi = wlss_qpi,
      wlsv_qpi = wlsv_qpi,
      shrink_qpi = shrink_qpi
    )
  }
  results_list <- vector("list", length(basef))
  
  # Process each pair of basef and resid
  for (i in seq_along(basef)) {
    results_list[[i]] <- process_pair(basef[[i]], resid[[i]])
  }
  
  ols_results <- lapply(results_list, `[[`, "ols_qpi")
  wlss_results <- lapply(results_list, `[[`, "wlss_qpi")
  wlsv_results <- lapply(results_list, `[[`, "wlsv_qpi")
  shrink_results <- lapply(results_list, `[[`, "shrink_qpi")
  basef <- lapply(results_list, `[[`, "basef")
  
  list(
    basef = basef,
    ols_I = ols_results,
    wlss_I = wlss_results,
    wlsv_I = wlsv_results,
    shrink_I = shrink_results
  )
}



### I-MinT ###

reconcile_f_imint <- function(basef_list, nodes, groups, resid_list, immu_set, smat) {
  
  ols_iMinT_results <- list()
  wlss_iMinT_results <- list()
  wlsv_iMinT_results <- list()
  shrink_iMinT_results <- list()
  
  # Iterate over each set of forecasts and residuals
  for (i in seq_along(basef_list)) {
    basef <- basef_list[[i]]
    resid <- resid_list[[i]]
    
    ols_iMinT_results[[i]] <- icombinef(fcasts = basef, nodes = NULL, groups = groups,
                                        weights = NULL, immute = immu_set,
                                        nonnegative = FALSE,
                                        algorithms = "lu",
                                        keep = "all")
    
    wlss_weight <- as.vector(smat %*% matrix(1, nrow = dim(smat)[2]))
    wlss_iMinT_results[[i]] <- icombinef(fcasts = basef, nodes = NULL, groups = groups,
                                         weights = 1 / wlss_weight, immute = immu_set,
                                         nonnegative = FALSE,
                                         algorithms = "lu",
                                         keep = "all")
    
    wlsv_weight <- diag(crossprod(resid) / nrow(resid))
    wlsv_iMinT_results[[i]] <- icombinef(fcasts = basef, nodes = NULL, groups = groups,
                                         weights = 1 / wlsv_weight, immute = immu_set,
                                         nonnegative = FALSE,
                                         algorithms = "lu",
                                         keep = "all")
    
    shrink_iMinT_results[[i]] <- iMinT(fcasts = basef, nodes = NULL, groups = groups,
                                       immute = immu_set, nonnegative = FALSE,
                                       residual = resid,
                                       covariance = "shr",
                                       algorithms = "lu",
                                       keep = "all")
  }
  reconciled_results <- list(
    basef = basef_list, 
    ols_iMinT = ols_iMinT_results,
    wlss_iMinT = wlss_iMinT_results,
    wlsv_iMinT = wlsv_iMinT_results,
    shrink_iMinT = shrink_iMinT_results
  )
  return(reconciled_results)
}


### Accuracy checks ###

# accuracy_check_v1 <- function(forecasts_list, test_sets, horizon_range = NULL) {
#   mse_results <- list()
#   
#   # Iterate over each forecast set
#   for (forecast_set_name in names(forecasts_list)) {
#     forecasts <- forecasts_list[[forecast_set_name]]
#     test_set_list <- test_sets[[forecast_set_name]]
#     
#     mse_matrix <- matrix(NA, nrow = length(test_set_list), ncol = 7)
#     colnames(mse_matrix) <- c("Total", "State", "Region", "Purpose", 
#                               "State_Purpose", "Region_Purpose", "All_series")
#     
#     # Iterate over each test set
#     for (i in seq_along(test_set_list)) {
#       actuals <- test_set_list[[i]]
#       
#       # Check if there are NAs in the relevant horizon range
#       if (any(is.na(actuals[horizon_range, ]))) {
#         next  # Skip MSE calculation if there are NAs in the specified horizon range
#       }
#       
#       num_cols <- ncol(actuals)
#       mse_values <- numeric(num_cols)
#       
#       # Calculate MSE for each series
#       for (j in 1:num_cols) {
#         # if (is.null(horizon_range)) {
#         #   mse_values[j] <- mean((actuals[, j] - forecasts[[i]][, j])^2, na.rm = TRUE)
#         # } else {
#           mse_values[j] <- mean((actuals[horizon_range, j] - forecasts[[i]][horizon_range, j])^2, na.rm = TRUE)
#         # }
#       }
#       
#       mse_categories <- c(
#         Total = mean(mse_values[1], na.rm = TRUE),
#         State = mean(mse_values[2:8], na.rm = TRUE),
#         Region = mean(mse_values[9:85], na.rm = TRUE),
#         Purpose = mean(mse_values[86:89], na.rm = TRUE),
#         State_Purpose = mean(mse_values[90:117], na.rm = TRUE),
#         Region_Purpose = mean(mse_values[118:425], na.rm = TRUE),
#         All_series = mean(mse_values[1:425], na.rm = TRUE)
#       )
#       mse_matrix[i, ] <- mse_categories
#     }
#     avg_mse <- colMeans(mse_matrix, na.rm = TRUE)  # Average MSEs across all test sets for this forecast set
#     mse_results[[forecast_set_name]] <- avg_mse
#   }
#   
#   result_df <- do.call(rbind, mse_results)
#   rownames(result_df) <- names(mse_results)
#   return(result_df)
# }
# 
# 
# mse_imint_ets1to6 <- accuracy_check(forecasts_list = etsf_imint, test_sets = ets_test_set, 
#                                     horizon_range = 1:6)

accuracy_check <- function(forecasts, test_sets, h = NULL) {
  mse_results <- list()
  # for base forecasts and reconciliation with different W estimators
  for (type in names(forecasts)) {
    forecasts_list <- forecasts[[type]]
    test_set_list <- test_sets[[type]]
    squared_errors <- list()
  
        # for each iterations (window)
        for (i in seq_along(test_set_list)) {
            actuals <- test_set_list[[i]]
            fcasts <- forecasts_list[[i]]
            squared_errors[[i]] <- (actuals - fcasts)^2
        }
  
    squared_errors <- simplify2array(squared_errors)
    mse <- apply(squared_errors, MARGIN = c(1, 2), mean, na.rm = TRUE)
    
    # average mse across each level
    mse_levels <- cbind(Total = rowMeans(as.matrix(mse[,1]), na.rm = TRUE),
      State = rowMeans(mse[,2:8], na.rm = TRUE),
      Region = rowMeans(mse[,9:85], na.rm = TRUE),
      Purpose = rowMeans(mse[,86:89], na.rm = TRUE),
      State_Purpose = rowMeans(mse[,90:117], na.rm = TRUE),
      Region_Purpose = rowMeans(mse[,118:425], na.rm = TRUE),
      All_series = rowMeans(mse[,1:425], na.rm = TRUE)
    )
    
    # average mse for the h-step-ahead
    mse_h <- colMeans(mse_levels[h, , drop = FALSE], na.rm = TRUE)
    
    mse_results[[type]] <- mse_h
  }
  result_df <- do.call(rbind, mse_results)
  rownames(result_df) <- names(mse_results)
  return(result_df)
}

# mse_test <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h = 1:6)

  
#################### Data importing and setting up #############################

header <- read_csv("tourism.csv", n_max = 3, col_names = FALSE)

header <- header %>%
  t() %>%
  as_tibble() %>%
  fill(V1:V2, .direction = "down") %>%
  unite(name, V1, V2, V3, na.rm = TRUE, sep = "/") %>%
  pull()

tourism <- read_csv("tourism.csv", skip = 3, col_names = header)

tourism <- tourism %>%
  fill(Year, .direction = "down") 

tourism <- tourism %>%
  pivot_longer(!Year:Month, names_to = "Variable", values_to = "Trips") %>%
  separate(Variable, c("State", "Region", "Purpose"), sep = "/")

tourism <- tourism %>%
  mutate(Year_month = yearmonth(ym(paste(Year, Month)))) %>%
  as_tsibble(key = c(State, Region, Purpose), index = Year_month) %>%
  dplyr::select(Year_month, State:Trips)

tourism <- tourism %>%
  mutate(State = recode(State,
                        #   `Australian Capital Territory` = "ACT",
                        `New South Wales` = "NSW",
                        `Northern Territory` = "N.T",
                        `Queensland` = "QLD",
                        `South Australia` = "S.A",
                        `Tasmania` = "TAS",
                        `Victoria` = "VIC",
                        `Western Australia` = "W.A"
  ),
  Purpose = recode(Purpose,
                   `Business` = "BUS",
                   `Visiting friends and relatives` = "VFR",
                   `Holiday` = "HOL",
                   `Other reason` = "OTH"))


# taking one period for ease to extract summing matrix and groups
tourjan1998 <- tourism %>% filter(Year_month == yearmonth("1998 Jan"))

tourjan1998 <- tourjan1998 %>%
  mutate(State_code = case_when(
    # State == "ACT" ~ 1,
    State == "NSW" ~ 1,
    State == "N.T" ~ 2,
    State == "QLD" ~ 3, 
    State == "S.A" ~ 4,
    State == "TAS" ~ 5,
    State == "VIC" ~ 6,
    State == "W.A" ~ 7,
  ))

state_counts <- tourjan1998 %>%
  group_by(State) %>%
  summarize(count = n())


gn <- rbind(tourjan1998$State_code, 
            rep(1:77, each = 4),
            rep(1:4, times = 77),
            cbind(t(rep(1:4, times = 56/4)), 
                  t(rep(5:8, times = 28/4)),
                  t(rep(9:12, times = 52/4)),
                  t(rep(13:16, times = 48/4)),
                  t(rep(17:20, times = 20/4)),
                  t(rep(21:24, times = 84/4)),
                  t(rep(25:28, times = 20/4))
            ))

tourism_bts <- tourism %>% 
  unite(StateRegionPurpose, State, Region, Purpose, sep = "_") %>% 
  pivot_wider(names_from = StateRegionPurpose, values_from = Trips) %>% 
  dplyr::select(-Year_month)


tourism_gts <- gts(tourism_bts, groups = gn, gnames = c("State", "Region", "Purpose", "StatePurpose"))


smat <- smatrix(tourism_gts)
groups <- get_groups(tourism_gts)

tourism_full <- t(smat %*% t(tourism_bts)) #computing full structure
immutable_set <- c(1,2,4,7) #immutable series 

######################## Base Forecasts ##############################

# Takes long, run only is new model fitting and base forecasts are needed.

### ETS forecasts
# start_time <- Sys.time()
# 
# ets_forecasts <- get_forecasts_ets(data=tourism_full, initial_train_size = 96, h = 12)
# 
# end_time <- Sys.time()
# time_taken_ets <- end_time-start_time
# 
# saveRDS(ets_forecasts, file="etsf.rds")


### ARIMA forecasts
# start_time <- Sys.time()
# 
# arima_forecasts <- get_forecasts_arima(data=tourism_full, initial_train_size = 96, h = 12)
# 
# end_time <- Sys.time()
# time_taken_arima <- end_time-start_time
# 
# saveRDS(arima_forecasts, file="arimaf.rds")


###################### Reconciliation ################################

ets_forecasts <- readRDS("etsf.rds")
basef_ets <- ets_forecasts$forecasts
resid_ets <- ets_forecasts$residuals
actuals_ets <- ets_forecasts$actuals


u_ets <- reconcile_f(basef_list=basef_ets, nodes=NULL, groups=groups, 
                      resid_list=resid_ets, smat=smat)
i_ets <- reconcile_f_qpi(basef=basef_ets, immu_set=immutable_set, 
                         smat=smat, resid=resid_ets)
iMinT_ets <- reconcile_f_imint(basef_list = basef_ets, nodes = NULL, 
                              groups = groups, resid_list = resid_ets, 
                              immu_set = immutable_set, smat = smat)


arima_forecasts <- readRDS("arimaf.rds")
basef_arima <- arima_forecasts$forecasts
resid_arima <- arima_forecasts$residuals
actuals_arima <- arima_forecasts$actuals

u_arima <- reconcile_f(basef_list=basef_arima, nodes=NULL, groups=groups, 
                       resid_list=resid_arima, smat=smat)
i_arima <- reconcile_f_qpi(basef=basef_arima, immu_set=immutable_set, 
                         smat=smat, resid=resid_arima)
iMinT_arima <- reconcile_f_imint(basef_list = basef_arima, nodes = NULL, 
                               groups = groups, resid_list = resid_arima, 
                               immu_set = immutable_set, smat = smat)




###################### Accuracy ################################

##### ETS #####
ets_test_set <- list(basef = actuals_ets, ols = actuals_ets,
                wlss = actuals_ets, wlsv = actuals_ets,
                shrink = actuals_ets)

etsf_U <- list(basef = basef_ets, ols = u_ets[["ols_U"]],
               wlss = u_ets[["wlss_U"]], wlsv = u_ets[["wlsv_U"]],
               shrink = u_ets[["shrink_U"]])

etsf_i <- list(basef = basef_ets, ols = i_ets[["ols_I"]],
               wlss = i_ets[["wlss_I"]], wlsv = i_ets[["wlsv_I"]],
               shrink = i_ets[["shrink_I"]])


etsf_imint <- list(basef = basef_ets, ols= iMinT_ets[["ols_iMinT"]],
                   wlss = iMinT_ets[["wlss_iMinT"]], wlsv = iMinT_ets[["wlsv_iMinT"]],
                   shrink = iMinT_ets[["shrink_iMinT"]])

row_names <- c("basef", "ols_U", "wlss_U", "wlsv_U", "shrink_U", 
               "ols_I", "wlss_I", "wlsv_I", "shrink_I", "ols_iMinT",
               "wlss_iMinT", "wlsv_iMinT", "shrink_iMinT")

### h=1
mse_u_ets1 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=1)
mse_i_ets1 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=1)
mse_imint_ets1 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=1)

mse_ets1 <- rbind(mse_u_ets1, mse_i_ets1[-1,], mse_imint_ets1[-1,])
rownames(mse_ets1) <- row_names

basef_mse_ets1 <- mse_ets1[1, ]
pc1 <- mse_ets1[-(1:5), ]
pc1 <- sweep(pc1, 2, basef_mse_ets1, FUN = function(x, y) ((x - y) / y) * 100)
pc1 <- rbind(basef_mse_ets1/(10^3), pc1)
rownames(pc1) <- rownames(mse_ets1[-(2:5),])
round(pc1,3)

as.data.frame(pc1) %>% 
  write.csv("pc/pc1.csv", row.names = TRUE)

### h=3
mse_u_ets3 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=3)
mse_i_ets3 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=3)
mse_imint_ets3 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=3)

mse_ets3 <- rbind(mse_u_ets3, mse_i_ets3[-1,], mse_imint_ets3[-1,])
rownames(mse_ets3) <- row_names

basef_mse_ets3 <- mse_ets3[1, ]
pc3 <- mse_ets3[-(1:5), ]
pc3 <- sweep(pc3, 2, basef_mse_ets3, FUN = function(x, y) ((x - y) / y) * 100)
pc3 <- rbind(basef_mse_ets3/(10^3), pc3)
rownames(pc3) <- rownames(mse_ets3[-(2:5),])
round(pc3,3)

as.data.frame(pc3) %>% 
  write.csv("pc/pc3.csv", row.names = TRUE)

### h=1-3
mse_u_ets1to3 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=1:3)
mse_i_ets1to3 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=1:3)
mse_imint_ets1to3 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=1:3)

mse_ets1to3 <- rbind(mse_u_ets1to3, mse_i_ets1to3[-1,], mse_imint_ets1to3[-1,])
rownames(mse_ets1to3) <- row_names

basef_mse_ets1to3 <- mse_ets1to3[1, ]
pc1to3 <- mse_ets1to3[-(1:5), ]
pc1to3 <- sweep(pc1to3, 2, basef_mse_ets1to3, FUN = function(x, y) ((x - y) / y) * 100)
pc1to3 <- rbind(basef_mse_ets1to3/(10^3), pc1to3)
rownames(pc1to3) <- rownames(mse_ets1to3[-(2:5),])
round(pc1to3,3)

as.data.frame(pc1to3) %>% 
  write.csv("pc/pc1to3.csv", row.names = TRUE)

### h=6
mse_u_ets6 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=6)
mse_i_ets6 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=6)
mse_imint_ets6 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=6)

mse_ets6 <- rbind(mse_u_ets6, mse_i_ets6[-1,], mse_imint_ets6[-1,])
rownames(mse_ets6) <- row_names

basef_mse_ets6 <- mse_ets6[1, ]
pc6 <- mse_ets6[-(1:5), ]
pc6 <- sweep(pc6, 2, basef_mse_ets6, FUN = function(x, y) ((x - y) / y) * 100)
pc6 <- rbind(basef_mse_ets6/(10^3), pc6)
rownames(pc6) <- rownames(mse_ets6[-(2:5),])
round(pc6,3)

as.data.frame(pc6) %>% 
  write.csv("pc/pc6.csv", row.names = TRUE)


### h=1-6
mse_u_ets1to6 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=1:6)
mse_i_ets1to6 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=1:6)
mse_imint_ets1to6 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=1:6)

mse_ets1to6 <- rbind(mse_u_ets1to6, mse_i_ets1to6[-1,], mse_imint_ets1to6[-1,])
rownames(mse_ets1to6) <- row_names

basef_mse_ets1to6 <- mse_ets1to6[1, ]
pc1to6 <- mse_ets1to6[-(1:5), ]
pc1to6 <- sweep(pc1to6, 2, basef_mse_ets1to6, FUN = function(x, y) ((x - y) / y) * 100)
pc1to6 <- rbind(basef_mse_ets1to6/(10^3), pc1to6)
rownames(pc1to6) <- rownames(mse_ets1to6[-(2:5),])
round(pc1to6,3)

as.data.frame(pc1to6) %>% 
  write.csv("pc/pc1to6.csv", row.names = TRUE)

### h=12
mse_u_ets12 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=12)
mse_i_ets12 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=12)
mse_imint_ets12 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=12)

mse_ets12 <- rbind(mse_u_ets12, mse_i_ets12[-1,], mse_imint_ets12[-1,])
rownames(mse_ets12) <- row_names

basef_mse_ets12 <- mse_ets12[1, ]
pc12 <- mse_ets12[-(1:5), ]
pc12 <- sweep(pc12, 2, basef_mse_ets12, FUN = function(x, y) ((x - y) / y) * 100)
pc12 <- rbind(basef_mse_ets12/(10^3), pc12)
rownames(pc12) <- rownames(mse_ets12[-(2:5),])
round(pc12,3)

as.data.frame(pc12) %>% 
  write.csv("pc/pc12.csv", row.names = TRUE)


### h=1-12
mse_u_ets1to12 <- accuracy_check(forecasts=etsf_U, test_sets=ets_test_set, h=1:12)
mse_i_ets1to12 <- accuracy_check(forecasts=etsf_i, test_sets=ets_test_set, h=1:12)
mse_imint_ets1to12 <- accuracy_check(forecasts=etsf_imint, test_sets=ets_test_set, h=1:12)

mse_ets1to12 <- rbind(mse_u_ets1to12, mse_i_ets1to12[-1,], mse_imint_ets1to12[-1,])
rownames(mse_ets1to12) <- row_names

basef_mse_ets1to12 <- mse_ets1to12[1, ]
pc1to12 <- mse_ets1to12[-(1:5), ]
pc1to12 <- sweep(pc1to12, 2, basef_mse_ets1to12, FUN = function(x, y) ((x - y) / y) * 100)
pc1to12 <- rbind(basef_mse_ets1to12/(10^3), pc1to12)
rownames(pc1to12) <- rownames(mse_ets1to12[-(2:5),])
round(pc1to12,3)

as.data.frame(pc1to12) %>% 
  write.csv("pc/pc1to12.csv", row.names = TRUE)




##### ARIMA #####

arima_test_set <- list(basef = actuals_arima, ols = actuals_arima,
                     wlss = actuals_arima, wlsv = actuals_arima,
                     shrink = actuals_arima)

arimaf_U <- list(basef = basef_arima, ols = u_arima[["ols_U"]],
               wlss = u_arima[["wlss_U"]], wlsv = u_arima[["wlsv_U"]],
               shrink = u_arima[["shrink_U"]])

arimaf_i <- list(basef = basef_arima, ols = i_arima[["ols_I"]],
               wlss = i_arima[["wlss_I"]], wlsv = i_arima[["wlsv_I"]],
               shrink = i_arima[["shrink_I"]])


arimaf_imint <- list(basef = basef_arima, ols= iMinT_arima[["ols_iMinT"]],
                   wlss = iMinT_arima[["wlss_iMinT"]], wlsv = iMinT_arima[["wlsv_iMinT"]],
                   shrink = iMinT_arima[["shrink_iMinT"]])


#### h=1
mse_u_arima1 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=1)
mse_i_arima1  <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=1)
mse_imint_arima1  <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=1)

mse_arima1 <- rbind(mse_u_arima1, mse_i_arima1[-1,], mse_imint_arima1[-1,])
rownames(mse_arima1) <- row_names

basef_mse_arima1 <- mse_arima1[1, ]
pc1_arima <- mse_arima1[-(1:5), ]
pc1_arima <- sweep(pc1_arima, 2, basef_mse_arima1, FUN = function(x, y) ((x - y) / y) * 100)
pc1_arima <- rbind(basef_mse_arima1/(10^3), pc1_arima)
rownames(pc1_arima) <- rownames(mse_arima1[-(2:5),])
round(pc1_arima,3)

as.data.frame(pc1_arima) %>% 
  write.csv("pc/pc1_arima.csv", row.names = TRUE)

### h=3
mse_u_arima3 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=3)
mse_i_arima3 <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=3)
mse_imint_arima3 <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=3)

mse_arima3 <- rbind(mse_u_arima3, mse_i_arima3[-1,], mse_imint_arima3[-1,])
rownames(mse_arima3) <- row_names

basef_mse_arima3 <- mse_arima3[1, ]
pc3_arima <- mse_arima3[-(1:5), ]
pc3_arima <- sweep(pc3_arima, 2, basef_mse_arima3, FUN = function(x, y) ((x - y) / y) * 100)
pc3_arima <- rbind(basef_mse_arima3/(10^3), pc3_arima)
rownames(pc3_arima) <- rownames(mse_arima3[-(2:5),])
round(pc3_arima,3)

as.data.frame(pc3_arima) %>% 
  write.csv("pc/pc3_arima.csv", row.names = TRUE)

### h=1-3
mse_u_arima1to3 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=1:3)
mse_i_arima1to3 <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=1:3)
mse_imint_arima1to3 <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=1:3)

mse_arima1to3 <- rbind(mse_u_arima1to3, mse_i_arima1to3[-1,], mse_imint_arima1to3[-1,])
rownames(mse_arima1to3) <- row_names

basef_mse_arima1to3 <- mse_arima1to3[1, ]
pc1to3_arima <- mse_arima1to3[-(1:5), ]
pc1to3_arima <- sweep(pc1to3_arima, 2, basef_mse_arima1to3, FUN = function(x, y) ((x - y) / y) * 100)
pc1to3_arima <- rbind(basef_mse_arima1to3/(10^3), pc1to3_arima)
rownames(pc1to3_arima) <- rownames(mse_arima1to3[-(2:5),])
round(pc1to3_arima,3)

as.data.frame(pc1to3_arima) %>% 
  write.csv("pc/pc1to3_arima.csv", row.names = TRUE)

### h=6
mse_u_arima6 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=6)
mse_i_arima6 <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=6)
mse_imint_arima6 <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=6)

mse_arima6 <- rbind(mse_u_arima6, mse_i_arima6[-1,], mse_imint_arima6[-1,])
rownames(mse_arima6) <- row_names

basef_mse_arima6 <- mse_arima6[1, ]
pc6_arima <- mse_arima6[-(1:5), ]
pc6_arima <- sweep(pc6_arima, 2, basef_mse_arima6, FUN = function(x, y) ((x - y) / y) * 100)
pc6_arima <- rbind(basef_mse_arima6/(10^3), pc6_arima)
rownames(pc6_arima) <- rownames(mse_arima6[-(2:5),])
round(pc6_arima,3)

as.data.frame(pc6_arima) %>% 
  write.csv("pc/pc6_arima.csv", row.names = TRUE)


### h=1-6
mse_u_arima1to6 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=1:6)
mse_i_arima1to6 <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=1:6)
mse_imint_arima1to6 <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=1:6)

mse_arima1to6 <- rbind(mse_u_arima1to6, mse_i_arima1to6[-1,], mse_imint_arima1to6[-1,])
rownames(mse_arima1to6) <- row_names

basef_mse_arima1to6 <- mse_arima1to6[1, ]
pc1to6_arima <- mse_arima1to6[-(1:5), ]
pc1to6_arima <- sweep(pc1to6_arima, 2, basef_mse_arima1to6, FUN = function(x, y) ((x - y) / y) * 100)
pc1to6_arima <- rbind(basef_mse_arima1to6/(10^3), pc1to6_arima)
rownames(pc1to6_arima) <- rownames(mse_arima1to6[-(2:5),])
round(pc1to6_arima,3)

as.data.frame(pc1to6_arima) %>% 
  write.csv("pc/pc1to6_arima.csv", row.names = TRUE)

### h=12
mse_u_arima12 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=12)
mse_i_arima12 <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=12)
mse_imint_arima12 <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=12)

mse_arima12 <- rbind(mse_u_arima12, mse_i_arima12[-1,], mse_imint_arima12[-1,])
rownames(mse_arima12) <- row_names

basef_mse_arima12 <- mse_arima12[1, ]
pc12_arima <- mse_arima12[-(1:5), ]
pc12_arima <- sweep(pc12_arima, 2, basef_mse_arima12, FUN = function(x, y) ((x - y) / y) * 100)
pc12_arima <- rbind(basef_mse_arima12/(10^3), pc12_arima)
rownames(pc12_arima) <- rownames(mse_arima12[-(2:5),])
round(pc12_arima,3)

as.data.frame(pc12_arima) %>% 
  write.csv("pc/pc12_arima.csv", row.names = TRUE)


### h=1-12
mse_u_arima1to12 <- accuracy_check(forecasts=arimaf_U, test_sets=arima_test_set, h=1:12)
mse_i_arima1to12 <- accuracy_check(forecasts=arimaf_i, test_sets=arima_test_set, h=1:12)
mse_imint_arima1to12 <- accuracy_check(forecasts=arimaf_imint, test_sets=arima_test_set, h=1:12)

mse_arima1to12 <- rbind(mse_u_arima1to12, mse_i_arima1to12[-1,], mse_imint_arima1to12[-1,])
rownames(mse_arima1to12) <- row_names

basef_mse_arima1to12 <- mse_arima1to12[1, ]
pc1to12_arima <- mse_arima1to12[-(1:5), ]
pc1to12_arima <- sweep(pc1to12_arima, 2, basef_mse_arima1to12, FUN = function(x, y) ((x - y) / y) * 100)
pc1to12_arima <- rbind(basef_mse_arima1to12/(10^3), pc1to12_arima)
rownames(pc1to12_arima) <- rownames(mse_arima1to12[-(2:5),])
round(pc1to12_arima,3)

as.data.frame(pc1to12_arima) %>% 
  write.csv("pc/pc1to12_arima.csv", row.names = TRUE)


