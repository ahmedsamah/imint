#source('chf.r')

library(tidyverse)
library(parallel)
library(ihts)
library(hts)

################## Functions #####################################

hts <- function(basis, sMat) {
  t(sMat %*% t(basis))
}

convert_to_list <- function(dataS1, index_col = "index") {
  data_list <- split(dataS1, dataS1[[index_col]])
  data_list <- lapply(data_list, function(df) df[ , !names(df) %in% index_col])
  return(data_list)
}

### Function for unconstrained reconciliation ###
reconcile_f <- function(basef_list, nodes, groups, data_list, smat) {
  
  basef_l <- list()
  ols_results <- list()
  wlss_results <- list()
  wlsv_results <- list()
  sam_results <- list()
  shrink_results <- list()
  n <- dim(smat)[1]
  
  # Iterate over each set of forecasts and residuals
  for (i in seq_along(basef_list)) {
    basef <- as.matrix(basef_list[[i]][301:324,])
    basef_l[[i]] <- basef
    fitted <- basef_list[[i]][1:300,]
    actual_basis <- data_list[[i]][1:300,]
    actual_full <- hts(actual_basis, smat)
    resid <- as.matrix(actual_full - fitted)
    cov_mat <- crossprod(resid)/nrow(resid)
    
    # Perform the reconciliation for each method
    ols_results[[i]] <- combinef(fcasts = basef, nodes = nodes, groups = groups,
                                 weights = NULL, nonnegative = FALSE,
                                 algorithms = "lu", keep = "all",
                                 parallel = TRUE)
    
    wlss_weight <- as.vector(smat %*% matrix(1, nrow = dim(smat)[2]))
    wlss_results[[i]] <- combinef(fcasts = basef, nodes = nodes, groups = groups,
                                  weights = 1 / wlss_weight, nonnegative = FALSE,
                                  algorithms = "lu", keep = "all",
                                  parallel = TRUE)
    
    wlsv_weight <- diag(crossprod(resid) / nrow(resid))
    wlsv_results[[i]] <- combinef(fcasts = basef, nodes = nodes, groups = groups,
                                  weights = 1 / wlsv_weight, nonnegative = FALSE,
                                  algorithms = "lu", keep = "all", 
                                  parallel = TRUE)
    
    sam_results[[i]] <- MinT(fcasts = basef, nodes = nodes, groups = groups,
                                nonnegative = FALSE, residual = resid,
                                covariance = "sam", algorithms = "lu",
                                keep = "all", parallel = TRUE)
    
    shrink_results[[i]] <- MinT(fcasts = basef, nodes = nodes, groups = groups,
                                nonnegative = FALSE, residual = resid,
                                covariance = "shr", algorithms = "lu",
                                keep = "all", parallel = TRUE)
  }
  reconciled_results <- list(
    basef = basef_l, 
    ols_U = ols_results,
    wlss_U = wlss_results,
    wlsv_U = wlsv_results,
    sam_U = sam_results,
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


reconcile_f_qpi <- function(basef_list, immu_set, smat, data_list) {
  
  # Function to process a single basef and resid pair
  process_pair <- function(basef, resid) {
    ols_weight <- diag(nrow(smat))
    ols_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = ols_weight)
    
    wlss_weight <- diag(as.vector(smat %*% rep(1, ncol(smat))))
    wlss_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = wlss_weight)
    
    wlsv_weight <- diag(diag(crossprod(resid) / nrow(resid)))
    wlsv_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = wlsv_weight)
    
    sam_weight <- crossprod(resid) / nrow(resid)
    sam_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = sam_weight)
    
    tar <- hts:::lowerD(resid)
    shrink <- hts:::shrink.estim(resid, tar)
    shrink_weight <- shrink[[1]]
    shrink_qpi <- qp_i(basef = basef, immu_set = immu_set, smat = smat, weights = shrink_weight)
    
    list(
      basef = basef,
      ols_qpi = ols_qpi,
      wlss_qpi = wlss_qpi,
      wlsv_qpi = wlsv_qpi,
      sam_qpi = sam_qpi,
      shrink_qpi = shrink_qpi
    )
  }
  
  results_list <- vector("list", length(basef_list))
  
  # Process each pair of basef and resid
  for (i in seq_along(basef_list)) {
    basef <- as.matrix(basef_list[[i]][301:324,])
    fitted <- basef_list[[i]][1:300,]
    actual_basis <- data_list[[i]][1:300,]
    actual_full <- hts(actual_basis, smat)
    resid <- as.matrix(actual_full - fitted)
    results_list[[i]] <- process_pair(basef, resid)
  }
  
  basef <- lapply(results_list, `[[`, "basef")
  ols_results <- lapply(results_list, `[[`, "ols_qpi")
  wlss_results <- lapply(results_list, `[[`, "wlss_qpi")
  wlsv_results <- lapply(results_list, `[[`, "wlsv_qpi")
  sam_results <- lapply(results_list, `[[`, "sam_qpi")
  shrink_results <- lapply(results_list, `[[`, "shrink_qpi")
  
  
  list(
    basef = basef,
    ols_I = ols_results,
    wlss_I = wlss_results,
    wlsv_I = wlsv_results,
    sam_I = sam_results,
    shrink_I = shrink_results
  )
}



# reconcile_f_i <- function(basef_list, data_list, immu_set=NULL, smat) {
#   
#   basef_l <- list()
#   ols_i_results <- list()
#   wlss_i_results <- list()
#   wlsv_i_results <- list()
#   sam_i_results <- list()
#   shrink_i_results <- list()
#   n <- dim(smat)[1]
#   
#   # Iterate over each set of forecasts and residuals
#   for (i in seq_along(basef_list)) {
#     basef <- as.matrix(basef_list[[i]][301:324,])
#     basef_l[[i]] <- basef
#     fitted <- basef_list[[i]][1:300,]
#     actual_basis <- data_list[[i]][1:300,]
#     actual_full <- hts(actual_basis, smat)
#     resid <- as.matrix(actual_full - fitted)
#     cov_mat <- crossprod(resid)/nrow(resid)
#     
#     # Perform the reconciliation for each method
#     ols_weight <- diag(n)
#     ols_weight <- methods::as(ols_weight, "sparseMatrix")
#     ols_i_results[[i]] <- forecast.reconcile(base_forecasts = basef, sMat = smat, 
#                                              weighting_matrix = ols_weight, immu_set=immu_set,
#                                              nonnegative=FALSE)
#     
#     wlss_weight <- diag(as.vector(smat %*% matrix(1, nrow = dim(smat)[2])))
#     wlss_weight <- methods::as(wlss_weight, "sparseMatrix")
#     wlss_i_results[[i]] <- forecast.reconcile(base_forecasts = basef, sMat = smat, 
#                                               weighting_matrix = wlss_weight, immu_set=immu_set,
#                                               nonnegative=FALSE)
#     
#     wlsv_weight <- diag(diag(cov_mat))
#     wlsv_weight <- methods::as(wlsv_weight, "sparseMatrix")
#     wlsv_i_results[[i]] <- forecast.reconcile(base_forecasts = basef, sMat = smat, 
#                                               weighting_matrix = wlsv_weight, immu_set=immu_set,
#                                               nonnegative=FALSE)
#     
#     sam_i_results[[i]] <- forecast.reconcile(base_forecasts = basef, sMat = smat, 
#                                               weighting_matrix = cov_mat, immu_set=immu_set,
#                                               nonnegative=FALSE)
#     
#     tar <- hts:::lowerD(resid)
#     shrink <- hts:::shrink.estim(resid, tar)
#     shrink_weight <- shrink[[1]]
#     lambda <- shrink[[2]]
#     shrink_weight <- methods::as(shrink_weight, "sparseMatrix")
#     shrink_i_results[[i]] <- forecast.reconcile(base_forecasts = basef, sMat = smat, 
#                                                 weighting_matrix = shrink_weight, immu_set=immu_set,
#                                                 nonnegative=FALSE)
#   }
#   reconciled_results <- list(
#     basef = basef_l, 
#     ols_I = ols_i_results,
#     wlss_I = wlss_i_results,
#     wlsv_I = wlsv_i_results,
#     sam_I = sam_i_results,
#     shrink_I = shrink_i_results
#   )
#   return(reconciled_results)
# }
# 





reconcile_f_imint <- function(basef_list, nodes, groups, data_list, immu_set, smat) {
  basef_l <- list()
  ols_iMinT_results <- list()
  wlss_iMinT_results <- list()
  wlsv_iMinT_results <- list()
  sam_iMinT_results <- list()
  shrink_iMinT_results <- list()

  # Iterate over each set of forecasts and residuals
  for (i in seq_along(basef_list)) {
    basef <- as.matrix(basef_list[[i]][301:324,])
    basef_l[[i]] <- basef
    fitted <- basef_list[[i]][1:300,]
    actual_basis <- data_list[[i]][1:300,]
    actual_full <- hts(actual_basis, smat)
    resid <- as.matrix(actual_full - fitted)

    ols_iMinT_results[[i]] <- icombinef(fcasts = basef, nodes = nodes, groups = NULL,
                                        weights = NULL, immute = immu_set,
                                        nonnegative = FALSE,
                                        algorithms = "lu",
                                        keep = "all")

    wlss_weight <- as.vector(smat %*% matrix(1, nrow = dim(smat)[2]))
    wlss_iMinT_results[[i]] <- icombinef(fcasts = basef, nodes = nodes, groups = NULL,
                                         weights = 1 / wlss_weight, immute = immu_set,
                                         nonnegative = FALSE,
                                         algorithms = "lu",
                                         keep = "all")

    wlsv_weight <- diag(crossprod(resid) / nrow(resid))
    wlsv_iMinT_results[[i]] <- icombinef(fcasts = basef, nodes = nodes, groups = NULL,
                                         weights = 1 / wlsv_weight, immute = immu_set,
                                         nonnegative = FALSE,
                                         algorithms = "lu",
                                         keep = "all")

    sam_iMinT_results[[i]] <- iMinT(fcasts = basef, nodes = nodes, groups = NULL,
                                       immute = immu_set, nonnegative = FALSE,
                                       residual = resid,
                                       covariance = "sam",
                                       algorithms = "lu",
                                       keep = "all")

    shrink_iMinT_results[[i]] <- iMinT(fcasts = basef, nodes = nodes, groups = NULL,
                                       immute = immu_set, nonnegative = FALSE,
                                       residual = resid,
                                       covariance = "shr",
                                       algorithms = "lu",
                                       keep = "all")
  }
  reconciled_results <- list(
    basef = basef_l,
    ols_iMinT = ols_iMinT_results,
    wlss_iMinT = wlss_iMinT_results,
    wlsv_iMinT = wlsv_iMinT_results,
    sam_iMinT_results = sam_iMinT_results,
    shrink_iMinT = shrink_iMinT_results
  )
  return(reconciled_results)
}



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
                        Middle = rowMeans(mse[,2:3], na.rm = TRUE),
                        Bottom = rowMeans(mse[,4:7], na.rm = TRUE),
                        All_series = rowMeans(mse[,1:7], na.rm = TRUE)
    )
    
    # average mse for the h-step-ahead
    mse_h <- colMeans(mse_levels[h, , drop = FALSE], na.rm = TRUE)
    
    mse_results[[type]] <- mse_h
  }
  result_df <- do.call(rbind, mse_results)
  rownames(result_df) <- names(mse_results)
  return(result_df)
}

########################## SCENARIO 1 ###########################

######################## Data import ##############################

smat <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
             diag(rep(1, 4)))
nodes <- list(2, rep(2, 2))

data1 <- read.csv('data/data_scenario1.csv', row.names = 1) %>% 
          select(!t)
dataS1 <- convert_to_list(data1, "index")

basef1 <- read.csv('basef/basef1.csv', row.names = 1) %>% 
  select(!t)

basef1_ets <- basef1 %>% filter(method=="ets")
basef1_arima <- basef1 %>% filter(method=="arima")
basef1_ea <- cbind(basef1_ets[,1], basef1_arima[,-c(1,8)])

basefS1_ets <- convert_to_list(basef1_ets[,-8], "index")
basefS1_ea <- convert_to_list(basef1_ea, "index")



################### Reconciliation ###############################

##### ETS+ETS

# u_ets <- reconcile_f_i(basef_list=basefS1_ets, data_list=dataS1, 
#                        immu_set = NULL, smat=smat)
# i_ets <- reconcile_f_i(basef_list=basefS1_ets, data_list=dataS1,
#                        immu_set=1, smat=smat) 
# imint_ets <- reconcile_f_imint(basef_list=basefS1_ets, nodes=nodes,
#                                data_list=dataS1, immu_set=1, smat=smat)

u_ets <- reconcile_f(basef_list=basefS1_ets, nodes=nodes, groups=NULL, 
                     data_list=dataS1, smat=smat)
i_ets <- reconcile_f_qpi(basef_list=basefS1_ets, immu_set=1, smat=smat, 
                         data_list=dataS1)
imint_ets <- reconcile_f_imint(basef_list=basefS1_ets, nodes=nodes,
                               data_list=dataS1, immu_set=1, smat=smat)

##### ETS+ARIMA

# u_ea <- reconcile_f_i(basef_list=basefS1_ea, data_list=dataS1, 
#                        immu_set = NULL, smat=smat) 
# i_ea <- reconcile_f_i(basef_list=basefS1_ea, data_list=dataS1,
#                        immu_set=1, smat=smat) 

u_ea <- reconcile_f(basef_list=basefS1_ea, nodes=nodes, groups=NULL, 
                     data_list=dataS1, smat=smat)
i_ea <- reconcile_f_qpi(basef_list=basefS1_ea, immu_set=1, smat=smat, 
                         data_list=dataS1)
imint_ea <- reconcile_f_imint(basef_list=basefS1_ea, nodes=nodes,
                              data_list=dataS1, immu_set=1, smat=smat)

###################### Accuracy ################################
test_set_d1 <- list()
for (i in seq_along(dataS1)) {
  test_basis <- dataS1[[i]][301:324,]
  test_full <- hts(test_basis, smat)
  test_set_d1[[i]] <- test_full
}

test_set_s1 <- list(basef = test_set_d1, ols = test_set_d1,
                 wlss = test_set_d1, wlsv = test_set_d1,
                 sam = test_set_d1, shrink = test_set_d1)

##### ETS + ETS #####

etsf_U <- list(basef = u_ets[["basef"]], ols = u_ets[["ols_U"]],
               wlss = u_ets[["wlss_U"]], wlsv = u_ets[["wlsv_U"]],
               sam = u_ets[["sam_U"]], shrink = u_ets[["shrink_U"]])

etsf_i <- list(basef = i_ets[["basef"]], ols = i_ets[["ols_I"]],
               wlss = i_ets[["wlss_I"]], wlsv = i_ets[["wlsv_I"]],
               sam= i_ets[["sam_I"]], shrink = i_ets[["shrink_I"]])

etsf_imint <- list(basef = imint_ets[["basef"]], ols= imint_ets[["ols_iMinT"]],
                   wlss = imint_ets[["wlss_iMinT"]], wlsv = imint_ets[["wlsv_iMinT"]],
                   sam = imint_ets[["sam_iMinT_results"]], shrink = imint_ets[["shrink_iMinT"]])

row_names <- c("basef", "ols_U", "wlss_U", "wlsv_U", "sam_U", "shrink_U", 
               "ols_I", "wlss_I", "wlsv_I", "sam_I", "shrink_I", 
               "ols_iMinT", "wlss_iMinT", "wlsv_iMinT", "sam_iMinT","shrink_iMinT")

### h=1
mse_u_ets1 <- accuracy_check(forecasts=etsf_U, test_sets=test_set_s1, h=1)
mse_i_ets1 <- accuracy_check(forecasts=etsf_i, test_sets=test_set_s1, h=1)
mse_imint_ets1 <- accuracy_check(forecasts=etsf_imint, test_sets=test_set_s1, h=1)

mse_ets1 <- rbind(mse_u_ets1, mse_i_ets1[-1,], mse_imint_ets1[-1,])
rownames(mse_ets1) <- row_names

basef_mse_ets1 <- mse_ets1[1, ]
pc1 <- mse_ets1[-(1:6), ]
pc1 <- sweep(pc1, 2, basef_mse_ets1, FUN = function(x, y) ((x - y) / y) * 100)
pc1 <- rbind(basef_mse_ets1, pc1)
rownames(pc1) <- rownames(mse_ets1[-(2:6),])
round(pc1,3)

as.data.frame(pc1) %>% 
  write.csv("pc/pc1.csv", row.names = TRUE)

### h=3
mse_u_ets3 <- accuracy_check(forecasts=etsf_U, test_sets=test_set_s1, h=3)
mse_i_ets3 <- accuracy_check(forecasts=etsf_i, test_sets=test_set_s1, h=3)
mse_imint_ets3 <- accuracy_check(forecasts=etsf_imint, test_sets=test_set_s1, h=3)

mse_ets3 <- rbind(mse_u_ets3, mse_i_ets3[-1,], mse_imint_ets3[-1,])
rownames(mse_ets3) <- row_names

basef_mse_ets3 <- mse_ets3[1, ]
pc3 <- mse_ets3[-(1:6), ]
pc3 <- sweep(pc3, 2, basef_mse_ets3, FUN = function(x, y) ((x - y) / y) * 100)
pc3 <- rbind(basef_mse_ets3, pc3)
rownames(pc3) <- rownames(mse_ets3[-(2:6),])
round(pc3,3)

as.data.frame(pc3) %>% 
  write.csv("pc/pc3.csv", row.names = TRUE)

### h=6
mse_u_ets6 <- accuracy_check(forecasts=etsf_U, test_sets=test_set_s1, h=6)
mse_i_ets6 <- accuracy_check(forecasts=etsf_i, test_sets=test_set_s1, h=6)
mse_imint_ets6 <- accuracy_check(forecasts=etsf_imint, test_sets=test_set_s1, h=6)

mse_ets6 <- rbind(mse_u_ets6, mse_i_ets6[-1,], mse_imint_ets6[-1,])
rownames(mse_ets6) <- row_names

basef_mse_ets6 <- mse_ets6[1, ]
pc6 <- mse_ets6[-(1:6), ]
pc6 <- sweep(pc6, 2, basef_mse_ets6, FUN = function(x, y) ((x - y) / y) * 100)
pc6 <- rbind(basef_mse_ets6, pc6)
rownames(pc6) <- rownames(mse_ets6[-(2:6),])
round(pc6,3)

as.data.frame(pc6) %>% 
  write.csv("pc/pc6.csv", row.names = TRUE)

### h=12
mse_u_ets12 <- accuracy_check(forecasts=etsf_U, test_sets=test_set_s1, h=12)
mse_i_ets12 <- accuracy_check(forecasts=etsf_i, test_sets=test_set_s1, h=12)
mse_imint_ets12 <- accuracy_check(forecasts=etsf_imint, test_sets=test_set_s1, h=12)

mse_ets12 <- rbind(mse_u_ets12, mse_i_ets12[-1,], mse_imint_ets12[-1,])
rownames(mse_ets12) <- row_names

basef_mse_ets12 <- mse_ets12[1, ]
pc12 <- mse_ets12[-(1:6), ]
pc12 <- sweep(pc12, 2, basef_mse_ets12, FUN = function(x, y) ((x - y) / y) * 100)
pc12 <- rbind(basef_mse_ets12, pc12)
rownames(pc12) <- rownames(mse_ets12[-(2:6),])
round(pc12,3)

as.data.frame(pc12) %>% 
  write.csv("pc/pc12.csv", row.names = TRUE)

### h=24
mse_u_ets24 <- accuracy_check(forecasts=etsf_U, test_sets=test_set_s1, h=24)
mse_i_ets24 <- accuracy_check(forecasts=etsf_i, test_sets=test_set_s1, h=24)
mse_imint_ets24 <- accuracy_check(forecasts=etsf_imint, test_sets=test_set_s1, h=24)

mse_ets24 <- rbind(mse_u_ets24, mse_i_ets24[-1,], mse_imint_ets24[-1,])
rownames(mse_ets24) <- row_names

basef_mse_ets24 <- mse_ets24[1, ]
pc24 <- mse_ets24[-(1:6), ]
pc24 <- sweep(pc24, 2, basef_mse_ets24, FUN = function(x, y) ((x - y) / y) * 100)
pc24 <- rbind(basef_mse_ets24, pc24)
rownames(pc24) <- rownames(mse_ets24[-(2:6),])
round(pc24,3)

as.data.frame(pc24) %>% 
  write.csv("pc/pc24.csv", row.names = TRUE)


##### ETS + ARIMA #####

eaf_U <- list(basef = u_ea[["basef"]], ols = u_ea[["ols_U"]],
               wlss = u_ea[["wlss_U"]], wlsv = u_ea[["wlsv_U"]],
               sam = u_ea[["sam_U"]], shrink = u_ea[["shrink_U"]])

eaf_i <- list(basef = i_ea[["basef"]], ols = i_ea[["ols_I"]],
               wlss = i_ea[["wlss_I"]], wlsv = i_ea[["wlsv_I"]],
               sam= i_ea[["sam_I"]], shrink = i_ea[["shrink_I"]])

eaf_imint <- list(basef = imint_ea[["basef"]], ols= imint_ea[["ols_iMinT"]],
                   wlss = imint_ea[["wlss_iMinT"]], wlsv = imint_ea[["wlsv_iMinT"]],
                   sam = imint_ea[["sam_iMinT_results"]], shrink = imint_ea[["shrink_iMinT"]])

row_names <- c("basef", "ols_U", "wlss_U", "wlsv_U", "sam_U", "shrink_U", 
               "ols_I", "wlss_I", "wlsv_I", "sam_I", "shrink_I", 
               "ols_iMinT", "wlss_iMinT", "wlsv_iMinT", "sam_iMinT","shrink_iMinT")

### h=1
mse_u_ea1 <- accuracy_check(forecasts=eaf_U, test_sets=test_set_s1, h=1)
mse_i_ea1 <- accuracy_check(forecasts=eaf_i, test_sets=test_set_s1, h=1)
mse_imint_ea1 <- accuracy_check(forecasts=eaf_imint, test_sets=test_set_s1, h=1)

mse_ea1 <- rbind(mse_u_ea1, mse_i_ea1[-1,], mse_imint_ea1[-1,])
rownames(mse_ea1) <- row_names

basef_mse_ea1 <- mse_ea1[1, ]
pc1_ea <- mse_ea1[-(1:6), ]
pc1_ea <- sweep(pc1_ea, 2, basef_mse_ea1, FUN = function(x, y) ((x - y) / y) * 100)
pc1_ea <- rbind(basef_mse_ea1, pc1_ea)
rownames(pc1_ea) <- rownames(mse_ea1[-(2:6),])
round(pc1_ea,3)

as.data.frame(pc1_ea) %>% 
  write.csv("pc/pc1_ea.csv", row.names = TRUE)

### h=3
mse_u_ea3 <- accuracy_check(forecasts=eaf_U, test_sets=test_set_s1, h=3)
mse_i_ea3 <- accuracy_check(forecasts=eaf_i, test_sets=test_set_s1, h=3)
mse_imint_ea3 <- accuracy_check(forecasts=eaf_imint, test_sets=test_set_s1, h=3)

mse_ea3 <- rbind(mse_u_ea3, mse_i_ea3[-1,], mse_imint_ea3[-1,])
rownames(mse_ea3) <- row_names

basef_mse_ea3 <- mse_ea3[1, ]
pc3_ea <- mse_ea3[-(1:6), ]
pc3_ea <- sweep(pc3_ea, 2, basef_mse_ea3, FUN = function(x, y) ((x - y) / y) * 100)
pc3_ea <- rbind(basef_mse_ea3, pc3_ea)
rownames(pc3_ea) <- rownames(mse_ea3[-(2:6),])
round(pc3_ea,3)

as.data.frame(pc3_ea) %>% 
  write.csv("pc/pc3_ea.csv", row.names = TRUE)

### h=6
mse_u_ea6 <- accuracy_check(forecasts=eaf_U, test_sets=test_set_s1, h=6)
mse_i_ea6 <- accuracy_check(forecasts=eaf_i, test_sets=test_set_s1, h=6)
mse_imint_ea6 <- accuracy_check(forecasts=eaf_imint, test_sets=test_set_s1, h=6)

mse_ea6 <- rbind(mse_u_ea6, mse_i_ea6[-1,], mse_imint_ea6[-1,])
rownames(mse_ea6) <- row_names

basef_mse_ea6 <- mse_ea6[1, ]
pc6_ea <- mse_ea6[-(1:6), ]
pc6_ea <- sweep(pc6_ea, 2, basef_mse_ea6, FUN = function(x, y) ((x - y) / y) * 100)
pc6_ea <- rbind(basef_mse_ea6, pc6_ea)
rownames(pc6_ea) <- rownames(mse_ea6[-(2:6),])
round(pc6_ea,3)

as.data.frame(pc6_ea) %>% 
  write.csv("pc/pc6_ea.csv", row.names = TRUE)


### h=12
mse_u_ea12 <- accuracy_check(forecasts=eaf_U, test_sets=test_set_s1, h=12)
mse_i_ea12 <- accuracy_check(forecasts=eaf_i, test_sets=test_set_s1, h=12)
mse_imint_ea12 <- accuracy_check(forecasts=eaf_imint, test_sets=test_set_s1, h=12)

mse_ea12 <- rbind(mse_u_ea12, mse_i_ea12[-1,], mse_imint_ea12[-1,])
rownames(mse_ea12) <- row_names

basef_mse_ea12 <- mse_ea12[1, ]
pc12_ea <- mse_ea12[-(1:6), ]
pc12_ea <- sweep(pc12_ea, 2, basef_mse_ea12, FUN = function(x, y) ((x - y) / y) * 100)
pc12_ea <- rbind(basef_mse_ea12, pc12_ea)
rownames(pc12_ea) <- rownames(mse_ea12[-(2:6),])
round(pc12_ea,3)

as.data.frame(pc12_ea) %>% 
  write.csv("pc/pc12_ea.csv", row.names = TRUE)

### h=24
mse_u_ea24 <- accuracy_check(forecasts=eaf_U, test_sets=test_set_s1, h=24)
mse_i_ea24 <- accuracy_check(forecasts=eaf_i, test_sets=test_set_s1, h=24)
mse_imint_ea24 <- accuracy_check(forecasts=eaf_imint, test_sets=test_set_s1, h=24)

mse_ea24 <- rbind(mse_u_ea24, mse_i_ea24[-1,], mse_imint_ea24[-1,])
rownames(mse_ea24) <- row_names

basef_mse_ea24 <- mse_ea24[1, ]
pc24_ea <- mse_ea24[-(1:6), ]
pc24_ea <- sweep(pc24_ea, 2, basef_mse_ea24, FUN = function(x, y) ((x - y) / y) * 100)
pc24_ea <- rbind(basef_mse_ea24, pc24_ea)
rownames(pc24_ea) <- rownames(mse_ea24[-(2:6),])
round(pc24_ea,3)

as.data.frame(pc24_ea) %>% 
  write.csv("pc/pc24_ea.csv", row.names = TRUE)




########################## SCENARIO 2 ###########################

######################## Data import ##############################

smat <- rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
              diag(rep(1, 4)))
nodes <- list(2, rep(2, 2))

data2 <- read.csv('data/data_scenario2.csv', row.names = 1) %>% 
  select(!t)
dataS2 <- convert_to_list(data2, "index")

basef2 <- read.csv('basef/basef2.csv', row.names = 1) %>% 
  select(!t)

basef2_ets <- basef2 %>% filter(method=="ets")
basef2_arima <- basef2 %>% filter(method=="arima")
basef2_ea <- cbind(basef2_ets[,1], basef2_arima[,-c(1,8)])

basefS2_ets <- convert_to_list(basef2_ets[,-8], "index")
basefS2_ea <- convert_to_list(basef2_ea, "index")



################### Reconciliation ###############################

##### ETS+ETS

# u_ets_s2 <- reconcile_f_i(basef_list=basefS2_ets, data_list=dataS2, 
#                        immu_set = NULL, smat=smat) 
# i_ets_s2 <- reconcile_f_i(basef_list=basefS2_ets, data_list=dataS2,
#                        immu_set=1, smat=smat) 

u_ets_s2 <- reconcile_f(basef_list=basefS2_ets, nodes=nodes, groups=NULL, 
                     data_list=dataS2, smat=smat)
i_ets_s2 <- reconcile_f_qpi(basef_list=basefS2_ets, immu_set=1, smat=smat, 
                         data_list=dataS2)
imint_ets_s2 <- reconcile_f_imint(basef_list=basefS2_ets, nodes=nodes,
                               data_list=dataS2, immu_set=1, smat=smat)

##### ETS+ARIMA

# u_ea <- reconcile_f_i(basef_list=basefS1_ea, data_list=dataS1, 
#                        immu_set = NULL, smat=smat) 
# i_ea <- reconcile_f_i(basef_list=basefS1_ea, data_list=dataS1,
#                        immu_set=1, smat=smat) 

u_ea_s2 <- reconcile_f(basef_list=basefS2_ea, nodes=nodes, groups=NULL, 
                    data_list=dataS2, smat=smat)
i_ea_s2 <- reconcile_f_qpi(basef_list=basefS2_ea, immu_set=1, smat=smat, 
                        data_list=dataS2)
imint_ea_s2 <- reconcile_f_imint(basef_list=basefS2_ea, nodes=nodes,
                              data_list=dataS2, immu_set=1, smat=smat)

###################### Accuracy ################################
test_set_d2 <- list()
for (i in seq_along(dataS2)) {
  test_basis <- dataS2[[i]][301:324,]
  test_full <- hts(test_basis, smat)
  test_set_d2[[i]] <- test_full
}

test_set_s2 <- list(basef = test_set_d2, ols = test_set_d2,
                    wlss = test_set_d2, wlsv = test_set_d2,
                    sam = test_set_d2, shrink = test_set_d2)

##### ETS + ETS #####

etsf_U_s2 <- list(basef = u_ets_s2[["basef"]], ols = u_ets_s2[["ols_U"]],
               wlss = u_ets_s2[["wlss_U"]], wlsv = u_ets_s2[["wlsv_U"]],
               sam = u_ets_s2[["sam_U"]], shrink = u_ets_s2[["shrink_U"]])

etsf_i_s2 <- list(basef = i_ets_s2[["basef"]], ols = i_ets_s2[["ols_I"]],
               wlss = i_ets_s2[["wlss_I"]], wlsv = i_ets_s2[["wlsv_I"]],
               sam= i_ets_s2[["sam_I"]], shrink = i_ets_s2[["shrink_I"]])

etsf_imint_s2 <- list(basef = imint_ets_s2[["basef"]], ols= imint_ets_s2[["ols_iMinT"]],
                   wlss = imint_ets_s2[["wlss_iMinT"]], wlsv = imint_ets_s2[["wlsv_iMinT"]],
                   sam = imint_ets_s2[["sam_iMinT_results"]], shrink = imint_ets_s2[["shrink_iMinT"]])

row_names <- c("basef", "ols_U", "wlss_U", "wlsv_U", "sam_U", "shrink_U", 
               "ols_I", "wlss_I", "wlsv_I", "sam_I", "shrink_I", 
               "ols_iMinT", "wlss_iMinT", "wlsv_iMinT", "sam_iMinT","shrink_iMinT")

### h=1
mse_u_ets1_s2 <- accuracy_check(forecasts=etsf_U_s2, test_sets=test_set_s2, h=1)
mse_i_ets1_s2 <- accuracy_check(forecasts=etsf_i_s2, test_sets=test_set_s2, h=1)
mse_imint_ets1_s2 <- accuracy_check(forecasts=etsf_imint_s2, test_sets=test_set_s2, h=1)

mse_ets1_s2 <- rbind(mse_u_ets1_s2, mse_i_ets1_s2[-1,], mse_imint_ets1_s2[-1,])
rownames(mse_ets1_s2) <- row_names

basef_mse_ets1_s2 <- mse_ets1_s2[1, ]
pc1_s2 <- mse_ets1_s2[-(1:6), ]
pc1_s2 <- sweep(pc1_s2, 2, basef_mse_ets1_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc1_s2 <- rbind(basef_mse_ets1_s2, pc1_s2)
rownames(pc1_s2) <- rownames(mse_ets1_s2[-(2:6),])
round(pc1_s2,3)

as.data.frame(pc1_s2) %>% 
  write.csv("pc/pc1_s2.csv", row.names = TRUE)

### h=3
mse_u_ets3_s2 <- accuracy_check(forecasts=etsf_U_s2, test_sets=test_set_s2, h=3)
mse_i_ets3_s2 <- accuracy_check(forecasts=etsf_i_s2, test_sets=test_set_s2, h=3)
mse_imint_ets3_s2 <- accuracy_check(forecasts=etsf_imint_s2, test_sets=test_set_s2, h=3)

mse_ets3_s2 <- rbind(mse_u_ets3_s2, mse_i_ets3_s2[-1,], mse_imint_ets3_s2[-1,])
rownames(mse_ets3_s2) <- row_names

basef_mse_ets3_s2 <- mse_ets3_s2[1, ]
pc3_s2 <- mse_ets3_s2[-(1:6), ]
pc3_s2 <- sweep(pc3_s2, 2, basef_mse_ets3_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc3_s2 <- rbind(basef_mse_ets3_s2, pc3_s2)
rownames(pc3_s2) <- rownames(mse_ets3_s2[-(2:6),])
round(pc3_s2,3)

as.data.frame(pc3_s2) %>% 
  write.csv("pc/pc3_s2.csv", row.names = TRUE)

### h=6
mse_u_ets6_s2 <- accuracy_check(forecasts=etsf_U_s2, test_sets=test_set_s2, h=6)
mse_i_ets6_s2 <- accuracy_check(forecasts=etsf_i_s2, test_sets=test_set_s2, h=6)
mse_imint_ets6_s2 <- accuracy_check(forecasts=etsf_imint_s2, test_sets=test_set_s2, h=6)

mse_ets6_s2 <- rbind(mse_u_ets6_s2, mse_i_ets6_s2[-1,], mse_imint_ets6_s2[-1,])
rownames(mse_ets6_s2) <- row_names

basef_mse_ets6_s2 <- mse_ets6_s2[1, ]
pc6_s2 <- mse_ets6_s2[-(1:6), ]
pc6_s2 <- sweep(pc6_s2, 2, basef_mse_ets6_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc6_s2 <- rbind(basef_mse_ets6_s2, pc6_s2)
rownames(pc6_s2) <- rownames(mse_ets6_s2[-(2:6),])
round(pc6_s2,3)

as.data.frame(pc6_s2) %>% 
  write.csv("pc/pc6_s2.csv", row.names = TRUE)

### h=12
mse_u_ets12_s2 <- accuracy_check(forecasts=etsf_U_s2, test_sets=test_set_s2, h=12)
mse_i_ets12_s2 <- accuracy_check(forecasts=etsf_i_s2, test_sets=test_set_s2, h=12)
mse_imint_ets12_s2 <- accuracy_check(forecasts=etsf_imint_s2, test_sets=test_set_s2, h=12)

mse_ets12_s2 <- rbind(mse_u_ets12_s2, mse_i_ets12_s2[-1,], mse_imint_ets12_s2[-1,])
rownames(mse_ets12_s2) <- row_names

basef_mse_ets12_s2 <- mse_ets12_s2[1, ]
pc12_s2 <- mse_ets12_s2[-(1:6), ]
pc12_s2 <- sweep(pc12_s2, 2, basef_mse_ets12_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc12_s2 <- rbind(basef_mse_ets12_s2, pc12_s2)
rownames(pc12_s2) <- rownames(mse_ets12_s2[-(2:6),])
round(pc12_s2,3)

as.data.frame(pc12_s2) %>% 
  write.csv("pc/pc12_s2.csv", row.names = TRUE)

### h=24
mse_u_ets24_s2 <- accuracy_check(forecasts=etsf_U_s2, test_sets=test_set_s2, h=24)
mse_i_ets24_s2 <- accuracy_check(forecasts=etsf_i_s2, test_sets=test_set_s2, h=24)
mse_imint_ets24_s2 <- accuracy_check(forecasts=etsf_imint_s2, test_sets=test_set_s2, h=24)

mse_ets24_s2 <- rbind(mse_u_ets24_s2, mse_i_ets24_s2[-1,], mse_imint_ets24_s2[-1,])
rownames(mse_ets24_s2) <- row_names

basef_mse_ets24_s2 <- mse_ets24_s2[1, ]
pc24_s2 <- mse_ets24_s2[-(1:6), ]
pc24_s2 <- sweep(pc24_s2, 2, basef_mse_ets24_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc24_s2 <- rbind(basef_mse_ets24_s2, pc24_s2)
rownames(pc24_s2) <- rownames(mse_ets24_s2[-(2:6),])
round(pc24_s2,3)

as.data.frame(pc24_s2) %>% 
  write.csv("pc/pc24_s2.csv", row.names = TRUE)


##### ETS + ARIMA #####

eaf_U_s2 <- list(basef = u_ea_s2[["basef"]], ols = u_ea_s2[["ols_U"]],
              wlss = u_ea_s2[["wlss_U"]], wlsv = u_ea_s2[["wlsv_U"]],
              sam = u_ea_s2[["sam_U"]], shrink = u_ea_s2[["shrink_U"]])

eaf_i_s2 <- list(basef = i_ea_s2[["basef"]], ols = i_ea_s2[["ols_I"]],
              wlss = i_ea_s2[["wlss_I"]], wlsv = i_ea_s2[["wlsv_I"]],
              sam= i_ea_s2[["sam_I"]], shrink = i_ea_s2[["shrink_I"]])

eaf_imint_s2 <- list(basef = imint_ea_s2[["basef"]], ols= imint_ea_s2[["ols_iMinT"]],
                  wlss = imint_ea_s2[["wlss_iMinT"]], wlsv = imint_ea_s2[["wlsv_iMinT"]],
                  sam = imint_ea_s2[["sam_iMinT_results"]], shrink = imint_ea_s2[["shrink_iMinT"]])

row_names <- c("basef", "ols_U", "wlss_U", "wlsv_U", "sam_U", "shrink_U", 
               "ols_I", "wlss_I", "wlsv_I", "sam_I", "shrink_I", 
               "ols_iMinT", "wlss_iMinT", "wlsv_iMinT", "sam_iMinT","shrink_iMinT")

### h=1
mse_u_ea1_s2 <- accuracy_check(forecasts=eaf_U_s2, test_sets=test_set_s2, h=1)
mse_i_ea1_s2 <- accuracy_check(forecasts=eaf_i_s2, test_sets=test_set_s2, h=1)
mse_imint_ea1_s2 <- accuracy_check(forecasts=eaf_imint_s2, test_sets=test_set_s2, h=1)

mse_ea1_s2 <- rbind(mse_u_ea1_s2, mse_i_ea1_s2[-1,], mse_imint_ea1_s2[-1,])
rownames(mse_ea1_s2) <- row_names

basef_mse_ea1_s2 <- mse_ea1_s2[1, ]
pc1_ea_s2 <- mse_ea1_s2[-(1:6), ]
pc1_ea_s2 <- sweep(pc1_ea_s2, 2, basef_mse_ea1_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc1_ea_s2 <- rbind(basef_mse_ea1_s2, pc1_ea_s2)
rownames(pc1_ea_s2) <- rownames(mse_ea1_s2[-(2:6),])
round(pc1_ea_s2,3)

as.data.frame(pc1_ea_s2) %>% 
  write.csv("pc/pc1_ea_s2.csv", row.names = TRUE)

### h=3
mse_u_ea3_s2 <- accuracy_check(forecasts=eaf_U_s2, test_sets=test_set_s2, h=3)
mse_i_ea3_s2 <- accuracy_check(forecasts=eaf_i_s2, test_sets=test_set_s2, h=3)
mse_imint_ea3_s2 <- accuracy_check(forecasts=eaf_imint_s2, test_sets=test_set_s2, h=3)

mse_ea3_s2 <- rbind(mse_u_ea3_s2, mse_i_ea3_s2[-1,], mse_imint_ea3_s2[-1,])
rownames(mse_ea3_s2) <- row_names

basef_mse_ea3_s2 <- mse_ea3_s2[1, ]
pc3_ea_s2 <- mse_ea3_s2[-(1:6), ]
pc3_ea_s2 <- sweep(pc3_ea_s2, 2, basef_mse_ea3_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc3_ea_s2 <- rbind(basef_mse_ea3_s2, pc3_ea_s2)
rownames(pc3_ea_s2) <- rownames(mse_ea3_s2[-(2:6),])
round(pc3_ea_s2,3)

as.data.frame(pc3_ea_s2) %>% 
  write.csv("pc/pc3_ea_s2.csv", row.names = TRUE)

### h=6
mse_u_ea6_s2 <- accuracy_check(forecasts=eaf_U_s2, test_sets=test_set_s2, h=6)
mse_i_ea6_s2 <- accuracy_check(forecasts=eaf_i_s2, test_sets=test_set_s2, h=6)
mse_imint_ea6_s2 <- accuracy_check(forecasts=eaf_imint_s2, test_sets=test_set_s2, h=6)

mse_ea6_s2 <- rbind(mse_u_ea6_s2, mse_i_ea6_s2[-1,], mse_imint_ea6_s2[-1,])
rownames(mse_ea6_s2) <- row_names

basef_mse_ea6_s2 <- mse_ea6_s2[1, ]
pc6_ea_s2 <- mse_ea6_s2[-(1:6), ]
pc6_ea_s2 <- sweep(pc6_ea_s2, 2, basef_mse_ea6_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc6_ea_s2 <- rbind(basef_mse_ea6_s2, pc6_ea_s2)
rownames(pc6_ea_s2) <- rownames(mse_ea6_s2[-(2:6),])
round(pc6_ea_s2,3)

as.data.frame(pc6_ea_s2) %>% 
  write.csv("pc/pc6_ea_s2.csv", row.names = TRUE)


### h=12
mse_u_ea12_s2 <- accuracy_check(forecasts=eaf_U_s2, test_sets=test_set_s2, h=12)
mse_i_ea12_s2 <- accuracy_check(forecasts=eaf_i_s2, test_sets=test_set_s2, h=12)
mse_imint_ea12_s2 <- accuracy_check(forecasts=eaf_imint_s2, test_sets=test_set_s2, h=12)

mse_ea12_s2 <- rbind(mse_u_ea12_s2, mse_i_ea12_s2[-1,], mse_imint_ea12_s2[-1,])
rownames(mse_ea12_s2) <- row_names

basef_mse_ea12_s2 <- mse_ea12_s2[1, ]
pc12_ea_s2 <- mse_ea12_s2[-(1:6), ]
pc12_ea_s2 <- sweep(pc12_ea_s2, 2, basef_mse_ea12_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc12_ea_s2 <- rbind(basef_mse_ea12_s2, pc12_ea_s2)
rownames(pc12_ea_s2) <- rownames(mse_ea12_s2[-(2:6),])
round(pc12_ea_s2,3)

as.data.frame(pc12_ea_s2) %>% 
  write.csv("pc/pc12_ea_s2.csv", row.names = TRUE)

### h=24
mse_u_ea24_s2 <- accuracy_check(forecasts=eaf_U_s2, test_sets=test_set_s2, h=24)
mse_i_ea24_s2 <- accuracy_check(forecasts=eaf_i_s2, test_sets=test_set_s2, h=24)
mse_imint_ea24_s2 <- accuracy_check(forecasts=eaf_imint_s2, test_sets=test_set_s2, h=24)

mse_ea24_s2 <- rbind(mse_u_ea24_s2, mse_i_ea24_s2[-1,], mse_imint_ea24_s2[-1,])
rownames(mse_ea24_s2) <- row_names

basef_mse_ea24_s2 <- mse_ea24_s2[1, ]
pc24_ea_s2 <- mse_ea24_s2[-(1:6), ]
pc24_ea_s2 <- sweep(pc24_ea_s2, 2, basef_mse_ea24_s2, FUN = function(x, y) ((x - y) / y) * 100)
pc24_ea_s2 <- rbind(basef_mse_ea24_s2, pc24_ea_s2)
rownames(pc24_ea_s2) <- rownames(mse_ea24_s2[-(2:6),])
round(pc24_ea_s2,3)

as.data.frame(pc24_ea_s2) %>% 
  write.csv("pc/pc24_ea_s2.csv", row.names = TRUE)

