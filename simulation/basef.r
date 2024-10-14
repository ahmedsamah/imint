# utility functions
library(forecast)
set.seed(42)
hts <- function(basis, sMat) {
  t(sMat %*% t(basis))
}

cal_basef <- function(series, method = 'arima'){
  series = hts(series, sMat)
  apply(series, 2, function(series){
    series = ts(series, frequency = 12)
    if (method == 'arima'){
      model = auto.arima(series)
    }
    if (method == 'ets'){
      model = ets(series)
    }
    c(fitted(model), forecast(model, h=24)$mean)
  })
}

# settings
sMat = rbind(matrix(c(1, 1, 0, 1, 1, 0, 1, 0, 1, 1, 0, 1), 3, 4),
             diag(rep(1, 4)))


# scenario
data1 = read.csv('data/data_scenario1.csv', row.names = 1)

# Function to run forecasting for a single index
run_forecast <- function(index) {
  series <- data1[data1$index == index,][1:300, 1:4]
  
  basef1 <- data.frame(cal_basef(series, method = 'arima'), 
                       method = 'arima', 
                       index = index, 
                       t = 1:324)
  
  basef2 <- data.frame(cal_basef(series, method = 'ets'), 
                       method = 'ets', 
                       index = index, 
                       t = 1:324)
  
  rbind(basef1, basef2)
}

library(parallel)
numCores <- detectCores() - 1 
cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(forecast) 
})

clusterExport(cl, list("cal_basef", "hts", "data1", "sMat", "run_forecast")) # Export necessary variables and functions to each worker

df_list <- parLapply(cl, 1:1000, run_forecast)  # Run in parallel
df <- do.call(rbind, df_list)  # Combine results
stopCluster(cl)

write.csv(df, 'basef/basef1.csv')

# scenario 2
data2 = read.csv('data/data_scenario2.csv', row.names = 1)

run_forecast <- function(index) {
  series <- data2[data2$index == index,][1:300, 1:4]
  
  basef1 <- data.frame(cal_basef(series, method = 'arima'), 
                       method = 'arima', 
                       index = index, 
                       t = 1:324)
  
  basef2 <- data.frame(cal_basef(series, method = 'ets'), 
                       method = 'ets', 
                       index = index, 
                       t = 1:324)
  
  rbind(basef1, basef2)
}
numCores <- detectCores() - 1 
cl <- makeCluster(numCores)

clusterEvalQ(cl, {
  library(forecast)
})

clusterExport(cl, list("cal_basef", "hts", "data2", "sMat", "run_forecast")) 

df_list <- parLapply(cl, 1:1000, run_forecast)  # Run in parallel
df <- do.call(rbind, df_list)  # Combine results
stopCluster(cl)

write.csv(df, 'basef/basef2.csv')

