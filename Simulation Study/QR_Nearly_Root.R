library(caret)
library(sn)
library(quantreg)
library(readxl)
library(forecast)
library(SuppDists)
library(readr)
library(dplyr)
library(parallel)
library(foreach)
library(doParallel)
library(progressr)
library(future)
library(future.apply)
library(progressr)
library(quantregForest)
library(openxlsx)

#Load correctly the file "functions.R", modifying the path
source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/To-Conformalise-or-not-to-Conformalise-in-Growth-at-Risk-Analysis/functions.R")

# File to save results of the simulation.
file_path <- "QR_Nearly_Root_Results.xlsx"

if (!file.exists(file_path)) {
  wb <- createWorkbook()
  addWorksheet(wb, "Data")
  saveWorkbook(wb, file_path, overwrite = TRUE)
  
} else {
  wb <- loadWorkbook(file_path)
}

run_simulation <- function(n, phi_ar2){
  
  n2 <- 100 # Number of test points
  
#------------------- Generate n + n2 points for the AR(1) Nearly Unit Root model, this is our DGP
  
  Y <- numeric(n+n2)
  Y[1] <- rnorm(1)
  for (i in 2:(n+n2)) {
    Y[i] <- phi_ar2 * Y[i-1] + rnorm(1)
  }
  # Create lagged variable
  Y_lag1 <- c(NA, Y[-length(Y)])

  # Prepare dataset for model excluding the first three NA values
  indices <- c((n+1):(n+n2))
  data <- data.frame(I = 1, Y = Y[-c(1:3,indices)], Y_lag1 = Y_lag1[-c(1:3,indices)])

  data_test <- data.frame(I = 1, Y = Y[indices], Y_lag1 = Y_lag1[indices])

  Y_test <- Y[c((n+1):(n+n2))]
  
  QQ <- c(seq(0.05, 0.95, by = 0.05),0.99) #Quantiles vector I want to estimate = alpha
  
#--------------- QUANTILE REGRESSION 
  
  QuantAR1_OOS <- matrix(0, n2, length(QQ))
  
  for (jq in 1:length(QQ)) {  
      
      QR1 <- rq(Y ~ Y_lag1, data = data, tau=QQ[jq])
      QuantAR1_OOS[,jq] <- as.matrix(data_test[,-2]) %*% coef(QR1)
      
  }

#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR1_OOS <- matrix(NA, n2, length(QQ))
  full_length <- nrow(data)
  test_length = floor(full_length*50/100)
  data_1 <- data[1:test_length,] #Train
  data_2 <- data[(test_length + 1) : full_length,] #Calibration

  x0 <- data_1[,-2]
  y0 <- data_1[,2]
  x1 <- data_2[,-2]
  y1 <- data_2[,2]
  x_test <- data_test[,-2]

  CQuantAR1_OOS <- qCQR_opposite(x0,y0,x1,y1,x_test,QQ)
  
  #---------------- CALIBRATION OF QR AND CQR
  
  coverageQRAR1 <- compute_coverage(Y_test, QuantAR1_OOS, QQ)
  coverageCQRAR1 <- compute_coverage(Y_test, CQuantAR1_OOS, QQ)
  
  ret <- list(coverageQRAR1, coverageCQRAR1)
  
  return(ret)
}


#------------------- Main Simulation Loop

vector_n <- c(101,201,1001)
vector_phi <- c(0.95, 1, 1.05)
QQ <- c(seq(0.05, 0.95, by = 0.05),0.99) #Quantiles vector I want to estimate = alpha
n2 <- 100  # Test set size
n3 <- 100  # Number of simulations
seeds <- 1:n3  
count <- 2  # Row counter for writing results to Excel


for (n in vector_n){
  for(phi in vector_phi){

  cl <- makeCluster(detectCores() - 1) 
  clusterExport(cl, varlist=c("run_simulation","n","phi")) 
  clusterEvalQ(cl, { 
    library(readxl)
    library(quantreg)
    library(quantregForest)
    source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/To-Conformalise-or-not-to-Conformalise-in-Growth-at-Risk-Analysis/functions.R")
  }) 

  # Run simulations in parallel
  results <- parLapply(cl, seeds, function(seed) {
    set.seed(seed)
    run_simulation(n, phi)
  })
  
  stopCluster(cl)

  # Extract results
  coveragetotQRAR1 <- matrix(NA,n3,length(QQ))
  coveragetotCQRAR1 <- matrix(NA,n3,length(QQ))

  index = 1
  for(res in results){
    coveragetotQRAR1[index,] <- res[[1]]
    coveragetotCQRAR1[index,] <- res[[2]]
    index <- index + 1
  }

  average_coverageQRAR1 <- compute_average_coverage(coveragetotQRAR1, QQ) 
  average_coverageCQRAR1 <- compute_average_coverage(coveragetotCQRAR1, QQ) 

  resultsPitSTQRAR1 <- compute_results(average_coverageQRAR1, n2*n3, QQ)
  resultsPitSTCQRAR1 <- compute_results(average_coverageCQRAR1, n2*n3, QQ)

  
#------------ Write Results in the excel file

  wb <- loadWorkbook(file_path)

  sheet <- "Data"

  #In the first column, write the parameters of the simulation:
  writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", phi = ", phi, ", QR AR(1)" ), startCol = 1, startRow = count, colNames = FALSE)
  writeData(wb, sheet = sheet, x = paste("n1 = ", n-3,", phi = ", phi, ", CQR AR(1)" ), startCol = 1, startRow = count+1, colNames = FALSE)

  #In the 2°,3°,4° columns, put how many inside, below and above CI
  writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$IsWithinCI) , startCol = 2, startRow = count, colNames = FALSE)
  writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$EmpiricalCoverage > resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI) , startCol = 3, startRow = count, colNames = FALSE)
  writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR1$EmpiricalCoverage < resultsPitSTQRAR1$Quantile & !resultsPitSTQRAR1$IsWithinCI) , startCol = 4, startRow = count, colNames = FALSE)

  writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$IsWithinCI) , startCol = 2, startRow = count + 1, colNames = FALSE)
  writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$EmpiricalCoverage > resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI) , startCol = 3, startRow = count + 1, colNames = FALSE)
  writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR1$EmpiricalCoverage < resultsPitSTCQRAR1$Quantile & !resultsPitSTCQRAR1$IsWithinCI) , startCol = 4, startRow = count + 1, colNames = FALSE)

  #In the 5° column, the MAE will be placed
  writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR1$Quantile-resultsPitSTQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
  writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR1$Quantile-resultsPitSTCQRAR1$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)

  saveWorkbook(wb, file_path, overwrite = TRUE)

  count <- count + 3 #So to leave an empty row

#------------- PLOT CALIBRATION CURVES OF QR and CQR
  
  # Create data frames for both datasets
  
  df1 <- data.frame(
    Quantile = rev(resultsPitSTCQRAR1$Quantile),
    EmpiricalCoverage = rev(resultsPitSTCQRAR1$EmpiricalCoverage),
    Group = "CQR QR"
  )
  
  df2 <- data.frame(
    Quantile = rev(resultsPitSTQRAR1$Quantile),
    EmpiricalCoverage = rev(resultsPitSTQRAR1$EmpiricalCoverage),
    Group = "QR"
  )
  
  df <- bind_rows(df1, df2)
  
  cqr_colors <- "#0000FF"
  qr_colors <- "#FF0000"
  
  legend_position <- if (n - 3 == 98) {
    c(0.85, 0.15)
  } else {
    "none"
  }
  
  p <- ggplot(df, aes(x = Quantile, y = EmpiricalCoverage, color = Group)) +
    geom_step(direction = "vh", linewidth = 1) + 
    geom_point(size = 3) + 
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) + 
    scale_color_manual(values = c("CQR QR" = cqr_colors, "QR" = qr_colors)) + 
    labs(title = paste("n =", n - 3, "phi =", phi), x = "Quantile Levels", y = "Empirical Coverage") + 
    theme_minimal() +
    theme(
      legend.position = legend_position, 
      legend.title = element_blank(),
      text = element_text(size = 15),
      plot.title = element_text(hjust = 0.5),
      axis.title.x = element_text(hjust = 0.5),
      axis.title.y = element_text(hjust = 0.5)
    )
  
  print(p)
  ggsave(filename = paste0("QR_Nearly_Root_n", n - 3,"_", phi, ".pdf"), plot = p, width = 7, height = 5)
  
  }
}