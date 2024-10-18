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
library(ggplot2)

#Load correctly the file "functions.R", modifying the path
source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")

file_path <- "QR_Cauchy_Results.xlsx"

# Check if the file exists
if (!file.exists(file_path)) {
  # Create a new workbook
  wb <- createWorkbook()
  
  # Add a worksheet named "Data"
  addWorksheet(wb, "Data")
  
  # Save the workbook (this creates the file)
  saveWorkbook(wb, file_path, overwrite = TRUE)
} else {
  # Load the existing workbook
  wb <- loadWorkbook(file_path)
}


run_simulation <- function(n){
  
  n2 <- 100 # Number of test points
  
#------------------- Generate n + n2 points for the AR(2) Cauchy error model, this is our DGP
  
  phi_ar2 <- c(0.5, -0.2)  # AR coefficients for AR(2)
  
  Y_ar2 <- numeric(n+n2)
  Y_ar2[1:2] <- rnorm(2) #Two random starting points
  for (i in 3:(n+n2)) {
    Y_ar2[i] <- phi_ar2[1] * Y_ar2[i-1] + phi_ar2[2] * Y_ar2[i-2] + rt(n = 1, df = 1) #Cauchy Errors!
  }
  # Create lagged variables
  Y_lag1_ar2 <- c(NA, Y_ar2[-length(Y_ar2)])
  Y_lag2_ar2 <- c(NA, NA, Y_ar2[1:(length(Y_ar2)-2)])

  # Prepare dataset, excluding the first three NA values
  indices <- c((n+1):(n+n2))
  data <- data.frame(I = 1, Y = Y_ar2[-c(1:3,indices)], Y_lag1 = Y_lag1_ar2[-c(1:3,indices)], Y_lag2 = Y_lag2_ar2[-c(1:3,indices)])

  data_test <- data.frame(I = 1, Y = Y_ar2[indices], Y_lag1 = Y_lag1_ar2[indices], Y_lag2 = Y_lag2_ar2[indices])

  Y_test <- Y_ar2[c((n+1):(n+n2))]
  
  QQ <- c(seq(0.05, 0.95, by = 0.05),0.99) #Quantiles vector I want to estimate = alpha


#--------------- QUANTILE REGRESSION 
   
  QuantAR2_OOS <- matrix(0, n2, length(QQ))
   
  for (jq in 1:length(QQ)) {  
   
      QR2 <- rq(Y ~ Y_lag1 + Y_lag2 , data = data, tau=(QQ[jq]))
      QuantAR2_OOS[,jq] <- as.matrix(data_test[,-c(2,5)]) %*% coef(QR2)
   
  }

#--------------- CONFORMALIZED QUANTILE REGRESSION
  
  CQuantAR2_OOS <- matrix(0, n2, length(QQ))

  full_length <- nrow(data)
  test_length = floor(full_length*50/100)
  data_1 <- data[1:test_length,] #Train
  data_2 <- data[(test_length + 1) : full_length,] #Calibration

  x0 <- data_1[,-c(2,5)]
  y0 <- data_1[,2]
  x1 <- data_2[,-c(2,5)]
  y1 <- data_2[,2]
  x_test <- data_test[,-c(2,5)]

  CQuantAR2_OOS <- qCQR_opposite(x0,y0,x1,y1,x_test,QQ)
         
  #---------------- CALIBRATION OF QR AND CQR

  coverageQRAR2 <- compute_coverage(Y_test, QuantAR2_OOS, QQ)
  coverageCQRAR2 <- compute_coverage(Y_test, CQuantAR2_OOS, QQ)
  
  ret <- list(coverageQRAR2, coverageCQRAR2)
  
  return(ret)
}


#------------------- Main Simulation Loop

vector_n <- c(101,201,1001)
QQ <- c(seq(0.05, 0.95, by = 0.05),0.99) #Quantiles vector I want to estimate = alpha
n2 <- 100  # Test set size
n3 <- 100  # Number of simulations
seeds <- 1:n3  # Vector of seeds, one per simulation
count <- 2  # Row counter for writing results to Excel

for (n in vector_n){

 # Setup parallel cluster
 cl <- makeCluster(detectCores() - 1) # Leave one core free to avoid freezing your system
 clusterExport(cl, varlist=c("run_simulation","n")) # Export the simulation function to each cluster node
 clusterEvalQ(cl, { 
   library(readxl)
   library(quantreg)
   library(quantregForest)
   source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/Thesis-Code/functions.R")
 }
 ) 

 # Run simulations in parallel
 results <- parLapply(cl, seeds, function(seed) {
  set.seed(seed)
  run_simulation(n)
 })


# Stop the cluster
stopCluster(cl)

# Extract results
coveragetotQRAR2 <- matrix(NA,n3,length(QQ))
coveragetotCQRAR2 <- matrix(NA,n3,length(QQ))

index = 1
for(res in results){
  coveragetotQRAR2[index,] <- res[[1]]
  coveragetotCQRAR2[index,] <- res[[2]]
  index <- index + 1
}


average_coverageQRAR2 <- compute_average_coverage(coveragetotQRAR2, QQ)
average_coverageCQRAR2 <- compute_average_coverage(coveragetotCQRAR2, QQ)

resultsPitSTQRAR2 <- compute_results(average_coverageQRAR2, n2*n3, QQ)
resultsPitSTCQRAR2 <- compute_results(average_coverageCQRAR2, n2*n3, QQ)

#------------ Write Results in the excel file

wb <- loadWorkbook(file_path)

# Access the worksheet 
sheet <- "Data"

#In the first column, write the parameters of the simulation:
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", QR AR(2)" ), startCol = 1, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = paste("n1 = ", n-3, ", CQR AR(2)" ), startCol = 1, startRow = count+1, colNames = FALSE)

#In the 2°,3°,4° columns, put how many inside, below and above CI

writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR2$IsWithinCI) , startCol = 2, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR2$EmpiricalCoverage > resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI) , startCol = 3, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTQRAR2$EmpiricalCoverage < resultsPitSTQRAR2$Quantile & !resultsPitSTQRAR2$IsWithinCI) , startCol = 4, startRow = count, colNames = FALSE)

writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR2$IsWithinCI) , startCol = 2, startRow = count+1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR2$EmpiricalCoverage > resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI) , startCol = 3, startRow = count + 1, colNames = FALSE)
writeData(wb, sheet = sheet, x = sum(resultsPitSTCQRAR2$EmpiricalCoverage < resultsPitSTCQRAR2$Quantile & !resultsPitSTCQRAR2$IsWithinCI) , startCol = 4, startRow = count + 1, colNames = FALSE)

#In the 5° column, the MAE will be placed

writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTQRAR2$Quantile-resultsPitSTQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count, colNames = FALSE)
writeData(wb, sheet = sheet, x = mean(abs(resultsPitSTCQRAR2$Quantile-resultsPitSTCQRAR2$EmpiricalCoverage)), startCol = 5, startRow = count+1, colNames = FALSE)

saveWorkbook(wb, file_path, overwrite = TRUE)

count <- count + 3 #So to leave an empty row


#------------- PLOT CALIBRATION CURVES OF QR and CQR 

# Create data frames for both datasets
df1 <- data.frame(
  Quantile = rev(resultsPitSTCQRAR2$Quantile),
  EmpiricalCoverage = rev(resultsPitSTCQRAR2$EmpiricalCoverage),
  Group = "CQR QR"
)

df2 <- data.frame(
  Quantile = rev(resultsPitSTQRAR2$Quantile),
  EmpiricalCoverage = rev(resultsPitSTQRAR2$EmpiricalCoverage),
  Group = "QR"
)

# Combine the data frames
df <- bind_rows(df1, df2)

# Define the colors
cqr_colors <- "#0000FF"
qr_colors <- "#FF0000"
# Define the legend position based on the value of n - 3
legend_position <- if ((n - 3) == 98) c(0.85, 0.15) else "none"

p <- ggplot(df, aes(x = Quantile, y = EmpiricalCoverage, color = Group)) +
  geom_step(direction = "vh", size = 1) + # Add staircase lines
  geom_point(size = 3) + # Add points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", size = 1) + # Add diagonal line
  scale_color_manual(values = c("CQR QR" = cqr_colors, "QR" = qr_colors)) + # Manual color scale
  labs(title = paste("n =", n - 3), x = "Quantile Levels", y = "Empirical Coverage") + # Add labels and title
  theme_minimal() +
  theme(
    legend.position = legend_position, # Set legend position based on condition
    legend.title = element_blank(),
    text = element_text(size = 15),
    plot.title = element_text(hjust = 0.5),
    axis.title.x = element_text(hjust = 0.5),
    axis.title.y = element_text(hjust = 0.5)
  )

# Print the plot
print(p)
# Save the plot as a PDF file
ggsave(filename = paste0("QR_Cauchy_n", n - 3, ".pdf"), plot = p, width = 7, height = 5)
}

