# This script replicates the analysis of Adrian et al., but implementing the CQR QR algorithm instead of QR

# NOTE: this script should be run separately for one quarter and four quarter ahead forecast horizons
# by manually changing the variable 'h'.

library(quantreg)
library(lubridate)
library(pracma)
library(readxl)
library(sn)
library(parallel)

# Clear workspace 
rm(list = ls())

# Set forecast horizon
h <- 1

loadsavedresults = FALSE; # TRUE If I have run the code already and I want to load results, stored in ResOOS_CQR_QR_H

# Load data and functions
source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/To-Conformalise-or-not-to-Conformalise-in-Growth-at-Risk-Analysis/functions.R")
file_path <- "DataVulnerabilityAppendix.xls"

# Read the file and  Filter data for 1973Q1-2015Q4
data <- read_excel(file_path)
data <- data[,1:3]
colnames(data)[1] <- "Time"
data$Time <- as.Date(data$Time)
data <- data[data$Time >= as.Date("1973-01-01") & data$Time <= as.Date("2015-10-01"), ]
X <- data[,2:3]
Time <- data$Time

# Set forecast settings
QQ <- seq(0.03, 0.99, by = 0.01) #We don't implement a smoothing, as in Adrian et al., but use a much finer quantile grid
deltaYY <- 0.1
YY <- seq(-20, 20, by = deltaYY)
jtFirstOOS <- which(year(data$Time) == 1993 & month(data$Time) == 1)

# Construct average growth rates
y <- X$A191RL1Q225SBEA
Yh <- matrix(0, nrow=length(y), ncol=4)

Yh <- stats::filter(y, rep(1/h, h), sides=1) #If h = 1, y = Yh
if (h>1){
  Yh[1:(h-1)] <- NA
}

#Construct matrices of regressors
Z <- as.matrix(cbind(1, X[,2], y))
ZGDPonly <- cbind(1, y)



if (loadsavedresults == FALSE) {
  
  {
  
  len_time <- length(data$Time)
  len_qq <- length(QQ)
  len_yy <- length(YY)

  # Raw quantiles
  YQ_NaNs <- matrix(NA, len_time, len_qq)

  YQ_low_OOSCO <- YQ_NaNs
  YQ_high_OOSCO <- YQ_NaNs
  YQ_OOSCO <- YQ_NaNs

  YQGDPonly_high_OOSCO <- YQ_NaNs
  YQGDPonly_low_OOSCO <- YQ_NaNs
  YQGDPonly_OOSCO <- YQ_NaNs

  # Probability integral transforms
  Pit_NaNs <- rep(NA, len_time)
  PitST_OOSCO <- Pit_NaNs
  PitSTGDPonly_OOSCO <- Pit_NaNs

  }

  for (jt in jtFirstOOS:(length(Time) - h)) { 
  
    YhRealized <- Yh[jt + h]

    if (lubridate::month(Time[jt]) == 1) {
      cat(sprintf("Now computing the real-time predictive densities in %d", year(Time[jt])), "\n")
    }
      
    #------- Conformalized Quantile Regression with both NFCI and GDP, out-of-sample
        
    # Split creating I1 and I2
    full_length <- length(Yh[(h + 1):jt])
    test_length = full_length*50/100
    y0 <- Yh[(h+1):(h+test_length)]
    y1 <- Yh[(h+1+test_length):jt] 
    
    x0 <- Z[1:test_length,]
    x1 <- Z[(test_length+1):(jt - h),] 
    x_test <- t(as.matrix(Z[jt,]))
    
    xGDPonly0 <- ZGDPonly[1:test_length,]
    xGDPonly1 <- ZGDPonly[(test_length+1):(jt - h),]
    xGDPonly_test <- t(as.matrix(ZGDPonly[jt,]))
        
    YQ_OOSCO[jt + h, ] <- qCQR_opposite(x0,y0,x1,y1,x_test,QQ)
    YQGDPonly_OOSCO[jt + h, ] <- qCQR_opposite(xGDPonly0,y0,xGDPonly1,y1,xGDPonly_test,QQ)
      
    PitST_OOSCO[jt + h] <- pst(YhRealized, QQ, YQ_OOSCO[jt + h, ]) # is the probability to observe a value < of YhRealized in this distribution 
    PitSTGDPonly_OOSCO[jt + h] <- pst(YhRealized, QQ, YQGDPonly_OOSCO[jt + h, ]) # is the probability to observe a value < of YhRealized in this distribution 

  }

}



if(loadsavedresults == FALSE) { # Save results
  
  filename <- paste("ResOOS_CQR_QR_H", h, ".RData", sep="")
  cat(paste("Saving results to file", filename, "\n"))
  save(
    YQ_OOSCO, YQGDPonly_OOSCO,PitST_OOSCO, PitSTGDPonly_OOSCO,
    file=filename
  )
  
} else { # Load results
  filename <- paste("ResOOS_CQR_QR_H", h, ".RData", sep="")
  cat(paste("Loading results from file", filename, "\n"))
  
  load(filename)
}
 

#------------------------------

# # (c)/(d): PITs ggplot (This plot is not present on the publication, it is just to get a visual output)

rvec <- seq(0, 1, by = 0.001)
zST_ecdf1 <- PITtest(PitST_OOSCO, rvec)
zSTGDPonly_ecdf1 <- PITtest(PitSTGDPonly_OOSCO, rvec)

# Create data frames for both datasets

df1 <- data.frame(
  Quantile = rvec,
  EmpiricalCoverage = zST_ecdf1,
  Group = "CQR QR NFCI"
)

df2 <- data.frame(
  Quantile = rvec,
  EmpiricalCoverage = zSTGDPonly_ecdf1,
  Group = "CQR QR GDPonly"
)

# Combine the data frames
df <- bind_rows(df1, df2)

# Define the colors
cqr_colors <- "#FFA500"
qr_colors <- "#008080"

# Define the legend position based on the value of n - 3
legend_position <- c(0.75, 0.15)

p <- ggplot(df, aes(x = Quantile, y = EmpiricalCoverage, color = Group)) +
  geom_step(direction = "vh", linewidth = 1) + # Use `linewidth` instead of `size` for lines
  geom_point(size = 0.5) + # Add points
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black", linewidth = 1) + # Use `linewidth` for the diagonal line
  scale_color_manual(values = c("CQR QR NFCI" = cqr_colors, "CQR QR GDPonly" = qr_colors)) + # Manual color scale
  labs(title = paste("h =", h), x = "Quantile Levels", y = "Empirical Coverage") + # Add labels and title
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

