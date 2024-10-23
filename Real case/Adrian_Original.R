# This R script allows to generate results obtained by Adrian et al.
# Their original code is in Matlab, and can be found here: https://www.openicpsr.org/openicpsr/project/113169/version/V1/view

# NOTE: this script should be run separately for one quarter and four quarter ahead forecast horizons
# by manually changing the variable 'h'.


# Load necessary packages
library(quantreg)
library(pracma)
library(readxl)
library(sn)
library(ggplot2)
library(quantreg)
library(dplyr)

# Clear workspace 
rm(list = ls())

# Set forecast horizon (run script separately for h = 1 and h = 4)
h <- 1

loadsavedresults = TRUE; # TRUE If I have run the code already and I want to load results, stored in ResOOS_H

# Load data and functions
source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/To-Conformalise-or-not-to-Conformalise-in-Growth-at-Risk-Analysis/functions.R")
file_path <- "Data_Original.xls" # Data containing GDP growth and NFCI

# Read the file and  Filter data for 1973Q1-2015Q4
data <- read_excel(file_path)
data <- data[,1:3]
colnames(data)[1] <- "Time"
data$Time <- as.Date(data$Time)
data <- data[data$Time >= as.Date("1973-01-01") & data$Time <= as.Date("2015-10-01"), ]
X <- data[,2:3]
Time <- data$Time

# Set forecast settings
QQ <- seq(0.05, 0.95, by = 0.05) # Vector of quantiles I want to predict
deltaYY <- 0.1
YY <- seq(-20, 20, by = deltaYY)
jtFirstOOS <- which(lubridate::year(data$Time) == 1993 & lubridate::month(data$Time) == 1)

# Construct average growth rates
y <- X$A191RL1Q225SBEA
Yh <- matrix(0, nrow=length(y), ncol=4)

Yh <- stats::filter(y, rep(1/h, h), sides=1) #If h = 1, y = Yh
if (h>1){
  Yh[1:(h-1)] <- NA
}
              
# Construct matrices of regressors
Z <- as.matrix(cbind(1, X[,2], y))
ZGDPonly <- cbind(1, y)

if (loadsavedresults == FALSE) {

{
len_time <- length(data$Time)
len_qq <- length(QQ)
len_yy <- length(YY)

# Raw quantiles
YQ_NaNs <- matrix(NA, len_time, len_qq)
YQ_OOS <- YQ_NaNs
YQGDPonly_OOS <- YQ_NaNs

# Probability integral transforms
Pit_NaNs <- rep(NA, len_time)
PitST_OOS <- Pit_NaNs
PitSTGDPonly_OOS <- Pit_NaNs
}

for (jt in jtFirstOOS:(length(Time) - h)) {

  YhRealized <- Yh[jt + h]
  
  if (lubridate::month(Time[jt]) == 1) {
    cat(sprintf("Now computing the real-time predictive densities in %d", lubridate::year(Time[jt])), "\n")
  }
    
  for (jq in 1:length(QQ)) {                                       
    # Quantile regression with both NFCI and GDP, out-of-sample
    b <- rq(Yh[(h + 1):jt] ~ Z[1:(jt - h),-1], QQ[jq])
    YQ_OOS[jt + h, jq] <- Z[jt, ] %*% coef(b)
      
    # Quantile regression with GDP only, out-of-sample
    bGDPonly <- rq(Yh[(h + 1):jt] ~ ZGDPonly[1:(jt- h),-1], tau=QQ[jq])
    YQGDPonly_OOS[jt + h, jq] <- ZGDPonly[jt, ] %*% coef(bGDPonly)

  }
    
  # Fit skewed-t distribution for quantile regression with NFCI and GDP, out-of-sample
  qqTarg <- YQ_OOS[jt + h, ]
  params <- QuantilesInterpolation(qqTarg, QQ)
  PitST_OOS[jt + h] <- sn::pst(YhRealized, params$lc, params$sc, params$sh, params$df) # is the probability to observe a value < of YhRealized in this distribution 
    
  # Fit skewed-t distribution for quantile regression with GDP only, out-of-sample
  qqTarg <- YQGDPonly_OOS[jt + h, ]
  params_GDPonly <- QuantilesInterpolation(qqTarg, QQ)
  PitSTGDPonly_OOS[jt + h] <- sn::pst(YhRealized, params_GDPonly$lc, params_GDPonly$sc, params_GDPonly$sh, params_GDPonly$df) # is the probability to observe a value < of YhRealized in this distribution 
    

  }
  
} 


if(loadsavedresults == FALSE) { # Save results
  
  filename <- paste("ResOOS_H", h, ".RData", sep="")
  cat(paste("Saving results to file", filename, "\n"))
  save(
    YQ_OOS, YQGDPonly_OOS,PitST_OOS, PitSTGDPonly_OOS,
    file=filename
  )
  
} else { # Load results
  filename <- paste("ResOOS_H", h, ".RData", sep="")
  cat(paste("Loading results from file", filename, "\n"))
  
  load(filename)
}

#------------------------------
 
# # (c)/(d): PITs ggplot

rvec <- seq(0, 1, by = 0.001)
zST_ecdf1 <- PITtest(PitST_OOS, rvec)
zSTGDPonly_ecdf1 <- PITtest(PitSTGDPonly_OOS, rvec)

# Create data frames for both datasets

df1 <- data.frame(
  Quantile = rvec,
  EmpiricalCoverage = zST_ecdf1,
  Group = "NFCI Adrian et al."
)
  
df2 <- data.frame(
  Quantile = rvec,
  EmpiricalCoverage = zSTGDPonly_ecdf1,
  Group = "GDPonly Adrian et al."
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
  scale_color_manual(values = c("NFCI Adrian et al." = cqr_colors, "GDPonly Adrian et al." = qr_colors)) + # Manual color scale
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
# Save the plot as a PDF file
#ggsave(filename = paste0("Adrian_Original_h", h,".pdf"), plot = p, width = 7, height = 5)




#------------------------------

# Plot quantile estimates and quantile levels 
len <- (jtFirstOOS+h):length(Time)
Time_Pred <- Time[len]
Quant_Pred <- YQ_OOS[len,]
Real <- Yh[len]
plot(Time_Pred, Quant_Pred[,1], type = 'l', col = 'blue', xlab = 'Time', ylab = 'QQ', xlim = range(Time_Pred),ylim = c(-20,20))
lines(Time_Pred, Quant_Pred[,10], type = 'l', col = 'blue')
lines(Time_Pred, Quant_Pred[,19], type = 'l', col = 'blue')

lines(Time_Pred, Real, type = 'l', col = 'red', lty = 2)

legend("bottomleft", 
       legend = c("Realization", "5th Percentile Adrian", "50th Percentile Adrian", "95th Percentile Adrian"), 
       col = c("red", "blue", "blue", "blue"), 
       lty = c(2, 1, 1, 1), 
       bty = "n",
       cex = 0.8)  # `bty = "n"` removes the box around the legend

ggsave(filename = paste0("Adrian_Estimates_Plot_H", h,".pdf"), plot = p, width = 7, height = 5)










library(ggplot2)
library(dplyr)

# Create a data frame for quantile predictions and real values
data <- data.frame(
  Time = Time_Pred,
  Real = Real,
  Q5 = Quant_Pred[, 1],
  Q50 = Quant_Pred[, 10],
  Q95 = Quant_Pred[, 19]
)

# Gather data into long format for ggplot2
data_long <- data %>%
  tidyr::gather(key = "Quantile", value = "Value", Q5, Q50, Q95)

# Plot using ggplot
plott <- ggplot(data = data_long, aes(x = Time)) +
  geom_line(aes(y = Value, color = Quantile, linetype = Quantile)) + # for quantile lines
  geom_line(aes(y = Real, color = "Realization"), linetype = "dashed") + # for real line
  scale_color_manual(values = c("Realization" = "red", "Q5" = "blue", "Q50" = "blue", "Q95" = "blue")) +
  scale_linetype_manual(values = c("Realization" = "dashed", "Q5" = "solid", "Q50" = "solid", "Q95" = "solid")) +
  labs(x = "Time", y = "QQ", color = "Legend", linetype = "Legend") +
  ylim(-20, 20) +
  theme_minimal() +
  theme(legend.position = "bottomleft", legend.title = element_blank())

ggsave(filename = paste0("Adrian_Estimates_Plot_H", h,".pdf"), plot = plott, width = 7, height = 5)


