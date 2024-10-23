#This script....

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
h <- 4

loadsavedresults = FALSE; # TRUE If I have run the code already and I want to load results, stored in ResOOS_CQR_QR_H

# Load data and functions
source("C:/Users/Pietro/Desktop/Pietro/Politecnico/Tesi/To-Conformalise-or-not-to-Conformalise-in-Growth-at-Risk-Analysis/functions.R")
file_path <- "Data_NFCI_div.xls"

# Read the file and  Filter data for 1973Q1-2015Q4
data <- read_excel(file_path)
data <- data[,-c(3:10)] # Making sure to have only components of NFCI i
colnames(data)[1] <- "Time"
data$Time <- as.Date(data$Time)
data <- data[data$Time >= as.Date("1973-01-01") & data$Time <= as.Date("2015-10-01"), ]
X <- data[,-1]
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

# Construct matrices of regressors
Z <- as.matrix(cbind(1, X[,-1], y))
ZGDPonly <- cbind(1, y)



if (loadsavedresults == FALSE) {
  
  {
    
    len_time <- length(data$Time)
    len_qq <- length(QQ)
    len_yy <- length(YY)
    
    # Raw quantiles
    YQ_NaNs <- matrix(NA, len_time, len_qq)
    
    YQ_low_OOSCOdiv <- YQ_NaNs
    YQ_high_OOSCOdiv <- YQ_NaNs
    YQ_OOSCOdiv <- YQ_NaNs
    
    YQGDPonly_high_OOSCOdiv <- YQ_NaNs
    YQGDPonly_low_OOSCOdiv <- YQ_NaNs
    YQGDPonly_OOSCOdiv <- YQ_NaNs
    
    # Probability integral transforms
    Pit_NaNs <- rep(NA, len_time)
    PitST_OOSCOdiv <- Pit_NaNs
    PitSTGDPonly_OOSCOdiv <- Pit_NaNs
    
  }
  
  for (jt in jtFirstOOS:(length(Time) - h)) { 
    
    YhRealized <- Yh[jt + h]
    
    if (lubridate::month(Time[jt]) == 1) {
      cat(sprintf("Now computing the real-time predictive densities in %d", year(Time[jt])), "\n")
    }
    
    #------- CQR QR with NFCI components and GDP, out-of-sample
    
    # Split creating I1 and I2
    full_length <- length(Yh[(h + 1):jt])
    test_length = full_length*50/100
    
    xGDPonly0 <- ZGDPonly[1:test_length,]
    xGDPonly1 <- ZGDPonly[(test_length + 1):(jt - h),]
    xGDPonly_test <- t(as.matrix(ZGDPonly[jt,]))
    
    y0 <- Yh[(h+1):(h+test_length)]
    y1 <- Yh[(h+1+test_length):jt]
    
    # I extract the "y" value (previous GDP value), to be added after PCA
    x0_y <- Z[1:test_length,"y"]
    x1_y <- Z[(test_length+1):(jt - h),"y"]
    Z_y <- Z[,"y"]
    
    # I extract the observations I need at this time step "jt" and I remove the columns that do not change
    Z_temp <- Z[1:jt, apply(Z[1:jt,], 2, function(x) var(x) > 1e-5) ]
    
    # I scale the data for PCA
    Z_temp <- scale(Z_temp) 
    
    
    x0 <- Z_temp[1:test_length,]
    x1 <- Z_temp[(test_length+1):(jt - h),] 


    # ---- Preparation for PCA

    pca_result <- prcomp(x0, center = TRUE, scale. = TRUE)
    cumulative_variance <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2)) 
    num_components <- which(cumulative_variance >= 0.90)[1] # Capture 90% of total variance
    selected_components <- pca_result$x[, 1:num_components]
     
    x1_projected <- predict(pca_result, newdata = x1)[, 1:num_components]
    Z_temp_projected <- predict(pca_result, newdata = Z_temp)[, 1:num_components]
     
    x0 <- as.matrix(cbind(cbind(selected_components),x0_y))
    x1 <- as.matrix(cbind(cbind(x1_projected),x1_y))
    Z_temp <- as.matrix(cbind(cbind(Z_temp_projected),Z_y))
    
    x_test <- t(as.matrix(Z_temp[jt,]))
    

    YQ_OOSCOdiv[jt + h, ] <- qCQR_opposite(x0,y0,x1,y1,x_test,QQ)
    YQGDPonly_OOSCOdiv[jt + h, ] <- qCQR_opposite(xGDPonly0,y0,xGDPonly1,y1,xGDPonly_test,QQ)
    
    PitST_OOSCOdiv[jt + h] <- pst(YhRealized, QQ, YQ_OOSCOdiv[jt + h, ]) # is the probability to observe a value < of YhRealized in this distribution 
    PitSTGDPonly_OOSCOdiv[jt + h] <- pst(YhRealized, QQ, YQGDPonly_OOSCOdiv[jt + h, ]) # is the probability to observe a value < of YhRealized in this distribution 
    
  }
  
}





if(loadsavedresults == FALSE) { # Save results
  
  filename <- paste("ResOOS_CQR_QR_NFCI_div_H", h, ".RData", sep="")
  cat(paste("Saving results to file", filename, "\n"))
  save(
    YQ_OOSCOdiv, YQGDPonly_OOSCOdiv,PitST_OOSCOdiv, PitSTGDPonly_OOSCOdiv,
    file=filename
  )
  
} else { # Load results
  filename <- paste("ResOOS_CQR_QR_NFCI_div_H", h, ".RData", sep="")
  cat(paste("Loading results from file", filename, "\n"))
  
  load(filename)
}


#------------------------------

# # (c)/(d): PITs ggplot (This plot is not present on the publication, it is just to get a visual output)

rvec <- seq(0, 1, by = 0.001)
zST_ecdf1 <- PITtest(PitST_OOSCOdiv, rvec)
zSTGDPonly_ecdf1 <- PITtest(PitSTGDPonly_OOSCOdiv, rvec)

# Create data frames for both datasets

df1 <- data.frame(
  Quantile = rvec,
  EmpiricalCoverage = zST_ecdf1,
  Group = "CQR QR NFCI div"
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
  scale_color_manual(values = c("CQR QR NFCI div" = cqr_colors, "CQR QR GDPonly" = qr_colors)) + # Manual color scale
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




#------------------------------

# Plot quantile estimates and quantile levels 
len <- (jtFirstOOS+h):length(Time)
Time_Pred <- Time[len]
Quant_Pred <- YQ_OOSCOdiv[len,]
Real <- Yh[len]
plot(Time_Pred, Quant_Pred[,3], type = 'l', col = 'blue', xlab = 'Time', ylab = 'QQ', xlim = range(Time_Pred),ylim = c(-30,20))
lines(Time_Pred, Quant_Pred[,48], type = 'l', col = 'blue')
lines(Time_Pred, Quant_Pred[,93], type = 'l', col = 'blue')
lines(Time_Pred, Real, type = 'l', col = 'red', lty = 2)

legend("bottomleft", 
       legend = c("Realization", "5th Percentile CQR QR NFCI", "50th Percentile  CQR QR NFCI", "95th Percentile  CQR QR NFCI"), 
       col = c("red", "blue", "blue", "blue"), 
       lty = c(2, 1, 1, 1), 
       bty = "n",
       cex = 0.8)  # `bty = "n"` removes the box around the legend

#ggsave(filename = paste0("CQR_QR_Estimates_Plot_H", h,".pdf"), plot = p, width = 7, height = 5)



























































    if (jt >= jtFirstOOS) { #Now the Out-of-Sample part
      if (month(Time[jt]) == 1) {
        cat(sprintf("Now computing the real-time predictive densities in %d", year(Time[jt])), "\n")
      }
      
      for (jq in 1:length(QQ)) {
        
        #------- Conformalized Quantile Regression with both NFCI and GDP, out-of-sample
        
        #Split creating I1 and I2
        full_length <- length(Yh[(h + 1):jt])
        test_length = full_length*50/100

        
        Yh1 <- Yh[(h+1):(h+test_length)]
        Yh2 <- Yh[(h+1+test_length):jt]
        Z1 <- Z[1:test_length,]
        Z2 <- Z[(test_length+1):(jt - h),]
        Z1_y <- Z1[,"y"]
        Z2_y <- Z2[,"y"]
        Z_y <- Z_temp[,"y"]
        ZGDPonly1 <- ZGDPonly[1:test_length,]
        ZGDPonly2 <- ZGDPonly[(test_length+1):(jt - h),]

        ## Step 2: Perform PCA
        Z1 <- Z1[,-c(1,ncol(Z1))] #devo rimuovere la colonna di 1!!!!
        
        pca_result <- prcomp(Z1, center = TRUE, scale. = TRUE)
        
        
        # Step 3: Calculate cumulative variance explained
        cumulative_variance <- cumsum(pca_result$sdev^2 / sum(pca_result$sdev^2))
        
        # Step 4: Determine the number of components to capture at least 80% of the variance
        num_components <- which(cumulative_variance >= 0.90)[1]
        
        # # Step 5: Extract the relevant principal components (scores)
        selected_components <- pca_result$x[, 1:num_components]

        Z1temp <- as.matrix(cbind(cbind(selected_components),Z1_y))
        Z2 <- Z2[, -c(1, ncol(Z2))] 
        Z2_projected <- predict(pca_result, newdata = Z2)[, 1:num_components]
        
        Z2temp <- as.matrix(cbind(cbind(Z2_projected),Z2_y))
        
        Ztemp <- Z[, -c(1, ncol(Z))]
        Z_projected <- predict(pca_result, newdata = Z)[, 1:num_components]
        
        Ztemp <- as.matrix(cbind(cbind(Z_projected),Z_y))
#         Z1temp <- Z1[, apply(Z1, 2, function(x) var(x) > 1e-5)]  
#         Z2temp <- Z2[, apply(Z1, 2, function(x) var(x) > 1e-5)]  
#         Ztemp <- Z[, apply(Z1, 2, function(x) var(x) > 1e-5)] 
#         
#          #If ncols > nrows, I remove some cols
#          ret <- remove_excess_columns(Z1temp)
#          Z1temp <- ret$modified_matrix
#          col_to_remove <- ret$removed_columns
#          
#          if (!is.null(col_to_remove) && length(col_to_remove) > 0) {
#            Z2temp <- Z2temp[,-col_to_remove]
#            Ztemp <- Ztemp[,-col_to_remove]
#          }
#          
#          
#          #If some cols are too autocorrelated, I remove them
#          ret <- remove_highly_correlated(Z1temp)
#          
#          Z1temp <- ret$cleaned_data
#          col_to_remove <- ret$removed_indices
#          
#          if (!is.null(col_to_remove) && length(col_to_remove) > 0) {
#            Z2temp <- Z2temp[,-col_to_remove]
#            Ztemp <- Ztemp[,-col_to_remove]
#          }
        

         YQ_high_OOSC[(h + 1):(h + test_length), jq] <- Inf 
         
         b_high <- rq(Yh1 ~ Z1temp, tau= (1-QQ[jq])) #Train on I1
         YQ_low_OOSC[(h + 1):(h + full_length-test_length), jq] <- as.vector(cbind(1, Z2temp) %*% coef(b_high)) #Evaluate on I2
         
        # Initialize a vector for errors
        E_i <- rep(NA, length(Yh2))
        
        # Calculate errors for each point in the test set I2
        for (i in 1:length(E_i)) {
          E_i[i] <- max(YQ_low_OOSC[h + i, jq] - Yh2[i], Yh2[i] - YQ_high_OOSC[h + i, jq])
        }
        
        # Compute Q(1-alpha)(E, I2)
        quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
        
        YQ_OOSC[jt + h, jq] <- c(1, Ztemp[jt,]) %*% coef(b_high) - quantile_E 
        

        
        #------- Quantile regression with GDP only, out-of-sample
        
        YQGDPonly_high_OOSC[(h + 1):(h + test_length), jq] <- Inf
        
        
        bGDPonly_high <- rq(Yh1 ~ ZGDPonly1[,-1], tau= (1-QQ[jq])) #Train on I1
        YQGDPonly_low_OOSC[(h + 1):(h + full_length-test_length), jq] <- as.vector(ZGDPonly2 %*% coef(bGDPonly_high)) #Evaluate on I2
        
        # Initialize a vector for errors
        E_i <- rep(NA, length(Yh2))
        
        # Calculate errors for each point in the test set I2
        for (i in 1:length(E_i)) {
          E_i[i] <- max(YQGDPonly_low_OOSC[h + i, jq] - Yh2[i], Yh2[i] - YQGDPonly_high_OOSC[h + i, jq])
        }
        
        # Compute Q(1-alpha)(E, I2)
        quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
        
        # YQ_low_adj_OOSC[jt + h, jq] <- Z[jt,] %*% coef(b_low) - quantile_E 
        YQGDPonly_OOSC[jt + h, jq] <- ZGDPonly[jt,] %*% coef(bGDPonly_high) - quantile_E 
        

        
        
        #------- Quantile regression with Unconditional QQ, out-of-sample
        
        
        YQunclow_OOSC[(h + 1):(h + test_length), jq] <- - Inf
        
        bunc_high <- rq(Yh1 ~ 1, tau=QQ[jq])
        YQunchigh_OOSC[(h + 1):(h + test_length), jq] <- rep(coef(bunc_high), length((h + 1):(h + test_length)))

        # Initialize a vector for errors
        E_i <- rep(NA, length(Yh2))
        
        # Calculate errors for each point in the test set I2
        for (i in 1:length(E_i)) {
          E_i[i] <- max(YQunclow_OOSC[h + i, jq] - Yh2[i], Yh2[i] - YQunchigh_OOSC[h + i, jq])
        }
        
        # Compute Q(1-??)(E, I2)
        quantile_E <- quantile(E_i, pmin(1, pmax(0, (QQ[jq]) * (1 + 1/length(Yh2)))))
        
        YQunc_OOSC[jt + h, jq] <- coef(bunc_high) + quantile_E
        
      }
      
      PST_OOSC[jt + h, ] <- dst(YY, QQ, rev(YQ_OOSC[jt + h, ]))
      QST_OOSC[jt + h, ] <- qst(QQ, rev(YQ_OOSC[jt + h, ]))
      CST_OOSC[jt + h, ] <- pst(YY, QQ, rev(YQ_OOSC[jt + h, ]))
      ScoreST_OOSC[jt + h] <- dst(YhRealized, QQ, rev(YQ_OOSC[jt + h, ]))
      PitST_OOSC[jt + h] <- pst(YhRealized, QQ, rev(YQ_OOSC[jt + h, ])) # is the probability to observe a value < of YhRealized in this distribution 
      
      PSTGDPonly_OOSC[jt + h, ] <- dst(YY, QQ, rev(YQGDPonly_OOSC[jt + h, ]))
      QSTGDPonly_OOSC[jt + h, ] <- qst(QQ, rev(YQGDPonly_OOSC[jt + h, ]))
      CSTGDPonly_OOSC[jt + h, ] <- pst(YY, QQ, rev(YQGDPonly_OOSC[jt + h, ]))
      ScoreSTGDPonly_OOSC[jt + h] <- dst(YhRealized, QQ, rev(YQGDPonly_OOSC[jt + h, ]))
      PitSTGDPonly_OOSC[jt + h] <- pst(YhRealized, QQ, rev(YQGDPonly_OOSC[jt + h, ])) # is the probability to observe a value < of YhRealized in this distribution 
      
      PSTunc_OOSC[jt + h, ] <- dst(YY, QQ, rev(YQunc_OOSC[jt + h, ]))
      QSTunc_OOSC[jt + h, ] <- qst(QQ, rev(YQunc_OOSC[jt + h, ]))
      CSTunc_OOSC[jt + h, ] <- pst(YY, QQ, rev(YQunc_OOSC[jt + h, ]))
      ScoreSTunc_OOSC[jt + h] <- dst(YhRealized, QQ, rev(YQunc_OOSC[jt + h, ]))
      PitSTunc_OOSC[jt + h] <- pst(YhRealized, QQ, rev(YQunc_OOSC[jt + h, ])) 
      
      # Compute out-of-sample entropy for the full model 
      Temp <- PST_OOSC[jt + h, ] * (YY < QST_OOSC[jt + h, jq50])

      non_zero_indexes <- (PSTunc_OOSC[jt + h, ] != 0) & (PST_OOSC[jt + h, ] != 0)
      
      PSTunc_OOSC_non_zero <- PSTunc_OOSC[jt + h, ][non_zero_indexes]
      PST_OOSC_non_zero <- PST_OOSC[jt + h, ][non_zero_indexes]
      Temp_non_zero <- Temp[non_zero_indexes]
      
      LeftEntropy_OOSC[jt + h] <- -sum((log(PSTunc_OOSC_non_zero) - log(PST_OOSC_non_zero)) * Temp_non_zero * deltaYY)
    }
    
  } 
  
}
 
 
     #  filename <- paste("ResOOSCOdivspacchettoPCA_H", h, ".RData",sep="")
     # cat(paste("Saving results to file", filename, "\n"))
     # 
     # # Save all the variables to the .RData file
     # save(
     #   YQ_ISC,      YQ_OOSC,      YQGDPonly_ISC,      YQGDPonly_OOSC,      YQunc_ISC,      YQunc_OOSC,
     #   PST_ISC,     PST_OOSC,     PSTGDPonly_ISC,     PSTGDPonly_OOSC,     PSTunc_ISC,     PSTunc_OOSC,
     #    QST_ISC,     QST_OOSC,     QSTGDPonly_ISC,     QSTGDPonly_OOSC,     QSTunc_ISC,     QSTunc_OOSC,
     #     CST_ISC,     CST_OOSC,     CSTGDPonly_ISC,     CSTGDPonly_OOSC,     CSTunc_ISC,     CSTunc_OOSC,
     #     STpar_ISC,   STpar_OOSC,   STparGDPonly_ISC,   STparGDPonly_OOSC,   STparunc_ISC,   STparunc_OOSC,
     #     ScoreST_ISC, ScoreST_OOSC, ScoreSTGDPonly_ISC, ScoreSTGDPonly_OOSC, ScoreSTunc_ISC, ScoreSTunc_OOSC,
     #    PitST_ISC,   PitST_OOSC,   PitSTGDPonly_ISC,   PitSTGDPonly_OOSC,   PitSTunc_ISC,   PitSTunc_OOSC,
     #     LeftEntropy_ISC, LeftEntropy_OOSC, 
     #     file=filename
     #   )
     # 

  
# #-------------------------------- per caricare dati salvati -----------------------
# 
# filename <- paste("ResOOS_H", h, "_50-50.RData", sep="")
# cat(paste("Loading results from file", filename, "\n"))
# 
# load(filename)
# 
# #-----------------------------------------------------------------------------------



PITtest_env <- new.env()
source("PITtest.r",local = PITtest_env)

rstestboot_env <- new.env()
source("rstestboot.r",local = rstestboot_env)







# Figure 10. Out-of-sample Predictions.             

# (a)/(b) QQ
par(mar = c(3, 3, 2, 1))  # Adjust the values as needed (bottom, left, top, right)

plot(Time, YQ_OOSC[, jq05], type = 'l', col = 'blue', xlab = 'Time', ylab = 'QQ', xlim = range(Time),ylim = c(-20,20))
lines(Time, YQ_OOSC[, jq50], type = 'l', col = 'blue', lty = 2)
lines(Time, YQ_OOSC[,jq95], type = 'l', col = 'blue', lty = 3)
lines(Time, Yh, type = 'l', col = 'red', lty = 1)
# lines(Time, YQ_ISC[, jq05], type = 'l', col = 'black', lty = 1)
# lines(Time, YQ_ISC[, jq50], type = 'l', col = 'black', lty = 2)
# lines(Time, YQ_ISC[, jq95], type = 'l', col = 'black', lty = 3)





# (c)/(d) Downside Entropy        
plot(Time, LeftEntropy_OOSC, type = 'l', col = 'blue', 
     xlab = 'Time', ylab = 'Entropy', xlim = range(Time))
lines(Time, LeftEntropy_ISC, type = 'l', col = 'black', lty = 2)



# Figure 11. Out-of-sample Accuracy.   
# (a)/(b) Predictive scores
plot(Time, ScoreST_OOSC, type = 'l', col = 'blue', xlab = 'Time', ylab = 'Scores')
lines(Time, ScoreSTGDPonly_OOSC, type = 'l', col = 'black', lty = 2)
legend('topleft', legend = c('GDP and NFCI', 'GDP only'))

# h = 4:> mean(ScoreST_OOSC - ScoreST_OOSCpaper, na.rm = TRUE) == 0.004654713. Predictive scores of my model are better. 
# h = 4:> sum(ScoreST_OOSC - ScoreST_OOSCpaper, na.rm = TRUE) == 0.4096147. Confirmed by using sum instead of mean
# h = 4:> mean(ScoreSTGDPonly_OOSC - ScoreSTGDPonly_OOSC1, na.rm = TRUE) == 0.04962592. Also the model with only GDP performs better
# h = 4:> sum(ScoreSTGDPonly_OOSC - ScoreSTGDPonly_OOSC1, na.rm = TRUE) == 4.367081. Confirmed by using sum instead of mean

# h = 1:> mean(ScoreST_OOSC - ScoreST_OOSC1, na.rm = TRUE) == 0.01211243 Predictive scores of my model are better. 
# h = 1:> sum(ScoreST_OOSC - ScoreST_OOSC1, na.rm = TRUE) == 1.102231 Confirmed by using sum instead of mean
# h = 1:> mean(ScoreSTGDPonly_OOSC - ScoreSTGDPonly_OOSC1, na.rm = TRUE) == 0.02495498 Also the model with only GDP performs better
# h = 1:> sum(ScoreSTGDPonly_OOSC - ScoreSTGDPonly_OOSC1, na.rm = TRUE) #== 2.270903 Confirmed by using sum instead of mean
plot(ecdf(ScoreSTGDPonly_OOSC), main="PS ECDF Comparison GDPonly", col="blue", lty=1, lwd=2)
lines(ecdf(ScoreSTGDPonly_OOSC1), col="red", lty=2, lwd=2)
legend("bottomright", legend=c("Mio", "Paper"), col=c("blue", "red"), lty=c(1, 2), lwd=2)





# (c)/(d): PITs
# The code below was modified from files provided by Barbara Rossi and
# Tatevik Sekhposyan implementing the specification tests for predictive
# densities described in Rossi and Sekhposyan (2017).
rvec <- seq(0.001, 1, by = 0.001)
zST_ecdf <- PITtest_env$PITtest(PitST_OOSC, rvec)
CnSS_model <- 1 - mean((zST_ecdf - seq(0, 1, length.out = length(zST_ecdf)))^2) / var(zST_ecdf) # h = 4: 0.67701012121 || h = 1: 0.9366  is the value. 1 is perfect calibrated
mean_squared_diff <- mean((zST_ecdf-rvec)^2) # h = 4: 0.02338473 || h = 1: 0.005109203

zSTGDPonly_ecdf <- PITtest_env$PITtest(PitSTGDPonly_OOSC, rvec)
CnSS_model_GDPonly <- 1 - mean((zSTGDPonly_ecdf - seq(0, 1, length.out = length(zSTGDPonly_ecdf)))^2) / var(zSTGDPonly_ecdf) # h = 4: 0.872747225 || h = 1: 0.95697 is the value. 1 is perfect calibrated
mean_squared_diff_GDPonly <- mean((zSTGDPonly_ecdf-rvec)^2) # h = 4: 0.009966178 || h = 1: 0.003439986

if (h == 1) {
  # Use asymptotic 5% critical value from Rossi and Sekhposyan (2017): 1.34
  kappa <- 1.34
  kappaGDPonly <- 1.34
  
} else if (h == 4) {
  # Compute bootstrapped 5% critical values
  PITs <- cbind(PitST_OOSC, PitSTGDPonly_OOSC)
  PITs <- PITs[(jtFirstOOS + h):nrow(PITs), , drop = FALSE]
  
  #testcritvalues <- matrix(NA, nrow = 1, ncol = 2, dimnames = list(NULL, c('GDP and NFCI', 'GDP only')))
  testcritvalues <- array(NA, dim = c(2, 3, 2))
  
  for (i in 1:2) {
    testcritvalues[,, i] <- round(rstestboot_env$rstestboot(PITs[, i])$critvalues[2] * 100) / 100
  }
  
  kappa <- testcritvalues[1, 2, 1] #different from Matlab due to seed in CVfinalbootstrapInoue
  kappaGDPonly <- testcritvalues[1, 2, 2] #different from Matlab due to seed in CVfinalbootstrapInoue
}

# Plot PIT for full quantile regression vs. quantile regression with GDP only
plot(rvec, zST_ecdf, type = 'l', col = 'blue', xlab = '??', ylab = 'Empirical CDF')
lines(rvec, zSTGDPonly_ecdf, type = 'l', col = 'red')
P <- sum(!is.na(PitST_OOSC)) #correct for both h = 1 and h = 4

#lines(rvec, rvec - (kappa / sqrt(P)), col = 'black',lty=2)
#lines(rvec, rvec + (kappa / sqrt(P)), col = 'black',lty=2)
lines(rvec,rvec , col = 'black',lty=2)

legend('bottomright', legend = c('GDP and NFCI', 'GDP only', 'Theoretical and 5% Critical Values'), cex = 0.5,fill = c('blue', 'red', 'black'))

#For h = 4, this PIT plot is different due to usage of a seed in CVfinalbootstrapInoue, even using the same as in Matlab, different results are obtained

