#----- Adapted CQR algorithm for Quantile Estimation

qCQR <- function(x0,y0,x1,y1,x_test,QQ) {
  
  #This is the function to apply CQR QR algorithm computing a Prediction interval of the form: (-Inf, A]
  
  # x0: A matrix with covariates for the training, and a first column full of 1 (I1)
  # y0: A vector of the target variable for the training (I1)
  # x1: A matrix with covariates for the calibration, and a first column full of 1 (I2)
  # y1: A vector of the target variable for the calibration (I2)
  # x_test: A matrix with covariates for the evaluation
  # QQ: A vector of quantile levels (e.g., 0.1, 0.5, 0.9) for which quantile predictions will be made (1-alpha)
  
  # CQuant_OOS: the output is a matrix with quantile predictions, one row for each test point and one column for each 
  #             quantile level
  
  Q_low <- matrix(NA, nrow(x1), length(QQ))
  Q_high <- matrix(NA, nrow(x1), length(QQ))
  CQuant_OOS <- matrix(0, nrow(x_test), length(QQ))
  
  for( jq in 1:length(QQ) ) {
    
    Q_low[, jq] <- -Inf 
    QR <- rq(formula, data = cbind(x0,Y = y0), tau=QQ[jq])    
    
    Q_high[1 : nrow(x1), jq] <- as.matrix(x1) %*% coef(QR) 
    
    # Initialize a vector for errors
    E_i <- rep(NA, nrow(x1))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(Q_low[i, jq] - y1[i], y1[i] - Q_high[i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
    quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(x1)))
    
    CQuant_OOS[,jq] <- as.matrix(x_test) %*% coef(QR) + quantile_E
    
  }
  
  return(CQuant_OOS)
  
}

qCQRF <- function(x0,y0,x1,y1,x_test,QQ) {
  
  #This is the function to apply CQR QRF algorithm computing a Prediction interval of the form: (-Inf, A]
  
  # x0: A matrix with covariates for the training, and a first column full of 1 (I1)
  # y0: A vector of the target variable for the training (I1)
  # x1: A matrix with covariates for the calibration, and a first column full of 1 (I2)
  # y1: A vector of the target variable for the calibration (I2)
  # x_test: A matrix with covariates for the evaluation
  # QQ: A vector of quantile levels (e.g., 0.1, 0.5, 0.9) for which quantile predictions will be made (1-alpha)
  
  # CQuant_OOS: the output is a matrix with quantile predictions, one row for each test point and one column for each 
  #             quantile level
  
  qrf_model <- quantregForest(x = x0, y = y0)
  CQuant_OOS <- matrix(0, nrow(x_test), length(QQ))
  Q_low <- matrix(NA, nrow(x1), length(QQ))
  Q_high <- matrix(NA, nrow(x1), length(QQ))
  
  for( jq in 1:length(QQ) ) {
    
    Q_low[1:nrow(x1), jq] <- -Inf 
    
    Q_high[1 : nrow(x1), jq] <- predict(qrf_model, newdata = x1, what = QQ[jq])
    
    # Initialize a vector for errors
    E_i <- rep(NA, nrow(x1))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(Q_low[i, jq] - y1[i], y1[i] - Q_high[i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
    quantile_E <- quantile(E_i, QQ[jq] * (1 + 1/nrow(x1)))
    
    CQuant_OOS[,jq] <- predict(qrf_model, newdata = x_test, what = QQ[jq]) + quantile_E
  }
  
  return(CQuant_OOS)
  
}

qCQR_opposite <- function(x0,y0,x1,y1,x_test,QQ) {
  
  #This is the function to apply CQR QR algorithm computing a Prediction interval of the form: [B, + Inf)
  
  # x0: A matrix with covariates for the training, and a first column full of 1 (I1)
  # y0: A vector of the target variable for the training (I1)
  # x1: A matrix with covariates for the calibration, and a first column full of 1 (I2)
  # y1: A vector of the target variable for the calibration (I2)
  # x_test: A matrix with covariates for the evaluation
  # QQ: A vector of quantile levels (e.g., 0.1, 0.5, 0.9) for which quantile predictions will be made (alpha)
  

  # CQuant_OOS: the output is a matrix with quantile predictions, one row for each test point and one column for each 
  #             quantile level
  
  
  Q_low <- matrix(NA, nrow(x1), length(QQ))
  Q_high <- matrix(NA, nrow(x1), length(QQ))
  CQuant_OOS <- matrix(0, nrow(x_test), length(QQ))
  
  for( jq in 1:length(QQ) ) {
    
    Q_high[, jq] <- Inf 
    QR <- rq(y0 ~ as.matrix(x0[,-1]), tau=(QQ[jq]))
    
    Q_low[1 : nrow(x1), jq] <- as.matrix(x1) %*% coef(QR) 
    
    # Initialize a vector for errors
    E_i <- rep(NA, nrow(x1))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(Q_low[i, jq] - y1[i], y1[i] - Q_high[i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
    quantile_E <- quantile(E_i, (1-QQ[jq]) * (1 + 1/nrow(x1)))
    
    CQuant_OOS[,jq] <- as.matrix(x_test) %*% coef(QR) - quantile_E
    
    
  }
  
  return(CQuant_OOS)
  
}

qCQRF_opposite <- function(x0,y0,x1,y1,x_test,QQ) {
  
  #This is the function to apply CQR QRF algorithm computing a Prediction interval of the form: [B, + Inf)
  
  # x0: A matrix with covariates for the training, and a first column full of 1 (I1)
  # y0: A vector of the target variable for the training (I1)
  # x1: A matrix with covariates for the calibration, and a first column full of 1 (I2)
  # y1: A vector of the target variable for the calibration (I2)
  # x_test: A matrix with covariates for the evaluation
  # QQ: A vector of quantile levels (e.g., 0.1, 0.5, 0.9) for which quantile predictions will be made (alpha)
  
  # CQuant_OOS: the output is a matrix with quantile predictions, one row for each test point and one column for each 
  #             quantile level
  
  x0 <- as.matrix(x0)
  x1 <- as.matrix(x1)
  x_test <- as.matrix(x_test)
  
  qrf_model <- quantregForest(x = x0, y = y0)
  CQuant_OOS <- matrix(0, nrow(x_test), length(QQ))
  Q_low <- matrix(NA, nrow(x1), length(QQ))
  Q_high <- matrix(NA, nrow(x1), length(QQ))
  
  for( jq in 1:length(QQ) ) {
    
    Q_low[1:nrow(x1), jq] <- predict(qrf_model, newdata = as.matrix(x1), what = (QQ[jq]))
    
    Q_high[1 : nrow(x1), jq] <- Inf
    
    # Initialize a vector for errors
    E_i <- rep(NA, nrow(x1))
    
    # Calculate errors for each point in the test set I2
    for (i in 1:length(E_i)) {
      E_i[i] <- max(Q_low[i, jq] - y1[i], y1[i] - Q_high[i, jq])
    }
    
    # Compute Q(QQ[jq])(E, I2) N.B 1 - alpha = QQ[jq]
    quantile_E <- quantile(E_i, (1-QQ[jq]) * (1 + 1/nrow(x1)))
    
    CQuant_OOS[,jq] <- predict(qrf_model, newdata = x_test, what = (QQ[jq])) - quantile_E
  }
  
  return(CQuant_OOS)
  
}


#----- General auxiliar functions

wilson_score_interval <- function(x, n, conf.level = 0.95) {
  # x: number of successes
  # n: number of trials
  # conf.level: confidence level (e.g., 0.95 for 95% confidence interval)
  
  # Calculate point estimate for proportion
  p_hat <- x / n
  
  # Find the Z value for the specified confidence level
  z <- qnorm(1 - (1 - conf.level) / 2)
  
  # Wilson score interval formula
  factor <- z^2 / (2 * n)
  denominator <- 1 + z^2 / n
  center_adjustment <- z * sqrt(p_hat * (1 - p_hat) / n + z^2 / (4 * n^2))
  
  lower_bound <- (p_hat + factor - center_adjustment) / denominator
  upper_bound <- (p_hat + factor + center_adjustment) / denominator
  
  # Adjust bounds to ensure they remain within [0, 1]
  lower_bound <- max(0, lower_bound)
  upper_bound <- min(1, upper_bound)
  
  return(c(lower = lower_bound, upper = upper_bound))
}

is_within_ci <- function(success_rate, ci_lower, ci_upper) {
  return(success_rate >= ci_lower & success_rate <= ci_upper)
}

make_formula <- function(num_lags, num_exog) {
  # Initialize an empty vector to store lag terms
  lag_terms <- c()
  
  # Add lag terms based on the value of num_lags
  if (num_lags >= 1) {
    lag_terms <- c(lag_terms, "Y_lag1")
  }
  if (num_lags >= 2) {
    lag_terms <- c(lag_terms, "Y_lag2")
  }
  if (num_lags >= 3) {
    lag_terms <- c(lag_terms, "Y_lag3")
  }
  
  # Create the lag terms part of the formula
  lag_terms <- paste(lag_terms, collapse = " + ")
  
  # Create the exogenous terms part of the formula only if num_exog is greater than 0
  if (num_exog > 0) {
    exog_terms <- paste(paste0("X", 1:num_exog), collapse = " + ")
    # Combine the parts to form the full formula
    formula <- as.formula(paste("Y ~", lag_terms, "+", exog_terms))
  } else {
    # If there are no exogenous terms, use only the lag terms
    formula <- as.formula(paste("Y ~", lag_terms))
  }
  
  return(formula)
}

compute_coverage <- function(Y_test, Quant_OOS, QQ) {
  

  coverage <- rep(NA,length(QQ))
  
  for (jq in 1:length(QQ)) {
    PitOOS <- rep(NA, length(Y_test))
    
    for (i in 1:length(Y_test)){
      if (Y_test[i] <= Quant_OOS[i,jq]) {
        PitOOS[i]  <- 1 
      }
      else 
        PitOOS[i]  <- 0
    }
    
  coverage[jq] <- sum(PitOOS) / length(PitOOS)
  }
  
  return(coverage)
    
}

compute_average_coverage <- function(coveragetot, QQ) {
  
  average_coverage <- rep(NA,length(QQ))
  for (jq in 1:length(QQ)) { 
    average_coverage[jq] <- mean(coveragetot[,jq])
  }
  return(average_coverage)
  
}

compute_results <- function(average_coverage, n, QQ){
  
  resultsPit <- data.frame(Quantile = numeric(), EmpiricalCoverage = numeric(), CI_Lower = numeric(), CI_Upper = numeric())
  
  
  for (jq in 1:length(QQ)) {
    
    success_rate <- average_coverage[jq]
    successes <- success_rate*n
    # Calculate Wilson score interval
    ci <- wilson_score_interval(successes, n)
    # Add to results data frame
    resultsPit <- rbind(resultsPit, data.frame(Quantile = QQ[jq], EmpiricalCoverage = success_rate, CI_Lower = ci[1], CI_Upper = ci[2]))
  
  }
  
  # Apply the function to each row and add the result as a new column
  resultsPit$IsWithinCI <- mapply(is_within_ci, resultsPit$Quantile, resultsPit$CI_Lower, resultsPit$CI_Upper)
  
  return (resultsPit)
}


#Da qui in giù ancora sistemare


pst <- function(X, QQ, qqtarg) { #Return the cumulated density, assuming it as a step funciton
  
  sapply(X, function(x) {
    if (length(qqtarg[qqtarg <= x]) == 0) {
      return(0) # If x is less than the smallest quantile value, the cumulative probability is 0
    }
    
    if (x > max(qqtarg, na.rm = TRUE)) {
      return(1) # If x is greater than all qqtarg, return the maximum cumulative probability
    }
    
    max_value_not_exceeding_x <- max(qqtarg[qqtarg <= x], na.rm = TRUE)
    quantile_index <- which(qqtarg == max_value_not_exceeding_x)
    
    
    return(QQ[min(quantile_index)])
  })
}

dst <- function(X, QQ, qqTarg) {
  # Function to find density for a single point
  find_density <- function(x) {
    density <- 0
    n <- length(qqTarg)
    
    for (i in 1:(n - 1)) {
      
      # Function to dynamically find an interval that is not too small
      find_wider_interval <- function(index) {
        left_index <- index
        right_index <- index + 1
        
        while (abs(qqTarg[right_index] - qqTarg[left_index]) < 0.01 && left_index > 1 && right_index < n) {
          if (left_index > 1) {
            left_index <- left_index - 1
          }
          if (right_index < n) {
            right_index <- right_index + 1
          }
        }
        
        return(c(left_index, right_index))
      }
      
      if (x >= qqTarg[i] && x <= qqTarg[i + 1]) {
        if (abs(qqTarg[i + 1] - qqTarg[i]) < 0.01) {
          indices <- find_wider_interval(i)
          density <- 1 / (qqTarg[indices[2]] - qqTarg[indices[1]]) * (QQ[indices[2]] - QQ[indices[1]])
        } else {
          density <- 1 / (qqTarg[i + 1] - qqTarg[i]) * (QQ[i + 1] - QQ[i])
        }
        
        if (density > 1) {
          cat("Interval:", qqTarg[i], "-", qqTarg[i + 1], "Density:", density, "\n") # Just to check if I get insanely high values
        }
        return(density)
      }
      
    }
    # If x is outside the range of given QQ, return 4 as the density
    return(0)
  }
  
  # Apply find_density to each element in X
  densities <- sapply(X, find_density)
  
  return(densities)
}

qst <- function(QQ, qqTarg) { # qst is always called giving in input QQ. Without smoothing, the estimated quantiles are exactly qqTarg!
  if (length(QQ)== length(qqTarg)) {
    return(qqTarg)
  }
  else
    cat("wrong dimensions")
}

remove_excess_columns <- function(Z) {
  # Check the dimensions of the matrix
  num_rows <- nrow(Z)
  num_cols <- ncol(Z)
  
  # Initialize a variable to store removed column indices
  removed_columns <- NULL
  
  # Check if the number of columns is greater than the number of rows
  if (num_cols >= num_rows) {
    # Calculate the difference
    diff <- num_cols - num_rows
    
    # Calculate the variance of each column
    column_variances <- apply(Z, 2, var)
    
    # Find the indices of columns with the least variance
    removed_columns <- order(column_variances)[1:(diff + 2)]
    
    # Remove the identified columns
    Z <- Z[, -removed_columns]
  }
  
  return(list(modified_matrix = Z, removed_columns = removed_columns))
}

remove_highly_correlated <- function(data, threshold = 0.99) {
  cor_matrix <- cor(data)
  highly_correlated <- which(abs(cor_matrix) > threshold & lower.tri(cor_matrix, diag = FALSE), arr.ind = TRUE)
  
  if (nrow(highly_correlated) == 0) {
    # No columns to remove
    cleaned_data <- data
    removed_indices <- integer(0)  # Empty vector for removed indices
  } else {
    # Get unique columns to remove
    cols_to_remove <- unique(highly_correlated[, 2])
    # Create the cleaned matrix
    cleaned_data <- data[, -cols_to_remove]
    removed_indices <- cols_to_remove
  }
  
  return(list(cleaned_data = cleaned_data, removed_indices = removed_indices))
}

calculate_percent_below <- function(PitST_OOSC,levels) {
  # Remove NA values from the vector
  PitST_OOSC <- na.omit(PitST_OOSC)
  
  # Initialize a vector to store the results
  percent_below <- numeric(length(levels))
  
  # Loop through each level and calculate the percentage of elements below that level
  for (i in seq_along(levels)) {
    level <- levels[i]
    percent_below[i] <- mean(PitST_OOSC < level) * 100
  }
  
  return(percent_below)
}