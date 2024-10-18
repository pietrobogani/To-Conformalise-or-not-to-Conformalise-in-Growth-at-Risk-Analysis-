PITtest <- function(PIT, rvec){
  z <- PIT[!is.na(PIT)]
  P <- length(z)
  cumcumz <- matrix(nrow = P, ncol = length(rvec))
  
  for(r in 1:length(rvec)){
    cumcumz[,r] <- (z < rvec[r]) - rvec[r]
  }
  
  v <- colSums(cumcumz) / sqrt(P)
  
  Qv <- v^2
  Qvabs <- abs(v)
  
  # These are calculated but not returned in the given MATLAB code
  # KvSTnaive <- max(Qvabs)
  # CVMvSTnaive <- mean(Qv)
  
  z_ecdf <- numeric(length(rvec))
  for(r in 1:length(rvec)){
    z_ecdf[r] <- mean(z < rvec[r])
  }
  
  return(z_ecdf)
}

# Usage:
# rvec <- seq(0, 1, 0.001) # for example
# PITvalues <- rnorm(100) # replace with your PIT values
# result <- PITtest(PITvalues, rvec)



