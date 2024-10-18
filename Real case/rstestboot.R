rstestboot <- function(z) {
  #function is correct! CVfinalbootstrapInoue though uses a seed, giving different numbers from Matlab equivalent
  rvec <- seq(0, 1, by = 0.001)
  P <- length(z)
  
  cumcumz <- matrix(0, nrow = P, ncol = length(rvec))
  for(r in 1:length(rvec)) {
    cumcumz[, r] <- (z < rvec[r]) - rvec[r]
  }
  
  v <- colSums(cumcumz) / sqrt(P)
  
  Qv <- v^2
  Qvabs <- abs(v)
  
  Kv <- max(Qvabs)
  CVMv <- mean(Qv)
  
  results <- c(Kv, CVMv)
  
  
  CVfinalbootstrapInoue_env <- new.env()
  source("CVfinalbootstrapInoue.r",local = CVfinalbootstrapInoue_env)
  
  
  # bootstrap the crit values
  el <- 12
  bootMC <- 200
  critvalues <- CVfinalbootstrapInoue_env$CVfinalbootstrapInoue(el, bootMC, z, rvec)
  #cat("z:", z, "\n")
  cat("Critical values:", critvalues, "\n")
  return(list(results = results, critvalues = critvalues))
}

# Assuming you have CVfinalbootstrapInoue function loaded in R (translated previously)
# Usage (as an example):
# z <- runif(1000)
# res <- rstestboot(z)
