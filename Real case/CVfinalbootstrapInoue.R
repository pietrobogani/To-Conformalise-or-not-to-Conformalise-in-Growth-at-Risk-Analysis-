CVfinalbootstrapInoue <- function(el, bootMC, pit, rvec) {
  set.seed(500) # Combines setting seeds for randn and rand
  
  KSv <- numeric(bootMC)
  CVMv <- numeric(bootMC)
  
  P <- length(pit)
  size_rvec <- length(rvec)
  
  for(bootrep in 1:bootMC) {
    z <- rnorm(n = P - el + 1) / sqrt(el)
    # if (bootrep == 1){
    #     cat("z:", z, "\n")            #Matlab with the same seed generates different numbers :(((
    # }
    # 
    CVfinalbootstrapInoue
    emp_cdf <- outer(pit, rvec, `<=`)
    
    K_star <- numeric(size_rvec)
    for(j in 1:(P - el + 1)) {
      K_star <- K_star + (z[j] / sqrt(P)) * colSums(emp_cdf[j:(j + el - 1), ] - matrix(rvec, nrow = el, ncol = size_rvec, byrow = TRUE))
    }
    
    KSv[bootrep] <- max(abs(K_star))
    CVMv[bootrep] <- mean(K_star^2)
  }

  KSv <- sort(KSv)
  CVMv <- sort(CVMv)
  cvKv1 <- KSv[round(bootMC * 0.90)]
  cvKv2 <- KSv[round(bootMC * 0.95)]
  cvKv3 <- KSv[round(bootMC * 0.99)]
  cvMv1 <- CVMv[round(bootMC * 0.90)]
  cvMv2 <- CVMv[round(bootMC * 0.95)]
  cvMv3 <- CVMv[round(bootMC * 0.99)]
  
  result <- matrix(c(cvKv1, cvKv2, cvKv3, cvMv1, cvMv2, cvMv3), nrow = 2, byrow = TRUE)
  
  return(result)
}

# Usage (as an example):
# el <- 10
# bootMC <- 100
# pit <- runif(1000)
# rvec <- seq(0, 1, 0.01)
# res <- CVfinalbootstrapInoue(el, bootMC, pit, rvec)
