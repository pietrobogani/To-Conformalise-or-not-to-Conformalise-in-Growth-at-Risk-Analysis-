QuantilesInterpolation <- function(qqTarg, QQ) {
  library(parallel)
  library(sn)
  
  # Set bounds for optimization
  LB <- c(-20, 0, -30)
  UB <- c(20, 50, 30)

  # Locate target quantiles
  jq50 <- which.min(abs(QQ - 0.50))
  jq25 <- which.min(abs(QQ - 0.25))
  jq75 <- which.min(abs(QQ - 0.75))
  jq05 <- which.min(abs(QQ - 0.05))
  jq95 <- which.min(abs(QQ - 0.95))
  
  # Set initial conditions for optimization (no option to provide them initially)

  iqn <- qnorm(0.75) - qnorm(0.25)
  lc0 <- qqTarg[jq50]
  sc0 <- (qqTarg[jq75] - qqTarg[jq25]) / iqn
  sh0 <- 0
  
  
  Select <- c(jq05, jq25, jq75, jq95)
  
  # Prepare for parallel execution
  cl <- makeCluster(detectCores() - 1) # Use one less than the total number of cores
  #clusterExport(cl, c("qqTarg", "QQ", "qnorm", "nlminb", "qst","lc0","sc0", "sh0","jq05", "jq25", "jq75", "jq95"))
  
  clusterEvalQ(cl, { 
    library(sn)
  }) 
  
  # Function to optimize in parallel for each degree of freedom
  optimize_fn_nlminb <- function(df,lc0, sc0, sh0, qqTarg, QQ, jq05, jq25, jq75, jq95) {
    objective_function <- function(p) {
      sum((qqTarg[c(jq05, jq25, jq75, jq95)] - qst(QQ[c(jq05, jq25, jq75, jq95)], xi=p[1], omega=p[2], alpha=p[3], nu=df))^2)
    }
    result <- nlminb(start = c(lc0, sc0, sh0), objective_function, lower = c(-20, 0, -30), upper = c(20, 50, 30))
    return(list(par = result$par, value = result$objective))
  }
  
  # Execute in parallel
  results <- parLapply(cl, 1:30, optimize_fn_nlminb,lc0, sc0, sh0, qqTarg, QQ, jq05, jq25, jq75, jq95)
  stopCluster(cl)
  
  # Extract results
  par <- matrix(unlist(lapply(results, function(x) x$par)), ncol = 3, byrow = TRUE)
  ssq <- sapply(results, function(x) x$value)
  
  # Find the best fit
  best_df <- which.min(ssq)
  
  list(lc = par[best_df, 1], sc = par[best_df, 2], sh = par[best_df, 3], df = best_df, ssq = ssq[best_df])
}
