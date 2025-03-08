# FUNCTIONS


# ===========================================================
# Functions from Nolde and Ziegel (2017) mainly for comparison
# ===========================================================

# Variance of a skewed t random variable with location zero and scale one
mean.st <- function(shape=4, skew=1)
{
  nu=shape; ga=skew
  k=gamma((nu+1)/2)/sqrt(pi*nu)/gamma(nu/2)
  m=2*k*nu/(nu-1)*(1-2/(1+ga^2))*(ga+1/ga)
  return(m)
}

var.st <- function(shape=4, skew=1)
{
  nu=shape; ga=skew
  
  # m=mean.st(shape=nu, skew=ga)
  # v=(ga+1/ga)^2*nu/(nu-2)*(1-3*ga^2/(1+ga^2)^2) - m^2
  
  k=gamma((nu+1)/2)/sqrt(pi*nu)/gamma(nu/2)
  v=nu/(nu-2)*(1-3*ga^2/(1+ga^2)^2)-4*k^2*(nu/(nu-1))^2*(1-2/(1+ga^2))^2
  v=v*(ga+1/ga)^2
  return(v)
}

# ===========================================================
# RISK MEASURE FORECASTS
# Risk measures considered: VaR_alpha, expectile_tau and (VaR_nu,Expected Shortfall_nu)
# Forecasting methods: 
# (1) model-based forecasts (using a functional of the assumed error distribution)
# (2) Filtered historical simulation:
# Nonparametric, sampling with replacement is used (uses "seed" for reproducibility)
# returns 1-period ahead risk measure forecasts based on the empirical distribution of the simulated returns
# (3) EVT (extreme value theory based forecasts)

# INPUTS:
# - fit: ugarchfit object "fit" (rugarch package) (after fitting negated log-returns = "log-losses")
# - alpha: VaR level, can be vectorized
# - tau: expectile level, >1/2 (upper tail), can be vectorized
# - nu: level for pair (VaR_nu, ES_nu), can be vectorized
# - qmodel: model-based quantiles corresponding to alpha level(s)
# - emodel: model-based expectile(s) corresponding to tau level(s)
# - esmodel: model-based expected shortfall corresponding to nu level(s)
# - seed: used for non-parameteric sampling for FHS
# - n: sample size used in FHS

# OUTPUTS: 
# - list of 1-period-ahead risk measure forecasts of negated log-return series under the considered estimation methods
# - Parameter estimates when fitting the GPD distribution to threshold excesses
# - Forecasts of conditional mean and volatility mu[t+1] and sigma[t+1]
# ===========================================================

RM.forecasts3 <- function(fit,qmodel,esmodel)
{
  frcst = ugarchforecast(fit,n.ahead=1)
  sigt1 = sigma(frcst)
  mut1 = fitted(frcst)
  
  VaRmodel=mut1 + sigt1*qmodel
  es.model= mut1 + sigt1*esmodel
  
  return(list(VaRmodel=VaRmodel,ESmodel=es.model, mut=mut1, sigt=sigt1))	
}



# ===========================================================
# CONDITIONAL CALIBRATION TESTS (CCT)
# ===========================================================

# ===================================================================
# ONE-SIDED TEST of super-calibration
# test statistic tau2 for average tests and tau4 for conditional tests
# ===================================================================

cct.1s.VaR <- function(x, r, lev=0.95) # super-calbration
{
  # browser()
  n = length(x) # out-of-sample size
  sf1 = .5
  sf2 = .2 
  sf3 = .1 # significance levels
  out.pv.cond = c() # output p-values
  
  # setting up appropriate identification functions for the risk measure
  v = (r>=x) - lev
  
  # conditional calibration tests	
  
  ## same test function as for the two-sided test
  i <- 1
  n.rej1 <- 0
  n.rej2 <- 0
  n.rej3 <- 0 # dates of rejection
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  pv.cond <- 1 # initial p-value
  N=n; vt=v
  while(i <= (N)){
    if(pv.cond >= sf1 && d1 == FALSE){
      n.rej1 <- n.rej1+1
    }
    if(pv.cond < sf1){
      d1 = TRUE
    }
    if(pv.cond >= sf2 && d2 == FALSE){
      n.rej2 <- n.rej2+1
    }
    if(pv.cond < sf2){
      d2 = TRUE
    }
    if(pv.cond >= sf3 && d3 == FALSE){
      n.rej3 <- n.rej3+1
    }
    if(pv.cond < sf3){
      d3 = TRUE
    }
    
    h=rbind(rep(1,i),abs(r[1:i]))
    
    q=dim(h)[1]
    
    zbar=(h %*% (vt[1:i]))/i # q x 1
    Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
    for (t in 1:i)
    {		
      Omega = Omega + vt[t]^2*(h[,t]%*%t(h[,t]))	# V is a scalar for k=1
    }
    Omega = Omega/i
    
    # Hommel's (1983) procedure
    tn <- sqrt(i) * diag(Omega)^(-1/2) * zbar 	#test statistic
    pi <- sort(pnorm(tn))
    cq <- sum(1/(1:q))
    pv.cond <- min(q*cq*min(pi/(1:q)),1)
    
    # pv.condB <- min(q*pi[1],1) # q*min(p_i) Bonferroni multiple test correction
    
    out.pv.cond <- cbind(out.pv.cond, pv.cond)
    i = i+1
  }
  
  return(list(out.pv.cond = out.pv.cond, n.rej1=n.rej1, n.rej2=n.rej2, n.rej3=n.rej3))
}


# One-sided Conditional Calibration Tests for k-variate risk measure with k>1
# The code is specific for (VaR,ES) pair
# Inputs:
# x: verifying observations
# r: risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# lev: risk measure level ("nu" for pair (VaR,ES))
# htype: type of the test function matrix h 

# Test functions
# 1: ht=diag(2) # identity matrix for the unconditional test
# 2: ht = cbind(c(1, abs(V1[t-1]), r1[t],0,0),c(0,0,0,1,sigt^(-1)) )
# 3: ht = cbind(c(1, abs(V1[t-1]), r1[t],0),c(0,0,0,sigt^(-1)) )

# Output:
# pv: p-value (minimum of the adj. p-values) for the cct test based on the Hommel's (1983) procedure
#     corresponding to one of the above mentioned choices of the test functions

cct.onesided2 <- function(x, r, lev=0.95, sigt) # sub-calibration
{
  # browser()
  n = length(x) # out-of-sample size
  sf1 = .5
  sf2 = .2 
  sf3 = .1 # significance levels
  out.pv.cond <- c() # output p-values
  
  # setting up appropriate identification functions for the risk measure
  ind = (x>r[,1])
  v1 = 1 - lev - ind
  v2 = r[,2] - r[,1] + 1/(1-lev)*ind*(r[,1]-x)
  v = rbind(v1,v2) # identification function matrix (k x n)
  
  
  # conditional calibration tests	
  # ht = cbind(c(1, abs(r1[t]),0,0),c(0,0,1,sigt^(-1)) )
  i <- 1
  n.rej1 <- 0
  n.rej2 <- 0
  n.rej3 <- 0 # dates of rejection
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  pv.cond <- 1 # initial p-value
  q=4; N=n
  while(i <= (N)){
    if(pv.cond >= sf1 && d1 == FALSE){
      n.rej1 <- n.rej1+1
    }
    if(pv.cond < sf1){
      d1 = TRUE
    }
    if(pv.cond >= sf2 && d2 == FALSE){
      n.rej2 <- n.rej2+1
    }
    if(pv.cond < sf2){
      d2 = TRUE
    }
    if(pv.cond >= sf3 && d3 == FALSE){
      n.rej3 <- n.rej3+1
    }
    if(pv.cond < sf3){
      d3 = TRUE
    }
    
    zsum = matrix(0, nrow=q,ncol=1) # storage variable
    Omega = matrix(0, nrow=q,ncol=q)   # covariance matrix
    
    for (t in 1:i)
    {
      ht = cbind(c(1, abs(r[t,1]),0,0),c(0,0,1,sigt[t]^(-1)) ) # q x k with q=q1+q2
      zt = ht %*% v[,t]	# q x 1 
      zsum = zsum + zt
      Omega = Omega + (zt %*% t(zt))	# q x q	
    }
    zbar = zsum/i
    Omega = Omega/i
    
    # Hommel's (1983) procedure
    
    tn <- sqrt(i) * diag(Omega)^(-1/2) * zbar 	#test statistic
    pi <- sort(pnorm(tn))
    cq <- sum(1/(1:q))
    pv.cond <- min(q*cq*min(pi/(1:q)),1)
    # pv.condB <- min(q*pi[1],1) # q*min(p_i) Bonferroni multiple test correction
    
    out.pv.cond <- cbind(out.pv.cond, pv.cond)
    i = i+1
  }
  
  return(list(out.pv.cond=out.pv.cond, n.rej1=n.rej1, n.rej2=n.rej2, n.rej3=n.rej3))
}


# ===================================================================
# Extracting and combining RM forecasts under different estimation methods for a given confidence level
# Input: index "k" for the RM level in "levels"
# ===================================================================


getVaRES_em <- function(k)
{
  load("Sim3norm3.RDATA")
  VaRout = out$VaR[,(2+k)]
  ESout = out$ES[,k]
  
  load("Sim3std3.RDATA")
  VaRout = rbind(VaRout,out$VaR[,(2+k)])
  ESout = rbind(ESout,out$ES[,k])
  
  load("Sim3sstd3.RDATA")
  VaRout = rbind(VaRout,out$VaR[,(2+k)])
  ESout = rbind(ESout,out$ES[,k])
  
  load("Empirical.RDATA")
  VaRout = rbind(VaRout,out$VaR[,(2+k)])
  ESout = rbind(ESout,out$ES[,k])
  
  return(list(VaR=VaRout,ES=ESout))
}





# ===========================================================
# Functional for e-backtest
# ===========================================================


# ===========================================================
# Get e-values for VaR
# Inputs:
# x: realized losses
# r: risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# z: auxiliary risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# lev: risk measure level ("nu" for pair (VaR,ES))
# w: time window
# ===================================================================

EVaR <- function(x,r,lev)
{
  (x > r) / (1-lev)
}

# E-backtesting VaR detecting structral change and comparison with Hoga and Demetrescu (2022)

evalue_VaR_mv <- function(x, r, lev=0.95, m)
{
  n = length(x) # out-of-sample size
  
  sf1 <- 1/3.2
  sf2 <- .1
  sf3 <- .05 # significance levels
  
  e.out.con <- c() # create a vector to store e-values for constant lambda
  e.out.akelly <- c(1) # GREE
  e.out.rkelly <- c(1) # GREL
  e.out.mix <- c() # GREM
  out.lambda <- c(0) # output parameters tuned for GREE
  out.lambda.rkelly <- c(0) # output parameters tuned for GREL
  
  # test martingale with lambda = lambda.c
  lambda.c <- 0.01
  
  # Initial e-value for constant lambda
  e.con <- 1
  
  # Initial e-value for GREE
  e.akelly <- 1
  
  # Initial e-value for GREL
  e.rkelly <- 1
  
  # Initial e-value for GREM
  e.mix <- 1
  
  # Loop calculating the test martingale for constant lambda
  i <- 1
  n.rej1.con <- 0
  n.rej2.con <- 0
  n.rej3.con <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.con <= (1/sf1) && d1 == FALSE){
      n.rej1.con <- n.rej1.con+1
    }
    if(e.con > (1/sf1)){
      d1 = TRUE
    }
    if(e.con <= (1/sf2) && d2 == FALSE){
      n.rej2.con <- n.rej2.con+1
    }
    if(e.con > (1/sf2)){
      d2 = TRUE
    }
    if(e.con <= (1/sf3) && d3 == FALSE){
      n.rej3.con <- n.rej3.con+1
    }
    if(e.con > (1/sf3)){
      d3 = TRUE
    }
    e.con <- e.con * (1 - lambda.c + lambda.c*EVaR(x=x[i],r=r[i],lev=lev))
    e.out.con <- cbind(e.out.con, e.con)
    i = i+1
  }
  
  # Loop calculating the test martingale for GREE
  i <- 2
  
  e.out.akelly <- 1
  n.rej1.akelly <- 1
  n.rej2.akelly <- 1
  n.rej3.akelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.akelly <= (1/sf1) && d1 == FALSE){
      n.rej1.akelly <- n.rej1.akelly+1
    }
    if(e.akelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.akelly <= (1/sf2) && d2 == FALSE){
      n.rej2.akelly <- n.rej2.akelly+1
    }
    if(e.akelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.akelly <= (1/sf3) && d3 == FALSE){
      n.rej3.akelly <- n.rej3.akelly+1
    }
    if(e.akelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value for GREE
    lambda <- (1 - lev) * (lev - mean((x[(max(1,i-m)):(i-1)]<=r[(max(1,i-m)):(i-1)]))) / (lev ^ 2 + (1 - 2*lev) * mean((x[(max(1,i-m)):(i-1)]<=r[(max(1,i-m)):(i-1)])))
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda <- cbind(out.lambda, lambda)
    
    e.akelly <- e.akelly * (1 - lambda + lambda * ((x[i]>r[i]) / (1-lev)))
    # e.akelly <- prod(1 - out.lambda + out.lambda * ((x[1:i+1]>r[1:i+1]) / (1-lev)), na.rm = FALSE)
    
    e.out.akelly <- cbind(e.out.akelly, e.akelly)
    
    i = i+1
  }
  
  
  # Loop calculating the test martingale for GREL
  i <- 2
  
  e.out.rkelly <- 1
  n.rej1.rkelly <- 1
  n.rej2.rkelly <- 1
  n.rej3.rkelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value for GREL
    lambda <- (1 - lev) * (lev - mean((x[(max(1,i-m)):(i-1)]<=r[i]))) / (lev ^ 2 + (1 - 2*lev) * mean((x[(max(1,i-m)):(i-1)]<=r[i])))
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * ((x[i]>r[i]) / (1-lev)))
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    i = i+1
  }
  
  
  # Loop calculating the test martingale for GREM
  i <- 1
  n.rej1.mix <- 0
  n.rej2.mix <- 0
  n.rej3.mix <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.mix <= (1/sf1) && d1 == FALSE){
      n.rej1.mix <- n.rej1.mix+1
    }
    if(e.mix > (1/sf1)){
      d1 = TRUE
    }
    if(e.mix <= (1/sf2) && d2 == FALSE){
      n.rej2.mix <- n.rej2.mix+1
    }
    if(e.mix > (1/sf2)){
      d2 = TRUE
    }
    if(e.mix <= (1/sf3) && d3 == FALSE){
      n.rej3.mix <- n.rej3.mix+1
    }
    if(e.mix > (1/sf3)){
      d3 = TRUE
    }
    
    
    
    e.mix <- (e.out.akelly[i] + e.out.rkelly[i])/2
    e.out.mix <- cbind(e.out.mix, e.mix)
    
    i = i+1
  }
  
  
  return(list(e.out.con = e.out.con, n.rej1.con = n.rej1.con, n.rej2.con = n.rej2.con, n.rej3.con = n.rej3.con,
              e.out.akelly = e.out.akelly, n.rej1.akelly = n.rej1.akelly, n.rej2.akelly = n.rej2.akelly, n.rej3.akelly = n.rej3.akelly, out.lambda = out.lambda,
              e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly,
              e.out.mix = e.out.mix, n.rej1.mix = n.rej1.mix, n.rej2.mix = n.rej2.mix, n.rej3.mix = n.rej3.mix))
}



# E-testing VaR for IID data and time series data


evalue.VaR <- function(x, r, lev=0.95)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels

  e.out.con <- c() # create a vector to store e-values for constant lambda
  e.out.akelly <- c(1) # GREE
  e.out.rkelly <- c(1) # GREL
  e.out.mix <- c() # GREM
  out.lambda <- c(0) # output parameters tuned for GREE
  out.lambda.rkelly <- c(0) # output parameters tuned for GREL
  
  # test martingale with lambda = lambda.c
  lambda.c <- 0.01 # lambda we use for time series data
  # lambda.c <- 0.025 # We use this lamabda for iid observations in Section 2.1 of "Simulation and data analysis for e-backtesting"
  
  # Initial e-value for constant lambda
  e.con <- 1
  
  # Initial e-value for GREE
  e.akelly <- 1
  
  # Initial e-value for GREL
  e.rkelly <- 1
  
  # Initial e-value for GREM
  e.mix <- 1

  # Loop calculating the test martingale for constant lambda
  i <- 1
  n.rej1.con <- 0
  n.rej2.con <- 0
  n.rej3.con <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.con <= (1/sf1) && d1 == FALSE){
      n.rej1.con <- n.rej1.con+1
    }
    if(e.con > (1/sf1)){
      d1 = TRUE
    }
    if(e.con <= (1/sf2) && d2 == FALSE){
      n.rej2.con <- n.rej2.con+1
    }
    if(e.con > (1/sf2)){
      d2 = TRUE
    }
    if(e.con <= (1/sf3) && d3 == FALSE){
      n.rej3.con <- n.rej3.con+1
    }
    if(e.con > (1/sf3)){
      d3 = TRUE
    }
    e.con <- e.con * (1 - lambda.c + lambda.c*EVaR(x=x[i],r=r[i],lev=lev))
    e.out.con <- cbind(e.out.con, e.con)
    i = i+1
  }

  # Loop calculating the test martingale for GREE
  i <- 2

  e.out.akelly <- 1
  n.rej1.akelly <- 1
  n.rej2.akelly <- 1
  n.rej3.akelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.akelly <= (1/sf1) && d1 == FALSE){
      n.rej1.akelly <- n.rej1.akelly+1
    }
    if(e.akelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.akelly <= (1/sf2) && d2 == FALSE){
      n.rej2.akelly <- n.rej2.akelly+1
    }
    if(e.akelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.akelly <= (1/sf3) && d3 == FALSE){
      n.rej3.akelly <- n.rej3.akelly+1
    }
    if(e.akelly > (1/sf3)){
      d3 = TRUE
    }

    # e-value for GREE
    lambda <- (1 - lev) * (lev - mean((x[1:(i-1)]<=r[1:(i-1)]))) / (lev ^ 2 + (1 - 2*lev) * mean((x[1:(i-1)]<=r[1:(i-1)])))
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda <- cbind(out.lambda, lambda)

    e.akelly <- e.akelly * (1 - lambda + lambda * ((x[i]>r[i]) / (1-lev)))
    # e.akelly <- prod(1 - out.lambda + out.lambda * ((x[1:i+1]>r[1:i+1]) / (1-lev)), na.rm = FALSE)

    e.out.akelly <- cbind(e.out.akelly, e.akelly)

    i = i+1
  }


  # Loop calculating the test martingale for GREL
  i <- 2

  e.out.rkelly <- 1
  n.rej1.rkelly <- 1
  n.rej2.rkelly <- 1
  n.rej3.rkelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }

    # e-value for GREL
    lambda <- (1 - lev) * (lev - mean((x[1:(i-1)]<=r[i]))) / (lev ^ 2 + (1 - 2*lev) * mean((x[1:(i-1)]<=r[i])))
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)


    e.rkelly <- e.rkelly * (1 - lambda + lambda * ((x[i]>r[i]) / (1-lev)))

    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)

    i = i+1
  }

  
  # Loop calculating the test martingale for GREM
  i <- 1
  n.rej1.mix <- 0
  n.rej2.mix <- 0
  n.rej3.mix <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.mix <= (1/sf1) && d1 == FALSE){
      n.rej1.mix <- n.rej1.mix+1
    }
    if(e.mix > (1/sf1)){
      d1 = TRUE
    }
    if(e.mix <= (1/sf2) && d2 == FALSE){
      n.rej2.mix <- n.rej2.mix+1
    }
    if(e.mix > (1/sf2)){
      d2 = TRUE
    }
    if(e.mix <= (1/sf3) && d3 == FALSE){
      n.rej3.mix <- n.rej3.mix+1
    }
    if(e.mix > (1/sf3)){
      d3 = TRUE
    }
    
    
    
    e.mix <- (e.out.akelly[i] + e.out.rkelly[i])/2
    e.out.mix <- cbind(e.out.mix, e.mix)
    
    i = i+1
  }
  
  
  return(list(e.out.con = e.out.con, n.rej1.con = n.rej1.con, n.rej2.con = n.rej2.con, n.rej3.con = n.rej3.con,
              e.out.akelly = e.out.akelly, n.rej1.akelly = n.rej1.akelly, n.rej2.akelly = n.rej2.akelly, n.rej3.akelly = n.rej3.akelly, out.lambda = out.lambda,
              e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly,
              e.out.mix = e.out.mix, n.rej1.mix = n.rej1.mix, n.rej2.mix = n.rej2.mix, n.rej3.mix = n.rej3.mix))
}

# e-test with constant lambda (for Type I error detection)

evalue.VaR.con <- function(x, r, lev=0.95)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.con <- c() # create a vector to store e-values for constant lambda
  
  # test martingale with lambda = lambda.c
  lambda.c <- 0.01 # lambda we use for time series data
  # lambda.c <- 0.025 # We use this lamabda for iid observations in Section 2.1 of "Simulation and data analysis for e-backtesting"
  
  # Initial e-value for constant lambda
  e.con <- 1
  
  # Loop calculating the test martingale for constant lambda
  i <- 1
  n.rej1.con <- 0
  n.rej2.con <- 0
  n.rej3.con <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.con <= (1/sf1) && d1 == FALSE){
      n.rej1.con <- n.rej1.con+1
    }
    if(e.con > (1/sf1)){
      d1 = TRUE
    }
    if(e.con <= (1/sf2) && d2 == FALSE){
      n.rej2.con <- n.rej2.con+1
    }
    if(e.con > (1/sf2)){
      d2 = TRUE
    }
    if(e.con <= (1/sf3) && d3 == FALSE){
      n.rej3.con <- n.rej3.con+1
    }
    if(e.con > (1/sf3)){
      d3 = TRUE
    }
    e.con <- e.con * (1 - lambda.c + lambda.c*EVaR(x=x[i],r=r[i],lev=lev))
    e.out.con <- cbind(e.out.con, e.con)
    i = i+1
  }
  
  
  return(list(e.out.con = e.out.con, n.rej1.con = n.rej1.con, n.rej2.con = n.rej2.con, n.rej3.con = n.rej3.con))
}

# e-test for VaR with known alternative Q normal (GRO method)

evalue.VaRQ <- function(x, r, lev=0.95)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.rkelly <- c() # e-process
  out.lambda.rkelly <- c() # output parameters
  
  # Initial e-value
  e.rkelly <- 1
  
  # Loop calculating the test martingale
  i <- 1
  n.rej1.rkelly <- 0
  n.rej2.rkelly <- 0
  n.rej3.rkelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value
    f = function(pm){
      function(x){
        log(1-pm+pm*(x > r[i])/(1-lev))*dnorm(x)
      }
    }
    fun = function(pm){
      integrate(f(pm), -Inf, Inf, rel.tol = 1e-3)$value
    }
    lambda = optimize(fun, interval = c(0,0.5), maximum = TRUE)$maximum
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda*EVaR(x=x[i],r=r[i],lev=lev))
    
    
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    
    i = i+1
  }
  
  
  return(list(e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly))
}


# ===========================================================
# Get e-values for ES
# Inputs:
# x: realized losses
# r: risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# z: auxiliary risk measure forecasts; n x k vector (n=out-of-sample size; k=number of risk measures)
# lev: risk measure level ("nu" for pair (VaR,ES))
# m: time window
# ===================================================================

EES <- function(x,r,z,lev)
{
  (pmax(x - z, 0)) / ((1 - lev) * (r - z))
}

# E-backtesting ES detecting structral change and comparison with Hoga and Demetrescu (2022)

evalue_ES_mv <- function(x, r, z, lev=0.875, m)
{
  n = length(x) # out-of-sample size
  
  sf1 <- 1/3.2
  sf2 <- .1
  sf3 <- .05 # significance levels
  
  e.out.con <- c() # create a vector to store e-values for constant lambda
  e.out.akelly <- c(1) # GREE
  e.out.rkelly <- c(1) # GREL
  e.out.mix <- c() # GREM
  out.lambda <- c(0) # output parameters tuned for GREE
  out.lambda.rkelly <- c(0) # output parameters tuned for GREM
  
  # test martingale with lambda = lambda.c
  lambda.c <- 0.01
  
  # Initial e-value for constant lambda
  e.con <- 1
  
  # Initial e-value for GREE
  e.akelly <- 1
  
  # Initial e-value for GREL
  e.rkelly <- 1
  
  # Initial e-value for GREM
  e.mix <- 1
  
  # Loop calculating the test martingale for constant lambda
  i <- 1
  n.rej1.con <- 0
  n.rej2.con <- 0
  n.rej3.con <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.con <= (1/sf1) && d1 == FALSE){
      n.rej1.con <- n.rej1.con+1
    }
    if(e.con > (1/sf1)){
      d1 = TRUE
    }
    if(e.con <= (1/sf2) && d2 == FALSE){
      n.rej2.con <- n.rej2.con+1
    }
    if(e.con > (1/sf2)){
      d2 = TRUE
    }
    if(e.con <= (1/sf3) && d3 == FALSE){
      n.rej3.con <- n.rej3.con+1
    }
    if(e.con > (1/sf3)){
      d3 = TRUE
    }
    e.con <- e.con * (1 - lambda.c + lambda.c * EES(x=x[i],r=r[i],z=z[i],lev=lev))
    e.out.con <- cbind(e.out.con, e.con)
    i = i+1
  }
  
  # Loop calculating the test martingale for GREE
  i <- 2
  n.rej1.akelly <- 1
  n.rej2.akelly <- 1
  n.rej3.akelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.akelly <= (1/sf1) && d1 == FALSE){
      n.rej1.akelly <- n.rej1.akelly+1
    }
    if(e.akelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.akelly <= (1/sf2) && d2 == FALSE){
      n.rej2.akelly <- n.rej2.akelly+1
    }
    if(e.akelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.akelly <= (1/sf3) && d3 == FALSE){
      n.rej3.akelly <- n.rej3.akelly+1
    }
    if(e.akelly > (1/sf3)){
      d3 = TRUE
    }
    
  
    lambda <- (mean(EES(x=x[(max(1,i-m)):(i-1)],r=r[(max(1,i-m)):(i-1)],z=z[(max(1,i-m)):(i-1)],lev=lev)) - 1) / 
      (mean((EES(x=x[(max(1,i-m)):(i-1)],r=r[(max(1,i-m)):(i-1)],z=z[(max(1,i-m)):(i-1)],lev=lev)
      - mean(EES(x=x[(max(1,i-m)):(i-1)],r=r[(max(1,i-m)):(i-1)],z=z[(max(1,i-m)):(i-1)],lev=lev)))^2) +
        (mean(EES(x=x[(max(1,i-m)):(i-1)],r=r[(max(1,i-m)):(i-1)],z=z[(max(1,i-m)):(i-1)],lev=lev)) - 1)^2)
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda <- cbind(out.lambda, lambda)
    
    e.akelly <- e.akelly * (1 - lambda + lambda * EES(x=x[i],r=r[i],z=z[i],lev=lev))
    
    e.out.akelly <- cbind(e.out.akelly, e.akelly)
    
    i = i+1
  }
  
  # Loop calculating the test martingale for GREL
  i <- 2
  n.rej1.rkelly <- 1
  n.rej2.rkelly <- 1
  n.rej3.rkelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    
    lambda <- (mean(EES(x=x[(max(1,i-m)):(i-1)],r=r[i],z=z[i],lev=lev)) - 1) / (mean((EES(x=x[(max(1,i-m)):(i-1)],r=r[i],z=z[i],lev=lev)
               - mean(EES(x=x[(max(1,i-m)):(i-1)],r=r[i],z=z[i],lev=lev)))^2) + (mean(EES(x=x[(max(1,i-m)):(i-1)],r=r[i],z=z[i],lev=lev)) - 1)^2)
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * EES(x=x[i],r=r[i],z=z[i],lev=lev))
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    i = i+1
  }
  
  
  # Loop calculating the test martingale for GREM
  i <- 1
  n.rej1.mix <- 0
  n.rej2.mix <- 0
  n.rej3.mix <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.mix <= (1/sf1) && d1 == FALSE){
      n.rej1.mix <- n.rej1.mix+1
    }
    if(e.mix > (1/sf1)){
      d1 = TRUE
    }
    if(e.mix <= (1/sf2) && d2 == FALSE){
      n.rej2.mix <- n.rej2.mix+1
    }
    if(e.mix > (1/sf2)){
      d2 = TRUE
    }
    if(e.mix <= (1/sf3) && d3 == FALSE){
      n.rej3.mix <- n.rej3.mix+1
    }
    if(e.mix > (1/sf3)){
      d3 = TRUE
    }
    
    
    
    e.mix <- (e.out.akelly[i] + e.out.rkelly[i])/2
    e.out.mix <- cbind(e.out.mix, e.mix)
    
    i = i+1
  }
  
  
  return(list(e.out.con = e.out.con, n.rej1.con = n.rej1.con, n.rej2.con = n.rej2.con, n.rej3.con = n.rej3.con,
              e.out.akelly = e.out.akelly, n.rej1.akelly = n.rej1.akelly, n.rej2.akelly = n.rej2.akelly, n.rej3.akelly = n.rej3.akelly,
              e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly,
              e.out.mix = e.out.mix, n.rej1.mix = n.rej1.mix, n.rej2.mix = n.rej2.mix, n.rej3.mix = n.rej3.mix,
              out.lambda = out.lambda, out.lambda.rkelly = out.lambda.rkelly))
}


# E-backtesting ES for financial data


evalue_ES_em <- function(x, r, z, lev=0.875, w)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels

  e.out.con <- c() # create a vector to store e-values for constant lambda
  e.out.akelly <- c() # GREE
  e.out.rkelly <- c() # GREL
  e.out.mix <- c() # GREM
  out.lambda <- c() # output parameters tuned for GREE
  out.lambda.rkelly <- c() # output parameters tuned for GREL
  # out.c.rkelly <- c() # output parameters tuned for GREM

  # test martingale with lambda = lambda.c
  lambda.c <- 0.01

  # Initial e-value for constant lambda
  e.con <- 1

  # Initial e-value for GREE
  e.akelly <- 1

  # Initial e-value for GREL
  e.rkelly <- 1
  
  # Initial e-value for GREM
  e.mix <- 1

  # Loop calculating the test martingale for constant lambda
  i <- 1
  n.rej1.con <- 0
  n.rej2.con <- 0
  n.rej3.con <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n-w)){
    if(e.con <= (1/sf1) && d1 == FALSE){
      n.rej1.con <- n.rej1.con+1
    }
    if(e.con > (1/sf1)){
      d1 = TRUE
    }
    if(e.con <= (1/sf2) && d2 == FALSE){
      n.rej2.con <- n.rej2.con+1
    }
    if(e.con > (1/sf2)){
      d2 = TRUE
    }
    if(e.con <= (1/sf3) && d3 == FALSE){
      n.rej3.con <- n.rej3.con+1
    }
    if(e.con > (1/sf3)){
      d3 = TRUE
    }
    e.con <- e.con * (1 - lambda.c + lambda.c * EES(x=x[i+w],r=r[i+w],z=z[i+w],lev=lev))
    e.out.con <- cbind(e.out.con, e.con)
    i = i+1
  }

  # Loop calculating the test martingale for GREE
  i <- 1
  n.rej1.akelly <- 0
  n.rej2.akelly <- 0
  n.rej3.akelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n-w)){
    if(e.akelly <= (1/sf1) && d1 == FALSE){
      n.rej1.akelly <- n.rej1.akelly+1
    }
    if(e.akelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.akelly <= (1/sf2) && d2 == FALSE){
      n.rej2.akelly <- n.rej2.akelly+1
    }
    if(e.akelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.akelly <= (1/sf3) && d3 == FALSE){
      n.rej3.akelly <- n.rej3.akelly+1
    }
    if(e.akelly > (1/sf3)){
      d3 = TRUE
    }

    # e-value for GREE
    
    lambda <- (mean(EES(x=x[i:(i+w-1)],r=r[i:(i+w-1)],z=z[i:(i+w-1)],lev=lev)) - 1) /
      (mean((EES(x=x[i:(i+w-1)],r=r[i:(i+w-1)],z=z[i:(i+w-1)],lev=lev) - mean(EES(x=x[i:(i+w-1)],r=r[i:(i+w-1)],z=z[i:(i+w-1)],lev=lev)))^2)
       + (mean(EES(x=x[i:(i+w-1)],r=r[i:(i+w-1)],z=z[i:(i+w-1)],lev=lev)) - 1)^2)
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda
    out.lambda <- cbind(out.lambda, lambda)
    
    
    e.akelly <- e.akelly * (1 - lambda + lambda * EES(x=x[i+w],r=r[i+w],z=z[i+w],lev=lev))
    e.out.akelly <- cbind(e.out.akelly, e.akelly)

    i = i+1
  }


  # Loop calculating the test martingale for GREL
  i <- 1
  n.rej1.rkelly <- 0
  n.rej2.rkelly <- 0
  n.rej3.rkelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n-w)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }

    # e-value for GREL
    
    lambda <- (mean(EES(x=x[i:(i+w-1)],r=r[i+w],z=z[i+w],lev=lev)) - 1) / (mean((EES(x=x[i:(i+w-1)],r=r[i+w],z=z[i+w],lev=lev)
              - mean(EES(x=x[i:(i+w-1)],r=r[i+w],z=z[i+w],lev=lev)))^2) + (mean(EES(x=x[i:(i+w-1)],r=r[i+w],z=z[i+w],lev=lev)) - 1)^2)

    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)


    e.rkelly <- e.rkelly * (1 - lambda + lambda * EES(x=x[i+w],r=r[i+w],z=z[i+w],lev=lev))
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)

    i = i+1
  }
  
  
  # Loop calculating the test martingale for GREM
  i <- 1
  n.rej1.mix <- 0
  n.rej2.mix <- 0
  n.rej3.mix <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n-w)){
    if(e.mix <= (1/sf1) && d1 == FALSE){
      n.rej1.mix <- n.rej1.mix+1
    }
    if(e.mix > (1/sf1)){
      d1 = TRUE
    }
    if(e.mix <= (1/sf2) && d2 == FALSE){
      n.rej2.mix <- n.rej2.mix+1
    }
    if(e.mix > (1/sf2)){
      d2 = TRUE
    }
    if(e.mix <= (1/sf3) && d3 == FALSE){
      n.rej3.mix <- n.rej3.mix+1
    }
    if(e.mix > (1/sf3)){
      d3 = TRUE
    }
    
    
    
    e.mix <- (e.out.akelly[i] + e.out.rkelly[i])/2
    e.out.mix <- cbind(e.out.mix, e.mix)
    
    i = i+1
  }


  return(list(e.out.con = e.out.con, n.rej1.con = n.rej1.con, n.rej2.con = n.rej2.con, n.rej3.con = n.rej3.con,
              e.out.akelly = e.out.akelly, n.rej1.akelly = n.rej1.akelly, n.rej2.akelly = n.rej2.akelly, n.rej3.akelly = n.rej3.akelly,
              e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly,
              e.out.mix = e.out.mix, n.rej1.mix = n.rej1.mix, n.rej2.mix = n.rej2.mix, n.rej3.mix = n.rej3.mix,
              out.lambda = out.lambda, out.lambda.rkelly = out.lambda.rkelly))
}


# E-testing ES for IID data and time series data


evalue.ES <- function(x, r, z, lev=0.875)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.con <- c() # create a vector to store e-values for constant lambda
  e.out.akelly <- c(1) # GREE
  e.out.rkelly <- c(1) # GREL
  e.out.mix <- c() # GREM
  out.lambda <- c(0) # output parameters tuned for GREE
  out.lambda.rkelly <- c(0) # output parameters tuned for GREL
  
  # test martingale with lambda = lambda.c
  lambda.c <- 0.01 # lambda we use for time series data
  # lambda.c <- 0.025 # We use this lamabda for iid observations in Section 2.1 of "Simulation and data analysis for e-backtesting"
  
  # Initial e-value for constant lambda
  e.con <- 1
  
  # Initial e-value for GREE
  e.akelly <- 1
  
  # Initial e-value for GREL
  e.rkelly <- 1
  
  # Initial e-value for GREM
  e.mix <- 1
  
  # Loop calculating the test martingale for constant lambda
  i <- 1
  n.rej1.con <- 0
  n.rej2.con <- 0
  n.rej3.con <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.con <= (1/sf1) && d1 == FALSE){
      n.rej1.con <- n.rej1.con+1
    }
    if(e.con > (1/sf1)){
      d1 = TRUE
    }
    if(e.con <= (1/sf2) && d2 == FALSE){
      n.rej2.con <- n.rej2.con+1
    }
    if(e.con > (1/sf2)){
      d2 = TRUE
    }
    if(e.con <= (1/sf3) && d3 == FALSE){
      n.rej3.con <- n.rej3.con+1
    }
    if(e.con > (1/sf3)){
      d3 = TRUE
    }
    e.con <- e.con * (1 - lambda.c + lambda.c * EES(x=x[i],r=r[i],z=z[i],lev=lev))
    e.out.con <- cbind(e.out.con, e.con)
    i = i+1
  }
  
  # Loop calculating the test martingale for GREE
  i <- 2
  n.rej1.akelly <- 1
  n.rej2.akelly <- 1
  n.rej3.akelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.akelly <= (1/sf1) && d1 == FALSE){
      n.rej1.akelly <- n.rej1.akelly+1
    }
    if(e.akelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.akelly <= (1/sf2) && d2 == FALSE){
      n.rej2.akelly <- n.rej2.akelly+1
    }
    if(e.akelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.akelly <= (1/sf3) && d3 == FALSE){
      n.rej3.akelly <- n.rej3.akelly+1
    }
    if(e.akelly > (1/sf3)){
      d3 = TRUE
    }
    
    
    # e-value for GREE
    lambda <- (mean(EES(x=x[1:(i-1)],r=r[1:(i-1)],z=z[1:(i-1)],lev=lev)) - 1) / (mean((EES(x=x[1:(i-1)],r=r[1:(i-1)],z=z[1:(i-1)],lev=lev)
              - mean(EES(x=x[1:(i-1)],r=r[1:(i-1)],z=z[1:(i-1)],lev=lev)))^2) + (mean(EES(x=x[1:(i-1)],r=r[1:(i-1)],z=z[1:(i-1)],lev=lev)) - 1)^2)
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda <- cbind(out.lambda, lambda)
    
    e.akelly <- e.akelly * (1 - lambda + lambda * EES(x=x[i],r=r[i],z=z[i],lev=lev))
    
    e.out.akelly <- cbind(e.out.akelly, e.akelly)
    
    i = i+1
  }
  
  # Loop calculating the test martingale for GREL
  i <- 2
  n.rej1.rkelly <- 1
  n.rej2.rkelly <- 1
  n.rej3.rkelly <- 1
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    
    # e-value for GREL
    
    lambda <- (mean(EES(x=x[1:(i-1)],r=r[i],z=z[i],lev=lev)) - 1) / (mean((EES(x=x[1:(i-1)],r=r[i],z=z[i],lev=lev)
              - mean(EES(x=x[1:(i-1)],r=r[i],z=z[i],lev=lev)))^2) + (mean(EES(x=x[1:(i-1)],r=r[i],z=z[i],lev=lev)) - 1)^2)
    lambda <- max(min(lambda, 1/2), 0) # tuning lambda, bounded by 0 and 1/2
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * EES(x=x[i],r=r[i],z=z[i],lev=lev))
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    i = i+1
  }
  
  
  # Loop calculating the test martingale for GREM
  i <- 1
  n.rej1.mix <- 0
  n.rej2.mix <- 0
  n.rej3.mix <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.mix <= (1/sf1) && d1 == FALSE){
      n.rej1.mix <- n.rej1.mix+1
    }
    if(e.mix > (1/sf1)){
      d1 = TRUE
    }
    if(e.mix <= (1/sf2) && d2 == FALSE){
      n.rej2.mix <- n.rej2.mix+1
    }
    if(e.mix > (1/sf2)){
      d2 = TRUE
    }
    if(e.mix <= (1/sf3) && d3 == FALSE){
      n.rej3.mix <- n.rej3.mix+1
    }
    if(e.mix > (1/sf3)){
      d3 = TRUE
    }
    
    
    
    e.mix <- (e.out.akelly[i] + e.out.rkelly[i])/2
    e.out.mix <- cbind(e.out.mix, e.mix)
    
    i = i+1
  }
  
  
  return(list(e.out.con = e.out.con, n.rej1.con = n.rej1.con, n.rej2.con = n.rej2.con, n.rej3.con = n.rej3.con,
              e.out.akelly = e.out.akelly, n.rej1.akelly = n.rej1.akelly, n.rej2.akelly = n.rej2.akelly, n.rej3.akelly = n.rej3.akelly,
              e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly,
              e.out.mix = e.out.mix, n.rej1.mix = n.rej1.mix, n.rej2.mix = n.rej2.mix, n.rej3.mix = n.rej3.mix,
              out.lambda = out.lambda, out.lambda.rkelly = out.lambda.rkelly))
}

# e-test for ES with known alternative Q (GRO)


# Function to compute E(X-z)_+ for X ~ N(0,1)
EtXmzp <- function(z)
{
  dnorm(z) - z*(1-pnorm(z))
}

# Function to compute E(X-z)_+^2 for X ~ N(0,1)
EtX2mzp <- function(z)
{
  (1+z^2)*(1-pnorm(z)) - dnorm(z)*z
}


evalue.ESQ <- function(x, r, z, lev=0.875)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.rkelly <- c() # e-process
  out.lambda.rkelly <- c() # output parameters
  
  # model-free e-statistic e_p
  e <- (pmax(x - z, 0)) / ((1 - lev) * (r - z))
  
  # Initial e-value
  e.rkelly <- 1
  
  # Loop calculating the test martingale
  i <- 1
  n.rej1.rkelly <- 0
  n.rej2.rkelly <- 0
  n.rej3.rkelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value
    f = function(pm){
      function(x){
        log(1-pm+pm*(pmax(x - z[i], 0))/((1-lev)*(r[i]-z[i])))*dnorm(x)
      }
    }
    fun = function(pm){
      integrate(f(pm), -Inf, Inf, rel.tol = 1e-2)$value
    }
    lambda = optimize(fun, interval = c(0,0.5), maximum = TRUE)$maximum
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * max(x[i] - z[i], 0) / ((1 - lev) * (r[i] - z[i])))
    
    
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    i = i+1
  }
  
  
  return(list(e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly))
}


# e-test for ES with known alternative Q AR-GARCH(1,1) (GRO method)

evalue.ES.GARCH <- function(x, r, z, mut, sigt, lev=0.875)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.rkelly <- c() # e-process
  out.lambda.rkelly <- c() # output parameters
  
  # Initial e-value
  e.rkelly <- 1
  
  # Loop calculating the test martingale
  i <- 1
  n.rej1.rkelly <- 0
  n.rej2.rkelly <- 0
  n.rej3.rkelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value
    
    # true parameters for skewed t rv
    nu=5; ga=1.5
    
    f = function(pm){
      function(x){
        log(1-pm+pm*(pmax(x - z[i], 0))/((1-lev)*(r[i]-z[i])))*dsgt(x,mu=mut[i],sigma=sigt[i],lambda=(ga^2-1)/(ga^2+1),p=2,q=nu/2)
      }
    }
    fun = function(pm){
      integrate(f(pm), -Inf, Inf, rel.tol = 1e-2)$value
    }
    lambda = optimize(fun, interval = c(0,0.5), maximum = TRUE)$maximum
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * max(x[i] - z[i], 0) / ((1 - lev) * (r[i] - z[i])))
    
    
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    
    i = i+1
  }
  
  
  return(list(e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly))
}


# e-test for GRO method in Example 7 (normal)

evalue.ESQ.eps <- function(x, r, z, eps, lev=0.875)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.rkelly <- c() # e-process
  out.lambda.rkelly <- c() # output parameters tuned
  
  # model-free e-statistic e_p
  e <- (pmax(x - z, 0)) / ((1 - lev) * (r - z))
  
  # Initial e-value
  e.rkelly <- 1
  
  # Loop calculating the test martingale
  i <- 1
  n.rej1.rkelly <- 0
  n.rej2.rkelly <- 0
  n.rej3.rkelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value
    f = function(pm){
      function(x){
        log(1-pm+pm*(pmax(x * (1 + eps[i]) - z[i], 0))/((1-lev)*(r[i]-z[i])))*dnorm(x)
      }
    }
    fun = function(pm){
      integrate(f(pm), -Inf, Inf, rel.tol = 1e-2)$value
    }
    lambda = optimize(fun, interval = c(0,0.5), maximum = TRUE)$maximum
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * max(x[i] - z[i], 0) / ((1 - lev) * (r[i] - z[i])))
    
    
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    i = i+1
  }
  
  
  return(list(e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly))
}


# e-test for GRO method in Example 7 (AR-GARCH(1,1))

evalue.ES.GARCH.eps <- function(x, r, z, mut, sigt, eps, lev=0.875)
{
  n = length(x) # out-of-sample size
  sf1 <- .5
  sf2 <- .2
  sf3 <- .1 # significance levels
  
  e.out.rkelly <- c() # e-process
  out.lambda.rkelly <- c() # output parameters tuned
  
  # model-free e-statistic e_p
  e <- (pmax(x - z, 0)) / ((1 - lev) * (r - z))
  
  # Initial e-value
  e.rkelly <- 1
  
  # Loop calculating the test martingale
  i <- 1
  n.rej1.rkelly <- 0
  n.rej2.rkelly <- 0
  n.rej3.rkelly <- 0
  d1 <- FALSE
  d2 <- FALSE
  d3 <- FALSE # dummmy variables
  while(i <= (n)){
    if(e.rkelly <= (1/sf1) && d1 == FALSE){
      n.rej1.rkelly <- n.rej1.rkelly+1
    }
    if(e.rkelly > (1/sf1)){
      d1 = TRUE
    }
    if(e.rkelly <= (1/sf2) && d2 == FALSE){
      n.rej2.rkelly <- n.rej2.rkelly+1
    }
    if(e.rkelly > (1/sf2)){
      d2 = TRUE
    }
    if(e.rkelly <= (1/sf3) && d3 == FALSE){
      n.rej3.rkelly <- n.rej3.rkelly+1
    }
    if(e.rkelly > (1/sf3)){
      d3 = TRUE
    }
    
    # e-value
    f = function(pm){
      function(x){
        log(1-pm+pm*(pmax(x * (1 + eps[i]) - z[i], 0))/((1-lev)*(r[i]-z[i])))*dsgt(x,mu=mut[i],sigma=sigt[i],lambda=(ga^2-1)/(ga^2+1),p=2,q=nu/2)
      }
    }
    fun = function(pm){
      integrate(f(pm), -Inf, Inf, rel.tol = 1e-2)$value
    }
    lambda = optimize(fun, interval = c(0,0.5), maximum = TRUE)$maximum
    
    e.rkelly <- e.rkelly * (1 - lambda + lambda * max(x[i] - z[i], 0) / ((1 - lev) * (r[i] - z[i])))
    
    
    out.lambda.rkelly <- cbind(out.lambda.rkelly, lambda)
    
    e.out.rkelly <- cbind(e.out.rkelly, e.rkelly)
    
    i = i+1
  }
  
  
  return(list(e.out.rkelly = e.out.rkelly, n.rej1.rkelly = n.rej1.rkelly, n.rej2.rkelly = n.rej2.rkelly, n.rej3.rkelly = n.rej3.rkelly, out.lambda.rkelly = out.lambda.rkelly))
}