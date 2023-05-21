# Monitoring method in Hoga and Demetrescu (2022)

# ===========================================================
# Calculate mu.uc.VaR for VaR MCS detectors
# ===========================================================

# Inputs: 
# k.max: maximum k for which mu.uc.VaR(k) is required
# alpha: risk level

mu.uc.VaR <- function(k.max, alpha){
  mu.alpha.VaR <- numeric( k.max ) 
  
  for(k in 1:k.max){
    for(i in 1:k){
      mu.alpha.VaR[k] <- mu.alpha.VaR[k] + abs(k*(1-alpha) - i) * choose(k, i) * (1-alpha)^i * alpha^(k-i)
    }
  }
  return( mu.alpha.VaR )
}


# ===========================================================
# Calculate mu.uc.ES for ES MCS detectors
# ===========================================================

# Inputs: 
# k.max: maximum k for which mu.uc.ES(k) is required
# alpha: risk level
# reps: number of Monte Carlo repetitions

mu.uc.ES <- function(k.max, alpha, reps){
  U.ges <- matrix(0, nrow=reps, ncol=k.max)
  for(i in 1:reps){
    U <- runif(k.max)                    # generate U[0,1] r.v.s
    U <- (U - (1-alpha)/2)*1*(U <= (1- alpha))    # censor U
    U.ges[i, ] <- abs(cumsum(U))         # sum censored U cumulatively
  }
  return( apply(U.ges, 2, mean) )
}


# ===========================================================
# Calculate mu.iid for MCS detectors
# ===========================================================

# Inputs: 
# k.max: maximum k for which mu.uc.ES(k) is required
# alpha: risk level
# reps: number of Monte Carlo repetitions

mu.iid.MC <- function(k.max, alpha, reps){
  I.ges <- matrix(0, nrow=reps, ncol=k.max)
  for(i in 1:reps){
    I <- rbinom(k.max, size=1, prob=1-alpha)
    for(k in 1:k.max){
      t.i <- which(I[1:k] == 1)
      
      if(length(t.i) > 0){
        I.ges[i, k]  <- sum( (diff( c(0, t.i, k+1) ))^2 )
        s            <- length(t.i)                       # number of VaR violations in expanding window
        P.S.eq.s     <- dbinom(x=s, size=k, prob=1-alpha)   # Prob. of S=s
        I.ges[i, k]  <- 1/P.S.eq.s * I.ges[i, k]
      }
      
    }
  }
  return( apply(I.ges, 2, mean) )
}


# ===========================================================
# Function to calculate critical values for VaR & ES MCS backtest (CUSUM & MOSUM & SBQ)
# ===========================================================

# Inputs: 
# m: Moving block length
# n: length of monitoring period
# alpha: risk level
# reps: number of simulation runs
# probs: significance levels
# mu.VaR: Values of mu_{uc}^{VaR}(k) for k=1,...n
# mu.ES: Values of mu_{uc}^{ES}(k) for k=1,...n
# a: weight for UC and IID

Simulated.cvs.MCS <- function(m, T, alpha, reps, probs, mu.VaR, mu.ES, mu.iid, a=0.5){    # Calculate sequence of finite-sample detectors
  # 0. Preliminaries
  sup.VaR.R   <- sup.VaR.M <- sup.VaR.S <- numeric(reps)
  sup.ES.R    <- sup.ES.M  <- sup.ES.S  <- numeric(reps)
  mu.alpha.ES <- mu.ES
  mu.alpha.VaR<- mu.VaR
  Delta.mk    <- c(1:m, rep.int(m, T-m))         # lengths of moving window 
  
  for(i in 1:reps){
    U      <- runif(T)
    I      <- 1*(U <= alpha)
    U      <- U*(U <= alpha)         # censor U's
    U.cens <- (U - alpha/2)*1*(U <= alpha)
    
    # 1. Moving Sum detector
    # 1.1 MCS.uc.ES
    
    sum.U.cens  <- cumsum(U.cens) - cumsum(c(rep(0, m), head(U.cens, -m )))
    MCS.uc.ES.M <- abs(sum.U.cens) / mu.alpha.ES[Delta.mk]
    
    # 1.2 MCS.uc.VaR
    sum.I        <- cumsum(I) - cumsum(c(rep(0, m), head(I, -m )))
    MCS.uc.VaR.M <- abs(sum.I - Delta.mk*alpha) / mu.alpha.VaR[Delta.mk]
    
    # 1.3 MCS.iid
    MCS.iid.M <- rep.int(0, times=T)
    for(k in 1:T){
      Ind <- I[(max(0,k-m)+1) : k]
      t.i <- which(Ind == 1)
      if(length(t.i) > 0){    # otherwise if length(t.i)==0, there are no violations s.t. MCS.iid=0
        MCS.iid.M[k] <- sum( ( diff(c( 0, t.i, Delta.mk[k]+1 )) )^2 )
        s            <- length(t.i)                                 # number of VaR violations in moving window
        P.S.eq.s     <- dbinom(x=s, size=Delta.mk[k], prob=alpha)   # Prob. of S=s
        MCS.iid.M[k] <- 1 / mu.iid[Delta.mk[k]] * 1/P.S.eq.s * MCS.iid.M[k]
      }
    }
    sup.ES.M[i]  <- max(a * (MCS.uc.ES.M + MCS.uc.VaR.M)  + (1-a) * MCS.iid.M)
    sup.VaR.M[i] <- max(a * MCS.uc.VaR.M                  + (1-a) * MCS.iid.M)
  }
  
  return( list(ES.M  = quantile(sup.ES.M, probs), VaR.M = quantile(sup.VaR.M, probs)) )
}

# ===========================================================
# Calculate MCS detectors for VaR & ES (CUSUM & MOSUM)
# ===========================================================

# Inputs: 
# U: sequence of PITs
# m: moving block length
# alpha: risk level
# mu.VaR: Null-hypothetical expectations
# mu.ES: Null-hypothetical expectations
# a: weight for UC and IID

MCS <- function(U, m, alpha, mu.VaR, mu.ES, mu.iid, a=0.5){    # Calculate sequence of finite-sample detectors
  # 0. Preliminaries
  I      <- 1*(U <= alpha)         # create hit sequence
  U      <- U*(U <= alpha)         # censor U's
  U.cens <- (U - alpha/2)*1*(U <= alpha)
  n      <- length(I)
  
  mu.alpha.ES  <- mu.ES           # save values of mu.uc.ES
  mu.alpha.VaR <- mu.VaR
  Delta.mk     <- c(1:m, rep.int(m, n-m))
  
  # 1. Moving Sum detector
  # 1.1 MCS.uc.ES
  sum.U.cens  <- cumsum(U.cens) - cumsum(c(rep(0, m), head(U.cens, -m )))
  MCS.uc.ES.M <- abs(sum.U.cens) / mu.alpha.ES[Delta.mk]
  
  # 1.2 MCS.uc.VaR
  sum.I        <- cumsum(I) - cumsum(c(rep(0, m), head(I, -m )))
  MCS.uc.VaR.M <- abs(sum.I - Delta.mk*alpha) / mu.alpha.VaR[Delta.mk]
  
  # 1.3 MCS.iid
  MCS.iid.M <- rep.int(0, times=n)
  for(k in 1:n){
    Ind <- I[(max(0,k-m)+1) : k]
    t.i <- which(Ind == 1)
    
    if(length(t.i) > 0){    # otherwise if length(t.i)==0, there are no violations s.t. MCS.iid=0
      MCS.iid.M[k] <- sum( ( diff(c( 0, t.i, Delta.mk[k]+1 )) )^2 )
      s            <- length(t.i)                                 # number of VaR violations in moving window
      P.S.eq.s     <- dbinom(x=s, size=Delta.mk[k], prob=alpha)   # Prob. of S=s
      MCS.iid.M[k] <- 1 / mu.iid[Delta.mk[k]] * 1/P.S.eq.s * MCS.iid.M[k]
    }
  }
  
  return(list( ES.M  = a*(MCS.uc.ES.M + MCS.uc.VaR.M) + (1-a)*MCS.iid.M,
               VaR.M = a*MCS.uc.VaR.M                 + (1-a)*MCS.iid.M) )
}


# VaR and ES monitoring

HDVaRES <- function(m, d, n, w, alpha, reps, sig.level, U, mu.VaR, mu.ES, mu.iid, cv.MCS){
  # Monte Carlo simulation detector
  res.MCS       <- MCS(U[1 : (n+d)], m, 1-alpha, mu.VaR, mu.ES, mu.iid)  # we use the same data for MCS test even though it does not include the first d lags
  sup.M.VaR.MCS <- max( res.MCS$VaR.M )
  sup.M.ES.MCS  <- max( res.MCS$ES.M )

  rej.VaR.MCS.M    <- (sup.M.VaR.MCS > cv.MCS$VaR.M)
  rej.ES.MCS.M     <- (sup.M.ES.MCS > cv.MCS$ES.M)
  
  
  # Number of passed backtests with monitoring
  # Calculate different (possibly infinite) break times
  break.VaR.MCS.M <- min(which(res.MCS$VaR.M > cv.MCS$VaR.M))
  break.VaR.MCS.M <- max( break.VaR.MCS.M - d, 0)
  
  break.ES.MCS.M   <- min(which(res.MCS$ES.M > cv.MCS$ES.M))
  break.ES.MCS.M   <- max( break.ES.MCS.M - d, 0)
  

return( list(sup.M.VaR.MCS = sup.M.VaR.MCS, sup.M.ES.MCS = sup.M.ES.MCS, cv.VaR = cv.MCS$VaR.M, cv.ES = cv.MCS$ES.M,
             rej.VaR.MCS.M = rej.VaR.MCS.M, rej.ES.MCS.M = rej.ES.MCS.M, break.VaR.MCS.M = break.VaR.MCS.M, break.ES.MCS.M = break.ES.MCS.M) )
}