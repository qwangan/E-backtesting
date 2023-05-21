library("spam")
library("rugarch")
library("fGarch")
library("MASS")
library("expectreg")
library("ismev")
library("lmom")
library("QRM")
library("skewt")
library("plotrix")
library("nloptr")
library("parallel")
library("sgt")

source("Rfns.R")
source("Rfns_HD22.R")


# Set number of cores
# NumCores <- detectCores()
NumCores <- 11

# Set number of runs
NumRuns <- 10000

# Set parameters
r.level <- 0.95 # levels of VaR and ES
n=250 # out-of-sample size to evaluate forecasts
w=n # presampling size
m     = 50          # moving-window length
d     = 5           # number of lags in LB-type test
reps  = 50000       # number of Monte Carlo runs
sig.level = .05    # significance level of monitoring procedures

set.seed(4123)

# Calculate mu.uc.ES
mu.VaR <- mu.uc.VaR(k.max=n+d, alpha=r.level)
mu.ES  <- mu.uc.ES( k.max=n+d, alpha=r.level, reps = 500000)
mu.iid <- mu.iid.MC(k.max=n+d, alpha=r.level, reps = 500000)

# Calculate critical values
cv.MCS    <- Simulated.cvs.MCS(m, n+d, alpha=1-r.level, reps, probs = 1-sig.level, mu.VaR, mu.ES, mu.iid)  # ... the simulated critical values

save(mu.VaR, mu.ES, mu.iid, cv.MCS, file = "CritVal.RData")
load("CritVal.RData")

# Function
rseed=sample(1:10^8,NumRuns*NumCores*2)
rseed = matrix(rseed, nrow  = NumRuns, byrow = TRUE)

simu <- function(iter){
  # "con" represents e-backtesting method with constant lambda,
  # "akelly" represents GREE method,
  # "rkelly" represents GREL method,
  # "mix" represents GREM method
  # "HD" represents the monitoring method in Hoga and Demetrescu (2022)
  
  # output values
  rejVaR1.con = 0
  rejVaR1.akelly = 0
  rejVaR1.rkelly = 0
  rejVaR1.mix = 0 # numbers of days where we reject for VaR (threshold 3.2)
  rejVaR2.con = 0
  rejVaR2.akelly = 0
  rejVaR2.rkelly = 0
  rejVaR2.mix = 0 # numbers of days where we reject for VaR (threshold 10)
  rejVaR3.con = 0
  rejVaR3.akelly = 0
  rejVaR3.rkelly = 0
  rejVaR3.mix = 0 # numbers of days where we reject for VaR (threshold 20)
  rejVaR.HD = 0

  rejnumVaR1.con <- rej.ARLVaR1.con <- 0
  rejnumVaR1.akelly <- rej.ARLVaR1.akelly <- 0
  rejnumVaR1.rkelly <- rej.ARLVaR1.rkelly <- 0
  rejnumVaR1.mix <- rej.ARLVaR1.mix <- 0 # used to count number of rejections
  rejnumVaR2.con <- rej.ARLVaR2.con <- 0
  rejnumVaR2.akelly <- rej.ARLVaR2.akelly <- 0
  rejnumVaR2.rkelly <- rej.ARLVaR2.rkelly <- 0
  rejnumVaR2.mix <- rej.ARLVaR2.mix <- 0
  rejnumVaR3.con <- rej.ARLVaR3.con <- 0
  rejnumVaR3.akelly <- rej.ARLVaR3.akelly <- 0
  rejnumVaR3.rkelly <- rej.ARLVaR3.rkelly <- 0
  rejnumVaR3.mix <- rej.ARLVaR3.mix <- 0
  rejnumVaR.HD <- rej.ARLVaR.HD <- 0

  rejES1.con = 0
  rejES1.akelly = 0
  rejES1.rkelly = 0
  rejES1.mix = 0 # numbers of days where we reject for ES (threshold 2)
  rejES2.con = 0
  rejES2.akelly = 0
  rejES2.rkelly = 0
  rejES2.mix = 0 # numbers of days where we reject for ES (threshold 5)
  rejES3.con = 0
  rejES3.akelly = 0
  rejES3.rkelly = 0
  rejES3.mix = 0 # numbers of days where we reject for ES (threshold 10)
  rejES.HD = 0

  rejnumES1.con <- rej.ARLES1.con <- 0
  rejnumES1.akelly <- rej.ARLES1.akelly <- 0
  rejnumES1.rkelly <- rej.ARLES1.rkelly <- 0
  rejnumES1.mix <- rej.ARLES1.mix <- 0
  rejnumES2.con <- rej.ARLES2.con <- 0
  rejnumES2.akelly <- rej.ARLES2.akelly <- 0
  rejnumES2.rkelly <- rej.ARLES2.rkelly <- 0
  rejnumES2.mix <- rej.ARLES2.mix <- 0
  rejnumES3.con <- rej.ARLES3.con <- 0
  rejnumES3.akelly <- rej.ARLES3.akelly <- 0
  rejnumES3.rkelly <- rej.ARLES3.rkelly <- 0
  rejnumES3.mix <- rej.ARLES3.mix <- 0
  rejnumES.HD <- rej.ARLES.HD <- 0

  # =======================================================
  # SIMULATION SET-UP
  # Synthetic dataset generated from an AR(1)-GARCH(1,1) process with skewed t innovations
  # AR-GARCH filter parameters:
  # mu= 0; ar1 = 0 # AR(1) part
  # omega=.00001; al=.04; be=.7 + .25 * (t > b*) # GARCH(1,1) parameters
  # Innovation distribution parameters:
  # nu=5 # shape parameter
  # ga=0.95 # skewness parameter
  # Burn-in period of 1000 points was used
  # =======================================================

  # b*
  b_star = 25 * (iter - 1)
  skip=1 # for backup seeds in case an error occurs

for(l in 1:NumRuns){

  if(l %% 100 == 0){
    cat("Core:", iter, "Run:", l, fill = TRUE)
  }

  seed = rseed[l,iter]
  cond=TRUE
  while (cond==TRUE){
  spec = garchSpec(model = list(mu = 0, ar = 0, omega = .00001, alpha = .04, beta = .7, skew = .95, shape = 5),
                   cond.dist = "sstd", rseed = seed)
  out = garchSim(spec, n = n+w, n.start = 1000, extended = TRUE)
  simdat = out$garch

  x1=-simdat[1:(w+b_star)] # time series to be used for fitting and forecasting (before structural change)

  spec = garchSpec(model = list(mu = 0, ar = 0, omega = .00001, alpha = .04, beta = .95, skew = .95, shape = 5),
                   cond.dist = "sstd", rseed = seed)
  out = garchSim(spec, n = n+w, n.start = 1000, extended = TRUE)
  simdat = out$garch

  x2=-simdat[(w+b_star+1):(w+n)] # time series to be used for fitting and forecasting (after structural change)

  x = c(x1,x2) # time series to be used for forecasting and testing
  y = x[(w+1):(w+n)] # time series to be used for testing

  # Estimation of parameters using presampling

    fit = garchFit(formula = ~ garch(1,1), data = x[1:w], cond.dist = "QMLE", include.mean = FALSE, trace = FALSE)
    fit.par = fit@fit$matcoef
    omega.hat = fit.par[1]
    alpha.hat = fit.par[2]
    beta.hat = fit.par[3]

    if(is.na(omega.hat)){omega.hat=1}
    if(is.na(alpha.hat)){alpha.hat=1}
    if(is.na(beta.hat)){beta.hat=1}

  if(is.character(fit)){cond=TRUE}
  else{
    cond=((beta.hat>=1)|(omega.hat<0) | (alpha.hat<0) |(beta.hat<0)|((alpha.hat + beta.hat)>=1))
  }

  if(cond){
    seed = rseed[skip,iter+NumCores]
    skip=skip+1
  }
  }

    tmp = omega.hat + alpha.hat * x[(w-d):(n+w-1)]^2
    sigt = filter(tmp, filter = c(beta.hat), method = "recursive", init = fit@h.t[w-d])

    resi = x[1:w] / fit@sigma.t
    F.n.resi       <- ecdf(resi) # Empirical cumulative distribution function of the residuals

  # Calculate VaR and ES predictions
    VaRa = quantile(resi, probs=r.level)
    VaRout = sqrt(sigt) * VaRa
    ESout = sqrt(sigt) * mean(resi[resi >= VaRa])

  # Calculate U as defined in Hoga and Demetrescu (2022)
    U  <- F.n.resi( x[(w - d + 1) : (n + w)] / sqrt(sigt) )

  # ===================================================
  # E-backtesting analysis of the simulation data
  # Data generating process: GARCH(1,1) with skewed t innovations
  # ===================================================


  ### summary for VaR

  tmp2=evalue_VaR_mv(x=y, r=VaRout,lev=r.level, m=50) # sequential test margtingale

  rejnumVaR1.con = rejnumVaR1.con + (tmp2$n.rej1.con < n)
  rejnumVaR1.akelly = rejnumVaR1.akelly + (tmp2$n.rej1.akelly < n)
  rejnumVaR1.rkelly = rejnumVaR1.rkelly + (tmp2$n.rej1.rkelly < n)
  rejnumVaR1.mix = rejnumVaR1.mix + (tmp2$n.rej1.mix < n)
  rejnumVaR2.con = rejnumVaR2.con + (tmp2$n.rej2.con < n)
  rejnumVaR2.akelly = rejnumVaR2.akelly + (tmp2$n.rej2.akelly < n)
  rejnumVaR2.rkelly = rejnumVaR2.rkelly + (tmp2$n.rej2.rkelly < n)
  rejnumVaR2.mix = rejnumVaR2.mix + (tmp2$n.rej2.mix < n)
  rejnumVaR3.con = rejnumVaR3.con + (tmp2$n.rej3.con < n)
  rejnumVaR3.akelly = rejnumVaR3.akelly + (tmp2$n.rej3.akelly < n)
  rejnumVaR3.rkelly = rejnumVaR3.rkelly + (tmp2$n.rej3.rkelly < n)
  rejnumVaR3.mix = rejnumVaR3.mix + (tmp2$n.rej3.mix < n)

  # numbers of days where we reject for VaR (threshold 3.2)
  if(((tmp2$n.rej1.con >= (b_star + 1)) & (tmp2$n.rej1.con != n))){
    rej.ARLVaR1.con = rej.ARLVaR1.con + 1 # only count rejection if detection occurred after actual break
    rejVaR1.con = rejVaR1.con + (tmp2$n.rej1.con - (b_star + 1))
  }
  if(((tmp2$n.rej1.akelly >= (b_star + 1)) & (tmp2$n.rej1.akelly != n))){
    rej.ARLVaR1.akelly = rej.ARLVaR1.akelly + 1
    rejVaR1.akelly = rejVaR1.akelly + (tmp2$n.rej1.akelly - (b_star + 1))
  }
  if(((tmp2$n.rej1.rkelly >= (b_star + 1)) & (tmp2$n.rej1.rkelly != n))){
    rej.ARLVaR1.rkelly = rej.ARLVaR1.rkelly + 1
    rejVaR1.rkelly = rejVaR1.rkelly + (tmp2$n.rej1.rkelly - (b_star + 1))
  }
  if(((tmp2$n.rej1.mix >= (b_star + 1)) & (tmp2$n.rej1.mix != n))){
    rej.ARLVaR1.mix = rej.ARLVaR1.mix + 1
    rejVaR1.mix = rejVaR1.mix + (tmp2$n.rej1.mix - (b_star + 1))
  }

  # numbers of days where we reject for VaR (threshold 10)
  if(((tmp2$n.rej2.con >= (b_star + 1)) & (tmp2$n.rej2.con != n))){
    rej.ARLVaR2.con = rej.ARLVaR2.con + 1 # only count rejection if detection occurred after actual break
    rejVaR2.con = rejVaR2.con + (tmp2$n.rej2.con - (b_star + 1))
  }
  if(((tmp2$n.rej2.akelly >= (b_star + 1)) & (tmp2$n.rej2.akelly != n))){
    rej.ARLVaR2.akelly = rej.ARLVaR2.akelly + 1
    rejVaR2.akelly = rejVaR2.akelly + (tmp2$n.rej2.akelly - (b_star + 1))
  }
  if(((tmp2$n.rej2.rkelly >= (b_star + 1)) & (tmp2$n.rej2.rkelly != n))){
    rej.ARLVaR2.rkelly = rej.ARLVaR2.rkelly + 1
    rejVaR2.rkelly = rejVaR2.rkelly + (tmp2$n.rej2.rkelly - (b_star + 1))
  }
  if(((tmp2$n.rej2.mix >= (b_star + 1)) & (tmp2$n.rej2.mix != n))){
    rej.ARLVaR2.mix = rej.ARLVaR2.mix + 1
    rejVaR2.mix = rejVaR2.mix + (tmp2$n.rej2.mix - (b_star + 1))
  }

  # numbers of days where we reject for VaR (threshold 20)
  if(((tmp2$n.rej3.con >= (b_star + 1)) & (tmp2$n.rej3.con != n))){
    rej.ARLVaR3.con = rej.ARLVaR3.con + 1 # only count rejection if detection occurred after actual break
    rejVaR3.con = rejVaR3.con + (tmp2$n.rej3.con - (b_star + 1))
  }
  if(((tmp2$n.rej3.akelly >= (b_star + 1)) & (tmp2$n.rej3.akelly != n))){
    rej.ARLVaR3.akelly = rej.ARLVaR3.akelly + 1
    rejVaR3.akelly = rejVaR3.akelly + (tmp2$n.rej3.akelly - (b_star + 1))
  }
  if(((tmp2$n.rej3.rkelly >= (b_star + 1)) & (tmp2$n.rej3.rkelly != n))){
    rej.ARLVaR3.rkelly = rej.ARLVaR3.rkelly + 1
    rejVaR3.rkelly = rejVaR3.rkelly + (tmp2$n.rej3.rkelly - (b_star + 1))
  }
  if(((tmp2$n.rej3.mix >= (b_star + 1)) & (tmp2$n.rej3.mix != n))){
    rej.ARLVaR3.mix = rej.ARLVaR3.mix + 1
    rejVaR3.mix = rejVaR3.mix + (tmp2$n.rej3.mix - (b_star + 1))
  }


  # ===================================================
  # Summary for ES

  tmp2=evalue_ES_mv(x=y, r=ESout, z=VaRout,lev=r.level, m=50)

  rejnumES1.con = rejnumES1.con + (tmp2$n.rej1.con < n)
  rejnumES1.akelly = rejnumES1.akelly + (tmp2$n.rej1.akelly < n)
  rejnumES1.rkelly = rejnumES1.rkelly + (tmp2$n.rej1.rkelly < n)
  rejnumES1.mix = rejnumES1.mix + (tmp2$n.rej1.mix < n)
  rejnumES2.con = rejnumES2.con + (tmp2$n.rej2.con < n)
  rejnumES2.akelly = rejnumES2.akelly + (tmp2$n.rej2.akelly < n)
  rejnumES2.rkelly = rejnumES2.rkelly + (tmp2$n.rej2.rkelly < n)
  rejnumES2.mix = rejnumES2.mix + (tmp2$n.rej2.mix < n)
  rejnumES3.con = rejnumES3.con + (tmp2$n.rej3.con < n)
  rejnumES3.akelly = rejnumES3.akelly + (tmp2$n.rej3.akelly < n)
  rejnumES3.rkelly = rejnumES3.rkelly + (tmp2$n.rej3.rkelly < n)
  rejnumES3.mix = rejnumES3.mix + (tmp2$n.rej3.mix < n)


  # numbers of days where we reject for ES (threshold 3.2)
  if(((tmp2$n.rej1.con >= (b_star + 1)) & (tmp2$n.rej1.con != n))){
    rej.ARLES1.con = rej.ARLES1.con + 1 # only count rejection if detection occurred after actual break
    rejES1.con = rejES1.con + (tmp2$n.rej1.con - (b_star + 1))
  }
  if(((tmp2$n.rej1.akelly >= (b_star + 1)) & (tmp2$n.rej1.akelly != n))){
    rej.ARLVaR1.akelly = rej.ARLVaR1.akelly + 1
    rejVaR1.akelly = rejVaR1.akelly + (tmp2$n.rej1.akelly - (b_star + 1))
  }
  if(((tmp2$n.rej1.rkelly >= (b_star + 1)) & (tmp2$n.rej1.rkelly != n))){
    rej.ARLVaR1.rkelly = rej.ARLVaR1.rkelly + 1
    rejVaR1.rkelly = rejVaR1.rkelly + (tmp2$n.rej1.rkelly - (b_star + 1))
  }
  if(((tmp2$n.rej1.mix >= (b_star + 1)) & (tmp2$n.rej1.mix != n))){
    rej.ARLVaR1.mix = rej.ARLVaR1.mix + 1
    rejVaR1.mix = rejVaR1.mix + (tmp2$n.rej1.mix - (b_star + 1))
  }

  # numbers of days where we reject for ES (threshold 10)
  if(((tmp2$n.rej2.con >= (b_star + 1)) & (tmp2$n.rej2.con != n))){
    rej.ARLES2.con = rej.ARLES2.con + 1 # only count rejection if detection occurred after actual break
    rejES2.con = rejES2.con + (tmp2$n.rej2.con - (b_star + 1))
  }
  if(((tmp2$n.rej2.akelly >= (b_star + 1)) & (tmp2$n.rej2.akelly != n))){
    rej.ARLES2.akelly = rej.ARLES2.akelly + 1
    rejES2.akelly = rejES2.akelly + (tmp2$n.rej2.akelly - (b_star + 1))
  }
  if(((tmp2$n.rej2.rkelly >= (b_star + 1)) & (tmp2$n.rej2.rkelly != n))){
    rej.ARLES2.rkelly = rej.ARLES2.rkelly + 1
    rejES2.rkelly = rejES2.rkelly + (tmp2$n.rej2.rkelly - (b_star + 1))
  }
  if(((tmp2$n.rej2.mix >= (b_star + 1)) & (tmp2$n.rej2.mix != n))){
    rej.ARLES2.mix = rej.ARLES2.mix + 1
    rejES2.mix = rejES2.mix + (tmp2$n.rej2.mix - (b_star + 1))
  }

  # numbers of days where we reject for ES (threshold 20)
  if(((tmp2$n.rej3.con >= (b_star + 1)) & (tmp2$n.rej3.con != n))){
    rej.ARLES3.con = rej.ARLES3.con + 1 # only count rejection if detection occurred after actual break
    rejES3.con = rejES3.con + (tmp2$n.rej3.con - (b_star + 1))
  }
  if(((tmp2$n.rej3.akelly >= (b_star + 1)) & (tmp2$n.rej3.akelly != n))){
    rej.ARLES3.akelly = rej.ARLES3.akelly + 1
    rejES3.akelly = rejES3.akelly + (tmp2$n.rej3.akelly - (b_star + 1))
  }
  if(((tmp2$n.rej3.rkelly >= (b_star + 1)) & (tmp2$n.rej3.rkelly != n))){
    rej.ARLES3.rkelly = rej.ARLES3.rkelly + 1
    rejES3.rkelly = rejES3.rkelly + (tmp2$n.rej3.rkelly - (b_star + 1))
  }
  if(((tmp2$n.rej3.mix >= (b_star + 1)) & (tmp2$n.rej3.mix != n))){
    rej.ARLES3.mix = rej.ARLES3.mix + 1
    rejES3.mix = rejES3.mix + (tmp2$n.rej3.mix - (b_star + 1))
  }



  # ===================================================
  # Monitoring method in Hoga and Demetrescu (2022)
  # ===================================================

  ### summary for VaR and ES
  tmp1 = HDVaRES(m, d, n, w, alpha=r.level, reps, sig.level, U, mu.VaR, mu.ES, mu.iid, cv.MCS)
  rejnumVaR.HD = rejnumVaR.HD + tmp1$rej.VaR.MCS.M
  rejnumES.HD = rejnumES.HD + tmp1$rej.ES.MCS.M

  if((is.na(tmp1$break.VaR.MCS.M) == FALSE) & (tmp1$break.VaR.MCS.M < Inf) & (tmp1$break.VaR.MCS.M >= (b_star + 1))){
    rej.ARLVaR.HD <- rej.ARLVaR.HD + 1  # only count rejection if detection occurred after actual break
    rejVaR.HD = rejVaR.HD + (tmp1$break.VaR.MCS.M - (b_star + 1))
  }
  if((is.na(tmp1$break.ES.MCS.M) == FALSE) & (tmp1$break.ES.MCS.M < Inf) & (tmp1$break.ES.MCS.M >= (b_star + 1))){
    rej.ARLES.HD <- rej.ARLES.HD + 1  # only count rejection if detection occurred after actual break
    rejES.HD = rejES.HD + (tmp1$break.ES.MCS.M - (b_star + 1))
  }
}

  rej.ARLVaR1.con = max(rej.ARLVaR1.con, 1)
  rej.ARLVaR1.akelly = max(rej.ARLVaR1.akelly, 1)
  rej.ARLVaR1.rkelly = max(rej.ARLVaR1.rkelly, 1)
  rej.ARLVaR1.mix = max(rej.ARLVaR1.mix, 1)
  rej.ARLVaR2.con = max(rej.ARLVaR2.con, 1)
  rej.ARLVaR2.akelly = max(rej.ARLVaR2.akelly, 1)
  rej.ARLVaR2.rkelly = max(rej.ARLVaR2.rkelly, 1)
  rej.ARLVaR2.mix = max(rej.ARLVaR2.mix, 1)
  rej.ARLVaR3.con = max(rej.ARLVaR3.con, 1)
  rej.ARLVaR3.akelly = max(rej.ARLVaR3.akelly, 1)
  rej.ARLVaR3.rkelly = max(rej.ARLVaR3.rkelly, 1)
  rej.ARLVaR3.mix = max(rej.ARLVaR3.mix, 1)
  rej.ARLVaR.HD = max(rej.ARLVaR.HD, 1)
  rej.ARLES1.con = max(rej.ARLES1.con, 1)
  rej.ARLES1.akelly = max(rej.ARLES1.akelly, 1)
  rej.ARLES1.rkelly = max(rej.ARLES1.rkelly, 1)
  rej.ARLES1.mix = max(rej.ARLES1.mix, 1)
  rej.ARLES2.con = max(rej.ARLES2.con, 1)
  rej.ARLES2.akelly = max(rej.ARLES2.akelly, 1)
  rej.ARLES2.rkelly = max(rej.ARLES2.rkelly, 1)
  rej.ARLES2.mix = max(rej.ARLES2.mix, 1)
  rej.ARLES3.con = max(rej.ARLES3.con, 1)
  rej.ARLES3.akelly = max(rej.ARLES3.akelly, 1)
  rej.ARLES3.rkelly = max(rej.ARLES3.rkelly, 1)
  rej.ARLES3.mix = max(rej.ARLES3.mix, 1)
  rej.ARLES.HD = max(rej.ARLES.HD, 1)



  return(c(rejVaR1.con/rej.ARLVaR1.con, rejVaR1.akelly/rej.ARLVaR1.akelly, rejVaR1.rkelly/rej.ARLVaR1.rkelly, rejVaR1.mix/rej.ARLVaR1.mix,
           rejVaR2.con/rej.ARLVaR2.con, rejVaR2.akelly/rej.ARLVaR2.con, rejVaR2.rkelly/rej.ARLVaR2.rkelly, rejVaR2.mix/rej.ARLVaR2.mix,
           rejVaR3.con/rej.ARLVaR3.con, rejVaR3.akelly/rej.ARLVaR3.akelly, rejVaR3.rkelly/rej.ARLVaR3.rkelly, rejVaR3.mix/rej.ARLVaR3.mix,
           rejnumVaR1.con/NumRuns, rejnumVaR1.akelly/NumRuns, rejnumVaR1.rkelly/NumRuns, rejnumVaR1.mix/NumRuns,
           rejnumVaR2.con/NumRuns, rejnumVaR2.akelly/NumRuns, rejnumVaR2.rkelly/NumRuns, rejnumVaR2.mix/NumRuns,
           rejnumVaR3.con/NumRuns, rejnumVaR3.akelly/NumRuns, rejnumVaR3.rkelly/NumRuns, rejnumVaR3.mix/NumRuns,
           rejES1.con/rej.ARLES1.con, rejES1.akelly/rej.ARLES1.akelly, rejES1.rkelly/rej.ARLES1.rkelly, rejES1.mix/rej.ARLES1.mix,
           rejES2.con/rej.ARLES2.con, rejES2.akelly/rej.ARLES2.akelly, rejES2.rkelly/rej.ARLES2.rkelly, rejES2.mix/rej.ARLES2.mix,
           rejES3.con/rej.ARLES3.con, rejES3.akelly/rej.ARLES3.akelly, rejES3.rkelly/rej.ARLES3.rkelly, rejES3.mix/rej.ARLES3.mix,
           rejnumES1.con/NumRuns, rejnumES1.akelly/NumRuns, rejnumES1.rkelly/NumRuns, rejnumES1.mix/NumRuns,
           rejnumES2.con/NumRuns, rejnumES2.akelly/NumRuns, rejnumES2.rkelly/NumRuns, rejnumES2.mix/NumRuns,
           rejnumES3.con/NumRuns, rejnumES3.akelly/NumRuns, rejnumES3.rkelly/NumRuns, rejnumES3.mix/NumRuns,
           rejVaR.HD/rej.ARLVaR.HD, rejES.HD/rej.ARLES.HD, rejnumVaR.HD/NumRuns, rejnumES.HD/NumRuns))
}


# Multiple-core computations
out <- do.call(cbind, mclapply(1:NumCores, simu, mc.cores = NumCores))
save(out, file="Result_compare.RDATA")


# Output values
load("Result_compare.RDATA")
avg.total = out
avg.rejVaR1.con = avg.total[1,]
avg.rejVaR1.akelly = avg.total[2,]
avg.rejVaR1.rkelly = avg.total[3,]
avg.rejVaR1.mix = avg.total[4,] # numbers of days where we reject for VaR (threshold 3.2)
avg.rejVaR2.con = avg.total[5,]
avg.rejVaR2.akelly = avg.total[6,]
avg.rejVaR2.rkelly = avg.total[7,]
avg.rejVaR2.mix = avg.total[8,] # numbers of days where we reject for VaR (threshold 10)
avg.rejVaR3.con = avg.total[9,]
avg.rejVaR3.akelly = avg.total[10,]
avg.rejVaR3.rkelly = avg.total[11,]
avg.rejVaR3.mix = avg.total[12,] # numbers of days where we reject for VaR (threshold 20)
avg.rejnumVaR1.con = avg.total[13,] * 100
avg.rejnumVaR1.akelly = avg.total[14,] * 100
avg.rejnumVaR1.rkelly = avg.total[15,] * 100
avg.rejnumVaR1.mix = avg.total[16,] * 100 # numbers of rejections for VaR (threshold 3.2)
avg.rejnumVaR2.con =avg.total[17,] * 100
avg.rejnumVaR2.akelly = avg.total[18,] * 100
avg.rejnumVaR2.rkelly = avg.total[19,] * 100
avg.rejnumVaR2.mix = avg.total[20,] * 100 # numbers of rejections for VaR (threshold 10)
avg.rejnumVaR3.con = avg.total[21,] * 100
avg.rejnumVaR3.akelly = avg.total[22,] * 100
avg.rejnumVaR3.rkelly = avg.total[23,] * 100
avg.rejnumVaR3.mix = avg.total[24,] * 100 # numbers of rejections for VaR (threshold 20)
avg.rejES1.con = avg.total[25,]
avg.rejES1.akelly = avg.total[26,]
avg.rejES1.rkelly = avg.total[27,]
avg.rejES1.mix = avg.total[28,] # numbers of days where we reject for ES (threshold 3.2)
avg.rejES2.con = avg.total[29,]
avg.rejES2.akelly = avg.total[30,]
avg.rejES2.rkelly = avg.total[31,]
avg.rejES2.mix = avg.total[32,] # numbers of days where we reject for ES (threshold 10)
avg.rejES3.con = avg.total[33,]
avg.rejES3.akelly = avg.total[34,]
avg.rejES3.rkelly = avg.total[35,]
avg.rejES3.mix = avg.total[36,] # numbers of days where we reject for ES (threshold 20)
avg.rejnumES1.con = avg.total[37,] * 100
avg.rejnumES1.akelly = avg.total[38,] * 100
avg.rejnumES1.rkelly = avg.total[39,] * 100
avg.rejnumES1.mix = avg.total[40,] * 100 # numbers of rejections for ES (threshold 3.2)
avg.rejnumES2.con = avg.total[41,] * 100
avg.rejnumES2.akelly = avg.total[42,] * 100
avg.rejnumES2.rkelly = avg.total[43,] * 100
avg.rejnumES2.mix = avg.total[44,] * 100 # numbers of rejections for ES (threshold 10)
avg.rejnumES3.con = avg.total[45,] * 100
avg.rejnumES3.akelly = avg.total[46,] * 100
avg.rejnumES3.rkelly = avg.total[47,] * 100
avg.rejnumES3.mix = avg.total[48,] * 100 # numbers of rejections for ES (threshold 20)
avg.rejVaR.HD = avg.total[49,]
avg.rejES.HD = avg.total[50,]
avg.rejnumVaR.HD = avg.total[51,] * 100
avg.rejnumES.HD = avg.total[52,] * 100 # Results in Hoga and Demetrescu (2022)


out.final = list(avg.rejVaR1.con = avg.rejVaR1.con, avg.rejVaR1.akelly = avg.rejVaR1.akelly, avg.rejVaR1.rkelly = avg.rejVaR1.rkelly, avg.rejVaR1.mix = avg.rejVaR1.mix,
                 avg.rejVaR2.con = avg.rejVaR2.con, avg.rejVaR2.akelly = avg.rejVaR2.akelly, avg.rejVaR2.rkelly = avg.rejVaR2.rkelly, avg.rejVaR2.mix = avg.rejVaR2.mix,
                 avg.rejVaR3.con = avg.rejVaR3.con, avg.rejVaR3.akelly = avg.rejVaR3.akelly, avg.rejVaR3.rkelly = avg.rejVaR3.rkelly, avg.rejVaR3.mix = avg.rejVaR3.mix,
                 avg.rejnumVaR1.con= avg.rejnumVaR1.con, avg.rejnumVaR1.akelly = avg.rejnumVaR1.akelly, avg.rejnumVaR1.rkelly = avg.rejnumVaR1.rkelly, avg.rejnumVaR1.mix = avg.rejnumVaR1.mix,
                 avg.rejnumVaR2.con= avg.rejnumVaR2.con, avg.rejnumVaR2.akelly = avg.rejnumVaR2.akelly, avg.rejnumVaR2.rkelly = avg.rejnumVaR2.rkelly, avg.rejnumVaR2.mix = avg.rejnumVaR2.mix,
                 avg.rejnumVaR3.con= avg.rejnumVaR3.con, avg.rejnumVaR3.akelly = avg.rejnumVaR3.akelly, avg.rejnumVaR3.rkelly = avg.rejnumVaR3.rkelly, avg.rejnumVaR3.mix = avg.rejnumVaR3.mix,
                 avg.rejES1.con = avg.rejES1.con, avg.rejES1.akelly = avg.rejES1.akelly, avg.rejES1.rkelly = avg.rejES1.rkelly, avg.rejES1.mix = avg.rejES1.mix,
                 avg.rejES2.con = avg.rejES2.con, avg.rejES2.akelly = avg.rejES2.akelly, avg.rejES2.rkelly = avg.rejES2.rkelly, avg.rejES2.mix = avg.rejES2.mix,
                 avg.rejES3.con = avg.rejES3.con, avg.rejES3.akelly = avg.rejES3.akelly, avg.rejES3.rkelly = avg.rejES3.rkelly, avg.rejES3.mix = avg.rejES3.mix,
                 avg.rejnumES1.con = avg.rejnumES1.con, avg.rejnumES1.akelly = avg.rejnumES1.akelly, avg.rejnumES1.rkelly = avg.rejnumES1.rkelly, avg.rejnumES1.mix = avg.rejnumES1.mix,
                 avg.rejnumES2.con = avg.rejnumES2.con, avg.rejnumES2.akelly = avg.rejnumES2.akelly, avg.rejnumES2.rkelly = avg.rejnumES2.rkelly, avg.rejnumES2.mix = avg.rejnumES2.mix,
                 avg.rejnumES3.con = avg.rejnumES3.con, avg.rejnumES3.akelly = avg.rejnumES3.akelly, avg.rejnumES3.rkelly = avg.rejnumES3.rkelly, avg.rejnumES3.mix = avg.rejnumES3.mix,
                 avg.rejVaR.HD = avg.rejVaR.HD, avg.rejES.HD = avg.rejES.HD, avg.rejnumVaR.HD = avg.rejnumVaR.HD, avg.rejnumES.HD = avg.rejnumES.HD)
save(out.final, file="Result_compare_final.RDATA")


avg.rejVaR1.con
avg.rejVaR1.akelly
avg.rejVaR1.rkelly
avg.rejVaR1.mix
avg.rejVaR2.con
avg.rejVaR2.akelly
avg.rejVaR2.rkelly
avg.rejVaR2.mix
avg.rejVaR3.con
avg.rejVaR3.akelly
avg.rejVaR3.rkelly
avg.rejVaR3.mix
avg.rejVaR.HD
avg.rejnumVaR1.con
avg.rejnumVaR1.akelly
avg.rejnumVaR1.rkelly
avg.rejnumVaR1.mix
avg.rejnumVaR2.con
avg.rejnumVaR2.akelly
avg.rejnumVaR2.rkelly
avg.rejnumVaR2.mix
avg.rejnumVaR3.con
avg.rejnumVaR3.akelly
avg.rejnumVaR3.rkelly
avg.rejnumVaR3.mix
avg.rejnumVaR.HD

avg.rejES1.con
avg.rejES1.akelly
avg.rejES1.rkelly
avg.rejES1.mix
avg.rejES2.con
avg.rejES2.akelly
avg.rejES2.rkelly
avg.rejES2.mix
avg.rejES3.con
avg.rejES3.akelly
avg.rejES3.rkelly
avg.rejES3.mix
avg.rejES.HD
avg.rejnumES1.con
avg.rejnumES1.akelly
avg.rejnumES1.rkelly
avg.rejnumES1.mix
avg.rejnumES2.con
avg.rejnumES2.akelly
avg.rejnumES2.rkelly
avg.rejnumES2.mix
avg.rejnumES3.con
avg.rejnumES3.akelly
avg.rejnumES3.rkelly
avg.rejnumES3.mix
avg.rejnumES.HD


# Plot the e-processes
load("Result_compare_final.RDATA")




setEPS()
postscript(file = "compareVaR.eps", width = 10, height = 8)
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
par(mar = c(0, 4.1, 4.1, 2.1))
matplot(seq(0,250,25), cbind(out.final$avg.rejnumVaR3.akelly[1:11], out.final$avg.rejnumVaR3.rkelly[1:11],
                             out.final$avg.rejnumVaR3.mix[1:11], out.final$avg.rejnumVaR.HD[1:11]),
        col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "", ylab = "detection percentage (%)", xaxt = "n", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1,1))
legend("topright", legend=c("GREE", "GREL", "GREM", "monitor"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)

par(mar = c(4.1, 4.1, 0, 2.1))
matplot(seq(0,250,25), cbind(out.final$avg.rejVaR3.akelly[1:11], out.final$avg.rejVaR3.rkelly[1:11],
                             out.final$avg.rejVaR3.mix[1:11], out.final$avg.rejVaR.HD[1:11]),
        col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "", ylab = "detection percentage (%)", xaxt = "n", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1,1))
legend("topright", legend=c("GREE", "GREL", "GREM", "monitor"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript(file = "compareES.eps", width = 10, height = 8)
layout(matrix(1:2, ncol = 1), widths = 1, heights = c(2,2), respect = FALSE)
par(mar = c(0, 4.1, 4.1, 2.1))
matplot(seq(0,250,25), cbind(out.final$avg.rejnumES3.akelly[1:11], out.final$avg.rejnumES3.rkelly[1:11],
                             out.final$avg.rejnumES3.mix[1:11], out.final$avg.rejnumES.HD[1:11]),
        col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "", ylab = "detection percentage (%)", xaxt = "n", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1,1))
legend("topright", legend=c("GREE", "GREL", "GREM", "monitor"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)

par(mar = c(4.1, 4.1, 0, 2.1))
matplot(seq(0,250,25), cbind(out.final$avg.rejES3.akelly[1:11], out.final$avg.rejES3.rkelly[1:11],
                             out.final$avg.rejES3.mix[1:11], out.final$avg.rejES.HD[1:11]),
        col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "", ylab = "detection percentage (%)", xaxt = "n", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1,1))
legend("topright", legend=c("GREE", "GREL", "GREM", "monitor"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)
dev.off()

