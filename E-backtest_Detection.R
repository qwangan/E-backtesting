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


# Set number of cores
NumCores <- 100

# Set number of runs
NumRuns <- 1000

# Set parameters
avec=.99 # vector of alpha levels for VaR
nvec=.975 # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR
n=2000 # out-of-sample size to evaluate forecasts
e.lim <- c(-1,5) # bounds for log e-values



# Function
set.seed(4123)
rseed=sample(1:10^8,NumRuns*2)
rseed = matrix(rseed, nrow  = NumCores, byrow = TRUE)

simu <- function(iter){
  # =======================================================
  # SIMULATION SET-UP
  # Synthetic dataset generated from an AR(1)-GARCH(1,1) process with skewed t innovations
  # AR-GARCH filter parameters:
  # mu=-.05; ar1 = .3 # AR(1) part
  # omega=.01; al=.1; be=.85 # GARCH(1,1) parameters
  # Innovation distribution parameters:
  # nu=5 # shape parameter
  # ga=1.5 # skewness parameter
  # Burn-in period of 1000 points was used
  # =======================================================
  
  
  spec = garchSpec(model = list(mu = -.05, ar = .3, omega = .01, alpha = .1, beta = .85, skew = 1.5, shape = 5),
                   cond.dist = "sstd", rseed = iter)
  out = garchSim(spec, n = 11000, n.start = 1000, extended = TRUE)
  simdat = out$garch
  
  
  n=2000 # out-of-sample size to evaluate forecasts
  
  w=500 # moving window size
  
  x=simdat[(w+1):(n+2*w)] # time series to be used for fitting and forecasting
  
  
  
  # =======================================================
  # Skewed-t innovations for first 1000 points and then change to normal innovations for last 1000 points
  # =======================================================
  
  # Skewed-t for first 1000 points
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")
  
  VaR.st <- matrix(nrow=n/2,ncol=length(VaR.levels))
  
  ES.st <- matrix(nrow=n/2,ncol=length(nvec))
  
  for(i in 1:(n/2))
    
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    
    # mean and std dev of skewed t rv
    nu=coef(fit)["shape"]; ga=coef(fit)["skew"]
    
    # calculate VaR of normalized skewed-t
    qmodel = qsgt(VaR.levels,mu=0,sigma=1,lambda=(ga^2-1)/(ga^2+1),p=2,q=nu/2)
    
    esmodel = NULL # expected shortfall for the assumed model/distribution
    
    # parameters of Hansen (1994)
    lam = -(ga^2-1)/(ga^2+1) # transfer from ga of Fernandez and Steel (1998) to lam of Hansen (1994)
    # minus sign is taken because we change from loss to return
    c = gamma((nu+1)/2) / (gamma(nu/2) * sqrt(pi * (nu-2)))
    a = 4*lam*c*(nu-2)/(nu-1)
    b = sqrt(1+3*lam^2-a^2)
    
    # Calculate explicit ES for skewed-t described in Patton et al. (2019)
    for(j in inu){
      # We make some transformation here because we are handling loss data
      if(qmodel[j] >= (a/b)){
        alpha_tilde1 = psgt(b/(1-lam)*(-qmodel[j]+a/b),mu=0,sigma=1,lambda=0,p=2,q=nu/2)
        es_t1 = sqrt((nu-2)/nu) * nu^(nu/2)/(2*(alpha_tilde1)*sqrt(pi))*gamma((nu-1)/2)/gamma(nu/2)*(qt(1-alpha_tilde1, df=nu)^2+nu)^((1-nu)/2) # ES for standardized t distribution
        esmodel = c(esmodel, -(alpha_tilde1)/(1-VaR.levels[j]) * (1-lam) * (-a/b - (1-lam)/b*es_t1))
      }else{
        lam2 = -lam
        ga2 = 1/ga
        nu2 = nu
        c2 = c
        a2 = 4*lam2*c2*(nu2-2)/(nu2-1)
        b2 = b
        alpha_tilde2 = psgt(b2/(1-lam2)*(qsgt(VaR.levels[j],mu=0,sigma=1,lambda=lam2,p=2,q=nu2/2)+a2/b2),mu=0,sigma=1,lambda=0,p=2,q=nu2/2)
        es_t2 = sqrt((nu2-2)/nu2) * nu2^(nu2/2)/(2*(alpha_tilde2)*sqrt(pi))*gamma((nu2-1)/2)/gamma(nu2/2)*(qt(1-alpha_tilde2, df=nu2)^2+nu2)^((1-nu2)/2)
        esmodel = c(esmodel, -(alpha_tilde2)/(1-VaR.levels[j]) * (1-lam2) * (-a2/b2 - (1-lam2)/b2*es_t2))
      }
    }
    
    
    
    foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)
    
    VaR.st[i,]=foc$VaRmodel
    
    ES.st[i,]=foc$ESmodel
  }
  
  # Normal for last 1000 points
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
  qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1) # VaR values
  esmodel = (ddist("norm", qdist("norm",p=nvec,mu=0,sigma=1), mu=0, sigma=1)) / (1 - nvec)
  
  VaR.norm <- matrix(nrow=n/2,ncol=length(VaR.levels))
  
  ES.norm <- matrix(nrow=n/2,ncol=length(nvec))
  
  
  # ----------------------------------------------------
  
  for(i in 1:n/2)
  {
    fit = ugarchfit(spec, x[(i+n/2):(i+n/2+w-1)], solver="hybrid")
    foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)
    
    VaR.norm[i,]=foc$VaRmodel
    
    ES.norm[i,]=foc$ESmodel
  }
  
  
  # Data preparation
  
  
  # --------------------------------------------------------
  # Verifying observations on which to assess forecasts
  y=simdat[(2*w+1):(n+2*w)] # 1000 is the in-sample to accomodate moving estimation windows of up to size 1000
  
  
  # --------------------------------------------------------
  # Forecasts
  
  #=====================================================
  # VaR_alpha
  #=====================================================
  
  VaRout = c(VaR.st[,1],VaR.norm[,1])
  VaRout = VaRout[1:(n)]
  
  
  #=====================================================
  # (VaR_nu, ES_nu)
  #=====================================================
  
  VaRoutb = c(VaR.st[,2], VaR.norm[,2])
  ESout = c(ES.st, ES.norm)
  
  
  
  VaRoutb = VaRoutb[1:(n)]; ESout = ESout[1:(n)]
  
  
  
  # ===================================================
  # E-backtesting analysis of the simulation data
  # Data generating process: AR(1)-GARCH(1,1) with skewed t innovations
  # "con" represents e-backtesting method with constant lambda,
  # "akelly" represents GREE method,
  # "rkelly" represents GREL method,
  # "mix" represents GREM method
  # ===================================================
  
  w <- 500 # time window
  nm = 1
  
  ### summary for VaR
  eVaR.con = matrix(nrow=nm, ncol=1) # VaR forecast
  eVaR.akelly = matrix(nrow=nm, ncol=1)
  eVaR.rkelly = matrix(nrow=nm, ncol=1)
  eVaR.mix = matrix(nrow=nm, ncol=1) # final e-values for VaR
  rejVaR1.con = matrix(nrow=nm, ncol=1)
  rejVaR1.akelly = matrix(nrow=nm, ncol=1)
  rejVaR1.rkelly = matrix(nrow=nm, ncol=1)
  rejVaR1.mix = matrix(nrow=nm, ncol=1) # numbers of days where we reject for VaR (threshold 2)
  rejVaR2.con = matrix(nrow=nm, ncol=1)
  rejVaR2.akelly = matrix(nrow=nm, ncol=1)
  rejVaR2.rkelly = matrix(nrow=nm, ncol=1)
  rejVaR2.mix = matrix(nrow=nm, ncol=1) # numbers of days where we reject for VaR (threshold 5)
  rejVaR3.con = matrix(nrow=nm, ncol=1)
  rejVaR3.akelly = matrix(nrow=nm, ncol=1)
  rejVaR3.rkelly = matrix(nrow=nm, ncol=1)
  rejVaR3.mix = matrix(nrow=nm, ncol=1) # numbers of days where we reject for VaR (threshold 10)
  tmVaR.con = matrix(nrow=nm, ncol=(length(y)))
  tmVaR.akelly = matrix(nrow=nm, ncol=(length(y)))
  tmVaR.rkelly = matrix(nrow=nm, ncol=(length(y)))
  tmVaR.mix = matrix(nrow=nm, ncol=(length(y))) # sequential test martingale for VaR
  
  for(i in 1:nm)
  {
    tmp=evalue_VaR_mv(x=y, r=VaRout, lev=avec, m=w) # sequential test margtingale
    eVaR.con[i] = log(tmp$e.out.con[length(tmp$e.out.con)])
    eVaR.akelly[i] = log(tmp$e.out.akelly[length(tmp$e.out.akelly)])
    eVaR.rkelly[i] = log(tmp$e.out.rkelly[length(tmp$e.out.rkelly)])
    eVaR.mix[i] = log(tmp$e.out.mix[length(tmp$e.out.mix)])
    rejVaR1.con[i] = tmp$n.rej1.con
    rejVaR1.akelly[i] = tmp$n.rej1.akelly
    rejVaR1.rkelly[i] = tmp$n.rej1.rkelly
    rejVaR1.mix[i] = tmp$n.rej1.mix
    rejVaR2.con[i] = tmp$n.rej2.con
    rejVaR2.akelly[i] = tmp$n.rej2.akelly
    rejVaR2.rkelly[i] = tmp$n.rej2.rkelly
    rejVaR2.mix[i] = tmp$n.rej2.mix
    rejVaR3.con[i] = tmp$n.rej3.con
    rejVaR3.akelly[i] = tmp$n.rej3.akelly
    rejVaR3.rkelly[i] = tmp$n.rej3.rkelly
    rejVaR3.mix[i] = tmp$n.rej3.mix
    tmVaR.con[i,1:length(tmp$e.out.con)] <- log(tmp$e.out.con)
    tmVaR.akelly[i,1:length(tmp$e.out.akelly)] <- log(tmp$e.out.akelly)
    tmVaR.rkelly[i,1:length(tmp$e.out.rkelly)] <- log(tmp$e.out.rkelly)
    tmVaR.mix[i,1:length(tmp$e.out.mix)] <- log(tmp$e.out.mix)
  }
  
  rejnumVaR1.con = (rejVaR1.con < n)
  rejnumVaR1.akelly = (rejVaR1.akelly < n)
  rejnumVaR1.rkelly = (rejVaR1.rkelly < n)
  rejnumVaR1.mix = (rejVaR1.mix < n)
  rejnumVaR2.con = (rejVaR2.con < n)
  rejnumVaR2.akelly = (rejVaR2.akelly < n)
  rejnumVaR2.rkelly = (rejVaR2.rkelly < n)
  rejnumVaR2.mix = (rejVaR2.mix < n)
  rejnumVaR3.con = (rejVaR3.con < n)
  rejnumVaR3.akelly = (rejVaR3.akelly < n)
  rejnumVaR3.rkelly = (rejVaR3.rkelly < n)
  rejnumVaR3.mix = (rejVaR3.mix < n)
  
  
  # ===================================================
  # Summary for ES
  eES.con = matrix(nrow=nm, ncol=1) # VaR-ES forecast
  eES.akelly = matrix(nrow=nm, ncol=1)
  eES.rkelly = matrix(nrow=nm, ncol=1)
  eES.mix = matrix(nrow=nm, ncol=1) # final e-values for ES
  rejES1.con = matrix(nrow=nm, ncol=1)
  rejES1.akelly = matrix(nrow=nm, ncol=1)
  rejES1.rkelly = matrix(nrow=nm, ncol=1)
  rejES1.mix = matrix(nrow=nm, ncol=1) # numbers of days where we reject for ES (threshold 2)
  rejES2.con = matrix(nrow=nm, ncol=1)
  rejES2.akelly = matrix(nrow=nm, ncol=1)
  rejES2.rkelly = matrix(nrow=nm, ncol=1)
  rejES2.mix = matrix(nrow=nm, ncol=1) # numbers of days where we reject for ES (threshold 5)
  rejES3.con = matrix(nrow=nm, ncol=1)
  rejES3.akelly = matrix(nrow=nm, ncol=1)
  rejES3.rkelly = matrix(nrow=nm, ncol=1)
  rejES3.mix = matrix(nrow=nm, ncol=1) # numbers of days where we reject for ES (threshold 10)
  tmES.con = matrix(nrow=nm, ncol=(length(y)))
  tmES.akelly = matrix(nrow=nm, ncol=(length(y)))
  tmES.rkelly = matrix(nrow=nm, ncol=(length(y)))
  tmES.mix = matrix(nrow=nm, ncol=(length(y))) # sequential test martingale for ES
  
  for(i in 1:nm)
  {
    tmp=evalue_ES_mv(x=y, r=ESout, z=VaRoutb, lev=nvec, m=w) # sequential test margtingale
    eES.con[i] = log(tmp$e.out.con[length(tmp$e.out.con)])
    eES.akelly[i] = log(tmp$e.out.akelly[length(tmp$e.out.akelly)])
    eES.rkelly[i] = log(tmp$e.out.rkelly[length(tmp$e.out.rkelly)])
    eES.mix[i] = log(tmp$e.out.mix[length(tmp$e.out.mix)])
    rejES1.con[i] = tmp$n.rej1.con
    rejES1.akelly[i] = tmp$n.rej1.akelly
    rejES1.rkelly[i] = tmp$n.rej1.rkelly
    rejES1.mix[i] = tmp$n.rej1.mix
    rejES2.con[i] = tmp$n.rej2.con
    rejES2.akelly[i] = tmp$n.rej2.akelly
    rejES2.rkelly[i] = tmp$n.rej2.rkelly
    rejES2.mix[i] = tmp$n.rej2.mix
    rejES3.con[i] = tmp$n.rej3.con
    rejES3.akelly[i] = tmp$n.rej3.akelly
    rejES3.rkelly[i] = tmp$n.rej3.rkelly
    rejES3.mix[i] = tmp$n.rej3.mix
    tmES.con[i,1:length(tmp$e.out.con)] <- log(tmp$e.out.con)
    tmES.akelly[i,1:length(tmp$e.out.akelly)] <- log(tmp$e.out.akelly)
    tmES.rkelly[i,1:length(tmp$e.out.rkelly)] <- log(tmp$e.out.rkelly)
    tmES.mix[i,1:length(tmp$e.out.mix)] <- log(tmp$e.out.mix)
  }
  
  rejnumES1.con = (rejES1.con < n)
  rejnumES1.akelly = (rejES1.akelly < n)
  rejnumES1.rkelly = (rejES1.rkelly < n)
  rejnumES1.mix = (rejES1.mix < n)
  rejnumES2.con = (rejES2.con < n)
  rejnumES2.akelly = (rejES2.akelly < n)
  rejnumES2.rkelly = (rejES2.rkelly < n)
  rejnumES2.mix = (rejES2.mix < n)
  rejnumES3.con = (rejES3.con < n)
  rejnumES3.akelly = (rejES3.akelly < n)
  rejnumES3.rkelly = (rejES3.rkelly < n)
  rejnumES3.mix = (rejES3.mix < n)
  
  return(c(eVaR.con, eVaR.akelly, eVaR.rkelly, eVaR.mix, rejVaR1.con, rejVaR1.akelly, rejVaR1.rkelly, rejVaR1.mix,
           rejVaR2.con, rejVaR2.akelly, rejVaR2.rkelly, rejVaR2.mix, rejVaR3.con, rejVaR3.akelly, rejVaR3.rkelly, rejVaR3.mix,
           rejnumVaR1.con, rejnumVaR1.akelly, rejnumVaR1.rkelly, rejnumVaR1.mix, rejnumVaR2.con, rejnumVaR2.akelly, rejnumVaR2.rkelly, rejnumVaR2.mix,
           rejnumVaR3.con, rejnumVaR3.akelly, rejnumVaR3.rkelly, rejnumVaR3.mix, tmVaR.con, tmVaR.akelly, tmVaR.rkelly, tmVaR.mix,
           eES.con, eES.akelly, eES.rkelly, eES.mix, rejES1.con, rejES1.akelly, rejES1.rkelly, rejES1.mix,
           rejES2.con, rejES2.akelly, rejES2.rkelly, rejES2.mix, rejES3.con, rejES3.akelly, rejES3.rkelly, rejES3.mix,
           rejnumES1.con, rejnumES1.akelly, rejnumES1.rkelly, rejnumES1.mix, rejnumES2.con, rejnumES2.akelly, rejnumES2.rkelly, rejnumES2.mix,
           rejnumES3.con, rejnumES3.akelly, rejnumES3.rkelly, rejnumES3.mix, tmES.con, tmES.akelly, tmES.rkelly, tmES.mix))
}

result = function(num){
  out = c()
  i = 0
  i.valid = 0
  while(i.valid <= (NumRuns/NumCores - 1) & i <= ncol(rseed)){
    i = i+1
    print(i)
    iter=rseed[num,i]
    val = tryCatch(simu(iter),error = function(e) NULL)
    if (is.null(val)){next}
    out = rbind(out,val)
    i.valid = i.valid+1
  }
  out
}

# Multiple-core computations
out <- do.call(rbind, mclapply(1:NumCores, result, mc.cores = NumCores))
save(out, file="Result_simulation.RDATA")


# Output values
nm=1
avg.total = c()
for(i in 1:((56+8*n)*nm)){
  avg = mean(out[,i][out[,i] != n])
  avg.total = cbind(avg.total, avg)
}
avg.eVaR.con = matrix(avg.total[1:(nm)],nrow=nm,ncol=1)
avg.eVaR.akelly = matrix(avg.total[(nm+1):(2*nm)],nrow=nm,ncol=1)
avg.eVaR.rkelly = matrix(avg.total[(2*nm+1):(3*nm)],nrow=nm,ncol=1)
avg.eVaR.mix = matrix(avg.total[(3*nm+1):(4*nm)],nrow=nm,ncol=1) # final e-values for VaR
avg.rejVaR1.con = matrix(avg.total[(4*nm+1):(5*nm)],nrow=nm,ncol=1)
avg.rejVaR1.akelly = matrix(avg.total[(5*nm+1):(6*nm)],nrow=nm,ncol=1)
avg.rejVaR1.rkelly = matrix(avg.total[(6*nm+1):(7*nm)],nrow=nm,ncol=1)
avg.rejVaR1.mix = matrix(avg.total[(7*nm+1):(8*nm)],nrow=nm,ncol=1) # numbers of days where we reject for VaR (threshold 2)
avg.rejVaR2.con = matrix(avg.total[(8*nm+1):(9*nm)],nrow=nm,ncol=1)
avg.rejVaR2.akelly = matrix(avg.total[(9*nm+1):(10*nm)],nrow=nm,ncol=1)
avg.rejVaR2.rkelly = matrix(avg.total[(10*nm+1):(11*nm)],nrow=nm,ncol=1)
avg.rejVaR2.mix = matrix(avg.total[(11*nm+1):(12*nm)],nrow=nm,ncol=1) # numbers of days where we reject for VaR (threshold 5)
avg.rejVaR3.con = matrix(avg.total[(12*nm+1):(13*nm)],nrow=nm,ncol=1)
avg.rejVaR3.akelly = matrix(avg.total[(13*nm+1):(14*nm)],nrow=nm,ncol=1)
avg.rejVaR3.rkelly = matrix(avg.total[(14*nm+1):(15*nm)],nrow=nm,ncol=1)
avg.rejVaR3.mix = matrix(avg.total[(15*nm+1):(16*nm)],nrow=nm,ncol=1) # numbers of days where we reject for VaR (threshold 10)
avg.rejnumVaR1.con = matrix(avg.total[(16*nm+1):(17*nm)],nrow=nm,ncol=1)
avg.rejnumVaR1.akelly = matrix(avg.total[(17*nm+1):(18*nm)],nrow=nm,ncol=1)
avg.rejnumVaR1.rkelly = matrix(avg.total[(18*nm+1):(19*nm)],nrow=nm,ncol=1)
avg.rejnumVaR1.mix = matrix(avg.total[(19*nm+1):(20*nm)],nrow=nm,ncol=1) # numbers of rejections for VaR (threshold 2)
avg.rejnumVaR2.con = matrix(avg.total[(20*nm+1):(21*nm)],nrow=nm,ncol=1)
avg.rejnumVaR2.akelly = matrix(avg.total[(21*nm+1):(22*nm)],nrow=nm,ncol=1)
avg.rejnumVaR2.rkelly = matrix(avg.total[(22*nm+1):(23*nm)],nrow=nm,ncol=1)
avg.rejnumVaR2.mix = matrix(avg.total[(23*nm+1):(24*nm)],nrow=nm,ncol=1) # numbers of rejections for VaR (threshold 5)
avg.rejnumVaR3.con = matrix(avg.total[(24*nm+1):(25*nm)],nrow=nm,ncol=1)
avg.rejnumVaR3.akelly = matrix(avg.total[(25*nm+1):(26*nm)],nrow=nm,ncol=1)
avg.rejnumVaR3.rkelly = matrix(avg.total[(26*nm+1):(27*nm)],nrow=nm,ncol=1)
avg.rejnumVaR3.mix = matrix(avg.total[(27*nm+1):(28*nm)],nrow=nm,ncol=1) # numbers of rejections for VaR (threshold 10)
avg.tmVaR.con = matrix(avg.total[(28*nm+1):((28+n)*nm)],nrow=nm,ncol=n)
avg.tmVaR.akelly = matrix(avg.total[((28+n)*nm+1):((28+2*n)*nm)],nrow=nm,ncol=n)
avg.tmVaR.rkelly = matrix(avg.total[((28+2*n)*nm+1):((28+3*n)*nm)],nrow=nm,ncol=n)
avg.tmVaR.mix = matrix(avg.total[((28+3*n)*nm+1):((28+4*n)*nm)],nrow=nm,ncol=n) # sequential test martingale for VaR
avg.eES.con = matrix(avg.total[((28+4*n)*nm+1):((29+4*n)*nm)],nrow=nm,ncol=1)
avg.eES.akelly = matrix(avg.total[((29+4*n)*nm+1):((30+4*n)*nm)],nrow=nm,ncol=1)
avg.eES.rkelly = matrix(avg.total[((30+4*n)*nm+1):((31+4*n)*nm)],nrow=nm,ncol=1)
avg.eES.mix = matrix(avg.total[((31+4*n)*nm+1):((32+4*n)*nm)],nrow=nm,ncol=1) # final e-values for ES
avg.rejES1.con = matrix(avg.total[((32+4*n)*nm+1):((33+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES1.akelly = matrix(avg.total[((33+4*n)*nm+1):((34+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES1.rkelly = matrix(avg.total[((34+4*n)*nm+1):((35+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES1.mix = matrix(avg.total[((35+4*n)*nm+1):((36+4*n)*nm)],nrow=nm,ncol=1) # numbers of days where we reject for ES (threshold 2)
avg.rejES2.con = matrix(avg.total[((36+4*n)*nm+1):((37+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES2.akelly = matrix(avg.total[((37+4*n)*nm+1):((38+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES2.rkelly = matrix(avg.total[((38+4*n)*nm+1):((39+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES2.mix = matrix(avg.total[((39+4*n)*nm+1):((40+4*n)*nm)],nrow=nm,ncol=1) # numbers of days where we reject for ES (threshold 5)
avg.rejES3.con = matrix(avg.total[((40+4*n)*nm+1):((41+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES3.akelly = matrix(avg.total[((41+4*n)*nm+1):((42+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES3.rkelly = matrix(avg.total[((42+4*n)*nm+1):((43+4*n)*nm)],nrow=nm,ncol=1)
avg.rejES3.mix = matrix(avg.total[((43+4*n)*nm+1):((44+4*n)*nm)],nrow=nm,ncol=1) # numbers of days where we reject for ES (threshold 10)
avg.rejnumES1.con = matrix(avg.total[((44+4*n)*nm+1):((45+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES1.akelly = matrix(avg.total[((45+4*n)*nm+1):((46+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES1.rkelly = matrix(avg.total[((46+4*n)*nm+1):((47+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES1.mix = matrix(avg.total[((47+4*n)*nm+1):((48+4*n)*nm)],nrow=nm,ncol=1) # numbers of rejections for ES (threshold 2)
avg.rejnumES2.con = matrix(avg.total[((48+4*n)*nm+1):((49+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES2.akelly = matrix(avg.total[((49+4*n)*nm+1):((50+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES2.rkelly = matrix(avg.total[((50+4*n)*nm+1):((51+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES2.mix = matrix(avg.total[((51+4*n)*nm+1):((52+4*n)*nm)],nrow=nm,ncol=1) # numbers of rejections for ES (threshold 5)
avg.rejnumES3.con = matrix(avg.total[((52+4*n)*nm+1):((53+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES3.akelly = matrix(avg.total[((53+4*n)*nm+1):((54+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES3.rkelly = matrix(avg.total[((54+4*n)*nm+1):((55+4*n)*nm)],nrow=nm,ncol=1)
avg.rejnumES3.mix = matrix(avg.total[((55+4*n)*nm+1):((56+4*n)*nm)],nrow=nm,ncol=1) # numbers of rejections for ES (threshold 10)
avg.tmES.con = matrix(avg.total[((56+4*n)*nm+1):((56+5*n)*nm)],nrow=nm,ncol=n)
avg.tmES.akelly = matrix(avg.total[((56+5*n)*nm+1):((56+6*n)*nm)],nrow=nm,ncol=n)
avg.tmES.rkelly = matrix(avg.total[((56+6*n)*nm+1):((56+7*n)*nm)],nrow=nm,ncol=n)
avg.tmES.mix = matrix(avg.total[((56+7*n)*nm+1):((56+8*n)*nm)],nrow=nm,ncol=n) # sequential test martingale for ES


out.final = list(avg.eVaR.con = avg.eVaR.con, avg.eVaR.akelly = avg.eVaR.akelly, avg.eVaR.rkelly = avg.eVaR.rkelly, avg.eVaR.mix = avg.eVaR.mix,
                 avg.rejVaR1.con = avg.rejVaR1.con, avg.rejVaR1.akelly = avg.rejVaR1.akelly, avg.rejVaR1.rkelly = avg.rejVaR1.rkelly, avg.rejVaR1.mix = avg.rejVaR1.mix,
                 avg.rejVaR2.con = avg.rejVaR2.con, avg.rejVaR2.akelly = avg.rejVaR2.akelly, avg.rejVaR2.rkelly = avg.rejVaR2.rkelly, avg.rejVaR2.mix = avg.rejVaR2.mix,
                 avg.rejVaR3.con = avg.rejVaR3.con, avg.rejVaR3.akelly = avg.rejVaR3.akelly, avg.rejVaR3.rkelly = avg.rejVaR3.rkelly, avg.rejVaR3.mix = avg.rejVaR3.mix,
                 avg.rejnumVaR1.con= avg.rejnumVaR1.con, avg.rejnumVaR1.akelly = avg.rejnumVaR1.akelly, avg.rejnumVaR1.rkelly = avg.rejnumVaR1.rkelly, avg.rejnumVaR1.mix = avg.rejnumVaR1.mix,
                 avg.rejnumVaR2.con= avg.rejnumVaR2.con, avg.rejnumVaR2.akelly = avg.rejnumVaR2.akelly, avg.rejnumVaR2.rkelly = avg.rejnumVaR2.rkelly, avg.rejnumVaR2.mix = avg.rejnumVaR2.mix,
                 avg.rejnumVaR3.con= avg.rejnumVaR3.con, avg.rejnumVaR3.akelly = avg.rejnumVaR3.akelly, avg.rejnumVaR3.rkelly = avg.rejnumVaR3.rkelly, avg.rejnumVaR3.mix = avg.rejnumVaR3.mix,
                 avg.tmVaR.con = avg.tmVaR.con, avg.tmVaR.akelly = avg.tmVaR.akelly, avg.tmVaR.rkelly = avg.tmVaR.rkelly, avg.tmVaR.mix = avg.tmVaR.mix,
                 avg.eES.con = avg.eES.con, avg.eES.akelly = avg.eES.akelly, avg.eES.rkelly = avg.eES.rkelly, avg.eES.mix = avg.eES.mix,
                 avg.rejES1.con = avg.rejES1.con, avg.rejES1.akelly = avg.rejES1.akelly, avg.rejES1.rkelly = avg.rejES1.rkelly, avg.rejES1.mix = avg.rejES1.mix,
                 avg.rejES2.con = avg.rejES2.con, avg.rejES2.akelly = avg.rejES2.akelly, avg.rejES2.rkelly = avg.rejES2.rkelly, avg.rejES2.mix = avg.rejES2.mix,
                 avg.rejES3.con = avg.rejES3.con, avg.rejES3.akelly = avg.rejES3.akelly, avg.rejES3.rkelly = avg.rejES3.rkelly, avg.rejES3.mix = avg.rejES3.mix,
                 avg.rejnumES1.con = avg.rejnumES1.con, avg.rejnumES1.akelly = avg.rejnumES1.akelly, avg.rejnumES1.rkelly = avg.rejnumES1.rkelly, avg.rejnumES1.mix = avg.rejnumES1.mix,
                 avg.rejnumES2.con = avg.rejnumES2.con, avg.rejnumES2.akelly = avg.rejnumES2.akelly, avg.rejnumES2.rkelly = avg.rejnumES2.rkelly, avg.rejnumES2.mix = avg.rejnumES2.mix,
                 avg.rejnumES3.con = avg.rejnumES3.con, avg.rejnumES3.akelly = avg.rejnumES3.akelly, avg.rejnumES3.rkelly = avg.rejnumES3.rkelly, avg.rejnumES3.mix = avg.rejnumES3.mix,
                 avg.tmES.con = avg.tmES.con, avg.tmES.akelly = avg.tmES.akelly, avg.tmES.rkelly = avg.tmES.rkelly, avg.tmES.mix = avg.tmES.mix)
save(out.final, file="Result_detection_final.RDATA")


avg.eVaR.con
avg.eVaR.akelly
avg.eVaR.rkelly
avg.eVaR.mix
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
avg.eES.con
avg.eES.akelly
avg.eES.rkelly
avg.eES.mix
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


# Plot the e-processes
load("Result_detection_final.RDATA")

setEPS()
postscript("VaRdetect.eps")
matplot(1:2000, cbind(out.final$avg.tmVaR.akelly[1,1:2000], out.final$avg.tmVaR.rkelly[1,1:2000], out.final$avg.tmVaR.mix[1,1:2000]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("GREE", "GREL", "GREM"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()



setEPS()
postscript("ESdetect.eps")
matplot(1:2000, cbind(out.final$avg.tmES.akelly[1,1:2000], out.final$avg.tmES.rkelly[1,1:2000], out.final$avg.tmES.mix[1,1:2000]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("GREE", "GREL", "GREM"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()




