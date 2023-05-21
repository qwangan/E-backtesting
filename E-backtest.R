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
avec=c(.95, .99) # vector of alpha levels for VaR
nvec=c(.875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR
n=500 # out-of-sample size to evaluate forecasts
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
  mut = (out$garch) - (out$sigma) * (out$eps)
  sigt = out$sigma
  simdat = out$garch
  
  
  n=500 # out-of-sample size to evaluate forecasts
  
  w=500 # moving window size
  
  x=simdat[501:(n+1000)] # time series to be used for fitting and forecasting
  
  
  # =======================================================
  # Normal innovations
  # =======================================================
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
  qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1) # VaR values
  esmodel = (ddist("norm", qdist("norm",p=nvec,mu=0,sigma=1), mu=0, sigma=1)) / (1 - nvec)
  
  VaR.norm <- matrix(nrow=n,ncol=length(VaR.levels))
  
  ES.norm <- matrix(nrow=n,ncol=length(nvec))
  
  
  # estimated parameters
  fit.par.norm <- matrix(nrow=n, ncol=5) # 5 model parameters
  mut.norm  <- vector(mode="numeric", length=n)
  sigt.norm  <- vector(mode="numeric", length=n)
  
  # ----------------------------------------------------
  
  for(i in 1:n)
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    fit.par.norm[i,] = coef(fit)
    foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)
    
    VaR.norm[i,]=foc$VaRmodel
    
    ES.norm[i,]=foc$ESmodel
    mut.norm[i] = foc$mut
    sigt.norm[i] = foc$sigt
  }
  
  
  # =======================================================
  # Student t innovations
  # =======================================================
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="std")
  
  VaR.t <- matrix(nrow=n,ncol=length(VaR.levels))
  
  ES.t <- matrix(nrow=n,ncol=length(nvec))
  
  fit.par.t <- matrix(nrow=n, ncol=6) # 6 model parameters
  mut.t  <- vector(mode="numeric", length=n)
  sigt.t  <- vector(mode="numeric", length=n)
  
  
  # ----------------------------------------------------
  for(i in 1:n)
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    fit.par.t[i,] = coef(fit)
    nu=coef(fit)["shape"]
    
    qmodel = qdist("std",p=VaR.levels,mu=0,sigma=1,shape=nu)
    
    # expected shortfall for the assumed model/distribution
    esmodel = (dt(qt(nvec, df=nu), df=nu)) / (1-nvec) * (nu + (qt(nvec, df=nu))^2) / (nu-1)
    esmodel = sqrt((nu-2)/nu)*esmodel
    
    foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)
    
    VaR.t[i,]=foc$VaRmodel
    
    ES.t[i,]=foc$ESmodel
    mut.t[i] = foc$mut
    sigt.t[i] = foc$sigt
  }
  
  
  # =======================================================
  # skewed Student t innovations
  # The version of Fernandez & Steel (1998)
  # =======================================================
  
  spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")
  
  VaR.st <- matrix(nrow=n,ncol=length(VaR.levels))
  
  ES.st <- matrix(nrow=n,ncol=length(nvec))
  fit.par.st <- matrix(nrow=n, ncol=7) # 7 model parameters
  mut.st  <- vector(mode="numeric", length=n)
  sigt.st  <- vector(mode="numeric", length=n)
  
  for(i in 1:n)
    
  {
    fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
    fit.par.st[i,] = coef(fit)
    
    # mean and std dev of skewed t rv
    nu=coef(fit)["shape"]; ga=coef(fit)["skew"]
    
    # m = mean.st(shape=nu,skew=ga)
    # s = sqrt(var.st(shape=nu,skew=ga))
    
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
    
    # for(j in inu)
    # esmodel = c(esmodel, integrate(function(x) x*dskt(x, df=nu, gamma=ga), qmodel[j]*s+m, Inf)$value/(1-VaR.levels[j]))
    # esmodel=(esmodel-m)/s
    
    
    foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)
    
    VaR.st[i,]=foc$VaRmodel
    
    ES.st[i,]=foc$ESmodel
    mut.st[i] = foc$mut
    sigt.st[i] = foc$sigt
  }
  
  
  # Data preparation
  
  # Parameters
  nu=5; xi=1.5
  mst = mean.st(shape=nu,skew=xi) # mean and sd for optimal forecast computations
  sst = sqrt(var.st(shape=nu,skew=xi))
  
  # recovering the mu[t] and sig[t] used in the data generation
  # keeping values relevant for the optimal forecast computations
  mut=mut[1001:(n+1000)]
  sigt=sigt[1001:(n+1000)]
  
  # --------------------------------------------------------
  # Verifying observations on which to assess forecasts
  y=simdat[1001:(n+1000)] # 1000 is the in-sample to accomodate moving estimation windows of up to size 1000
  
  
  # --------------------------------------------------------
  # Forecasts
  
  #=====================================================
  # VaR_alpha
  #=====================================================
  
  levels = avec
  VaRa=qsstd(p=levels[1],nu=nu,xi=xi)  # a-VaR for the innovation distribution
  VaRopt=mut + sigt*VaRa
  VaRout1 = rbind(VaR.norm[,1], VaR.t[,1], VaR.st[,1], VaRopt)
  
  VaRa=qsstd(p=levels[2],nu=nu,xi=xi)  # a-VaR for the innovation distribution
  VaRopt=mut + sigt*VaRa
  VaRout2 = rbind(VaR.norm[,2], VaR.t[,2], VaR.st[,2], VaRopt)
  
  VaRout1 = VaRout1[1:4,1:(n)]
  VaRout2 = VaRout2[1:4,1:(n)]
  
  
  #=====================================================
  # (VaR_nu, ES_nu)
  #=====================================================
  
  levels = nvec
  VaRa = qskt(p=levels[1],df=nu,gamma=xi) # for skew t variable in standard form of the density
  ESa=integrate(function(x) x*dskt(x, df=nu, gamma=xi), VaRa, Inf)$value/(1-levels[1])
  
  VaRa = (VaRa-mst)/sst
  ESa=(ESa-mst)/sst # adjustment for skew t distribution with mean zero and sd=1
  VaRopt=mut + sigt*VaRa
  ESopt = mut + sigt*ESa
  
  VaRout1b = rbind(VaR.norm[,3], VaR.t[,3], VaR.st[,3], VaRopt)
  ESout1 = rbind(ES.norm[,1], ES.t[,1], ES.st[,1], ESopt)
  
  
  VaRa = qskt(p=levels[2],df=nu,gamma=xi) # for skew t variable in standard form of the density
  ESa=integrate(function(x) x*dskt(x, df=nu, gamma=xi), VaRa, Inf)$value/(1-levels[2])
  
  VaRa = (VaRa-mst)/sst
  ESa=(ESa-mst)/sst # adjustment for skew t distribution with mean zero and sd=1
  VaRopt=mut + sigt*VaRa
  ESopt = mut + sigt*ESa
  
  VaRout2b = rbind(VaR.norm[,4], VaR.t[,4], VaR.st[,4], VaRopt)
  ESout2 = rbind(ES.norm[,2], ES.t[,2], ES.st[,2], ESopt)
  
  
  VaRout1b = VaRout1b[1:4,1:(n)]; ESout1 = ESout1[1:4,1:(n)]
  VaRout2b = VaRout2b[1:4,1:(n)]; ESout2 = ESout2[1:4,1:(n)]
  
  
  
  # ===================================================
  # E-backtesting analysis of the simulation data
  # Data generating process: AR(1)-GARCH(1,1) with skewed t innovations
  # "con" represents e-backtesting method with constant lambda,
  # "akelly" represents GREE method,
  # "rkelly" represents GREL method,
  # "mix" represents GREM method
  # ===================================================
  
  err <- .1 # error for under- and over-reporting
  w <- 0 # time window
  nm = 3
  
  ### summary for VaR
  fcVaR = c(rowMeans(VaRout1), rowMeans(VaRout2)) # VaR forecast
  eVaR.con = matrix(nrow=2*nm, ncol=3)
  eVaR.akelly = matrix(nrow=2*nm, ncol=3)
  eVaR.rkelly = matrix(nrow=2*nm, ncol=3)
  eVaR.mix = matrix(nrow=2*nm, ncol=3) # final e-values for VaR
  rejVaR1.con = matrix(nrow=2*nm, ncol=3)
  rejVaR1.akelly = matrix(nrow=2*nm, ncol=3)
  rejVaR1.rkelly = matrix(nrow=2*nm, ncol=3)
  rejVaR1.mix = matrix(nrow=2*nm, ncol=3) # numbers of days where we reject for VaR (threshold 2)
  rejVaR2.con = matrix(nrow=2*nm, ncol=3)
  rejVaR2.akelly = matrix(nrow=2*nm, ncol=3)
  rejVaR2.rkelly = matrix(nrow=2*nm, ncol=3)
  rejVaR2.mix = matrix(nrow=2*nm, ncol=3) # numbers of days where we reject for VaR (threshold 5)
  rejVaR3.con = matrix(nrow=2*nm, ncol=3)
  rejVaR3.akelly = matrix(nrow=2*nm, ncol=3)
  rejVaR3.rkelly = matrix(nrow=2*nm, ncol=3)
  rejVaR3.mix = matrix(nrow=2*nm, ncol=3) # numbers of days where we reject for VaR (threshold 10)
  tmVaR.con = matrix(nrow=2*nm, ncol=(length(y)-w))
  tmVaR.akelly = matrix(nrow=2*nm, ncol=(length(y)-w))
  tmVaR.rkelly = matrix(nrow=2*nm, ncol=(length(y)-w))
  tmVaR.mix = matrix(nrow=2*nm, ncol=(length(y)-w)) # sequential test martingale for VaR
  
  for(i in 1:nm)
  {
    tmp2=evalue.VaR(x=y, r=VaRout1[i,],lev=avec[1]) # sequential test margtingale
    tmp3=evalue.VaR(x=y, r=VaRout1[i,] * (1-err),lev=avec[1]) # sequential test margtingale for under-reporting
    tmp4=evalue.VaR(x=y, r=VaRout1[i,] * (1+err),lev=avec[1]) # sequential test margtingale for over-reporting
    eVaR.con[i,] = log(c(tmp2$e.out.con[length(tmp2$e.out.con)], tmp3$e.out.con[length(tmp3$e.out.con)],
                         tmp4$e.out.con[length(tmp4$e.out.con)]))
    eVaR.akelly[i,] = log(c(tmp2$e.out.akelly[length(tmp2$e.out.akelly)], tmp3$e.out.akelly[length(tmp3$e.out.akelly)],
                            tmp4$e.out.akelly[length(tmp4$e.out.akelly)]))
    eVaR.rkelly[i,] = log(c(tmp2$e.out.rkelly[length(tmp2$e.out.rkelly)], tmp3$e.out.rkelly[length(tmp3$e.out.rkelly)],
                            tmp4$e.out.rkelly[length(tmp4$e.out.rkelly)]))
    eVaR.mix[i,] = log(c(tmp2$e.out.mix[length(tmp2$e.out.mix)], tmp3$e.out.mix[length(tmp3$e.out.mix)],
                            tmp4$e.out.mix[length(tmp4$e.out.mix)]))
    rejVaR1.con[i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp4$n.rej1.con)
    rejVaR1.akelly[i,] = cbind(tmp2$n.rej1.akelly, tmp3$n.rej1.akelly, tmp4$n.rej1.akelly)
    rejVaR1.rkelly[i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp4$n.rej1.rkelly)
    rejVaR1.mix[i,] = cbind(tmp2$n.rej1.mix, tmp3$n.rej1.mix, tmp4$n.rej1.mix)
    rejVaR2.con[i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp4$n.rej2.con)
    rejVaR2.akelly[i,] = cbind(tmp2$n.rej2.akelly, tmp3$n.rej2.akelly, tmp4$n.rej2.akelly)
    rejVaR2.rkelly[i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp4$n.rej2.rkelly)
    rejVaR2.mix[i,] = cbind(tmp2$n.rej2.mix, tmp3$n.rej2.mix, tmp4$n.rej2.mix)
    rejVaR3.con[i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp4$n.rej3.con)
    rejVaR3.akelly[i,] = cbind(tmp2$n.rej3.akelly, tmp3$n.rej3.akelly, tmp4$n.rej3.akelly)
    rejVaR3.rkelly[i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp4$n.rej3.rkelly)
    rejVaR3.mix[i,] = cbind(tmp2$n.rej3.mix, tmp3$n.rej3.mix, tmp4$n.rej3.mix)
    tmVaR.con[i,1:length(tmp2$e.out.con)] <- log(tmp2$e.out.con)
    tmVaR.akelly[i,1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
    tmVaR.rkelly[i,1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
    tmVaR.mix[i,1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
    
    
    tmp2=evalue.VaR(x=y, r=VaRout2[i,],lev=avec[2])
    tmp3=evalue.VaR(x=y, r=VaRout2[i,] * (1-err),lev=avec[2])
    tmp4=evalue.VaR(x=y, r=VaRout2[i,] * (1+err),lev=avec[2])
    eVaR.con[nm+i,] = log(c(tmp2$e.out.con[length(tmp2$e.out.con)], tmp3$e.out.con[length(tmp3$e.out.con)],
                            tmp4$e.out.con[length(tmp4$e.out.con)]))
    eVaR.akelly[nm+i,] = log(c(tmp2$e.out.akelly[length(tmp2$e.out.akelly)], tmp3$e.out.akelly[length(tmp3$e.out.akelly)],
                               tmp4$e.out.akelly[length(tmp4$e.out.akelly)]))
    eVaR.rkelly[nm+i,] = log(c(tmp2$e.out.rkelly[length(tmp2$e.out.rkelly)], tmp3$e.out.rkelly[length(tmp3$e.out.rkelly)],
                               tmp4$e.out.rkelly[length(tmp4$e.out.rkelly)]))
    eVaR.mix[nm+i,] = log(c(tmp2$e.out.mix[length(tmp2$e.out.mix)], tmp3$e.out.mix[length(tmp3$e.out.mix)],
                               tmp4$e.out.mix[length(tmp4$e.out.mix)]))
    rejVaR1.con[nm+i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp4$n.rej1.con)
    rejVaR1.akelly[nm+i,] = cbind(tmp2$n.rej1.akelly, tmp3$n.rej1.akelly, tmp4$n.rej1.akelly)
    rejVaR1.rkelly[nm+i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp4$n.rej1.rkelly)
    rejVaR1.mix[nm+i,] = cbind(tmp2$n.rej1.mix, tmp3$n.rej1.mix, tmp4$n.rej1.mix)
    rejVaR2.con[nm+i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp4$n.rej2.con)
    rejVaR2.akelly[nm+i,] = cbind(tmp2$n.rej2.akelly, tmp3$n.rej2.akelly, tmp4$n.rej2.akelly)
    rejVaR2.rkelly[nm+i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp4$n.rej2.rkelly)
    rejVaR2.mix[nm+i,] = cbind(tmp2$n.rej2.mix, tmp3$n.rej2.mix, tmp4$n.rej2.mix)
    rejVaR3.con[nm+i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp4$n.rej3.con)
    rejVaR3.akelly[nm+i,] = cbind(tmp2$n.rej3.akelly, tmp3$n.rej3.akelly, tmp4$n.rej3.akelly)
    rejVaR3.rkelly[nm+i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp4$n.rej3.rkelly)
    rejVaR3.mix[nm+i,] = cbind(tmp2$n.rej3.mix, tmp3$n.rej3.mix, tmp4$n.rej3.mix)
    tmVaR.con[nm+i,1:length(tmp2$e.out.con)] <- log(tmp2$e.out.con)
    tmVaR.akelly[nm+i,1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
    tmVaR.rkelly[nm+i,1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
    tmVaR.mix[nm+i,1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
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
  fcES = rbind(cbind(rowMeans(VaRout1b), rowMeans(ESout1)), cbind(rowMeans(VaRout2b), rowMeans(ESout2))) # VaR-ES forecast
  eES.con = matrix(nrow=2*nm, ncol=5) 
  eES.akelly = matrix(nrow=2*nm, ncol=5)
  eES.rkelly = matrix(nrow=2*nm, ncol=5)
  eES.mix = matrix(nrow=2*nm, ncol=5) # final e-values for ES
  rejES1.con = matrix(nrow=2*nm, ncol=5)
  rejES1.akelly = matrix(nrow=2*nm, ncol=5)
  rejES1.rkelly = matrix(nrow=2*nm, ncol=5)
  rejES1.mix = matrix(nrow=2*nm, ncol=5) # numbers of days where we reject for ES (threshold 2)
  rejES2.con = matrix(nrow=2*nm, ncol=5)
  rejES2.akelly = matrix(nrow=2*nm, ncol=5)
  rejES2.rkelly = matrix(nrow=2*nm, ncol=5)
  rejES2.mix = matrix(nrow=2*nm, ncol=5) # numbers of days where we reject for ES (threshold 5)
  rejES3.con = matrix(nrow=2*nm, ncol=5)
  rejES3.akelly = matrix(nrow=2*nm, ncol=5)
  rejES3.rkelly = matrix(nrow=2*nm, ncol=5)
  rejES3.mix = matrix(nrow=2*nm, ncol=5) # numbers of days where we reject for ES (threshold 10)
  tmES.con = matrix(nrow=2*nm, ncol=(length(y)-w))
  tmES.akelly = matrix(nrow=2*nm, ncol=(length(y)-w))
  tmES.rkelly = matrix(nrow=2*nm, ncol=(length(y)-w))
  tmES.mix = matrix(nrow=2*nm, ncol=(length(y)-w)) # sequential test martingale for ES
  
  for(i in 1:nm)
  {
    tmp2=evalue.ES(x=y, r=ESout1[i,], z=VaRout1b[i,],lev=nvec[1]) # sequential test margtingale
    tmp3=evalue.ES(x=y, r=ESout1[i,] * (1-err) + err * ESout1[i,] * ((ESout1[i,] * (1-err)) <= VaRout1b[i,]), z=VaRout1b[i,], lev=nvec[1]) # sequential test margtingale for under-reporting ES
    tmp4=evalue.ES(x=y, r=ESout1[i,] * (1+err), z=VaRout1b[i,],lev=nvec[1]) # sequential test margtingale for over-reporting ES
    tmp5=evalue.ES(x=y, r=ESout1[i,] * (1-err), z=VaRout1b[i,] * (1-err), lev=nvec[1]) # sequential test margtingale for under-reporting both VaR and ES
    tmp6=evalue.ES(x=y, r=ESout1[i,] * (1+err), z=VaRout1b[i,] * (1+err),lev=nvec[1]) # sequential test margtingale for over-reporting both VaR and ES
    eES.con[i,] = log(c(tmp2$e.out.con[length(tmp2$e.out.con)], tmp3$e.out.con[length(tmp3$e.out.con)],
                        tmp4$e.out.con[length(tmp4$e.out.con)], tmp5$e.out.con[length(tmp5$e.out.con)],
                        tmp6$e.out.con[length(tmp6$e.out.con)]))
    eES.akelly[i,] = log(c(tmp2$e.out.akelly[length(tmp2$e.out.akelly)], tmp3$e.out.akelly[length(tmp3$e.out.akelly)],
                           tmp4$e.out.akelly[length(tmp4$e.out.akelly)], tmp5$e.out.akelly[length(tmp5$e.out.akelly)],
                           tmp6$e.out.akelly[length(tmp6$e.out.akelly)]))
    eES.rkelly[i,] = log(c(tmp2$e.out.rkelly[length(tmp2$e.out.rkelly)], tmp3$e.out.rkelly[length(tmp3$e.out.rkelly)],
                           tmp4$e.out.rkelly[length(tmp4$e.out.rkelly)], tmp5$e.out.rkelly[length(tmp5$e.out.rkelly)],
                           tmp6$e.out.rkelly[length(tmp6$e.out.rkelly)]))
    eES.mix[i,] = log(c(tmp2$e.out.mix[length(tmp2$e.out.mix)], tmp3$e.out.mix[length(tmp3$e.out.mix)],
                           tmp4$e.out.mix[length(tmp4$e.out.mix)], tmp5$e.out.mix[length(tmp5$e.out.mix)],
                           tmp6$e.out.mix[length(tmp6$e.out.mix)]))
    rejES1.con[i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp4$n.rej1.con, tmp5$n.rej1.con, tmp6$n.rej1.con)
    rejES1.akelly[i,] = cbind(tmp2$n.rej1.akelly, tmp3$n.rej1.akelly, tmp4$n.rej1.akelly, tmp5$n.rej1.akelly, tmp6$n.rej1.akelly)
    rejES1.rkelly[i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp4$n.rej1.rkelly, tmp5$n.rej1.rkelly, tmp6$n.rej1.rkelly)
    rejES1.mix[i,] = cbind(tmp2$n.rej1.mix, tmp3$n.rej1.mix, tmp4$n.rej1.mix, tmp5$n.rej1.mix, tmp6$n.rej1.mix)
    rejES2.con[i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp4$n.rej2.con, tmp5$n.rej2.con, tmp6$n.rej2.con)
    rejES2.akelly[i,] = cbind(tmp2$n.rej2.akelly, tmp3$n.rej2.akelly, tmp4$n.rej2.akelly, tmp5$n.rej2.akelly, tmp6$n.rej2.akelly)
    rejES2.rkelly[i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp4$n.rej2.rkelly, tmp5$n.rej2.rkelly, tmp6$n.rej2.rkelly)
    rejES2.mix[i,] = cbind(tmp2$n.rej2.mix, tmp3$n.rej2.mix, tmp4$n.rej2.mix, tmp5$n.rej2.mix, tmp6$n.rej2.mix)
    rejES3.con[i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp4$n.rej3.con, tmp5$n.rej3.con, tmp6$n.rej3.con)
    rejES3.akelly[i,] = cbind(tmp2$n.rej3.akelly, tmp3$n.rej3.akelly, tmp4$n.rej3.akelly, tmp5$n.rej3.akelly, tmp6$n.rej3.akelly)
    rejES3.rkelly[i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp4$n.rej3.rkelly, tmp5$n.rej3.rkelly, tmp6$n.rej3.rkelly)
    rejES3.mix[i,] = cbind(tmp2$n.rej3.mix, tmp3$n.rej3.mix, tmp4$n.rej3.mix, tmp5$n.rej3.mix, tmp6$n.rej3.mix)
    tmES.con[i,1:length(tmp2$e.out.con)] <- log(tmp2$e.out.con)
    tmES.akelly[i,1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
    tmES.rkelly[i,1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
    tmES.mix[i,1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
    
    
    tmp2=evalue.ES(x=y,  r=ESout2[i,], z=VaRout2b[i,],lev=nvec[2])
    tmp3=evalue.ES(x=y, r=ESout2[i,] * (1-err) + err * ESout2[i,] * ((ESout2[i,] * (1-err)) <= VaRout2b[i,]), z=VaRout2b[i,],lev=nvec[2])
    tmp4=evalue.ES(x=y, r=ESout2[i,] * (1+err), z=VaRout2b[i,],lev=nvec[2])
    tmp5=evalue.ES(x=y, r=ESout2[i,] * (1-err), z=VaRout2b[i,] * (1-err),lev=nvec[2])
    tmp6=evalue.ES(x=y, r=ESout2[i,] * (1+err), z=VaRout2b[i,] * (1+err),lev=nvec[2])
    eES.con[nm+i,] = log(c(tmp2$e.out.con[length(tmp2$e.out.con)], tmp3$e.out.con[length(tmp2$e.out.con)],
                           tmp4$e.out.con[length(tmp2$e.out.con)], tmp5$e.out.con[length(tmp2$e.out.con)],
                           tmp6$e.out.con[length(tmp2$e.out.con)]))
    eES.akelly[nm+i,] = log(c(tmp2$e.out.akelly[length(tmp2$e.out.akelly)], tmp3$e.out.akelly[length(tmp2$e.out.akelly)],
                              tmp4$e.out.akelly[length(tmp2$e.out.akelly)], tmp5$e.out.akelly[length(tmp2$e.out.akelly)],
                              tmp6$e.out.akelly[length(tmp2$e.out.akelly)]))
    eES.rkelly[nm+i,] = log(c(tmp2$e.out.rkelly[length(tmp2$e.out.rkelly)], tmp3$e.out.rkelly[length(tmp3$e.out.rkelly)],
                              tmp4$e.out.rkelly[length(tmp4$e.out.rkelly)], tmp5$e.out.rkelly[length(tmp5$e.out.rkelly)],
                              tmp6$e.out.rkelly[length(tmp6$e.out.rkelly)]))
    eES.mix[nm+i,] = log(c(tmp2$e.out.mix[length(tmp2$e.out.mix)], tmp3$e.out.mix[length(tmp3$e.out.mix)],
                              tmp4$e.out.mix[length(tmp4$e.out.mix)], tmp5$e.out.mix[length(tmp5$e.out.mix)],
                              tmp6$e.out.mix[length(tmp6$e.out.mix)]))
    rejES1.con[nm+i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp4$n.rej1.con, tmp5$n.rej1.con, tmp6$n.rej1.con)
    rejES1.akelly[nm+i,] = cbind(tmp2$n.rej1.akelly, tmp3$n.rej1.akelly, tmp4$n.rej1.akelly, tmp5$n.rej1.akelly, tmp6$n.rej1.akelly)
    rejES1.rkelly[nm+i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp4$n.rej1.rkelly, tmp5$n.rej1.rkelly, tmp6$n.rej1.rkelly)
    rejES1.mix[nm+i,] = cbind(tmp2$n.rej1.mix, tmp3$n.rej1.mix, tmp4$n.rej1.mix, tmp5$n.rej1.mix, tmp6$n.rej1.mix)
    rejES2.con[nm+i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp4$n.rej2.con, tmp5$n.rej2.con, tmp6$n.rej2.con)
    rejES2.akelly[nm+i,] = cbind(tmp2$n.rej2.akelly, tmp3$n.rej2.akelly, tmp4$n.rej2.akelly, tmp5$n.rej2.akelly, tmp6$n.rej2.akelly)
    rejES2.rkelly[nm+i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp4$n.rej2.rkelly, tmp5$n.rej2.rkelly, tmp6$n.rej2.rkelly)
    rejES2.mix[nm+i,] = cbind(tmp2$n.rej2.mix, tmp3$n.rej2.mix, tmp4$n.rej2.mix, tmp5$n.rej2.mix, tmp6$n.rej2.mix)
    rejES3.con[nm+i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp4$n.rej3.con, tmp5$n.rej3.con, tmp6$n.rej3.con)
    rejES3.akelly[nm+i,] = cbind(tmp2$n.rej3.akelly, tmp3$n.rej3.akelly, tmp4$n.rej3.akelly, tmp5$n.rej3.akelly, tmp6$n.rej3.akelly)
    rejES3.rkelly[nm+i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp4$n.rej3.rkelly, tmp5$n.rej3.rkelly, tmp6$n.rej3.rkelly)
    rejES3.mix[nm+i,] = cbind(tmp2$n.rej3.mix, tmp3$n.rej3.mix, tmp4$n.rej3.mix, tmp5$n.rej3.mix, tmp6$n.rej3.mix)
    tmES.con[nm+i,1:length(tmp2$e.out.con)] <- log(tmp2$e.out.con)
    tmES.akelly[nm+i,1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
    tmES.rkelly[nm+i,1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
    tmES.mix[nm+i,1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
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
  
  return(c(fcVaR, eVaR.con, eVaR.akelly, eVaR.rkelly, eVaR.mix, rejVaR1.con, rejVaR1.akelly, rejVaR1.rkelly, rejVaR1.mix,
           rejVaR2.con, rejVaR2.akelly, rejVaR2.rkelly, rejVaR2.mix, rejVaR3.con, rejVaR3.akelly, rejVaR3.rkelly, rejVaR3.mix,
           rejnumVaR1.con, rejnumVaR1.akelly, rejnumVaR1.rkelly, rejnumVaR1.mix, rejnumVaR2.con, rejnumVaR2.akelly, rejnumVaR2.rkelly, rejnumVaR2.mix,
           rejnumVaR3.con, rejnumVaR3.akelly, rejnumVaR3.rkelly, rejnumVaR3.mix, tmVaR.con, tmVaR.akelly, tmVaR.rkelly, tmVaR.mix,
           fcES, eES.con, eES.akelly, eES.rkelly, eES.mix, rejES1.con, rejES1.akelly, rejES1.rkelly, rejES1.mix,
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
avg.total = c()
for(i in 1:(1368+48*n)){
  avg = mean(out[,i][out[,i] != n])
  avg.total = cbind(avg.total, avg)
}
avg.fcVaR = matrix(avg.total[1:8],nrow=8,ncol=1)  # VaR forecast
avg.eVaR.con = matrix(avg.total[9:26],nrow=6,ncol=3)
avg.eVaR.akelly = matrix(avg.total[27:44],nrow=6,ncol=3)
avg.eVaR.rkelly = matrix(avg.total[45:62],nrow=6,ncol=3)
avg.eVaR.mix = matrix(avg.total[63:80],nrow=6,ncol=3) # final e-values for VaR
avg.rejVaR1.con = matrix(avg.total[81:98],nrow=6,ncol=3)
avg.rejVaR1.akelly = matrix(avg.total[99:116],nrow=6,ncol=3)
avg.rejVaR1.rkelly = matrix(avg.total[117:134],nrow=6,ncol=3)
avg.rejVaR1.mix = matrix(avg.total[135:152],nrow=6,ncol=3) # numbers of days where we reject for VaR (threshold 2)
avg.rejVaR2.con = matrix(avg.total[153:170],nrow=6,ncol=3)
avg.rejVaR2.akelly = matrix(avg.total[171:188],nrow=6,ncol=3)
avg.rejVaR2.rkelly = matrix(avg.total[189:206],nrow=6,ncol=3)
avg.rejVaR2.mix = matrix(avg.total[207:224],nrow=6,ncol=3) # numbers of days where we reject for VaR (threshold 5)
avg.rejVaR3.con = matrix(avg.total[225:242],nrow=6,ncol=3)
avg.rejVaR3.akelly = matrix(avg.total[243:260],nrow=6,ncol=3)
avg.rejVaR3.rkelly = matrix(avg.total[261:278],nrow=6,ncol=3)
avg.rejVaR3.mix = matrix(avg.total[279:296],nrow=6,ncol=3) # numbers of days where we reject for VaR (threshold 10)
avg.rejnumVaR1.con = matrix(avg.total[297:314],nrow=6,ncol=3)
avg.rejnumVaR1.akelly = matrix(avg.total[315:332],nrow=6,ncol=3)
avg.rejnumVaR1.rkelly = matrix(avg.total[333:350],nrow=6,ncol=3)
avg.rejnumVaR1.mix = matrix(avg.total[351:368],nrow=6,ncol=3) # numbers of rejections for VaR (threshold 2)
avg.rejnumVaR2.con = matrix(avg.total[369:386],nrow=6,ncol=3)
avg.rejnumVaR2.akelly = matrix(avg.total[387:404],nrow=6,ncol=3)
avg.rejnumVaR2.rkelly = matrix(avg.total[405:422],nrow=6,ncol=3)
avg.rejnumVaR2.mix = matrix(avg.total[423:440],nrow=6,ncol=3) # numbers of rejections for VaR (threshold 5)
avg.rejnumVaR3.con = matrix(avg.total[441:458],nrow=6,ncol=3)
avg.rejnumVaR3.akelly = matrix(avg.total[459:476],nrow=6,ncol=3)
avg.rejnumVaR3.rkelly = matrix(avg.total[477:494],nrow=6,ncol=3)
avg.rejnumVaR3.mix = matrix(avg.total[495:512],nrow=6,ncol=3) # numbers of rejections for VaR (threshold 10)
avg.tmVaR.con = matrix(avg.total[513:(512+6*n)],nrow=6,ncol=n)
avg.tmVaR.akelly = matrix(avg.total[(513+6*n):(512+12*n)],nrow=6,ncol=n)
avg.tmVaR.rkelly = matrix(avg.total[(513+12*n):(512+18*n)],nrow=6,ncol=n)
avg.tmVaR.mix = matrix(avg.total[(513+18*n):(512+24*n)],nrow=6,ncol=n) # sequential test martingale for VaR
avg.fcES = matrix(avg.total[(513+24*n):(528+24*n)],nrow=8,ncol=2) # VaR-ES forecast
avg.eES.con = matrix(avg.total[(529+24*n):(558+24*n)],nrow=6,ncol=5)
avg.eES.akelly = matrix(avg.total[(559+24*n):(588+24*n)],nrow=6,ncol=5)
avg.eES.rkelly = matrix(avg.total[(589+24*n):(618+24*n)],nrow=6,ncol=5)
avg.eES.mix = matrix(avg.total[(619+24*n):(648+24*n)],nrow=6,ncol=5) # final e-values for ES
avg.rejES1.con = matrix(avg.total[(649+24*n):(678+24*n)],nrow=6,ncol=5)
avg.rejES1.akelly = matrix(avg.total[(679+24*n):(708+24*n)],nrow=6,ncol=5)
avg.rejES1.rkelly = matrix(avg.total[(709+24*n):(738+24*n)],nrow=6,ncol=5)
avg.rejES1.mix = matrix(avg.total[(739+24*n):(768+24*n)],nrow=6,ncol=5) # numbers of days where we reject for ES (threshold 2)
avg.rejES2.con = matrix(avg.total[(769+24*n):(798+24*n)],nrow=6,ncol=5)
avg.rejES2.akelly = matrix(avg.total[(799+24*n):(828+24*n)],nrow=6,ncol=5)
avg.rejES2.rkelly = matrix(avg.total[(829+24*n):(858+24*n)],nrow=6,ncol=5)
avg.rejES2.mix = matrix(avg.total[(859+24*n):(888+24*n)],nrow=6,ncol=5) # numbers of days where we reject for ES (threshold 5)
avg.rejES3.con = matrix(avg.total[(889+24*n):(918+24*n)],nrow=6,ncol=5)
avg.rejES3.akelly = matrix(avg.total[(919+24*n):(948+24*n)],nrow=6,ncol=5)
avg.rejES3.rkelly = matrix(avg.total[(949+24*n):(978+24*n)],nrow=6,ncol=5)
avg.rejES3.mix = matrix(avg.total[(979+24*n):(1008+24*n)],nrow=6,ncol=5) # numbers of days where we reject for ES (threshold 10)
avg.rejnumES1.con = matrix(avg.total[(1009+24*n):(1038+24*n)],nrow=6,ncol=5)
avg.rejnumES1.akelly = matrix(avg.total[(1039+24*n):(1068+24*n)],nrow=6,ncol=5)
avg.rejnumES1.rkelly = matrix(avg.total[(1069+24*n):(1098+24*n)],nrow=6,ncol=5)
avg.rejnumES1.mix = matrix(avg.total[(1099+24*n):(1128+24*n)],nrow=6,ncol=5) # numbers of rejections for ES (threshold 2)
avg.rejnumES2.con = matrix(avg.total[(1129+24*n):(1158+24*n)],nrow=6,ncol=5)
avg.rejnumES2.akelly = matrix(avg.total[(1159+24*n):(1188+24*n)],nrow=6,ncol=5)
avg.rejnumES2.rkelly = matrix(avg.total[(1189+24*n):(1218+24*n)],nrow=6,ncol=5)
avg.rejnumES2.mix = matrix(avg.total[(1219+24*n):(1248+24*n)],nrow=6,ncol=5) # numbers of rejections for ES (threshold 5)
avg.rejnumES3.con = matrix(avg.total[(1249+24*n):(1278+24*n)],nrow=6,ncol=5)
avg.rejnumES3.akelly = matrix(avg.total[(1279+24*n):(1308+24*n)],nrow=6,ncol=5)
avg.rejnumES3.rkelly = matrix(avg.total[(1309+24*n):(1338+24*n)],nrow=6,ncol=5)
avg.rejnumES3.mix = matrix(avg.total[(1339+24*n):(1368+24*n)],nrow=6,ncol=5) # numbers of rejections for ES (threshold 10)
avg.tmES.con = matrix(avg.total[(1369+24*n):(1368+30*n)],nrow=6,ncol=n)
avg.tmES.akelly = matrix(avg.total[(1369+30*n):(1368+36*n)],nrow=6,ncol=n)
avg.tmES.rkelly = matrix(avg.total[(1369+36*n):(1368+42*n)],nrow=6,ncol=n)
avg.tmES.mix = matrix(avg.total[(1369+42*n):(1368+48*n)],nrow=6,ncol=n) # sequential test martingale for ES


out.final = list(avg.fcVaR = avg.fcVaR, avg.eVaR.con = avg.eVaR.con, avg.eVaR.akelly = avg.eVaR.akelly, avg.eVaR.rkelly = avg.eVaR.rkelly, avg.eVaR.mix = avg.eVaR.mix,
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
save(out.final, file="Result_simulation_final.RDATA")


avg.fcVaR
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
avg.fcES
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
load("Result_simulation_final.RDATA")

setEPS()
postscript("VaR95con.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.con[1,1:500], out.final$avg.tmVaR.con[2,1:500], out.final$avg.tmVaR.con[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("VaR99con.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.con[4,1:500], out.final$avg.tmVaR.con[5,1:500], out.final$avg.tmVaR.con[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"),  type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()

setEPS()
postscript("VaR95akelly.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.akelly[1,1:500], out.final$avg.tmVaR.akelly[2,1:500], out.final$avg.tmVaR.akelly[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("VaR99akelly.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.akelly[4,1:500], out.final$avg.tmVaR.akelly[5,1:500], out.final$avg.tmVaR.akelly[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()

setEPS()
postscript("VaR95rkelly.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.rkelly[1,1:500], out.final$avg.tmVaR.rkelly[2,1:500], out.final$avg.tmVaR.rkelly[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("VaR99rkelly.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.rkelly[4,1:500], out.final$avg.tmVaR.rkelly[5,1:500], out.final$avg.tmVaR.rkelly[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()

setEPS()
postscript("VaR95mix.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.mix[1,1:500], out.final$avg.tmVaR.mix[2,1:500], out.final$avg.tmVaR.mix[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("VaR99mix.eps")
matplot(1:500, cbind(out.final$avg.tmVaR.mix[4,1:500], out.final$avg.tmVaR.mix[5,1:500], out.final$avg.tmVaR.mix[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()



setEPS()
postscript("ES875con.eps")
matplot(1:500, cbind(out.final$avg.tmES.con[1,1:500], out.final$avg.tmES.con[2,1:500], out.final$avg.tmES.con[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975con.eps")
matplot(1:500, cbind(out.final$avg.tmES.con[4,1:500], out.final$avg.tmES.con[5,1:500], out.final$avg.tmES.con[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"),  type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()

setEPS()
postscript("ES875akelly.eps")
matplot(1:500, cbind(out.final$avg.tmES.akelly[1,1:500], out.final$avg.tmES.akelly[2,1:500], out.final$avg.tmES.akelly[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975akelly.eps")
matplot(1:500, cbind(out.final$avg.tmES.akelly[4,1:500], out.final$avg.tmES.akelly[5,1:500], out.final$avg.tmES.akelly[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()


setEPS()
postscript("ES875rkelly.eps")
matplot(1:500, cbind(out.final$avg.tmES.rkelly[1,1:500], out.final$avg.tmES.rkelly[2,1:500], out.final$avg.tmES.rkelly[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975rkelly.eps")
matplot(1:500, cbind(out.final$avg.tmES.rkelly[4,1:500], out.final$avg.tmES.rkelly[5,1:500], out.final$avg.tmES.rkelly[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()

setEPS()
postscript("ES875mix.eps")
matplot(1:500, cbind(out.final$avg.tmES.mix[1,1:500], out.final$avg.tmES.mix[2,1:500], out.final$avg.tmES.mix[3,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975mix.eps")
matplot(1:500, cbind(out.final$avg.tmES.mix[4,1:500], out.final$avg.tmES.mix[5,1:500], out.final$avg.tmES.mix[6,1:500]),
        ylim=e.lim, col=c("blue", "red", "black"), type = "l", lty = 1,
        xlab = "number of days", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = c(1,1,1))
legend("topleft", legend=c("normal", "t", "skewed-t"),
       col=c("blue", "red", "black"), lwd = 1, cex=1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
dev.off()
