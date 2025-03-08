# install.packages("rugarch")
# install.packages("fGarch")
# install.packages("MASS")
# install.packages("expectreg")
# install.packages("ismev")
# install.packages("lmom")
# install.packages("QRM")
# install.packages("skewt")
# install.packages("spam")
# install.packages("plotrix")


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



# NASDAQ real data
dat = rev(read.csv("nasdaq2021.csv")$Close) # NASDAQ Composite index (Jan 16,1996 - Dec 31,2021)
dates = as.Date(rev(read.csv("nasdaq2021.csv")$Date))[-1] # dates for returns
x = - log(dat[-1]/dat[-length(dat)])*100 # negated percentage log returns


N=length(x) # sample size
w=500 # moving window size (same as in simulations)
(n=N-w) # out-of-sample size


nvec=.975 # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- nvec # VaR levels for VaR
inu = 1:(length(nvec)) #index set for nu levels of VaR


# =======================================================
# Normal innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1) # VaR values
esmodel = (ddist("norm", qdist("norm",p=nvec,mu=0,sigma=1), mu=0, sigma=1)) / (1 - nvec)

VaR <- matrix(nrow=n,ncol=length(VaR.levels))

ES <- matrix(nrow=n,ncol=length(nvec))


# estimated parameters
fit.par <- matrix(nrow=n, ncol=5) # 5 model parameters
mut  <- vector(mode="numeric", length=n)
sigt  <- vector(mode="numeric", length=n)

# ----------------------------------------------------

for(i in 1:n)
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)

  VaR[i]=foc$VaRmodel

  ES[i]=foc$ESmodel
  mut[i] = foc$mut
  sigt[i] = foc$sigt

  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,ES=ES,par=fit.par,mut=mut,sigt=sigt)
save(out, file="Sim3norm3.RDATA")

# =======================================================
# Student t innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="std")

VaR <- matrix(nrow=n,ncol=length(VaR.levels))

ES <- matrix(nrow=n,ncol=length(nvec))
fit.par <- matrix(nrow=n, ncol=6) # 6 model parameters
mut  <- vector(mode="numeric", length=n)
sigt  <- vector(mode="numeric", length=n)


# ----------------------------------------------------
for(i in 1:n)
{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)
  nu=coef(fit)["shape"]

  qmodel = qdist("std",p=VaR.levels,mu=0,sigma=1,shape=nu)

  # expected shortfall for the assumed model/distribution
  esmodel = (dt(qt(nvec, df=nu), df=nu)) / (1-nvec) * (nu + (qt(nvec, df=nu))^2) / (nu-1)
  esmodel = sqrt((nu-2)/nu)*esmodel

  foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)

  VaR[i]=foc$VaRmodel

  ES[i]=foc$ESmodel
  mut[i] = foc$mut
  sigt[i] = foc$sigt

  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,ES=ES,par=fit.par,mut=mut,sigt=sigt)

save(out, file="Sim3std3.RDATA")


# =======================================================
# skewed Student t innovations
# The version of Fernandez & Steel (1998)
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")

VaR <- matrix(nrow=n,ncol=length(VaR.levels))

ES <- matrix(nrow=n,ncol=length(nvec))
fit.par <- matrix(nrow=n, ncol=7) # 7 model parameters
mut  <- vector(mode="numeric", length=n)
sigt  <- vector(mode="numeric", length=n)

for(i in 1:n)

{
  fit = ugarchfit(spec, x[i:(i+w-1)], solver="hybrid")
  fit.par[i,] = coef(fit)

  # mean and std dev of skewed t rv
  nu=coef(fit)["shape"]; ga=coef(fit)["skew"]

  # m = mean.st(shape=nu,skew=ga)
  # s = sqrt(var.st(shape=nu,skew=ga))

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

  # emodel = (est(asy=tvec,shape=nu,skew=ga)-m)/s

  foc=RM.forecasts3(fit=fit,qmodel=qmodel,esmodel=esmodel)

  VaR[i]=foc$VaRmodel

  ES[i]=foc$ESmodel
  mut[i] = foc$mut
  sigt[i] = foc$sigt

  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,ES=ES,par=fit.par,mut=mut,sigt=sigt)

save(out, file="Sim3sstd3.RDATA")


# # =======================================================
# # Empirical estimations
# # =======================================================
# 
# VaR <- matrix(nrow=n,ncol=length(VaR.levels))
# 
# ES <- matrix(nrow=n,ncol=length(nvec))
# 
# for(i in 1:n)
# {
#   VaR[i] = quantile(x[i:(i+w-1)], probs=VaR.levels)
#   ES[i] = mean(x[i:(i+w-1)][x[i:(i+w-1)] >= VaR[i]])
# 
#   if(i %% 100 == 0) print(i)
# }
# 
# out=list(VaR=VaR,ES=ES)
# save(out, file="Empirical.RDATA")


# NASDAQ real data
dates = as.Date(rev(read.csv("nasdaq2021.csv")$Date))[-1] # dates for returns
N=length(x) # sample size

# Length of data you want to choose (choose one and comment out others)

# m=N-2*w+2 # NASDAQ Composite index (Jan 3,2000 - Dec 31,2021)
m=4781-w # NASDAQ Composite index (Jan 3,2005 - Dec 31,2021)
# m=3522-w # NASDAQ Composite index (Jan 4,2010 - Dec 31,2021)
# m=2060-w # NASDAQ Composite index (Jan 2,2015 - Dec 31,2021)


n=m-2
y=tail(x,n)
dates=tail(dates,n)
dates = as.Date(dates[-c(1:(w))])


# --------------------------------------------------------
# Forecasts

#=====================================================
# (VaR_nu, ES_nu)
#=====================================================

load("Sim3norm3.RDATA")
VaRoutb = out$VaR
ESout = out$ES
sigt = out$sigt

load("Sim3std3.RDATA")
VaRoutb = cbind(VaRoutb,out$VaR)
ESout = cbind(ESout,out$ES)
sigt = cbind(sigt, out$sigt)

load("Sim3sstd3.RDATA")
VaRoutb = cbind(VaRoutb,out$VaR)
ESout = cbind(ESout,out$ES)
sigt = cbind(sigt, out$sigt)

# load("Empirical.RDATA")
# VaRoutb = rbind(VaRoutb,out$VaR)
# ESout = rbind(ESout,out$ES)
# sigt = rbind(sigt, out$sigt)



VaRoutb = t(VaRoutb)[1:3,(N-w-n+1):(N-w)]; ESout = t(ESout)[1:3,(N-w-n+1):(N-w)]
sigt = t(sigt)[1:3,(N-w-n+1):(N-w)]



# ===================================================
# E-backtesting analysis of the simulation data
# Data generating process: AR(1)-GARCH(1,1) with skewed t innovations
# "con" represents e-backtesting method with constant lambda,
# "akelly" represents GREE method,
# "rkelly" represents GREL method,
# "mix" represents GREM method
# ===================================================

nm = 3
err <- .1 # error for under- and over-reporting
e.lim <- c(-1,5) # bounds for e-values
w <- 500 # time window


# ===================================================
# Summary for ES
pES = matrix(nrow=nm, ncol=2) # final p-values for ES
pprocES = matrix(nrow=2*nm+2, ncol=n) # p-process
rejES1 = matrix(nrow=nm, ncol=2) # numbers of days where we reject for ES (threshold 0.5)
rejES2 = matrix(nrow=nm, ncol=2) # numbers of days where we reject for ES (threshold 0.2)
rejES3 = matrix(nrow=nm, ncol=2) # numbers of days where we reject for ES (threshold 0.1)



for(i in 1:nm)
{
  tmp=cct.onesided2(x=y, r=cbind(VaRoutb[i,], ESout[i,]),lev=nvec,sigt=sigt[i,])  # sequential p-values
  # tmp3=cct.onesided2(x=y, r=cbind(VaRoutb[i,], ESout[i,]*(1-err)),lev=nvec[i],sigt=sigt) # sequential p-values for under-reporting ES
  tmp5=cct.onesided2(x=y, r=cbind(VaRoutb[i,], ESout[i,]*(1+err)),lev=nvec,sigt=sigt[i,]) # sequential p-values for over-reporting ES
  pES[i,] = c(tmp$out.pv.cond[length(tmp$out.pv.cond)], tmp5$out.pv.cond[length(tmp5$out.pv.cond)])
  pprocES[i,1:length(tmp$out.pv.cond)] = tmp$out.pv.cond
  pprocES[nm+i,1:length(tmp5$out.pv.cond)] = tmp5$out.pv.cond
  rejES1[i,] = cbind(tmp$n.rej1, tmp5$n.rej1)
  rejES2[i,] = cbind(tmp$n.rej2, tmp5$n.rej2)
  rejES3[i,] = cbind(tmp$n.rej3, tmp5$n.rej3)
}



pES
rejES1
rejES2
rejES3


# plots for data from different years (choose one and comment out others)


# plots from 2005

setEPS()
postscript("traditional_empirical2005.eps")
plot(dates[1:3000],pprocES[1,1:3000],type="l",xlab="dates", ylab="p-values", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=c(-0.1,1.1))
lines(dates[1:3000],pprocES[2,1:3000],type="l",col="green")
lines(dates[1:3000],pprocES[3,1:3000],type="l",col="blue")
lines(dates[1:3000],pprocES[4,1:3000],type="l",col="black")
abline(h = 0.5, lwd = 1, lty = 2)
abline(h = 0.2, lwd = 1, lty = 2)
abline(h = 0.1, lwd = 1, lty = 2)
legend("topright", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("orange", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()