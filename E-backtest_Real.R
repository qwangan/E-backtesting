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


avec=c(.95, .99) # vector of alpha levels for VaR
nvec=c(.875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR


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

  VaR[i,]=foc$VaRmodel

  ES[i,]=foc$ESmodel
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

  VaR[i,]=foc$VaRmodel

  ES[i,]=foc$ESmodel
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

  VaR[i,]=foc$VaRmodel

  ES[i,]=foc$ESmodel
  mut[i] = foc$mut
  sigt[i] = foc$sigt

  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,ES=ES,par=fit.par,mut=mut,sigt=sigt)

save(out, file="Sim3sstd3.RDATA")


# =======================================================
# Empirical estimations
# =======================================================

VaR <- matrix(nrow=n,ncol=length(VaR.levels))

ES <- matrix(nrow=n,ncol=length(nvec))

for(i in 1:n)
{
  VaR[i,] = quantile(x[i:(i+w-1)], probs=VaR.levels)
  ES[i,] = cbind(mean(x[i:(i+w-1)][x[i:(i+w-1)] >= VaR[i,3]]),mean(x[i:(i+w-1)][x[i:(i+w-1)] >= VaR[i,4]]))

  if(i %% 100 == 0) print(i)
}

out=list(VaR=VaR,ES=ES)
save(out, file="Empirical.RDATA")


# NASDAQ real data
dates = as.Date(rev(read.csv("nasdaq2021.csv")$Date))[-1] # dates for returns
N=length(x) # sample size

# Length of data you want to choose (choose one and comment out others)

m=N-w+2 # NASDAQ Composite index (Jan 3,2000 - Dec 31,2021)
# m=4781 # NASDAQ Composite index (Jan 3,2005 - Dec 31,2021)
# m=3522 # NASDAQ Composite index (Jan 4,2010 - Dec 31,2021)
# m=2060 # NASDAQ Composite index (Jan 2,2015 - Dec 31,2021)


n=m-2
y=tail(x,n)
dates=tail(dates,n)
dates = as.Date(dates[-c(1:(w))])


# --------------------------------------------------------
# Forecasts

#=====================================================
# (VaR_nu, ES_nu)
#=====================================================

out=getVaRES_em(k=1)
VaRout1b = out$VaR; ESout1 = out$ES

out=getVaRES_em(k=2)
VaRout2b = out$VaR; ESout2 = out$ES


VaRout1b = VaRout1b[1:4,(N-w-n+1):(N-w)]; ESout1 = ESout1[1:4,(N-w-n+1):(N-w)]
VaRout2b = VaRout2b[1:4,(N-w-n+1):(N-w)]; ESout2 = ESout2[1:4,(N-w-n+1):(N-w)]



# time series plot
setEPS()
postscript("NASDAQ.eps")
plot(dates[-c(1:(w))],y[-c(1:(w))],type="l",xlab="dates", ylab="negated percentage log returns")
dev.off()


# risk forecasts plot

setEPS()
postscript("NASDAQ_forecast.eps")
plot(dates,ESout2[1,-c(1:(w))],type="l",xlab="dates", ylab="ES forecast",col="orange", ylim=c(1.2,10), lty=1)
lines(dates,ESout2[2,-c(1:(w))],type="l",col=3)
lines(dates,ESout2[3,-c(1:(w))],type="l",col=4, lty = 3)
lines(dates,ESout2[4,-c(1:(w))],type="l",col=2, lwd = 2)
legend("topleft", inset=c(0.12,0), legend=c("normal", "t", "skewed-t", "empirical"),
       col=c("orange",3,4,2), lwd = c(1,1,1,2), lty = c(1,1,3,1), cex=1)
dev.off()



# ===================================================
# E-backtesting analysis of the simulation data
# Data generating process: AR(1)-GARCH(1,1) with skewed t innovations
# "con" represents e-backtesting method with constant lambda,
# "akelly" represents GREE method,
# "rkelly" represents GREL method,
# "mix" represents GREM method
# ===================================================

nm = 4
err <- .1 # error for under- and over-reporting
e.lim <- c(-1,5) # bounds for e-values
w <- 500 # time window


# ===================================================
# Summary for ES
ES.lambda.akelly = matrix(nrow=2*nm, ncol=(length(y)-w)) # lambda for GREE
ES.lambda.rkelly = matrix(nrow=2*nm, ncol=(length(y)-w)) # lambda for GREL
fcES = matrix(nrow=2*nm, ncol=2) # VaR-ES forecast
eES.con = matrix(nrow=2*nm, ncol=5) 
eES.akelly = matrix(nrow=2*nm, ncol=5)
eES.rkelly = matrix(nrow=2*nm, ncol=5)
eES.mix = matrix(nrow=2*nm, ncol=5) # final e-values for ES
tmES.con = matrix(nrow=2*(nm+1), ncol=(length(y)-w))
tmES.akelly = matrix(nrow=2*(nm+1), ncol=(length(y)-w))
tmES.rkelly = matrix(nrow=2*(nm+1), ncol=(length(y)-w))
tmES.mix = matrix(nrow=2*(nm+1), ncol=(length(y)-w)) # sequential test martingale for ES
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

for(i in 1:nm)
{
  tmp2=evalue_ES_em(x=y, r=ESout1[i,], z=VaRout1b[i,],lev=nvec[1], w=500) # sequential test margtingale
  tmp3=evalue_ES_em(x=y, r=ESout1[i,] * (1-err), z=VaRout1b[i,], lev=nvec[1], w=500) # sequential test margtingale for under-reporting ES
  tmp4=evalue_ES_em(x=y, r=ESout1[i,] * (1+err), z=VaRout1b[i,],lev=nvec[1], w=500) # sequential test margtingale for over-reporting ES
  tmp5=evalue_ES_em(x=y, r=ESout1[i,] * (1-err), z=VaRout1b[i,] * (1-err), lev=nvec[1], w=500) # sequential test margtingale for under-reporting both VaR and ES
  tmp6=evalue_ES_em(x=y, r=ESout1[i,] * (1+err), z=VaRout1b[i,] * (1+err),lev=nvec[1], w=500) # sequential test margtingale for over-reporting both VaR and ES
  fcES[i,] = c(mean(VaRout1b[i,]), mean(ESout1[i,]))
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
  tmES.con[i,1:length(tmp2$e.out.con)] <- log(tmp2$e.out.con)
  rejES1.con[i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp4$n.rej1.con, tmp5$n.rej1.con, tmp6$n.rej1.con)
  rejES2.con[i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp4$n.rej2.con, tmp5$n.rej2.con, tmp6$n.rej2.con)
  rejES3.con[i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp4$n.rej3.con, tmp5$n.rej3.con, tmp6$n.rej3.con)
  tmES.akelly[i,1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
  rejES1.akelly[i,] = cbind(tmp2$n.rej1.akelly, tmp3$n.rej1.akelly, tmp4$n.rej1.akelly, tmp5$n.rej1.akelly, tmp6$n.rej1.akelly)
  rejES2.akelly[i,] = cbind(tmp2$n.rej2.akelly, tmp3$n.rej2.akelly, tmp4$n.rej2.akelly, tmp5$n.rej2.akelly, tmp6$n.rej2.akelly)
  rejES3.akelly[i,] = cbind(tmp2$n.rej3.akelly, tmp3$n.rej3.akelly, tmp4$n.rej3.akelly, tmp5$n.rej3.akelly, tmp6$n.rej3.akelly)
  tmES.rkelly[i,1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
  rejES1.rkelly[i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp4$n.rej1.rkelly, tmp5$n.rej1.rkelly, tmp6$n.rej1.rkelly)
  rejES2.rkelly[i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp4$n.rej2.rkelly, tmp5$n.rej2.rkelly, tmp6$n.rej2.rkelly)
  rejES3.rkelly[i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp4$n.rej3.rkelly, tmp5$n.rej3.rkelly, tmp6$n.rej3.rkelly)
  tmES.mix[i,1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
  rejES1.mix[i,] = cbind(tmp2$n.rej1.mix, tmp3$n.rej1.mix, tmp4$n.rej1.mix, tmp5$n.rej1.mix, tmp6$n.rej1.mix)
  rejES2.mix[i,] = cbind(tmp2$n.rej2.mix, tmp3$n.rej2.mix, tmp4$n.rej2.mix, tmp5$n.rej2.mix, tmp6$n.rej2.mix)
  rejES3.mix[i,] = cbind(tmp2$n.rej3.mix, tmp3$n.rej3.mix, tmp4$n.rej3.mix, tmp5$n.rej3.mix, tmp6$n.rej3.mix)
  ES.lambda.akelly[i,1:length(tmp2$out.lambda)] = tmp2$out.lambda
  ES.lambda.rkelly[i,1:length(tmp2$out.lambda.rkelly)] = tmp2$out.lambda.rkelly
  if(i == 3){
    tmES.con[nm+1,1:length(tmp4$e.out.con)] <- log(tmp4$e.out.con)
    tmES.akelly[nm+1,1:length(tmp4$e.out.akelly)] <- log(tmp4$e.out.akelly)
    tmES.rkelly[nm+1,1:length(tmp4$e.out.rkelly)] <- log(tmp4$e.out.rkelly)
    tmES.mix[nm+1,1:length(tmp4$e.out.mix)] <- log(tmp4$e.out.mix)
  }
  
  
  tmp2=evalue_ES_em(x=y,  r=ESout2[i,], z=VaRout2b[i,],lev=nvec[2], w=500)
  tmp3=evalue_ES_em(x=y, r=ESout2[i,] * (1-err), z=VaRout2b[i,],lev=nvec[2], w=500)
  tmp4=evalue_ES_em(x=y, r=ESout2[i,] * (1+err), z=VaRout2b[i,],lev=nvec[2], w=500)
  tmp5=evalue_ES_em(x=y, r=ESout2[i,] * (1-err), z=VaRout2b[i,] * (1-err),lev=nvec[2], w=500)
  tmp6=evalue_ES_em(x=y, r=ESout2[i,] * (1+err), z=VaRout2b[i,] * (1+err),lev=nvec[2], w=500)
  fcES[nm+i,] = c(mean(VaRout2b[i,]), mean(ESout2[i,]))
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
  tmES.con[nm+i+1,1:length(tmp2$e.out.con)] <- log(tmp2$e.out.con)
  rejES1.con[nm+i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp4$n.rej1.con, tmp5$n.rej1.con, tmp6$n.rej1.con)
  rejES2.con[nm+i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp4$n.rej2.con, tmp5$n.rej2.con, tmp6$n.rej2.con)
  rejES3.con[nm+i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp4$n.rej3.con, tmp5$n.rej3.con, tmp6$n.rej3.con)
  tmES.akelly[nm+i+1,1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
  rejES1.akelly[nm+i,] = cbind(tmp2$n.rej1.akelly, tmp3$n.rej1.akelly, tmp4$n.rej1.akelly, tmp5$n.rej1.akelly, tmp6$n.rej1.akelly)
  rejES2.akelly[nm+i,] = cbind(tmp2$n.rej2.akelly, tmp3$n.rej2.akelly, tmp4$n.rej2.akelly, tmp5$n.rej2.akelly, tmp6$n.rej2.akelly)
  rejES3.akelly[nm+i,] = cbind(tmp2$n.rej3.akelly, tmp3$n.rej3.akelly, tmp4$n.rej3.akelly, tmp5$n.rej3.akelly, tmp6$n.rej3.akelly)
  tmES.rkelly[nm+i+1,1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
  rejES1.rkelly[nm+i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp4$n.rej1.rkelly, tmp5$n.rej1.rkelly, tmp6$n.rej1.rkelly)
  rejES2.rkelly[nm+i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp4$n.rej2.rkelly, tmp5$n.rej2.rkelly, tmp6$n.rej2.rkelly)
  rejES3.rkelly[nm+i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp4$n.rej3.rkelly, tmp5$n.rej3.rkelly, tmp6$n.rej3.rkelly)
  tmES.mix[nm+i+1,1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
  rejES1.mix[nm+i,] = cbind(tmp2$n.rej1.mix, tmp3$n.rej1.mix, tmp4$n.rej1.mix, tmp5$n.rej1.mix, tmp6$n.rej1.mix)
  rejES2.mix[nm+i,] = cbind(tmp2$n.rej2.mix, tmp3$n.rej2.mix, tmp4$n.rej2.mix, tmp5$n.rej2.mix, tmp6$n.rej2.mix)
  rejES3.mix[nm+i,] = cbind(tmp2$n.rej3.mix, tmp3$n.rej3.mix, tmp4$n.rej3.mix, tmp5$n.rej3.mix, tmp6$n.rej3.mix)
  ES.lambda.akelly[nm+i,1:length(tmp2$out.lambda)] = tmp2$out.lambda
  ES.lambda.rkelly[nm+i,1:length(tmp2$out.lambda.rkelly)] = tmp2$out.lambda.rkelly
  if(i == 3){
    tmES.con[2*nm+2,1:length(tmp4$e.out.con)] <- log(tmp4$e.out.con)
    tmES.akelly[2*nm+2,1:length(tmp4$e.out.akelly)] <- log(tmp4$e.out.akelly)
    tmES.rkelly[2*nm+2,1:length(tmp4$e.out.rkelly)] <- log(tmp4$e.out.rkelly)
    tmES.mix[2*nm+2,1:length(tmp4$e.out.mix)] <- log(tmp4$e.out.mix)
  }
}

fcES
eES.con
eES.akelly
eES.rkelly
eES.mix
rejES1.con
rejES2.con
rejES3.con
rejES1.akelly
rejES2.akelly
rejES3.akelly
rejES1.rkelly
rejES2.rkelly
rejES3.rkelly
rejES1.mix
rejES2.mix
rejES3.mix


# plots for data from different years (choose one and comment out others)

# plots from 2000

setEPS()
postscript("ES875con_empirical_all.eps")
plot(dates[1:2000],tmES.con[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:2000],tmES.con[2,1:2000],type="l",col="green")
lines(dates[1:2000],tmES.con[3,1:2000],type="l",col="blue")
lines(dates[1:2000],tmES.con[4,1:2000],type="l",col="red")
lines(dates[1:2000],tmES.con[5,1:2000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875akelly_empirical_all.eps")
plot(dates[1:2000],tmES.akelly[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:2000],tmES.akelly[2,1:2000],type="l",col="green")
lines(dates[1:2000],tmES.akelly[3,1:2000],type="l",col="blue")
lines(dates[1:2000],tmES.akelly[4,1:2000],type="l",col="red")
lines(dates[1:2000],tmES.akelly[5,1:2000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875rkelly_empirical_all.eps")
plot(dates[1:2000],tmES.rkelly[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:2000],tmES.rkelly[2,1:2000],type="l",col="green")
lines(dates[1:2000],tmES.rkelly[3,1:2000],type="l",col="blue")
lines(dates[1:2000],tmES.rkelly[4,1:2000],type="l",col="red")
lines(dates[1:2000],tmES.rkelly[5,1:2000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875mix_empirical_all.eps")
plot(dates[1:2000],tmES.mix[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:2000],tmES.mix[2,1:2000],type="l",col="green")
lines(dates[1:2000],tmES.mix[3,1:2000],type="l",col="blue")
lines(dates[1:2000],tmES.mix[4,1:2000],type="l",col="red")
lines(dates[1:2000],tmES.mix[5,1:2000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()


setEPS()
postscript("ES975con_empirical_all.eps")
plot(dates[1:4000],tmES.con[6,1:4000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:4000],tmES.con[7,1:4000],type="l",col="green")
lines(dates[1:4000],tmES.con[8,1:4000],type="l",col="blue")
lines(dates[1:4000],tmES.con[9,1:4000],type="l",col="red")
lines(dates[1:4000],tmES.con[10,1:4000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975akelly_empirical_all.eps")
plot(dates[1:4000],tmES.akelly[6,1:4000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:4000],tmES.akelly[7,1:4000],type="l",col="green")
lines(dates[1:4000],tmES.akelly[8,1:4000],type="l",col="blue")
lines(dates[1:4000],tmES.akelly[9,1:4000],type="l",col="red")
lines(dates[1:4000],tmES.akelly[10,1:4000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()


setEPS()
postscript("ES975rkelly_empirical_all.eps")
plot(dates[1:4000],tmES.rkelly[6,1:4000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:4000],tmES.rkelly[7,1:4000],type="l",col="green")
lines(dates[1:4000],tmES.rkelly[8,1:4000],type="l",col="blue")
lines(dates[1:4000],tmES.rkelly[9,1:4000],type="l",col="red")
lines(dates[1:4000],tmES.rkelly[10,1:4000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975mix_empirical_all.eps")
plot(dates[1:4000],tmES.mix[6,1:4000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:4000],tmES.mix[7,1:4000],type="l",col="green")
lines(dates[1:4000],tmES.mix[8,1:4000],type="l",col="blue")
lines(dates[1:4000],tmES.mix[9,1:4000],type="l",col="red")
lines(dates[1:4000],tmES.mix[10,1:4000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
       col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
dev.off()


# plots from 2005

# setEPS()
# postscript("ES875con_empirical2005.eps")
# plot(dates[1:3000],tmES.con[1,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.con[2,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.con[3,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.con[4,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.con[5,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875akelly_empirical2005.eps")
# plot(dates[1:3000],tmES.akelly[1,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.akelly[2,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.akelly[3,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.akelly[4,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.akelly[5,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875rkelly_empirical2005.eps")
# plot(dates[1:3000],tmES.rkelly[1,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.rkelly[2,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.rkelly[3,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.rkelly[4,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.rkelly[5,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875mix_empirical2005.eps")
# plot(dates[1:3000],tmES.mix[1,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.mix[2,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.mix[3,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.mix[4,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.mix[5,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# 
# setEPS()
# postscript("ES975con_empirical2005.eps")
# plot(dates[1:3500],tmES.con[6,1:3500],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3500],tmES.con[7,1:3500],type="l",col="green")
# lines(dates[1:3500],tmES.con[8,1:3500],type="l",col="blue")
# lines(dates[1:3500],tmES.con[9,1:3500],type="l",col="red")
# lines(dates[1:3500],tmES.con[10,1:3500],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES975akelly_empirical2005.eps")
# plot(dates[1:3500],tmES.akelly[6,1:3500],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3500],tmES.akelly[7,1:3500],type="l",col="green")
# lines(dates[1:3500],tmES.akelly[8,1:3500],type="l",col="blue")
# lines(dates[1:3500],tmES.akelly[9,1:3500],type="l",col="red")
# lines(dates[1:3500],tmES.akelly[10,1:3500],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# 
# setEPS()
# postscript("ES975rkelly_empirical2005.eps")
# plot(dates[1:3500],tmES.rkelly[6,1:3500],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3500],tmES.rkelly[7,1:3500],type="l",col="green")
# lines(dates[1:3500],tmES.rkelly[8,1:3500],type="l",col="blue")
# lines(dates[1:3500],tmES.rkelly[9,1:3500],type="l",col="red")
# lines(dates[1:3500],tmES.rkelly[10,1:3500],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES975mix_empirical2005.eps")
# plot(dates[1:3500],tmES.mix[6,1:3500],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3500],tmES.mix[7,1:3500],type="l",col="green")
# lines(dates[1:3500],tmES.mix[8,1:3500],type="l",col="blue")
# lines(dates[1:3500],tmES.mix[9,1:3500],type="l",col="red")
# lines(dates[1:3500],tmES.mix[10,1:3500],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()


# plots from 2010


# setEPS()
# postscript("ES875con_empirical2010.eps")
# plot(dates[1:2000],tmES.con[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:2000],tmES.con[2,1:2000],type="l",col="green")
# lines(dates[1:2000],tmES.con[3,1:2000],type="l",col="blue")
# lines(dates[1:2000],tmES.con[4,1:2000],type="l",col="red")
# lines(dates[1:2000],tmES.con[5,1:2000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875akelly_empirical2010.eps")
# plot(dates[1:2000],tmES.akelly[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:2000],tmES.akelly[2,1:2000],type="l",col="green")
# lines(dates[1:2000],tmES.akelly[3,1:2000],type="l",col="blue")
# lines(dates[1:2000],tmES.akelly[4,1:2000],type="l",col="red")
# lines(dates[1:2000],tmES.akelly[5,1:2000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875rkelly_empirical2010.eps")
# plot(dates[1:2000],tmES.rkelly[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:2000],tmES.rkelly[2,1:2000],type="l",col="green")
# lines(dates[1:2000],tmES.rkelly[3,1:2000],type="l",col="blue")
# lines(dates[1:2000],tmES.rkelly[4,1:2000],type="l",col="red")
# lines(dates[1:2000],tmES.rkelly[5,1:2000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875mix_empirical2010.eps")
# plot(dates[1:2000],tmES.mix[1,1:2000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:2000],tmES.mix[2,1:2000],type="l",col="green")
# lines(dates[1:2000],tmES.mix[3,1:2000],type="l",col="blue")
# lines(dates[1:2000],tmES.mix[4,1:2000],type="l",col="red")
# lines(dates[1:2000],tmES.mix[5,1:2000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# 
# setEPS()
# postscript("ES975con_empirical2010.eps")
# plot(dates[1:3000],tmES.con[6,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.con[7,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.con[8,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.con[9,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.con[10,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES975akelly_empirical2010.eps")
# plot(dates[1:3000],tmES.akelly[6,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.akelly[7,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.akelly[8,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.akelly[9,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.akelly[10,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# 
# setEPS()
# postscript("ES975rkelly_empirical2010.eps")
# plot(dates[1:3000],tmES.rkelly[6,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.rkelly[7,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.rkelly[8,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.rkelly[9,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.rkelly[10,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES975mix_empirical2010.eps")
# plot(dates[1:3000],tmES.mix[6,1:3000],type="l",xlab="dates", ylab="e-process", col="orange", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
# lines(dates[1:3000],tmES.mix[7,1:3000],type="l",col="green")
# lines(dates[1:3000],tmES.mix[8,1:3000],type="l",col="blue")
# lines(dates[1:3000],tmES.mix[9,1:3000],type="l",col="red")
# lines(dates[1:3000],tmES.mix[10,1:3000],type="l",col="black")
# abline(h = log(2), lwd = 1, lty = 2)
# abline(h = log(5), lwd = 1, lty = 2)
# abline(h = log(10), lwd = 1, lty = 2)
# legend("topleft", legend=c("normal", "t", "skewed-t", "empirical", "st +10% ES"),
#        col=c("orange", "green", "blue", "red", "black"), lwd = 1, cex=1)
# dev.off()