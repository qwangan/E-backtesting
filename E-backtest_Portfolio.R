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
library("quadprog")
library("nloptr")
library("parallel")
library("sgt")

source("Rfns.R")


# Portfolio
n_stock = 22 # number of stocks
VZ = rev(read.csv("portfolio.csv")$VZ_R)[-1]*100 # Portfolio negated percentage log returns (Jan 5, 2001 - Dec 31, 2021)
AT = rev(read.csv("portfolio.csv")$AT_R)[-1]*100
TWX = rev(read.csv("portfolio.csv")$TWX_R)[-1]*100
HD = rev(read.csv("portfolio.csv")$HD_R)[-1]*100
WMT = rev(read.csv("portfolio.csv")$WMT_R)[-1]*100
PG = rev(read.csv("portfolio.csv")$PG_R)[-1]*100
XOM = rev(read.csv("portfolio.csv")$XOM_R)[-1]*100
CVX = rev(read.csv("portfolio.csv")$CVX_R)[-1]*100
C = rev(read.csv("portfolio.csv")$C_R)[-1]*100
BAC = rev(read.csv("portfolio.csv")$BAC_R)[-1]*100
PFE = rev(read.csv("portfolio.csv")$PFE_R)[-1]*100
JNJ = rev(read.csv("portfolio.csv")$JNJ_R)[-1]*100
GE = rev(read.csv("portfolio.csv")$GE_R)[-1]*100
UPS = rev(read.csv("portfolio.csv")$UPS_R)[-1]*100
DD = rev(read.csv("portfolio.csv")$DD_R)[-1]*100
DOW = rev(read.csv("portfolio.csv")$DOW_R)[-1]*100
WY = rev(read.csv("portfolio.csv")$WY_R)[-1]*100
SPG = rev(read.csv("portfolio.csv")$SPG_R)[-1]*100
MSFT = rev(read.csv("portfolio.csv")$MSFT_R)[-1]*100
IBM = rev(read.csv("portfolio.csv")$IBM_R)[-1]*100
EXC = rev(read.csv("portfolio.csv")$EXC_R)[-1]*100
SO = rev(read.csv("portfolio.csv")$SO_R)[-1]*100
portfolio = cbind(VZ, AT, TWX, HD, WMT, PG, XOM, CVX, C, BAC, PFE, JNJ, GE, UPS, DD, DOW, WY, SPG, MSFT, IBM, EXC, SO)
x = rev(read.csv("portfolio.csv")$NegLogReturn)[-1]*100



# Parameters

N=length(x) # sample size
w=500 # moving window size (same as in simulations)
(n=N-w) # out-of-sample size
y = matrix(nrow = n, ncol = 3) # returns observed by the regulator


avec=c(.95, .99) # vector of alpha levels for VaR
nvec=c(.875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR

# =======================================================
# The bank uses mean-variance model to optimize portfolio
# maximize E[wX] - ra / 2 * Var(wX)
# =======================================================

ra = 1 # level of risk aversion

# =======================================================
# Normal innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="norm")
qmodel = qdist("norm",p=VaR.levels,mu=0,sigma=1) # VaR values
esmodel = (ddist("norm", qdist("norm",p=nvec,mu=0,sigma=1), mu=0, sigma=1)) / (1 - nvec)

VaR_norm <- matrix(nrow=n,ncol=length(VaR.levels))

ES_norm <- matrix(nrow=n,ncol=length(nvec))

weight_norm = matrix(0, nrow = n, ncol = n_stock) # record weights


# ----------------------------------------------------

for(i in 1:n)
{
  # Fit each stock to AR(1)-GARCH(1,1)
  dummy_portfolio <- 1 - is.na(portfolio[i,])
  n_stock_tmp = sum(dummy_portfolio) # num of valid stocks
  portfolio_tmp = portfolio[,dummy_portfolio]

  mut  <- vector(mode="numeric", length=n_stock_tmp)
  sigt  <- vector(mode="numeric", length=n_stock_tmp)
  res = matrix(nrow = w, ncol = n_stock_tmp) # matrix of the residuals

  for(j in 1:n_stock_tmp){
    fit = ugarchfit(spec, portfolio_tmp[i:(i+w-1), j], solver="hybrid")
    res[,j] = residuals(fit)
    frcst = ugarchforecast(fit,n.ahead=1)
    mut[j] = fitted(frcst)
    sigt[j] = sigma(frcst)
  }
  Corr = cor(res)

  # Solve optimal weights with constrained optimization
  eval_f = function(w){
    - w %*% mut - ra/2 * (w*sigt) %*% Corr %*% (w*sigt)
  }
  eval_g_eq = function(w){
    sum(w)-1
  }
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  weight = nloptr(x0=rep(1/n_stock_tmp,n_stock_tmp), eval_f=eval_f, lb=rep(0,n_stock_tmp), ub=rep(1,n_stock_tmp),
                  eval_g_eq=eval_g_eq, opts=opts)$solution
  weight_norm[i,dummy_portfolio] = weight

  # Calculate VaR and ES for the portfolio
  VaR_norm[i,]=weight %*% mut + sqrt((weight*sigt) %*% Corr %*% (weight*sigt)) * qmodel
  ES_norm[i,]=weight %*% mut + sqrt((weight*sigt) %*% Corr %*% (weight*sigt)) * esmodel

  # Report return
  y[i,1] = weight %*% portfolio_tmp[(i+w),]

  if(i %% 100 == 0) print(i)
}

out=list(VaR_norm=VaR_norm, ES_norm=ES_norm)
save(out, file="VaRES_norm.RDATA")
out=y[,1]
save(out, file="observed_norm.RDATA")
save(weight_norm, file="weight_norm.RDATA")

# =======================================================
# Student t innovations
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="std")

VaR_t <- matrix(nrow=n,ncol=length(VaR.levels))

ES_t <- matrix(nrow=n,ncol=length(nvec))

weight_t = matrix(0, nrow = n, ncol = n_stock) # record weights

# ----------------------------------------------------
for(i in 1:n)
{
  # Fit each stock to AR(1)-GARCH(1,1)
  dummy_portfolio <- 1 - is.na(portfolio[i,])
  n_stock_tmp = sum(dummy_portfolio) # num of valid stocks
  portfolio_tmp = portfolio[,dummy_portfolio]

  mut  <- vector(mode="numeric", length=n_stock_tmp)
  sigt  <- vector(mode="numeric", length=n_stock_tmp)
  res = matrix(nrow = w, ncol = n_stock_tmp) # matrix of the residuals

  for(j in 1:n_stock_tmp){
    fit = ugarchfit(spec, portfolio_tmp[i:(i+w-1), j], solver="hybrid")
    res[,j] = residuals(fit)
    frcst = ugarchforecast(fit,n.ahead=1)
    mut[j] = fitted(frcst)
    sigt[j] = sigma(frcst)
  }
  Corr = cor(res)

  # Solve optimal weights with constrained optimization
  eval_f = function(w){
    - w %*% mut - ra/2 * (w*sigt) %*% Corr %*% (w*sigt)
  }
  eval_g_eq = function(w){
    sum(w)-1
  }
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  weight = nloptr(x0=rep(1/n_stock_tmp,n_stock_tmp), eval_f=eval_f, lb=rep(0,n_stock_tmp), ub=rep(1,n_stock_tmp),
                  eval_g_eq=eval_g_eq, opts=opts)$solution
  weight_t[i,dummy_portfolio] = weight

  # Calculate VaR and ES for the portfolio
  NLR_portfolio = portfolio_tmp %*% weight

  fit = ugarchfit(spec, NLR_portfolio[i:(i+w-1)], solver="hybrid")
  nu=coef(fit)["shape"]

  qmodel = qdist("std",p=VaR.levels,mu=0,sigma=1,shape=nu)

  # expected shortfall for the assumed model/distribution
  esmodel = (dt(qt(nvec, df=nu), df=nu)) / (1-nvec) * (nu + (qt(nvec, df=nu))^2) / (nu-1)
  esmodel = sqrt((nu-2)/nu)*esmodel


  VaR_t[i,]=weight %*% mut + sqrt((weight*sigt) %*% Corr %*% (weight*sigt)) * qmodel
  ES_t[i,]=weight %*% mut + sqrt((weight*sigt) %*% Corr %*% (weight*sigt)) * esmodel

  # Report return
  y[i,2] = NLR_portfolio[(i+w)]

  if(i %% 100 == 0) print(i)
}

out=list(VaR_t=VaR_t, ES_t=ES_t)
save(out, file="VaRES_t.RDATA")
out=y[,2]
save(out, file="observed_t.RDATA")
save(weight_t, file="weight_t.RDATA")


# =======================================================
# skewed Student t innovations
# The version of Fernandez & Steel (1998)
# =======================================================

spec = ugarchspec(mean.model=list(armaOrder=c(1,0),include.mean=T), distribution.model="sstd")

VaR_st <- matrix(nrow=n,ncol=length(VaR.levels))

ES_st <- matrix(nrow=n,ncol=length(nvec))

weight_st = matrix(0, nrow = n, ncol = n_stock) # record weights

for(i in 1:n)

{
  # Fit each stock to AR(1)-GARCH(1,1)
  dummy_portfolio <- 1 - is.na(portfolio[i,])
  n_stock_tmp = sum(dummy_portfolio) # num of valid stocks
  portfolio_tmp = portfolio[,dummy_portfolio]

  mut  <- vector(mode="numeric", length=n_stock_tmp)
  sigt  <- vector(mode="numeric", length=n_stock_tmp)
  res = matrix(nrow = w, ncol = n_stock_tmp) # matrix of the residuals

  for(j in 1:n_stock_tmp){
    fit = ugarchfit(spec, portfolio_tmp[i:(i+w-1), j], solver="hybrid")
    res[,j] = residuals(fit)
    frcst = ugarchforecast(fit,n.ahead=1)
    mut[j] = fitted(frcst)
    sigt[j] = sigma(frcst)
  }
  Corr = cor(res)

  # Solve optimal weights with constrained optimization
  eval_f = function(w){
    - w %*% mut - ra/2 * (w*sigt) %*% Corr %*% (w*sigt)
  }
  eval_g_eq = function(w){
    sum(w)-1
  }
  local_opts <- list( "algorithm" = "NLOPT_LD_MMA", "xtol_rel" = 1.0e-15 )
  opts <- list( "algorithm"= "NLOPT_GN_ISRES",
                "xtol_rel"= 1.0e-15,
                "maxeval"= 160000,
                "local_opts" = local_opts,
                "print_level" = 0 )
  weight = nloptr(x0=rep(1/n_stock_tmp,n_stock_tmp), eval_f=eval_f, lb=rep(0,n_stock_tmp), ub=rep(1,n_stock_tmp),
                  eval_g_eq=eval_g_eq, opts=opts)$solution
  weight_st[i,dummy_portfolio] = weight

  # Calculate VaR and ES for the portfolio
  NLR_portfolio = portfolio_tmp %*% weight

  fit = ugarchfit(spec, NLR_portfolio[i:(i+w-1)], solver="hybrid")

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


  VaR_st[i,]=weight %*% mut + sqrt((weight*sigt) %*% Corr %*% (weight*sigt)) * qmodel
  ES_st[i,]=weight %*% mut + sqrt((weight*sigt) %*% Corr %*% (weight*sigt)) * esmodel


  # Report return
  y[i,3] = NLR_portfolio[(i+w)]

  if(i %% 100 == 0) print(i)
}

out=list(VaR_st=VaR_st, ES_st=ES_st)
save(out, file="VaRES_st.RDATA")
out=y[,3]
save(out, file="observed_st.RDATA")
save(weight_st, file="weight_st.RDATA")



out=list(VaR_norm=VaR_norm, ES_norm=ES_norm, VaR_t=VaR_t, ES_t=ES_t, VaR_st=VaR_st, ES_st=ES_st)
save(out, file="VaRES_portfolio.RDATA")
save(y, file="observed_return.RDATA")



# Portfolio
dates = as.Date(rev(read.csv("portfolio.csv")$Date))[-1] # dates for returns

# Length of data
m = N-w+2 # Jan 3, 2005 - Dec 31, 2021
n=m-2

load("VaRES_portfolio.RDATA")
load("observed_return.RDATA")
y=tail(y,n)
dates=tail(dates,n)
dates = as.Date(dates[-c(1:(w))])


# --------------------------------------------------------
# Forecasts

#=====================================================
# (VaR_nu, ES_nu)
#=====================================================

VaRout1b = rbind(out$VaR_norm[,3], out$VaR_t[,3], out$VaR_st[,3])
ESout1 = rbind(out$ES_norm[,1], out$ES_t[,1], out$ES_st[,1])

VaRout2b = rbind(out$VaR_norm[,4], out$VaR_t[,4], out$VaR_st[,4])
ESout2 = rbind(out$ES_norm[,2], out$ES_t[,2], out$ES_st[,2])

VaRout1b = VaRout1b[1:3,(N-w-n+1):(N-w)]; ESout1 = ESout1[1:3,(N-w-n+1):(N-w)]
VaRout2b = VaRout2b[1:3,(N-w-n+1):(N-w)]; ESout2 = ESout2[1:3,(N-w-n+1):(N-w)]



# Portfolio plots

setEPS()
postscript("Portfolio.eps")
plot(dates,y[-c(1:(w)),1],type="l",xlab="dates", ylab="negated percentage log returns", ylim=c(-8,8), col="red", lty=1, lwd=2)
lines(dates,y[-c(1:(w)),2],type="l",col="green",lwd=1.5)
lines(dates,y[-c(1:(w)),3],type="l",col="blue", lty = 3,lwd=0.8)
legend("topleft", inset=c(0,0), legend=c("normal", "t", "skewed-t"),
       col=c("red","green","blue"), lwd = c(2,1.5,0.8), lty = c(1,1,3), cex=1)
dev.off()


setEPS()
postscript("Portfolio_ES975.eps")
plot(dates,ESout2[1,-c(1:(w))],type="l",xlab="dates", ylab="ES forecast",col="red", ylim=c(1.2,7), lty=1)
lines(dates,ESout2[2,-c(1:(w))],type="l",col="green")
lines(dates,ESout2[3,-c(1:(w))],type="l",col="blue", lty = 3)
legend("topleft", inset=c(0,0), legend=c("normal", "t", "skewed-t"),
       col=c("red","green","blue"), lwd = c(1,1,1), lty = c(1,1,3), cex=1)
dev.off()



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
ES.lambda.akelly = matrix(nrow=2*nm, ncol=length(dates)) # lambda for GREE
ES.lambda.rkelly = matrix(nrow=2*nm, ncol=length(dates)) # lambda for GREL
fcES = matrix(nrow=2*nm, ncol=2) # VaR-ES forecast
eES.con = matrix(nrow=2*nm, ncol=5) 
eES.akelly = matrix(nrow=2*nm, ncol=5)
eES.rkelly = matrix(nrow=2*nm, ncol=5)
eES.mix = matrix(nrow=2*nm, ncol=5) # final e-values for ES
tmES.con = matrix(nrow=2*(nm+1), ncol=length(dates))
tmES.akelly = matrix(nrow=2*(nm+1), ncol=length(dates))
tmES.rkelly = matrix(nrow=2*(nm+1), ncol=length(dates))
tmES.mix = matrix(nrow=2*(nm+1), ncol=length(dates)) # sequential test martingale for ES
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
  tmp2=evalue_ES_em(x=y[,i], r=ESout1[i,], z=VaRout1b[i,],lev=nvec[1], w=500) # sequential test margtingale
  tmp3=evalue_ES_em(x=y[,i], r=ESout1[i,] * (1-err), z=VaRout1b[i,], lev=nvec[1], w=500) # sequential test margtingale for under-reporting ES
  tmp4=evalue_ES_em(x=y[,i], r=ESout1[i,] * (1+err), z=VaRout1b[i,],lev=nvec[1], w=500) # sequential test margtingale for over-reporting ES
  tmp5=evalue_ES_em(x=y[,i], r=ESout1[i,] * (1-err), z=VaRout1b[i,] * (1-err), lev=nvec[1], w=500) # sequential test margtingale for under-reporting both VaR and ES
  tmp6=evalue_ES_em(x=y[,i], r=ESout1[i,] * (1+err), z=VaRout1b[i,] * (1+err),lev=nvec[1], w=500) # sequential test margtingale for over-reporting both VaR and ES
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
  
  
  tmp2=evalue_ES_em(x=y[,i],  r=ESout2[i,], z=VaRout2b[i,],lev=nvec[2], w=500)
  tmp3=evalue_ES_em(x=y[,i], r=ESout2[i,] * (1-err), z=VaRout2b[i,],lev=nvec[2], w=500)
  tmp4=evalue_ES_em(x=y[,i], r=ESout2[i,] * (1+err), z=VaRout2b[i,],lev=nvec[2], w=500)
  tmp5=evalue_ES_em(x=y[,i], r=ESout2[i,] * (1-err), z=VaRout2b[i,] * (1-err),lev=nvec[2], w=500)
  tmp6=evalue_ES_em(x=y[,i], r=ESout2[i,] * (1+err), z=VaRout2b[i,] * (1+err),lev=nvec[2], w=500)
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

out=list(fcES=fcES, eES.con=eES.con, eES.akelly=eES.akelly, eES.rkelly=eES.rkelly, eES.mix=eES.mix,
         rejES1.con=rejES1.con, rejES2.con=rejES2.con, rejES3.con=rejES3.con,
         rejES1.akelly=rejES1.akelly, rejES2.akelly=rejES2.akelly, rejES3.akelly=rejES3.akelly,
         rejES1.rkelly=rejES1.rkelly, rejES2.rkelly=rejES2.rkelly, rejES3.rkelly=rejES3.rkelly,
         rejES1.mix=rejES1.mix, rejES2.mix=rejES2.mix, rejES3.mix=rejES3.mix,
         tmES.con=tmES.con, tmES.akelly=tmES.akelly, tmES.rkelly=tmES.rkelly, tmES.mix=tmES.mix)
save(out, file="Result_portfolio.RDATA")


# Portfolio plots
load("Result_portfolio.RDATA")

setEPS()
postscript("ES875con_portfolio2005.eps")
plot(dates[1:3000],out$tmES.con[1,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.con[2,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.con[3,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.con[4,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875akelly_portfolio2005.eps")
plot(dates[1:3000],out$tmES.akelly[1,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.akelly[2,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.akelly[3,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.akelly[4,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875rkelly_portfolio2005.eps")
plot(dates[1:3000],out$tmES.rkelly[1,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.rkelly[2,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.rkelly[3,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.rkelly[4,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875mix_portfolio2005.eps")
plot(dates[1:3000],out$tmES.mix[1,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.mix[2,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.mix[3,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.mix[4,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()


setEPS()
postscript("ES975con_portfolio2005.eps")
plot(dates[1:3000],out$tmES.con[5,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.con[6,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.con[7,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.con[8,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975akelly_portfolio2005.eps")
plot(dates[1:3000],out$tmES.akelly[5,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.akelly[6,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.akelly[7,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.akelly[8,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()


setEPS()
postscript("ES975rkelly_portfolio2005.eps")
plot(dates[1:3000],out$tmES.rkelly[5,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.rkelly[6,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.rkelly[7,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.rkelly[8,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()


setEPS()
postscript("ES975mix_portfolio2005.eps")
plot(dates[1:3000],out$tmES.mix[5,1:3000],type="l",xlab="dates", ylab="e-process", col="red", lty=1, cex.lab = 1.4, cex.axis = 1.4, pch = 16, ylim=e.lim)
lines(dates[1:3000],out$tmES.mix[6,1:3000],type="l",col="green")
lines(dates[1:3000],out$tmES.mix[7,1:3000],type="l",col="blue")
lines(dates[1:3000],out$tmES.mix[8,1:3000],type="l",col="black")
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("normal", "t", "skewed-t", "st +10% ES"),
       col=c("red", "green", "blue", "black"), lwd = 1, cex=1)
dev.off()
