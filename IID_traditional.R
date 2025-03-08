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
# library("foreach")
# library("doParallel")


source("Rfns.R")


# Set number of cores
NumCores <- detectCores()

# Set number of runs
NumRuns <- 1000


# Set parameters
avec=c(.95, .99) # vector of alpha levels for VaR
nvec=c(.875, .975) # vector of nu levels for (VaR_nu, ES_nu)
VaR.levels <- c(avec,nvec) # VaR levels for VaR on its own (1:3) and in pair with ES (4:6)
inu = (length(avec)+1):(length(avec)+length(nvec)) #index set for nu levels of VaR
n = 1000 # sample size


# Function

simu <- function(iter) {
  print(iter)
  
  # Simulate normal distribution
  set.seed(iter)
  simnorm = rnorm(n)
  simnorm.forecast = rnorm(500) # data in estimation window
  
  y = simnorm
  
  # calculate VaR and ES
  normVaR = qnorm(VaR.levels, 0, 1)
  normES = (dnorm(qnorm(nvec, 0, 1), 0, 1)) / (1 - nvec)
  VaRout = rep(normVaR[1], length(y))
  VaRout = rbind(VaRout, rep(normVaR[2], length(y)))
  VaRoutb = rep(normVaR[3], length(y))
  VaRoutb = rbind(VaRoutb, rep(normVaR[4], length(y)))
  ESout = rep(normES[1], length(y))
  ESout = rbind(ESout, rep(normES[2], length(y)))
  
  
  # testing procedure
  
  err <- .1 # error for under- and over-reporting
  w <- 0 # time window
  
  # VaR
  
  pVaR = matrix(nrow=2, ncol=2) # final p-values for VaR
  rejVaR1 = matrix(nrow=2, ncol=2) # numbers of days where we reject for VaR (threshold 0.5)
  rejVaR2 = matrix(nrow=2, ncol=2) # numbers of days where we reject for VaR (threshold 0.2)
  rejVaR3 = matrix(nrow=2, ncol=2) # numbers of days where we reject for VaR (threshold 0.1)
  
  for(i in 1:2)
  {
    tmp2=cct.1s.VaR(x=y, r=VaRout[i,],lev=avec[i]) # sequential p-values
    tmp3=cct.1s.VaR(x=y, r=VaRout[i,] * (1-err),lev=avec[i]) # sequential p-values for under-reporting
    pVaR[i,] = c(tmp2$out.pv.cond[length(tmp2$out.pv.cond)], tmp3$out.pv.cond[length(tmp3$out.pv.cond)])
    rejVaR1[i,] = cbind(tmp2$n.rej1, tmp3$n.rej1)
    rejVaR2[i,] = cbind(tmp2$n.rej2, tmp3$n.rej2)
    rejVaR3[i,] = cbind(tmp2$n.rej3, tmp3$n.rej3)
  }
  
  rejnumVaR1 = (rejVaR1 < length(y))
  rejnumVaR2 = (rejVaR2 < length(y))
  rejnumVaR3 = (rejVaR3 < length(y))
  
  # ES
  
  pES = matrix(nrow=2, ncol=3) # final p-values for ES
  rejES1 = matrix(nrow=2, ncol=3) # numbers of days where we reject for ES (threshold 0.5)
  rejES2 = matrix(nrow=2, ncol=3) # numbers of days where we reject for ES (threshold 0.2)
  rejES3 = matrix(nrow=2, ncol=3) # numbers of days where we reject for ES (threshold 0.1)
  
  # estimate sigma
  y.extend = c(y, simnorm.forecast)
  sigt = c()
  for(i in 1:length(y)){
    sigt= cbind(sigt, sqrt(var(y.extend[i:(i+499)])))
  }
  
  
  for(i in 1:2)
  {
    tmp2=cct.onesided2(x=y, r=cbind(VaRoutb[i,], ESout[i,]),lev=nvec[i],sigt=sigt)  # sequential p-values
    tmp3=cct.onesided2(x=y, r=cbind(VaRoutb[i,], ESout[i,]*(1-err)),lev=nvec[i],sigt=sigt) # sequential p-values for under-reporting ES
    tmp5=cct.onesided2(x=y, r=cbind(VaRoutb[i,], ESout[i,])*(1-err),lev=nvec[i],sigt=sigt) # sequential p-values for under-reporting both VaR and ES
    pES[i,] = c(tmp2$out.pv.cond[length(tmp2$out.pv.cond)], tmp3$out.pv.cond[length(tmp3$out.pv.cond)], tmp5$out.pv.cond[length(tmp5$out.pv.cond)])
    rejES1[i,] = cbind(tmp2$n.rej1, tmp3$n.rej1, tmp5$n.rej1)
    rejES2[i,] = cbind(tmp2$n.rej2, tmp3$n.rej2, tmp5$n.rej2)
    rejES3[i,] = cbind(tmp2$n.rej3, tmp3$n.rej3, tmp5$n.rej3)
  }
  
  rejnumES1 = (rejES1 < length(y))
  rejnumES2 = (rejES2 < length(y))
  rejnumES3 = (rejES3 < length(y))
  
  return(c(pVaR, rejVaR1, rejVaR2, rejVaR3, rejnumVaR1, rejnumVaR2, rejnumVaR3,
           pES, rejES1, rejES2, rejES3, rejnumES1, rejnumES2, rejnumES3))
}


# Multiple-core computations
out <- do.call(rbind, mclapply(1:NumRuns, simu, mc.cores = NumCores))
save(out, file="Result_traditional.RDATA")


# Output values
avg.total = c()
for(i in 1:70){
  avg = mean(out[,i][out[,i] != n])
  avg.total = cbind(avg.total, avg)
}
avg.pVaR <- matrix(avg.total[1:4],nrow=2,ncol=2)
avg.rejVaR1 <- matrix(avg.total[5:8],nrow=2,ncol=2)
avg.rejVaR2 <- matrix(avg.total[9:12],nrow=2,ncol=2)
avg.rejVaR3 <- matrix(avg.total[13:16],nrow=2,ncol=2)
avg.rejnumVaR1 <- matrix(avg.total[17:20],nrow=2,ncol=2)
avg.rejnumVaR2 <- matrix(avg.total[21:24],nrow=2,ncol=2)
avg.rejnumVaR3 <- matrix(avg.total[25:28],nrow=2,ncol=2)
avg.pES <- matrix(avg.total[29:34],nrow=2,ncol=3)
avg.rejES1 <- matrix(avg.total[35:40],nrow=2,ncol=3)
avg.rejES2 <- matrix(avg.total[41:46],nrow=2,ncol=3)
avg.rejES3 <- matrix(avg.total[47:52],nrow=2,ncol=3)
avg.rejnumES1 <- matrix(avg.total[53:58],nrow=2,ncol=3)
avg.rejnumES2 <- matrix(avg.total[59:64],nrow=2,ncol=3)
avg.rejnumES3 <- matrix(avg.total[65:70],nrow=2,ncol=3)


out.final = list(avg.pVaR = avg.pVaR, avg.rejVaR1 = avg.rejVaR1, avg.rejVaR2 = avg.rejVaR2, avg.rejVaR3 = avg.rejVaR3,
                 avg.rejnumVaR1= avg.rejnumVaR1, avg.rejnumVaR2= avg.rejnumVaR2, avg.rejnumVaR3= avg.rejnumVaR3,
                 avg.pES = avg.pES, avg.rejES1 = avg.rejES1, avg.rejES2 = avg.rejES2, avg.rejES3 = avg.rejES3,
                 avg.rejnumES1 = avg.rejnumES1, avg.rejnumES2 = avg.rejnumES2, avg.rejnumES3 = avg.rejnumES3)
save(out.final, file="Result_traditional_final.RDATA")


avg.pVaR
avg.rejVaR1
avg.rejVaR2
avg.rejVaR3
avg.rejnumVaR1
avg.rejnumVaR2
avg.rejnumVaR3
avg.pES
avg.rejES1
avg.rejES2
avg.rejES3
avg.rejnumES1
avg.rejnumES2
avg.rejnumES3