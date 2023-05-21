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
n = 250 # sample size
e.lim <- c(-1,5) # bounds for e-values


# Function

simu <- function(iter) {
  print(iter)
  
  # Simulate normal distribution
  set.seed(iter)
  simnorm = rnorm(n)
  # simnorm.forecast = rnorm(500) # data in estimation window
  
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
  
  
  # e-testing
  
  err <- .1 # error for under- and over-reporting
  w <- 0 # time window
  
  ### summary for VaR
  # "con" represents e-backtesting method with constant lambda,
  # "rkelly" represents GREE/GREL method,
  # "rkellyQ" represents GRO method
  
  eVaR.con = matrix(nrow=2, ncol=2)
  eVaR.rkelly = matrix(nrow=2, ncol=2)
  eVaR.rkellyQ = matrix(nrow=2, ncol=2) # final e-values for VaR
  rejVaR1.con = matrix(nrow=2, ncol=2)
  rejVaR1.rkelly = matrix(nrow=2, ncol=2)
  rejVaR1.rkellyQ = matrix(nrow=2, ncol=2) # numbers of days where we reject for VaR (threshold 2)
  rejVaR2.con = matrix(nrow=2, ncol=2)
  rejVaR2.rkelly = matrix(nrow=2, ncol=2)
  rejVaR2.rkellyQ = matrix(nrow=2, ncol=2) # numbers of days where we reject for VaR (threshold 5)
  rejVaR3.con = matrix(nrow=2, ncol=2)
  rejVaR3.rkelly = matrix(nrow=2, ncol=2)
  rejVaR3.rkellyQ = matrix(nrow=2, ncol=2) # numbers of days where we reject for VaR (threshold 10)
  tmVaR.con = matrix(nrow=2, ncol=(length(y)-w))
  tmVaR.rkelly = matrix(nrow=2, ncol=(length(y)-w))
  tmVaR.rkellyQ = matrix(nrow=2, ncol=(length(y)-w)) # sequential test martingale for VaR
  VaR.lambda.rkelly = matrix(nrow=2, ncol=(length(y)-w)) # lambda for rkelly (unknown distribution)
  VaR.lambda.rkellyQ = matrix(nrow=2, ncol=(length(y)-w)) # lambda for rkelly (known distribution)
  
  for(i in 1:2)
  {
    tmp2=evalue.VaR(x=y, r=VaRout[i,],lev=avec[i]) # sequential test margtingale
    tmp3=evalue.VaR(x=y, r=VaRout[i,] * (1-err),lev=avec[i]) # sequential test margtingale for under-reporting
    tmp4=evalue.VaRQ(x=y, r=VaRout[i,],lev=avec[i]) # sequential test margtingale (known Q)
    tmp5=evalue.VaRQ(x=y, r=VaRout[i,] * (1-err),lev=avec[i]) # sequential test margtingale for under-reporting (known Q)
    eVaR.con[i,] = log(c(tmp2$e.out.con[length(tmp2$e.out.con)], tmp3$e.out.con[length(tmp3$e.out.con)]))
    eVaR.rkelly[i,] = log(c(tmp2$e.out.rkelly[length(tmp2$e.out.rkelly)], tmp3$e.out.rkelly[length(tmp3$e.out.rkelly)]))
    eVaR.rkellyQ[i,] = log(c(tmp4$e.out.rkelly[length(tmp4$e.out.rkelly)], tmp5$e.out.rkelly[length(tmp5$e.out.rkelly)]))
    rejVaR1.con[i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con)
    rejVaR1.rkelly[i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly)
    rejVaR1.rkellyQ[i,] = cbind(tmp4$n.rej1.rkelly, tmp5$n.rej1.rkelly)
    rejVaR2.con[i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con)
    rejVaR2.rkelly[i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly)
    rejVaR2.rkellyQ[i,] = cbind(tmp4$n.rej2.rkelly, tmp5$n.rej2.rkelly)
    rejVaR3.con[i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con)
    rejVaR3.rkelly[i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly)
    rejVaR3.rkellyQ[i,] = cbind(tmp4$n.rej3.rkelly, tmp5$n.rej3.rkelly)
    tmVaR.con[i,1:length(tmp3$e.out.con)] <- log(tmp3$e.out.con)
    tmVaR.rkelly[i,1:length(tmp3$e.out.rkelly)] <- log(tmp3$e.out.rkelly)
    tmVaR.rkellyQ[i,1:length(tmp5$e.out.rkelly)] <- log(tmp5$e.out.rkelly)
    VaR.lambda.rkelly[i,1:length(tmp3$e.out.rkelly)] <- tmp3$out.lambda.rkelly
    VaR.lambda.rkellyQ[i,1:length(tmp5$e.out.rkelly)] <- tmp5$out.lambda.rkelly
  }
  
  rejnumVaR1.con = (rejVaR1.con < length(y)) # whether we reject or not
  rejnumVaR1.rkelly = (rejVaR1.rkelly < length(y))
  rejnumVaR1.rkellyQ = (rejVaR1.rkellyQ < length(y))
  rejnumVaR2.con = (rejVaR2.con < length(y)) # whether we reject or not
  rejnumVaR2.rkelly = (rejVaR2.rkelly < length(y))
  rejnumVaR2.rkellyQ = (rejVaR2.rkellyQ < length(y))
  rejnumVaR3.con = (rejVaR3.con < length(y)) # whether we reject or not
  rejnumVaR3.rkelly = (rejVaR3.rkelly < length(y))
  rejnumVaR3.rkellyQ = (rejVaR3.rkellyQ < length(y))
  
  
  
  # ===================================================
  # Summary for ES
  eES.con = matrix(nrow=2, ncol=3)
  eES.rkelly = matrix(nrow=2, ncol=3)
  eES.rkellyQ = matrix(nrow=2, ncol=3) # final e-values for ES
  rejES1.con = matrix(nrow=2, ncol=3)
  rejES1.rkelly = matrix(nrow=2, ncol=3)
  rejES1.rkellyQ = matrix(nrow=2, ncol=3) # numbers of days where we reject for ES (threshold 2)
  rejES2.con = matrix(nrow=2, ncol=3)
  rejES2.rkelly = matrix(nrow=2, ncol=3)
  rejES2.rkellyQ = matrix(nrow=2, ncol=3) # numbers of days where we reject for ES (threshold 5)
  rejES3.con = matrix(nrow=2, ncol=3)
  rejES3.rkelly = matrix(nrow=2, ncol=3)
  rejES3.rkellyQ = matrix(nrow=2, ncol=3) # numbers of days where we reject for ES (threshold 10)
  tmES.con = matrix(nrow=2, ncol=(length(y)-w))
  tmES.rkelly = matrix(nrow=2, ncol=(length(y)-w))
  tmES.rkellyQ = matrix(nrow=2, ncol=(length(y)-w)) # sequential test martingale for ES
  ES.lambda.rkelly = matrix(nrow=2, ncol=(length(y)-w)) # lambda for rkelly (unknown distribution)
  ES.lambda.rkellyQ = matrix(nrow=2, ncol=(length(y)-w)) # lambda for rkelly (known distribution)
  
  for(i in 1:2)
  {
    tmp2=evalue.ES(x=y, r=ESout[i,], z=VaRoutb[i,],lev=nvec[i]) # sequential test margtingale
    tmp3=evalue.ES(x=y, r=ESout[i,] * (1-err), z=VaRoutb[i,], lev=nvec[i]) # sequential test margtingale for under-reporting ES
    tmp5=evalue.ES(x=y, r=ESout[i,] * (1-err), z=VaRoutb[i,] * (1-err), lev=nvec[i]) # sequential test margtingale for under-reporting both VaR and ES
    tmp4=evalue.ESQ(x=y, r=ESout[i,], z=VaRoutb[i,],lev=nvec[i]) # sequential test margtingale (known Q)
    tmp6=evalue.ESQ(x=y, r=ESout[i,] * (1-err), z=VaRoutb[i,], lev=nvec[i]) # sequential test margtingale for under-reporting ES (known Q)
    tmp7=evalue.ESQ(x=y, r=ESout[i,] * (1-err), z=VaRoutb[i,] * (1-err), lev=nvec[i]) # sequential test margtingale for under-reporting both VaR and ES (known Q)
    eES.con[i,] = log(c(tmp2$e.out.con[length(tmp2$e.out.con)], tmp3$e.out.con[length(tmp3$e.out.con)], tmp5$e.out.con[length(tmp5$e.out.con)]))
    eES.rkelly[i,] = log(c(tmp2$e.out.rkelly[length(tmp2$e.out.rkelly)], tmp3$e.out.rkelly[length(tmp3$e.out.rkelly)],
                           tmp5$e.out.rkelly[length(tmp5$e.out.rkelly)]))
    eES.rkellyQ[i,] = log(c(tmp4$e.out.rkelly[length(tmp4$e.out.rkelly)], tmp6$e.out.rkelly[length(tmp6$e.out.rkelly)],
                            tmp7$e.out.rkelly[length(tmp7$e.out.rkelly)]))
    rejES1.con[i,] = cbind(tmp2$n.rej1.con, tmp3$n.rej1.con, tmp5$n.rej1.con)
    rejES1.rkelly[i,] = cbind(tmp2$n.rej1.rkelly, tmp3$n.rej1.rkelly, tmp5$n.rej1.rkelly)
    rejES1.rkellyQ[i,] = cbind(tmp4$n.rej1.rkelly, tmp6$n.rej1.rkelly, tmp7$n.rej1.rkelly)
    rejES2.con[i,] = cbind(tmp2$n.rej2.con, tmp3$n.rej2.con, tmp5$n.rej2.con)
    rejES2.rkelly[i,] = cbind(tmp2$n.rej2.rkelly, tmp3$n.rej2.rkelly, tmp5$n.rej2.rkelly)
    rejES2.rkellyQ[i,] = cbind(tmp4$n.rej2.rkelly, tmp6$n.rej2.rkelly, tmp7$n.rej2.rkelly)
    rejES3.con[i,] = cbind(tmp2$n.rej3.con, tmp3$n.rej3.con, tmp5$n.rej3.con)
    rejES3.rkelly[i,] = cbind(tmp2$n.rej3.rkelly, tmp3$n.rej3.rkelly, tmp5$n.rej3.rkelly)
    rejES3.rkellyQ[i,] = cbind(tmp4$n.rej3.rkelly, tmp6$n.rej3.rkelly, tmp7$n.rej3.rkelly)
    tmES.con[i,1:length(tmp3$e.out.con)] <- log(tmp3$e.out.con)
    tmES.rkelly[i,1:length(tmp3$e.out.rkelly)] <- log(tmp3$e.out.rkelly)
    tmES.rkellyQ[i,1:length(tmp6$e.out.rkelly)] <- log(tmp6$e.out.rkelly)
    ES.lambda.rkelly[i,1:length(tmp3$e.out.rkelly)] <- tmp3$out.lambda.rkelly
    ES.lambda.rkellyQ[i,1:length(tmp6$e.out.rkelly)] <- tmp6$out.lambda.rkelly
  }
  
  rejnumES1.con = (rejES1.con < length(y)) # whether we reject or not
  rejnumES1.rkelly = (rejES1.rkelly < length(y))
  rejnumES1.rkellyQ = (rejES1.rkellyQ < length(y))
  rejnumES2.con = (rejES2.con < length(y)) # whether we reject or not
  rejnumES2.rkelly = (rejES2.rkelly < length(y))
  rejnumES2.rkellyQ = (rejES2.rkellyQ < length(y))
  rejnumES3.con = (rejES3.con < length(y)) # whether we reject or not
  rejnumES3.rkelly = (rejES3.rkelly < length(y))
  rejnumES3.rkellyQ = (rejES3.rkellyQ < length(y))
  
  
  return(c(eVaR.con, eVaR.rkelly, eVaR.rkellyQ, rejVaR1.con, rejVaR1.rkelly, rejVaR1.rkellyQ,
           rejVaR2.con, rejVaR2.rkelly, rejVaR2.rkellyQ, rejVaR3.con, rejVaR3.rkelly, rejVaR3.rkellyQ,
           rejnumVaR1.con, rejnumVaR1.rkelly, rejnumVaR1.rkellyQ, rejnumVaR2.con, rejnumVaR2.rkelly, rejnumVaR2.rkellyQ,
           rejnumVaR3.con, rejnumVaR3.rkelly, rejnumVaR3.rkellyQ, tmVaR.con, tmVaR.rkelly, tmVaR.rkellyQ,
           eES.con, eES.rkelly, eES.rkellyQ, rejES1.con, rejES1.rkelly, rejES1.rkellyQ,
           rejES2.con, rejES2.rkelly, rejES2.rkellyQ, rejES3.con, rejES3.rkelly, rejES3.rkellyQ,
           rejnumES1.con, rejnumES1.rkelly, rejnumES1.rkellyQ, rejnumES2.con, rejnumES2.rkelly, rejnumES2.rkellyQ,
           rejnumES3.con, rejnumES3.rkelly, rejnumES3.rkellyQ, tmES.con, tmES.rkelly, tmES.rkellyQ,
           VaR.lambda.rkelly, VaR.lambda.rkellyQ, ES.lambda.rkelly, ES.lambda.rkellyQ))
}


# Multiple-core computations
out <- do.call(rbind, mclapply(1:NumRuns, simu, mc.cores = NumCores))
save(out, file="Result_iid.RDATA")


# Output values
load("Result_iid.RDATA")
avg.total = c()
for(i in 1:(210+20*n)){
  avg = mean(out[,i][out[,i] != n])
  avg.total = cbind(avg.total, avg)
}
avg.eVaR.con = matrix(avg.total[1:4],nrow=2,ncol=2)
avg.eVaR.rkelly = matrix(avg.total[5:8],nrow=2,ncol=2)
avg.eVaR.rkellyQ = matrix(avg.total[9:12],nrow=2,ncol=2) # final e-values for VaR
avg.rejVaR1.con = matrix(avg.total[13:16],nrow=2,ncol=2)
avg.rejVaR1.rkelly = matrix(avg.total[17:20],nrow=2,ncol=2)
avg.rejVaR1.rkellyQ = matrix(avg.total[21:24],nrow=2,ncol=2) # numbers of days where we reject for VaR (threshold 2)
avg.rejVaR2.con = matrix(avg.total[25:28],nrow=2,ncol=2)
avg.rejVaR2.rkelly = matrix(avg.total[29:32],nrow=2,ncol=2)
avg.rejVaR2.rkellyQ = matrix(avg.total[33:36],nrow=2,ncol=2) # numbers of days where we reject for VaR (threshold 5)
avg.rejVaR3.con = matrix(avg.total[37:40],nrow=2,ncol=2)
avg.rejVaR3.rkelly = matrix(avg.total[41:44],nrow=2,ncol=2)
avg.rejVaR3.rkellyQ = matrix(avg.total[45:48],nrow=2,ncol=2) # numbers of days where we reject for VaR (threshold 10)
avg.rejnumVaR1.con = matrix(avg.total[49:52],nrow=2,ncol=2)
avg.rejnumVaR1.rkelly = matrix(avg.total[53:56],nrow=2,ncol=2)
avg.rejnumVaR1.rkellyQ = matrix(avg.total[57:60],nrow=2,ncol=2) # numbers of rejections for VaR (threshold 2)
avg.rejnumVaR2.con = matrix(avg.total[61:64],nrow=2,ncol=2)
avg.rejnumVaR2.rkelly = matrix(avg.total[65:68],nrow=2,ncol=2)
avg.rejnumVaR2.rkellyQ = matrix(avg.total[69:72],nrow=2,ncol=2) # numbers of rejections for VaR (threshold 5)
avg.rejnumVaR3.con = matrix(avg.total[73:76],nrow=2,ncol=2)
avg.rejnumVaR3.rkelly = matrix(avg.total[77:80],nrow=2,ncol=2)
avg.rejnumVaR3.rkellyQ = matrix(avg.total[81:84],nrow=2,ncol=2) # numbers of rejections for VaR (threshold 10)
avg.tmVaR.con = matrix(avg.total[85:(84+2*n)],nrow=2,ncol=n)
avg.tmVaR.rkelly = matrix(avg.total[(85+2*n):(84+4*n)],nrow=2,ncol=n)
avg.tmVaR.rkellyQ = matrix(avg.total[(85+4*n):(84+6*n)],nrow=2,ncol=n) # sequential test martingale for VaR
avg.eES.con = matrix(avg.total[(85+6*n):(90+6*n)],nrow=2,ncol=3)
avg.eES.rkelly = matrix(avg.total[(91+6*n):(96+6*n)],nrow=2,ncol=3)
avg.eES.rkellyQ = matrix(avg.total[(97+6*n):(102+6*n)],nrow=2,ncol=3) # final e-values for ES
avg.rejES1.con = matrix(avg.total[(103+6*n):(108+6*n)],nrow=2,ncol=3)
avg.rejES1.rkelly = matrix(avg.total[(109+6*n):(114+6*n)],nrow=2,ncol=3)
avg.rejES1.rkellyQ = matrix(avg.total[(115+6*n):(120+6*n)],nrow=2,ncol=3) # numbers of days where we reject for ES (threshold 2)
avg.rejES2.con = matrix(avg.total[(121+6*n):(126+6*n)],nrow=2,ncol=3)
avg.rejES2.rkelly = matrix(avg.total[(127+6*n):(132+6*n)],nrow=2,ncol=3)
avg.rejES2.rkellyQ = matrix(avg.total[(133+6*n):(138+6*n)],nrow=2,ncol=3) # numbers of days where we reject for ES (threshold 5)
avg.rejES3.con = matrix(avg.total[(139+6*n):(144+6*n)],nrow=2,ncol=3)
avg.rejES3.rkelly = matrix(avg.total[(145+6*n):(150+6*n)],nrow=2,ncol=3)
avg.rejES3.rkellyQ = matrix(avg.total[(151+6*n):(156+6*n)],nrow=2,ncol=3) # numbers of days where we reject for ES (threshold 10)
avg.rejnumES1.con = matrix(avg.total[(157+6*n):(162+6*n)],nrow=2,ncol=3)
avg.rejnumES1.rkelly = matrix(avg.total[(163+6*n):(168+6*n)],nrow=2,ncol=3)
avg.rejnumES1.rkellyQ = matrix(avg.total[(169+6*n):(174+6*n)],nrow=2,ncol=3) # numbers of rejections for ES (threshold 2)
avg.rejnumES2.con = matrix(avg.total[(175+6*n):(180+6*n)],nrow=2,ncol=3)
avg.rejnumES2.rkelly = matrix(avg.total[(181+6*n):(186+6*n)],nrow=2,ncol=3)
avg.rejnumES2.rkellyQ = matrix(avg.total[(187+6*n):(192+6*n)],nrow=2,ncol=3) # numbers of rejections for ES (threshold 5)
avg.rejnumES3.con = matrix(avg.total[(193+6*n):(198+6*n)],nrow=2,ncol=3)
avg.rejnumES3.rkelly = matrix(avg.total[(199+6*n):(204+6*n)],nrow=2,ncol=3)
avg.rejnumES3.rkellyQ = matrix(avg.total[(205+6*n):(210+6*n)],nrow=2,ncol=3) # numbers of rejections for ES (threshold 10)
avg.tmES.con = matrix(avg.total[(211+6*n):(210+8*n)],nrow=2,ncol=n)
avg.tmES.rkelly = matrix(avg.total[(211+8*n):(210+10*n)],nrow=2,ncol=n)
avg.tmES.rkellyQ = matrix(avg.total[(211+10*n):(210+12*n)],nrow=2,ncol=n) # sequential test martingale for ES
avg.VaR.lambda.rkelly = matrix(avg.total[(211+12*n):(210+14*n)],nrow=2,ncol=n)
avg.VaR.lambda.rkellyQ = matrix(avg.total[(211+14*n):(210+16*n)],nrow=2,ncol=n)
avg.ES.lambda.rkelly = matrix(avg.total[(211+16*n):(210+18*n)],nrow=2,ncol=n)
avg.ES.lambda.rkellyQ = matrix(avg.total[(211+18*n):(210+20*n)],nrow=2,ncol=n)


out.final = list(avg.eVaR.con = avg.eVaR.con, avg.eVaR.rkelly = avg.eVaR.rkelly, avg.eVaR.rkellyQ = avg.eVaR.rkellyQ,
                 avg.rejVaR1.con = avg.rejVaR1.con, avg.rejVaR1.rkelly = avg.rejVaR1.rkelly, avg.rejVaR1.rkellyQ = avg.rejVaR1.rkellyQ,
                 avg.rejVaR2.con = avg.rejVaR2.con, avg.rejVaR2.rkelly = avg.rejVaR2.rkelly, avg.rejVaR2.rkellyQ = avg.rejVaR2.rkellyQ,
                 avg.rejVaR3.con = avg.rejVaR3.con, avg.rejVaR3.rkelly = avg.rejVaR3.rkelly, avg.rejVaR3.rkellyQ = avg.rejVaR3.rkellyQ,
                 avg.rejnumVaR1.con= avg.rejnumVaR1.con, avg.rejnumVaR1.rkelly = avg.rejnumVaR1.rkelly, avg.rejnumVaR1.rkellyQ = avg.rejnumVaR1.rkellyQ,
                 avg.rejnumVaR2.con= avg.rejnumVaR2.con, avg.rejnumVaR2.rkelly = avg.rejnumVaR2.rkelly, avg.rejnumVaR2.rkellyQ = avg.rejnumVaR2.rkellyQ,
                 avg.rejnumVaR3.con= avg.rejnumVaR3.con, avg.rejnumVaR3.rkelly = avg.rejnumVaR3.rkelly, avg.rejnumVaR3.rkellyQ = avg.rejnumVaR3.rkellyQ,
                 avg.tmVaR.con = avg.tmVaR.con, avg.tmVaR.rkelly = avg.tmVaR.rkelly, avg.tmVaR.rkellyQ = avg.tmVaR.rkellyQ,
                 avg.eES.con = avg.eES.con, avg.eES.rkelly = avg.eES.rkelly, avg.eES.rkellyQ = avg.eES.rkellyQ,
                 avg.rejES1.con = avg.rejES1.con, avg.rejES1.rkelly = avg.rejES1.rkelly, avg.rejES1.rkellyQ = avg.rejES1.rkellyQ,
                 avg.rejES2.con = avg.rejES2.con, avg.rejES2.rkelly = avg.rejES2.rkelly, avg.rejES2.rkellyQ = avg.rejES2.rkellyQ,
                 avg.rejES3.con = avg.rejES3.con, avg.rejES3.rkelly = avg.rejES3.rkelly, avg.rejES3.rkellyQ = avg.rejES3.rkellyQ,
                 avg.rejnumES1.con = avg.rejnumES1.con, avg.rejnumES1.rkelly = avg.rejnumES1.rkelly, avg.rejnumES1.rkellyQ = avg.rejnumES1.rkellyQ,
                 avg.rejnumES2.con = avg.rejnumES2.con, avg.rejnumES2.rkelly = avg.rejnumES2.rkelly, avg.rejnumES2.rkellyQ = avg.rejnumES2.rkellyQ,
                 avg.rejnumES3.con = avg.rejnumES3.con, avg.rejnumES3.rkelly = avg.rejnumES3.rkelly, avg.rejnumES3.rkellyQ = avg.rejnumES3.rkellyQ,
                 avg.tmES.con = avg.tmES.con, avg.tmES.rkelly = avg.tmES.rkelly, avg.tmES.rkellyQ = avg.tmES.rkellyQ,
                 avg.VaR.lambda.rkelly = avg.VaR.lambda.rkelly, avg.VaR.lambda.rkellyQ = avg.VaR.lambda.rkellyQ,
                 avg.ES.lambda.rkelly = avg.ES.lambda.rkelly, avg.ES.lambda.rkellyQ = avg.ES.lambda.rkellyQ)
save(out.final, file="Result_iid_final.RDATA")


avg.eVaR.con
avg.eVaR.rkelly
avg.eVaR.rkellyQ
avg.rejVaR1.con
avg.rejVaR1.rkelly
avg.rejVaR1.rkellyQ
avg.rejVaR2.con
avg.rejVaR2.rkelly
avg.rejVaR2.rkellyQ
avg.rejVaR3.con
avg.rejVaR3.rkelly
avg.rejVaR3.rkellyQ
avg.rejnumVaR1.con
avg.rejnumVaR1.rkelly
avg.rejnumVaR1.rkellyQ
avg.rejnumVaR2.con
avg.rejnumVaR2.rkelly
avg.rejnumVaR2.rkellyQ
avg.rejnumVaR3.con
avg.rejnumVaR3.rkelly
avg.rejnumVaR3.rkellyQ
avg.eES.con
avg.eES.rkelly
avg.eES.rkellyQ
avg.rejES1.con
avg.rejES1.rkelly
avg.rejES1.rkellyQ
avg.rejES2.con
avg.rejES2.rkelly
avg.rejES2.rkellyQ
avg.rejES3.con
avg.rejES3.rkelly
avg.rejES3.rkellyQ
avg.rejnumES1.con
avg.rejnumES1.rkelly
avg.rejnumES1.rkellyQ
avg.rejnumES2.con
avg.rejnumES2.rkelly
avg.rejnumES2.rkellyQ
avg.rejnumES3.con
avg.rejnumES3.rkelly
avg.rejnumES3.rkellyQ


# Plot the e-processes

load("Result_iid_final.RDATA")

setEPS()
postscript("VaR95_iid.eps")
matplot(1:1000, cbind(out.final$avg.tmVaR.con[1,1:1000], out.final$avg.tmVaR.rkelly[1,1:1000], out.final$avg.tmVaR.rkellyQ[1,1:1000]),
        ylim=e.lim, col=c("red", "blue", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("lambda=.025", "GREE/GREL", "GRO"),
       col=c("red", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("VaR99_iid.eps")
matplot(1:1000, cbind(out.final$avg.tmVaR.con[2,1:1000], out.final$avg.tmVaR.rkelly[2,1:1000], out.final$avg.tmVaR.rkellyQ[2,1:1000]),
        ylim=e.lim, col=c("red", "blue", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("lambda=.025", "GREE/GREL", "GRO"),
       col=c("red", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES875_iid.eps")
matplot(1:1000, cbind(out.final$avg.tmES.con[1,1:1000], out.final$avg.tmES.rkelly[1,1:1000], out.final$avg.tmES.rkellyQ[1,1:1000]),
        ylim=e.lim, col=c("red", "blue", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("lambda=.025", "GREE/GREL", "GRO"),
       col=c("red", "blue", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("ES975_iid.eps")
matplot(1:1000, cbind(out.final$avg.tmES.con[2,1:1000], out.final$avg.tmES.rkelly[2,1:1000], out.final$avg.tmES.rkellyQ[2,1:1000]),
        ylim=e.lim, col=c("red", "blue", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("lambda=.025", "GREL", "GRO"),
       col=c("red", "blue", "black"), lwd = 1, cex=1)
dev.off()

# setEPS()
# postscript("VaR95_iid_lambda.eps")
# matplot(1:1000, cbind(out.final$avg.VaR.lambda.rkelly[1,1:1000], out.final$avg.VaR.lambda.rkellyQ[1,1:1000]),
#         col=c("red", "black"), type = "l", lty = 1,
#         xlab = "number of data", ylab = "lambda", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
# legend("topleft", legend=c("functional Kelly (unknown Q)", "functional Kelly (known Q)"),
#        col=c("red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("VaR99_iid_lambda.eps")
# matplot(1:1000, cbind(out.final$avg.VaR.lambda.rkelly[2,1:1000], out.final$avg.VaR.lambda.rkellyQ[2,1:1000]),
#         col=c("red", "black"), type = "l", lty = 1,
#         xlab = "number of data", ylab = "lambda", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
# legend("topleft", legend=c("functional Kelly (unknown Q)", "functional Kelly (known Q)"),
#        col=c("red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES875_iid_lambda.eps")
# matplot(1:1000, cbind(out.final$avg.ES.lambda.rkelly[1,1:1000], out.final$avg.ES.lambda.rkellyQ[1,1:1000]),
#         col=c("red", "black"), type = "l", lty = 1,
#         xlab = "number of data", ylab = "lambda", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
# legend("topleft", legend=c("functional Kelly (unknown Q)", "functional Kelly (known Q)"),
#        col=c("red", "black"), lwd = 1, cex=1)
# dev.off()
# 
# setEPS()
# postscript("ES975_iid_lambda.eps")
# matplot(1:1000, cbind(out.final$avg.ES.lambda.rkelly[2,1:1000], out.final$avg.ES.lambda.rkellyQ[2,1:1000]),
#         col=c("red", "black"), type = "l", lty = 1,
#         xlab = "number of data", ylab = "lambda", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
# legend("topleft", legend=c("functional Kelly (unknown Q)", "functional Kelly (known Q)"),
#        col=c("red", "black"), lwd = 1, cex=1)
# dev.off()