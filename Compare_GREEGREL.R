library("purrr")
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

# Set parameters
n = 1000 # sample size
delay = 10
e.lim <- c(-2,6) # bounds for e-values

# Simulate normal distribution
set.seed(233)
simnorm = rnorm(n+delay)

# Compare methods in Example 7 of the paper

# plot the loss and forecast figures

# Example where GREE outperforms GREL
eps = (1:(n+delay)) / (n+delay)
y = simnorm * (1+eps)

# calculate VaR and ES
VaRoutb = 1.48 * (1 + eps)
ESout = 1.86 * (1 + eps)
EStrue = (ddist("norm", qdist("norm",p=0.95,mu=0,sigma=1), mu=0, sigma=1)) / (1 - 0.95) * (1+eps)


# Plot
setEPS()
postscript("GREE_loss_forecast.eps")
matplot((delay+1):(n+delay), cbind(y[(delay+1):(n+delay)], ESout[(delay+1):(n+delay)], EStrue[(delay+1):(n+delay)]), type = "l",
        lty = c(1,1,2), lwd = c(1,1,2), col = c("black", "red", "blue"), xlab="number of data", ylab="loss and ES forecast")
legend("topleft", legend=c("loss", "ES forecast", "true ES"), col = c("black", "red", "blue"), lty = c(1,1,2), lwd = c(1,1,2), cex=1)
dev.off()


# Example where GREL outperforms GREE
eps = rdunif(n+delay, 5, -5) / 10
y = simnorm

# calculate VaR and ES
VaRoutb = 1.64 + eps
ESout = 2.06 + eps
EStrue = rep((ddist("norm", qdist("norm",p=0.95,mu=0,sigma=1), mu=0, sigma=1)) / (1 - 0.95), n+delay)

# Plot
setEPS()
postscript("GREL_loss_forecast.eps")
matplot((delay+1):(n+delay), cbind(y[(delay+1):(n+delay)], ESout[(delay+1):(n+delay)], EStrue[(delay+1):(n+delay)]), type = "l",
        lty = c(1,1,2), lwd = c(1,1,2), col = c("black", "red", "blue"), xlab="number of data", ylab="loss and ES forecast", ylim = c(-3.5,4))
legend("topleft", legend=c("loss", "ES forecast", "true ES"), col = c("black", "red", "blue"), lty = c(1,1,2), lwd = c(1,1,2), cex=1)
dev.off()


# Example of non-linear trend
eps = sin((1:(n+delay))*0.01)
y = simnorm * (1+eps)

# calculate VaR and ES
VaRoutb = 1.48 * (1 + eps)
ESout = 1.86 * (1 + eps)
EStrue = (ddist("norm", qdist("norm",p=0.95,mu=0,sigma=1), mu=0, sigma=1)) / (1 - 0.95) * (1+eps)

# Plot
setEPS()
postscript("nonlinear_loss_forecast.eps")
matplot((delay+1):(n+delay), cbind(y[(delay+1):(n+delay)], ESout[(delay+1):(n+delay)], EStrue[(delay+1):(n+delay)]), type = "l",
        lty = c(1,1,2), lwd = c(1,1,2), col = c("black", "red", "blue"), xlab="number of data", ylab="loss and ES forecast", ylim = c(-3.5,4))
legend("topleft", legend=c("loss", "ES forecast", "true ES"), inset=c(0.35,0), col = c("black", "red", "blue"), lty = c(1,1,2), lwd = c(1,1,2), cex=1)
dev.off()


# Set number of cores
# NumCores <- detectCores()
NumCores <- 25

# Set number of runs
NumRuns <- 1000



# Function
rseed=sample(1:10^8,NumRuns*2)
rseed = matrix(rseed, nrow  = NumCores, byrow = TRUE)

simu <- function(iter){
  set.seed(iter)
  simnorm = rnorm(n+delay)
  
  # e-testing
  
  # ===================================================
  # Summary for ES
  
  # Output of Example 1
  tmES1.akelly = vector(mode="numeric", length=n+delay) # e-process for GREE
  tmES1.rkelly = vector(mode="numeric", length=n+delay) # e-process for GREL
  tmES1.mix = vector(mode="numeric", length=n+delay) # e-process for GREM
  tmES1.rkellyQ = vector(mode="numeric", length=n+delay) # e-process for GRO
  
  # Output of Example 2
  tmES2.akelly = vector(mode="numeric", length=n+delay)
  tmES2.rkelly = vector(mode="numeric", length=n+delay)
  tmES2.mix = vector(mode="numeric", length=n+delay)
  tmES2.rkellyQ = vector(mode="numeric", length=n+delay)
  
  # Output of Example 3
  tmES3.akelly = vector(mode="numeric", length=n+delay)
  tmES3.rkelly = vector(mode="numeric", length=n+delay)
  tmES3.mix = vector(mode="numeric", length=n+delay)
  tmES3.rkellyQ = vector(mode="numeric", length=n+delay)
  
  # Example 1: GREE outperforms GREL
  eps = (1:(n+delay)) / (n+delay)
  y = simnorm * (1+eps)
  
  # calculate VaR and ES
  VaRoutb = 1.48 * (1 + eps)
  ESout = 1.86 * (1 + eps)
  EStrue = (ddist("norm", qdist("norm",p=0.95,mu=0,sigma=1), mu=0, sigma=1)) / (1 - 0.95) * (1+eps)
  
  tmp2=evalue.ES(x=y, r=ESout, z=VaRoutb,lev=0.95)
  tmp3=evalue.ESQ.eps(x=y, r=ESout, z=VaRoutb, eps = eps, lev=0.95)
  tmES1.akelly[1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
  tmES1.rkelly[1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
  tmES1.mix[1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
  tmES1.rkellyQ[1:length(tmp3$e.out.rkelly)] <- log(tmp3$e.out.rkelly)
  
  # Example 2: GREL outperforms GREE
  eps = rdunif(n+delay, 5, -5) / 10
  y = simnorm
  
  # calculate VaR and ES
  VaRoutb = 1.64 + eps
  ESout = 2.06 + eps
  EStrue = rep((ddist("norm", qdist("norm",p=0.95,mu=0,sigma=1), mu=0, sigma=1)) / (1 - 0.95), n+delay)
  
  tmp2=evalue.ES(x=y, r=ESout, z=VaRoutb,lev=0.95)
  tmp3=evalue.ESQ(x=y, r=ESout, z=VaRoutb, lev=0.95)
  tmES2.akelly[1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
  tmES2.rkelly[1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
  tmES2.mix[1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
  tmES2.rkellyQ[1:length(tmp3$e.out.rkelly)] <- log(tmp3$e.out.rkelly)
  
  # Example 3: non-linear trend
  eps = sin((1:(n+delay))*0.01)
  y = simnorm * (1+eps)
  
  # calculate VaR and ES
  VaRoutb = 1.48 * (1 + eps)
  ESout = 1.86 * (1 + eps)
  EStrue = (ddist("norm", qdist("norm",p=0.95,mu=0,sigma=1), mu=0, sigma=1)) / (1 - 0.95) * (1+eps)
  
  tmp2=evalue.ES(x=y, r=ESout, z=VaRoutb,lev=0.95)
  tmp3=evalue.ESQ.eps(x=y, r=ESout, z=VaRoutb, eps = eps, lev=0.95)
  tmES3.akelly[1:length(tmp2$e.out.akelly)] <- log(tmp2$e.out.akelly)
  tmES3.rkelly[1:length(tmp2$e.out.rkelly)] <- log(tmp2$e.out.rkelly)
  tmES3.mix[1:length(tmp2$e.out.mix)] <- log(tmp2$e.out.mix)
  tmES3.rkellyQ[1:length(tmp3$e.out.rkelly)] <- log(tmp3$e.out.rkelly)
  
  tmES1.akelly <- tail(tmES1.akelly, n)
  tmES1.rkelly <- tail(tmES1.rkelly, n)
  tmES1.mix <- tail(tmES1.mix, n)
  tmES1.rkellyQ <- tail(tmES1.rkellyQ, n)
  tmES2.akelly <- tail(tmES2.akelly, n)
  tmES2.rkelly <- tail(tmES2.rkelly, n)
  tmES2.mix <- tail(tmES2.mix, n)
  tmES2.rkellyQ <- tail(tmES2.rkellyQ, n)
  tmES3.akelly <- tail(tmES3.akelly, n)
  tmES3.rkelly <- tail(tmES3.rkelly, n)
  tmES3.mix <- tail(tmES3.mix, n)
  tmES3.rkellyQ <- tail(tmES3.rkellyQ, n)
  
  return(c(tmES1.akelly, tmES1.rkelly, tmES1.mix, tmES1.rkellyQ, tmES2.akelly, tmES2.rkelly, tmES2.mix, tmES2.rkellyQ, tmES3.akelly, tmES3.rkelly, tmES3.mix, tmES3.rkellyQ))
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
save(out, file="Result_compareGREEGREL.RDATA")


# Output values
load("Result_compareGREEGREL.RDATA")
avg.total = colMeans(out)
avg.tmES1.akelly = avg.total[1:n]
avg.tmES1.rkelly = avg.total[(n+1):(2*n)]
avg.tmES1.mix = avg.total[(2*n+1):(3*n)]
avg.tmES1.rkellyQ = avg.total[(3*n+1):(4*n)]
avg.tmES2.akelly = avg.total[(4*n+1):(5*n)]
avg.tmES2.rkelly = avg.total[(5*n+1):(6*n)]
avg.tmES2.mix = avg.total[(6*n+1):(7*n)]
avg.tmES2.rkellyQ = avg.total[(7*n+1):(8*n)]
avg.tmES3.akelly = avg.total[(8*n+1):(9*n)]
avg.tmES3.rkelly = avg.total[(9*n+1):(10*n)]
avg.tmES3.mix = avg.total[(10*n+1):(11*n)]
avg.tmES3.rkellyQ = avg.total[(11*n+1):(12*n)]

out.final = list(avg.tmES1.akelly = avg.tmES1.akelly, avg.tmES1.rkelly = avg.tmES1.rkelly, avg.tmES1.mix = avg.tmES1.mix, avg.tmES1.rkellyQ = avg.tmES1.rkellyQ,
                 avg.tmES2.akelly = avg.tmES2.akelly, avg.tmES2.rkelly = avg.tmES2.rkelly, avg.tmES2.mix = avg.tmES2.mix, avg.tmES2.rkellyQ = avg.tmES2.rkellyQ,
                 avg.tmES3.akelly = avg.tmES3.akelly, avg.tmES3.rkelly = avg.tmES3.rkelly, avg.tmES3.mix = avg.tmES3.mix, avg.tmES3.rkellyQ = avg.tmES3.rkellyQ)
save(out.final, file="Result_compareGREEGREL_final.RDATA")


load("Result_compareGREEGREL_final.RDATA")
setEPS()
postscript("GREE.eps")
matplot(1:n, cbind(out.final$avg.tmES1.akelly, out.final$avg.tmES1.rkelly, out.final$avg.tmES1.mix, out.final$avg.tmES1.rkellyQ),
        ylim=e.lim, col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("GREE", "GREL", "GREM", "GRO"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("GREL.eps")
matplot(1:n, cbind(out.final$avg.tmES2.akelly, out.final$avg.tmES2.rkelly, out.final$avg.tmES2.mix, out.final$avg.tmES2.rkellyQ),
        ylim=e.lim, col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("GREE", "GREL", "GREM", "GRO"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)
dev.off()

setEPS()
postscript("GREE_nonlinear.eps")
matplot(1:n, cbind(out.final$avg.tmES3.akelly, out.final$avg.tmES3.rkelly, out.final$avg.tmES3.mix, out.final$avg.tmES3.rkellyQ),
        ylim=e.lim, col=c("red", "blue", "darkgreen", "black"), type = "l", lty = 1,
        xlab = "number of data", ylab = "e-process", cex.lab = 1.4, cex.axis = 1.4, pch = 16, lwd = 1)
abline(h = log(2), lwd = 1, lty = 2)
abline(h = log(5), lwd = 1, lty = 2)
abline(h = log(10), lwd = 1, lty = 2)
legend("topleft", legend=c("GREE", "GREL", "GREM", "GRO"),
       col=c("red", "blue", "darkgreen", "black"), lwd = 1, cex=1)
dev.off()