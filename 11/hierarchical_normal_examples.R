library(mvtnorm) # for dmvnorm

# our data
y <-  c(62, 60, 63, 59,                   # group A
        63, 67, 71, 64, 65, 66,           # group B
        68, 66, 71, 67, 68, 68,           # group C
        56, 62, 60, 61, 63, 64, 63, 59)   # group D
numParams <- 7
numThetas <- 4

#################################
## Metropolis Algorithm
################################

# note that the parameter space we sample in is for
#  thetas, mu, logSigma, logTau

# unnormalized posterior
log_unnorm_post <- function(paramPack, yDat){
  tau <- exp(paramPack$logTau)
  sigma <- exp(paramPack$logSigma)
  return(
    -.0002 * paramPack$logSigma - .0001/sigma/sigma +
      -.0002 * paramPack$logTau - .0001/tau/tau +
      dmvnorm(paramPack$thetas, rep(paramPack$mu, each = numThetas),  sigma = tau^2 * diag(numThetas), log = T) +
      dmvnorm(x = yDat, mean = rep(paramPack$thetas, times = c(4,6,6,8)), sigma = sigma * sigma * diag(length(yDat)), log = T)
  )
}

# function that samples new parameters
# we won't hardcode tuning parameters

# try 1 first
qSigmaCovDiagSqrt <- .5  ## PRIMARY TUNING PARAMETER

sampleQ <- function(oldParams){
  newLogTau <- rnorm(n = 1, mean = oldParams$logTau, sd = qSigmaCovDiagSqrt)
  newLogSigma <- rnorm(n = 1, mean = oldParams$logSigma, sd = qSigmaCovDiagSqrt)
  newMu <- rnorm(n = 1, mean = oldParams$mu, sd = qSigmaCovDiagSqrt)
  newThetas <- rnorm(n = numThetas, mean = oldParams$thetas, sd = rep(qSigmaCovDiagSqrt, numThetas))
  return(list(thetas = newThetas, mu = newMu, logSigma = newLogSigma, logTau = newLogTau))
}

# start sampling
startParams <- list()
startParams$thetas <- c(60, 60, 60, 60)
startParams$mu <- 60
startParams$logSigma <- 0 
startParams$logTau <- 0
numIters <- 50000
storeEveryMultiple <- 100
samples <- vector(mode="list", numIters/storeEveryMultiple)
samples[[1]] <- startParams
currentParams <- startParams
numAcceptances <- 0
for(i in 2:numIters){
  proposal <- sampleQ(currentParams)
  logAcceptRatio <- log_unnorm_post(proposal, y) - log_unnorm_post(currentParams, y)
  if( log(runif(n=1)) < logAcceptRatio ){
    currentParams <- proposal # change currentParams for next iter
    numAcceptances <- numAcceptances + 1
  }else{
    # don't change currentParams
  }
  
  # only store every tenth draw
  if(i %% storeEveryMultiple == 0){
    cat("iter: ", i, " out of ", numIters, ". ")
    cat("current log unnorm post dens: ", log_unnorm_post(currentParams, y), "\n")
    samples[[i/storeEveryMultiple + 1]] <- currentParams
  } 
}


cat("Completed. Final acceptance rate: ", numAcceptances/numIters, "\n")
# extract things and plot them
thetaSamples <- t(sapply(samples, "[[", 1))
plot.ts(thetaSamples)
library(GGally)
ggpairs(as.data.frame(thetaSamples))

muSamples <- sapply(samples, "[[", 2)
plot.ts(muSamples)
hist(muSamples)

logSigmaSamples <- sapply(samples, "[[", 3)
plot.ts(logSigmaSamples)
hist(exp(2*logSigmaSamples), xlab = "sigma squared")

logTauSamples <- sapply(samples, "[[", 4)
plot.ts(logTauSamples)
hist(exp(logTauSamples), xlab = "tau")

#
pairs(cbind(logTauSamples, logSigmaSamples, muSamples))



#################################
## Gibbs
################################

library(invgamma)
rm(samples, thetaSamples)

# assume the paramPack elements are
# thetas, mu, tau, sigma


numIters <- 1000
J <- numThetas 
N <- length(y)
ybars <- c(mean(y[1:4]), mean(y[5:10]), mean(y[11:16]), mean(y[17:24]))
njs <- c(4,6,6,8)

# functions that sample one thing at a time
# no tuning required!
sampleThetas <- function(paramPack, whichTheta, y){
  nj <- njs[whichTheta]
  ybardotj <- ybars[whichTheta]
  precision <- nj/paramPack$sigma^2 + 1/paramPack$tau^2
  thisMean <- (ybardotj*nj/paramPack$sigma^2 + paramPack$mu/paramPack$tau^2)/precision
  returnPack <- paramPack
  returnPack$thetas[whichTheta] <- rnorm(n = 1, mean = thisMean, sd = sqrt(1/precision))
  return(returnPack)
}
sampleMu <- function(paramPack, y){
  thetaBar <- mean(paramPack$thetas)
  precision <- J/paramPack$tau^2 + 1/100
  thisMean <- (thetaBar*J/paramPack$tau^2 + 60/100)/precision
  returnPack <- paramPack
  returnPack$mu <- rnorm(n = 1, mean = thisMean, sd = sqrt(1/precision) )
  return(returnPack)
}
sampleTau <- function(paramPack, y){
  tauSq <- rinvgamma(n = 1, shape = .0001 + J/2, scale = .0001 + sum((paramPack$thetas - paramPack$mu)^2)/2)
  returnPack <- paramPack
  returnPack$tau <- sqrt(tauSq)
  return(returnPack)
}
sampleSigma <- function(paramPack, y){
  sigSq <- rinvgamma(n = 1, shape = .0001 + N/2, scale = .0001 + sum((y - rep(paramPack$thetas, times = c(4,6,6,8)))^2)/2)
  returnPack <- paramPack
  returnPack$sigma <- sqrt(sigSq)
  return(returnPack)
}

# commence sampling
currentParam <- list()
currentParam$thetas <- c(60, 60, 60, 60)
currentParam$mu <- 60
currentParam$sigma <- 1
currentParam$tau <- 1
loggedParams <- vector(mode = "list", length = numIters)

for(i in 1:numIters){

  cat("iteration ", i, " out of ", numIters, "\n")
  
  # step 1: change thetas
  for(j in 1:numThetas){
    currentParam <- sampleThetas(currentParam, j, y)
  }
  
  # step 2: change mu
  currentParam <- sampleMu(currentParam, y)
  
  # step 3: change tauSquared
  currentParam <- sampleTau(currentParam, y)
  
  # step 4: change sigma
  currentParam <- sampleSigma(currentParam, y)
  
  # store every complete sweep through the parameters
  loggedParams[[i]] <- currentParam
  
}


# plots
thetaSamples <- t(sapply(loggedParams, "[[", 1))
plot.ts(thetaSamples)

muSamples <- sapply(loggedParams, "[[", 2)
plot.ts(muSamples)

sigmaSquaredSamples <- sapply(loggedParams, "[[", 3)^2
plot.ts(sigmaSquaredSamples)

tauSamples <- sapply(loggedParams, "[[", 4)
plot.ts(tauSamples)
