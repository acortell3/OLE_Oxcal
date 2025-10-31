######## SCRIPT TO DO A PHASE MODEL USING NIMBLE CARBON. WE WILL USE TRAPEZOIDAL AND UNIFORM MODEL PRIORS

## Load dates
dates <- readRDS("../Simu_data/Uncal_YearsBP_100dates.rds")
colnames(dates) <- c("Year", "Sd")

## Load required packages
library(rcarbon)
library(nimbleCarbon)
library(coda)

## UNIFORM
# Simulate Observed Data
n <- nrow(dates)
cra <- round(dates$Year)
cra.error <- dates$Sd

### UNIFORM

# Define NIMBLE Model
phasemodel <- nimbleCode({
  for (i in 1:N){
    #  Likelihood
    theta[i] ~ dunif(alpha[1],alpha[2]);
    # Calibration
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sd[i] <- (sigma[i]^2+sigmaCurve[i]^2)^(1/2);
    X[i] ~ dnorm(mean=mu[i],sd=sd[i]);
  }
  # Prior
  alpha[1] ~ dunif(0,50000);
  alpha[2] ~ T(dunif(0,50000),alpha[1],50000)
})  

#define constant, data, and inits:
data("intcal20") 
constants <- list(N=n,calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma)
data <- list(X=cra,sigma=cra.error)
m.dates = medCal(calibrate(cra,cra.error,verbose = FALSE))
inits <- list(alpha=c(4000,8000),theta=m.dates)

#Run MCMC
mcmc.samples<- nimbleMCMC(code = phasemodel,constants = constants,data = data,niter = 3000000, nchains = 3, thin=3000, nburnin = 600000, progressBar = FALSE, monitors=c('alpha','theta'), inits=inits, samplesAsCodaMCMC=TRUE,setSeed = c(1,2,3), WAIC = TRUE)

statistics_uniform <- list("rhat" = gelman.diag(mcmc.samples$samples),
			   "ESS" = effectiveSize(mcmc.samples$samples))

saveRDS(mcmc.samples, "../Results/NC_uniform.rds")
saveRDS(statistics_uniform, "../Results/NC_uni_stats.rds")

#Plot Posteriors
#par(mfrow=c(1,2))
#postHPDplot(mcmc.samples[,'alpha[2]'],xlim=c(6400,5800),xlab='Cal BP',ylab='Posterior Probability',main='Start of Phase')
#abline(v=a,lty=2)
#postHPDplot(mcmc.samples[,'alpha[1]'],xlim=c(5600,5300),xlab='Cal BP',ylab='Posterior Probability',main='End of Phase')
#abline(v=b,lty=2)


## TRAPEZOID
# Simulate Observed Data
n <- nrow(dates)
cra <- round(dates$Year)
cra.error <- dates$Sd

# Define NIMBLE Model
phasemodel <- nimbleCode({
  for (i in 1:N){
    #  Likelihood
    theta[i] ~ dTrapezoidal(a=a,m1=m1,m2=m2,b=b);
    # Calibration
    mu[i] <- interpLin(z=theta[i], x=calBP[], y=C14BP[]);
    sigmaCurve[i] <- interpLin(z=theta[i], x=calBP[], y=C14err[]);
    sd[i] <- (sigma[i]^2+sigmaCurve[i]^2)^(1/2);
    X[i] ~ dnorm(mean=mu[i],sd=sd[i]);
  }
  
  # Prior
  a ~ dunif(0,10000);
  m1 ~ dunif(0,10000);
  m2 ~ dunif(0,10000);
  b ~ dunif(0,10000);
  cons ~ dconstraint(a > m1 & m1 > m2 & m2 > b);

})  

#define constant, data, and inits:
data("intcal20") 
constants <- list(N=n,calBP=intcal20$CalBP,C14BP=intcal20$C14Age,C14err=intcal20$C14Age.sigma,cons=1)
data <- list(X=cra,sigma=cra.error)
m.dates = medCal(calibrate(cra,cra.error,verbose = FALSE))
inits <- list(a = 6500, m1 = 5500, m2 = 4500, b = 3500,theta=m.dates)

#Run MCMC
mcmc.samples<- nimbleMCMC(code = phasemodel,constants = constants,data = data,niter = 1000000, nchains = 3, thin=1000, nburnin = 200000, progressBar = FALSE, monitors=c('a','m1','m2','b','theta'), inits=inits, samplesAsCodaMCMC=TRUE,setSeed= c(1,2,3), WAIC = TRUE)

statistics_trapezoid <- list("rhat" = gelman.diag(mcmc.samples$samples),
			   "ESS" = effectiveSize(mcmc.samples$samples))

saveRDS(mcmc.samples, "../Results/NC_trapezoid.rds")
saveRDS(statistics_trapezoid, "../Results/NC_tra_stats.rds")
#Plot Posteriors
#par(mfrow=c(1,2))
#postHPDplot(mcmc.samples[,'a'],xlim=c(7000,3000),xlab='Cal BP',ylab='Posterior Probability',main='Start of Phase')
#abline(v=a,lty=2)
#postHPDplot(mcmc.samples[,'b'],xlim=c(5600,5300),xlab='Cal BP',ylab='Posterior Probability',main='End of Phase')
#abline(v=b,lty=2)

