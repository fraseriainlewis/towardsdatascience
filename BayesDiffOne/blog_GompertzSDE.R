
#########################################################################
############ GENERAL SETUP - get the data ###############################
#########################################################################
rm(list=ls()); # reset workspace

# clone this repo
# git clone https://github.com/fraseriainlewis/towardsdatascience.git
# 
# 
setwd("towardsdatascience/BayesDiffOne");

# this loads in the data in a format needed for the modeling
# this provides data.frame "data1" 
source("data_mgt.R"); # provide data.frame data1

# We now have the data
#########################################################################
############ MODEL DEFINITION + MCMC setup in R-nimble  #################
#########################################################################
#
# Computational Note: 
#
# As coded below this fits 8 chains, each of 5 000 000 iterations and thin=500 
# parallelized with one chain per cpu - so assumes 8 cpus. 
#
# This can take 60 mins across 8 cores. To get a quick example working 
# then reduce length of chains, and also set cores to how many
# are physically available. 
#
#########################################################################
# the general structure of the use of parallel is taken from an example
# on the R-nimble site. Basic idea is that the model + mcmc config is 
# wrapped into a function which is then distributed across nodes. 
# 
print(Sys.time())

library(parallel)

Nnodes<-8; # set to <= number of physical cores

this_cluster <- makeCluster(Nnodes) # uses parallel library

run_MCMC_allcode <- function(seeds,dat) {
###########################################
library(nimble, warn.conflicts = FALSE)
###########################################

  ## set up custom density - this is the slice density
  ## how to setup custom densities is documented in the R-nimble manual
  ## note the formula is for the log of the density
  ## this function does not need registered
  ## dXXXX = d is for density
  dmysde <- nimbleFunction(
    run = function(x = double(0),                  # random variable value
                   x0 = double(0, default = 0),
                   t = double(0, default = 0),
                   alpha = double(0, default = 0),    
                   K = double(0, default = 0),
                   sigma = double(0, default = 1), # standard deviation
                   
                   log = integer(0, default = 0)) {
                         returnType(double(0))
                         # log pdf
                         logProb <- (1/2)*(
                                            -log(2*pi*t) - 2*log(x) - ((
                                                  ((t*sigma^2)/2) - log(K) + log(x) + 
                                                  (exp(-t*alpha))*(log(K) - log(x0)))^2)/(t*sigma^2) - 
                                          2*log(sigma)
                                          )
        if(log) return(logProb)
        else return(exp(logProb)) 
    })

  ## this is not strictly needed but here as a check - the mcmc sampler should generate the variates
  ## and so a stop() is used here to catch if this is used, in which case we need to switch sampler type, to RW
  rmysde <- nimbleFunction(
    run = function(n = integer(0),
                   x0 = double(0, default = 0),
                    t = double(0, default = 0),
                   alpha = double(0, default = 0),    
                   K = double(0, default = 0),
                   sigma = double(0, default = 1)) {
                                      stop('rmynorm should never be run')
                                      returnType(double(0))
                             return(0)
    })
  
  ## if using parallel then this essential to ensure each compute node has the custom functions
  ## in scope. If not using parallel this is unnecessary
  assign('dmysde', dmysde, envir = .GlobalEnv)
  assign('rmysde', rmysde, envir = .GlobalEnv)


  ## This is the model definition - a BUGS model, format similar to JAGS, OpenBUGS  
  code <- nimbleCode({
  
    sigma1 ~ dunif(0, 1000) # flat prior for sd
    sigma2 ~ dunif(0, 1000) # flat prior for sd
    sigma3 ~ dunif(0, 1000) # flat prior for sd
  
    a ~ dunif(0,1)  # flat prior for Gompertz parameter - must be non-negative 
    K ~ dnorm(0,sd=1000) # diffuse Gaussian prior, must be non-negative but is a large value
  
    mu_t ~ dnorm(0,sd=1000) # diffuse Gaussian prior for mean
    sigma_t ~ dunif(0, 1000) # flat prior for sd
 
      for(i in 1:n) { # for each observation - each transition from X(t) to X(t+1) given X(t)
     
         # the model includes iid Normal errors - measurement error
         # without measurement error then only the third line below is needed
        
         Observation1_latent[i]  ~ dnorm(mu_t,sd=sigma_t); # a prior as Observation1_latent is a random variable 
         Observation1[i] ~ dnorm(Observation1_latent[i],sd=sigma1); # Observation1 = data likelihood with sd=sigma1
         
         # this is the SDE part, and note that Observation2_latent is a random variable - not data as it's latent due to 
         # measurement error. Sigma3 = the sd for the Brownian motion diffusion
         Observation2_latent[i] ~ dmysde(x0=Observation1_latent[i],t=Time2[i]-Time1[i],alpha=a,K=K,sigma=sigma3);  
      
         # this is data likelihood with sd=sigma2
         Observation2[i]~dnorm(Observation2_latent[i],sd=sigma2); 
     
      }
  
   # we want to store forecasts - 20 time points, each with a prediction posterior density
   # we create these as separate vectors to avoid monitoring every data point (which are static in any case)
   yPred2_latent[1:20] <- Observation2_latent[231:250] # these indexes are hard-coded based on dataset data1
   yPred2[1:20] <- Observation2[231:250]
  
   # as above but we also want to estimate 10 fitted values - as opposed to forecasts
   # not all fitted values, just enough to give a flavour 
   yfit2_latent[1:10] <- Observation2_latent[251:260] # these indexes are hard-coded based on dataset data1
   yfit2[1:10] <- Observation2[251:260]
   
  })

  ## this line is IMPORTANT - we pass a list to the outer MCMC function and then extract out the data
  data1=dat[[1]];
  
  constants <- list(n = nrow(data1))#
  data <- list(Observation1 = data1$Observation1,Observation2 = data1$Observation2,
             Time1 = data1$Time1,Time2=data1$Time2)
  inits <- list(a=0.1,K=1000,sigma2=10,sigma3=0.1,sigma1=0.1,mu_t=0.1,sigma_t=1.0,
              Observation1_latent=rep(100,nrow(data1)),
              Observation2_latent=rep(100,nrow(data1)),
              Observation2=rep(100,nrow(data1))
              )
              

  model <- nimbleModel(code, constants = constants, data = data1, inits = inits) # build model

  myMCMCconf <- configureMCMC(model) 

  myMCMCconf$addMonitors('yPred2[]')
  myMCMCconf$addMonitors('yPred2_latent[]')
  myMCMCconf$addMonitors('yfit2[]')
  myMCMCconf$addMonitors('yfit2_latent[]')

  myMCMC<-buildMCMC(myMCMCconf)

  Cmodel <- compileNimble(model,resetFunctions = TRUE,showCompilerOutput = FALSE)
  CmodelMCMC <- compileNimble(myMCMC, project = Cmodel,resetFunctions = TRUE,showCompilerOutput = FALSE)

  ## note seeds is passed as argument to outer function
  ## use coda format since easier for multiple chains
  MCMCsamples <- runMCMC(CmodelMCMC, niter = 5000000,thin=500,nchains=1,setSeed=seeds,
                   samplesAsCodaMCMC = TRUE)

return(MCMCsamples);
}


#########################################################################
############ Run the actual MCMC sampler                #################
#########################################################################

chain_output <- parLapply(cl = this_cluster, X = c(1:Nnodes)+1000, ## NOTE this is the seeds argument!
                          fun = run_MCMC_allcode, 
                          dat=list(data1) # this is one extra argument
                          )

# It's good practice to close the cluster when you're done with it.
stopCluster(this_cluster)

print(Sys.time());

library(coda);
myresGp<-mcmc.list(chain_output[[1]],
                   chain_output[[2]],
                   chain_output[[3]],
                   chain_output[[4]],
                   chain_output[[5]],
                   chain_output[[6]],
                   chain_output[[7]],
                   chain_output[[8]])

# load("mcmcSDE_set2.RData"); #provides myresGp
save(myresGp,file="mcmcSDE_set2.RData");

#########################################################################
### Analysis and summary of MCMC results  ###############################
#########################################################################
# some basic good practice MCMC checks

# look for burn-in
plot(window(myresGp[,c("K","a"),],start=1))

# 2K burnin seems fine (on 500 thin)
myresGp2<-window(myresGp,start=2000); # overwrite dropping burn-in

# check for convergence
gelman.diag(myresGp2); # looks pretty good
# fine.

summary(myresGp2)

# now compute posterior prediction intervals for forecasts
myFcast<-summary(myresGp2[,c(8:27),]);# 8:27, actual, 28:47 = latent, 48:57, fits
lower<-myFcast$quantiles[,"2.5%"]
mid<-myFcast$quantiles[,"50%"]
upper<-myFcast$quantiles[,"97.5%"]

# now compute posterior prediction intervals for fitted values
myFits<-summary(myresGp2[,c(48:57),]);# 8:27, actual, 28:47 = latent, 48:57, fits
midFit<-myFits$quantiles[,"50%"]


#########################################################################
### A summary plot                        ###############################
#########################################################################
############# plot
interactive<-FALSE;
if(interactive){png("plotNLSDE1.png",width=480*1.5,height=480*1.5,pointsize=14);}

# Setup plotting environment - uses classic plot() and par()
mycol<-"black";
par(mar=c(5,6,5,5));
par(cex.axis=2.0);par(cex.lab=2.5);par(bg="white");par(fg=mycol);par(col.axis=mycol);par(col=mycol);par(col.main=mycol);
par(cex.main=2.0);par(col.lab=mycol);par(las=1);par(xaxs="r");par(yaxs="r");
par(mfrow=c(1,1));
# end setup


plot(data1$Time2[1:231],data1$Observation2[1:231],pch=21,col="green",bg="darkgreen",type="n",axes=FALSE,ylab="",xlab="",xlim=c(0,25))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
grid(col="grey60")

lines(data1$Time2[1:231],data1$Observation2[1:231],lwd=2,lty=1,col="green");
points(data1$Time2[1:231],data1$Observation2[1:231],pch=23,col="green",bg="darkgreen");


axis(1,seq(0,25,5),padj=0.3,cex.axis=1.5);
mtext("Time (days)",1,line=3.5,cex=2);par(las=1);
axis(2,cex.axis=1.5);par(las=3);
mtext(expression('Volume mm' ^3), 2, line=4,cex=2);par(las=1);
axis(4,padj=0.3,cex.axis=1.5);par(las=3);
mtext("",4,line=3,cex=2);
title("Bayesian diffusion forecasting",line=2);par(las=1);


# now plot the posterior predictive intervals 95%
x<-data1$Time2[231:250];#seq(23.1,25.0,by=0.1); # time points
for(i in 1:length(lower)){
lines(c(x[i],x[i]),c(lower[i],upper[i]),pch=22,col="grey60",bg="grey60")
  points(c(x[i]),c(mid[i]),pch="+",col="black",bg="black")

  }

x<-data1$Time1[251:260]
#points(c(x[1],x[1]),c(lower[1],upper[1]),pch=22,col="red",bg="black")
for(i in 1:length(lower)){
  points(c(x[i]),c(midFit[i]),pch=24,col="blue",bg="magenta",cex=1.0)
    }
box();

if(interactive){dev.off();}

#############  The END ###################################################


