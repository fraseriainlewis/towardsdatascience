###############################################################################################################
# This file has data hard coded from https://zenodo.org/record/3574531#.YG_n3jHisuV
# It them fits an SDE model to this and 
#
###############################################################################################################
# PARALLEL VERSION
#
###############################################################################################################
rm(list=ls())
############## DATA
data1 <-
structure(list(ID = c(38L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 
38L, 38L, 38L, 38L, 38L, 38L), Time1 = c(11, 12, 13, 14, 17, 
18, 19, 20, 25, 26, 27, 28, 32, 33), Observation1 = c(109.08, 
161.68000000000001, 275.70999999999998, 231.21000000000001, 299.63, 
455.69, 540, 713.75999999999999, 723.82000000000005, 819.04999999999995, 
791.57000000000005, 681.34000000000003, 999.44000000000005, 1067.6500000000001
), ID = c(38L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 38L, 
38L, 38L, 38L), Time2 = c(12, 13, 14, 17, 18, 19, 20, 25, 26, 
27, 28, 32, 33, 34), Observation2 = c(161.68000000000001, 275.70999999999998, 
231.21000000000001, 299.63, 455.69, 540, 713.75999999999999, 
723.82000000000005, 819.04999999999995, 791.57000000000005, 681.34000000000003, 
999.44000000000005, 1067.6500000000001, 1132.1600000000001)), class = "data.frame", row.names = 331:344)

library(parallel)

Nnodes<-4;

this_cluster <- makeCluster(Nnodes)

run_MCMC_allcode <- function(seeds,dat) {
###########################################
library(nimble, warn.conflicts = FALSE)
###########################################

## custom density dXXXX
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

## this is not strictly needed but here as a check - the mcmc sampler should generate the variate
## and so a stop() is used here to catch if this is used, in which case we need to switch sampler type
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
  
if(FALSE){ # note - not actually necessary to register 
  deregisterDistributions('dmysde')

registerDistributions(list(
    dmysde = list(
        BUGSdist = "dmysde(x0,t,alpha,K,sigma)",
        types = c('value = double(0)', 'x0=double(0)','t=double(0)','alpha=double(0)','K=double(0)','sigma=double(0)')
        )
    ))
}
# this is essential to ensure parallel works in scope  
assign('dmysde', dmysde, envir = .GlobalEnv)
assign('rmysde', rmysde, envir = .GlobalEnv)



code <- nimbleCode({
  
  sigma1 ~ dunif(0, 1000)
  sigma2 ~ dunif(0, 1000)
  sigma3 ~ dunif(0, 1000)
  
  a ~ dunif(0,1)
  K ~ dnorm(0,sd=1000)
  
   mu_t ~ dnorm(0,sd=1000)
  sigma_t ~ dunif(0, 1000)
 
  for(i in 1:n) {
     
      Observation1_latent[i]  ~ dnorm(mu_t,sd=sigma_t); # note = precision
      Observation1[i] ~ dnorm(Observation1_latent[i],sd=sigma1); # note = precision
 
      Observation2_latent[i] ~ dmysde(x0=Observation1_latent[i],t=Time2[i]-Time1[i],alpha=a,K=K,sigma=sigma3); # note sigma3=sd 
      
      Observation2[i]~dnorm(Observation2_latent[i],sd=sigma2); # note precision
     
       }
  
})

data1=dat[[1]];
constants <- list(n = nrow(data1))
data <- list(Observation1 = data1$Observation1,Observation2 = data1$Observation2,
             Time1 = data1$Time1,Time2=data1$Time2)
inits <- list(a=0.1,K=1000,sigma2=10,sigma3=0.1,sigma1=0.1,mu_t=0.1,sigma_t=1.0,
              Observation1_latent=rep(0.1,nrow(data1)),
              Observation2_latent=rep(0.1,nrow(data1))
              )
              

model <- nimbleModel(code, constants = constants, data = data1, inits = inits) # build model

myMCMCconf <- configureMCMC(model) 

print(myMCMCconf);
myMCMC<-buildMCMC(myMCMCconf)

Cmodel <- compileNimble(model,resetFunctions = TRUE,showCompilerOutput = FALSE)
CmodelMCMC <- compileNimble(myMCMC, project = Cmodel,resetFunctions = TRUE,showCompilerOutput = FALSE)

MCMCsamples <- runMCMC(CmodelMCMC, niter = 5000000,thin=500,nchains=1,setSeed=seeds,
                   samplesAsCodaMCMC = TRUE)

return(MCMCsamples);
}


chain_output <- parLapply(cl = this_cluster, X = 1:Nnodes, ## NOTE this is the seeds argument!
                          fun = run_MCMC_allcode, 
                          dat=list(data1) # this is one extra argument
                          )

# It's good practice to close the cluster when you're done with it.
stopCluster(this_cluster)



myres<-mcmc.list(chain_output[[1]],chain_output[[2]],chain_output[[3]],chain_output[[4]])

#load("mcmcOrigFitData.RData");# provides myres mcmc.list

save(myres,file="mcmcOrigFitData.RData");

# assume 2000 burn-in n.b. at 500 thin
plot(window(myres,start=2000),ask=TRUE)

gelman.diag(window(myres,start=2000))
# fine.

summary(window(myres,start=2000))




# \[Sigma] -> 0.09, \[Alpha] -> 0.11, K -> 1107 

####################
library(ggplot2);

mydat<-read.csv("data/sdedata4.csv",header=FALSE);# 4*231*2, 

mydat.m<-rbind(cbind(1,matrix(data=mydat[1:(231*2),1],ncol=2,byrow=TRUE)),
               cbind(2,matrix(data=mydat[((231*2)+1):(231*4),1],ncol=2,byrow=TRUE)),
               cbind(3,matrix(data=mydat[((231*4)+1):(231*6),1],ncol=2,byrow=TRUE)),
               cbind(4,matrix(data=mydat[((231*6)+1):(231*8),1],ncol=2,byrow=TRUE))
               )
mydat.df<-as.data.frame(mydat.m);
names(mydat.df)<-c("set","t","X")
mydat.df$set<-as.factor(mydat.df$set)

save(mydat.df,file="simsData.RData")

#mydat.l<-split(mydat.df,list(mydat.df[,1]))

#######---- Plot of four realizations - from gompertz_notes_v2.nb

p1<-ggplot(mydat.df,aes(t,X,color=set));
p1+geom_line(lwd=1.5)+
  labs(x="Time (days)",y=expression('Volume mm' ^3),
       color="Series")+theme(axis.line = element_line(size = 1, colour = "grey80"))+
  theme(axis.title.x = element_text(color = "grey20", size = 20),
        axis.title.y = element_text(color = "grey20", size = 20))
  
ggsave("sims.png",width=30,height=20,units="cm")

  #theme(axis.line = element_line(size = 1, colour = "grey80"))+

# +geom_line(color="darkblue",lwd=1)+geom_point(color="skyblue",size=3)+
# + +
  
  