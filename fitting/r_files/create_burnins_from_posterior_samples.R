library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

# load equations (define functions)
source("equations.R")


#Intervention parameters
MDAi<-1 # I: MDAi=1, M: MDAi=0
MDAon<-0 # 0:no MDA, 1: MDA (starts off at zero for the burnin)
MDAf <-1 #Time between MDAs (years)
NumberMDAs <- 20
Dburnin<-300
DpostMDA<-10
# Coverage
my_coverage<-c(65,65,65,65,65,
               65,65,65,65,65,
               65,81.7,80.1,80,84.10,
               84.10,74.70,79.40,87.10,86.20)/100

# load fixed parameters (define functions)
source("fixed-parameters.R")

# From: "fitting_burnin_65000_UU.Rdata"
O1<-l[[1]]
O2<-l[[2]]

indv2<-which(l[[2]][,1]-l[[2]][,3]>0 & l[[2]][,3]-l[[2]][,2]>0 & l[[1]][,1]< 0.15 & l[[1]][,1] > 0.05)
num_samples<-length(indv2)

# run burn in for each sample and store

# Burn in details
yinit <- c(10, 10, 10, 10, 10, 10, 10, 10, 20, 20)
times <- seq(0, Dburnin, by = 0.1)
t0<-(0:100)*my_p/100-my_p
lg0<-0*t0+10
ls0<-1*lg0
# Not relevant for burn-in but set them anyway:
my_dom_c_SS<-1
my_dom_c_SG<-3
my_psiSOR<-my_psiGR*(10^0.5) 
my_r<-0

#setup parallel backend to use many processors
parallel::detectCores()
n.cores <- parallel::detectCores() - 1
#create the cluster
my.cluster <- parallel::makeCluster(
  n.cores, 
  type = "FORK" # Alternative is "PSOCK". FORK works for a Mac
)

#check cluster definition (optional)
print(my.cluster)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)
#check if it is registered (optional)
foreach::getDoParRegistered()
#how many workers are available? (optional)
foreach::getDoParWorkers()

burnins <-list()

# Run burnin for each sample, determine equilibrium proportion of worms
OutputB <- foreach(
  j=1:num_samples#, 
  #.combine=cbind
) %dopar% {
  
  samp<-indv2[j]
  my_CGR_GG<-l[[2]][samp,1]
  my_CGR_SS<-l[[2]][samp,2]
  my_CGR_SG<-l[[2]][samp,3]
  
  # Run burnin
  parameters <- c(m=my_m, beta=my_beta, Omega=my_Omega,
                  muH=my_muH, p=my_p, deltaH0=my_deltaH0,
                  deltaHinfinity=my_deltaHinfinity, cH=my_cH,
                  r=my_r, lambda0=my_lambda0, 
                  lambdaMAX=my_lambdaMAX, psiGR=my_psiGR,
                  nu=my_nu, omega=my_omega, psiSOR=my_psiSOR, 
                  sigmaW=my_sigmaW, deltaV0=my_deltaV0, cV=my_cV, 
                  CGR_GG=my_CGR_GG, CGR_SS=my_CGR_SS, CGR_SG=my_CGR_SG,
                  dom_c_SS=my_dom_c_SS, dom_c_SG=my_dom_c_SG,
                  kW=my_kW, zeta=my_zeta, sigmaM0=my_sigmaM0, 
                  alphaH=my_alphaH, g=my_g, sigmaL=my_sigmaL,
                  muL=my_muL, alphaV=my_alphaV, epsilon=my_epsilon, 
                  coverage=my_coverage)
  
  youtb <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, 
                MDAon=0, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=0, postMDA=0,
                fMDA=MDAf,icw=0)
  
  tempOutputB <-youtb
  
}

parallel::stopCluster(cl = my.cluster)

## Save fitting data
save(OutputB, file = filename) 










