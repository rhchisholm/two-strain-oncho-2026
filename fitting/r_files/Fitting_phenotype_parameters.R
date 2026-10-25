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

# Not relevant for burn-in but set them anyway:
my_dom_c_SS<-1
my_dom_c_SG<-3
my_psiSOR<-my_psiGR*(10^0.5) 
my_r<-0

# Burn in details
yinit <- c(10, 10, 10, 10, 10, 10, 10, 10, 20, 20)
times <- seq(0, Dburnin, by = 0.1)
t0<-(0:100)*my_p/100-my_p
lg0<-0*t0+10
ls0<-1*lg0

# Specify LHS distributions
if (is.null(seednumber)) {
  seednumber <- sample.int(.Machine$integer.max, 1)
}
set.seed(seednumber)
LHSsamples<- 1000 # Number of LHS samples
LHSparameters<- 3 # number of LHS parameters 
uniformsamples <- randomLHS(LHSsamples,LHSparameters) # LHS samples from U[0,1]
minmaxvalues <- matrix(0, nrow = LHSparameters, ncol = 3) 
# Order of parameters:
# 1 CGG
# 2 CSS
# 3 CSG
minmaxvalues[,1] <- c(0.9,0.01,0.01) # minimum values of each of the parameters
minmaxvalues[,2] <- c(1,1,1) # maximum values of each of the parameters
minmaxvalues[,3] <- minmaxvalues[,2] - minmaxvalues[,1] # width of distribution of each parameter

# transform LHS samples from U[0,1] to U[min,max]
for(j in 1:LHSparameters){
  uniformsamples[,j] <- uniformsamples[,j]*minmaxvalues[j,3]+minmaxvalues[j,1]
}

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

# Run burnin, determine equilibrium proportion of worms
Output <- foreach(
  i=1:LHSsamples #, 
  #.combine=cbind
  ) %dopar% {
  
  # Specify parameter values from sample
  my_CGR_GG<-uniformsamples[i,1] 
  my_CGR_SS<-uniformsamples[i,2] 
  my_CGR_SG<-uniformsamples[i,3]
  
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

  RecSOR_proportion_worms <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
  RecSOR_proportion_mf <- 1 - (youtb[,5])/(youtb[,5]+youtb[,9])
  
  RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
  RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
  
  RecPrevalenceGR<-prevalence(youtb[,5],my_k0,my_k1)
  RecPrevalenceSOR<-prevalence(youtb[,9],my_k0,my_k1)
  RecPrevalence<- prevalence(youtb[,9]+youtb[,5],my_k0,my_k1)#1-(1-RecPrevalenceSOR)*(1-RecPrevalenceGR)
  
  RecFinalPrevalence<-tail(RecPrevalence,n=1)
  RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
  RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
  
  tempOutputM = c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                  RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence,
                  uniformsamples[i,1],uniformsamples[i,2],uniformsamples[i,3])
  
  tempOutputM
  
  }

parallel::stopCluster(cl = my.cluster)

## Save fitting data
Output <- do.call(cbind, Output)
l <- vector("list", 2)
l[[1]]<-t(Output)
l[[2]]<-uniformsamples

save(l, file = filename) 



