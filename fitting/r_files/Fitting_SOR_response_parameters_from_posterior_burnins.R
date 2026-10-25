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
my_zeta<-my_zetaX

indv2<-which(l[[2]][,1]-l[[2]][,3]>0 & l[[2]][,3]-l[[2]][,2]>0 & l[[1]][,1]< 0.15 & l[[1]][,1] > 0.05)
num_samples<-length(indv2)

# Specify LHS distributions
LHSsamples<- 100 # Number of LHS samples per burnin
LHSparameters<- 4 # number of LHS parameters 
uniformsamples <- randomLHS(LHSsamples,LHSparameters) # LHS samples from U[0,1]
minmaxvalues <- matrix(0, nrow = LHSparameters, ncol = 4) 
# Order of parameters:
# 1 dom_css = s_SOR,SOR
# 2 dom_csg = s_SOR,GR
# 3 ratef = multiplier on psiGR to get psiSOR
# 4 r = multiplier on lambda_max and omega_max to get new model transition rates

minmaxvalues[,1] <- c(0,0,1,0) # minimum values of each of the parameters
minmaxvalues[,2] <- c(20,20,10,1) # maximum values of each of the parameters
minmaxvalues[,3] <- minmaxvalues[,2] - minmaxvalues[,1] # width of distribution of each parameter

# transform LHS samples from U[0,1] to U[min,max]
for(k in 1:LHSparameters){
  uniformsamples[,k] <- uniformsamples[,k]*minmaxvalues[k,3]+minmaxvalues[k,1]
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

ll <- vector("list", num_samples)

for(n in 1:num_samples){
  
  youtb<-OutputB[[n]]
  SOR_W_ProportionBI <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
  Final_SOR_proportion_wormsBI<-tail(SOR_W_ProportionBI,n=1)
  tailyALLBI<-tail(youtb,n=1)
  tailyBI<-tailyALLBI[,2:9]
  samp<-indv2[n]
  my_CGR_GG<-l[[2]][samp,1]
  my_CGR_SS<-l[[2]][samp,2]
  my_CGR_SG<-l[[2]][samp,3]
  
  ## define initial conditions for MDA scenario from burn in
  cBI<-my_coverage[1]
  if(my_zeta!=0){
  yinitBI <- c(tailyBI,0*cBI+(1-cBI)*tailyALLBI[,10], 0*cBI+(1-cBI)*tailyALLBI[,11], cBI*tailyALLBI[,10:11])
  }else{
    yinitBI <-tailyALLBI[,2:11]
  }
  timesBI <- seq(0, MDAf, by = 0.01)
  
  ## Calculate lags
  t0BI<-(0:100)*my_p/100-my_p
  lg0BI<-0*t0BI+tailyBI[1]
  ls0BI<-0*t0BI+tailyBI[5]
  
  # Run burnin, determine equilibrium proportion of worms
  OutputM <- foreach(
    j=1:LHSsamples, 
    .combine=cbind
  ) %dopar% {
    
    # Specify parameter values from sample
    my_dom_c_SS<-uniformsamples[j,1] 
    my_dom_c_SG<-uniformsamples[j,2] 
    my_psiSOR<-my_psiGR*uniformsamples[j,3]
    my_r<-uniformsamples[j,4] 
    
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
    
    ## run MDA
    MDAon<-1
    if(my_zeta!=0){
      system.time(
        yout <- dede(y = yinitBI, times = timesBI, func = derivs_MDA, parms = parameters, 
                   MDAon=MDAon, t0=t0BI, lg0=lg0BI, ls0=ls0BI, expMIN=0, expMAX=1, 
                   postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI))
    }else{
      system.time(
        yout <- dede(y = yinitBI, times = timesBI, func = derivs_MDA_noc, parms = parameters, 
                     MDAon=MDAon, t0=t0BI, lg0=lg0BI, ls0=ls0BI, expMIN=0, expMAX=1, 
                     postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI))
    }
    
    ## combine time series from burn in to end of first MDA
    yout[,1]<-yout[,1]+tail(youtb[,1],n=1)
    if(length(yout[1,])>length(youtb[1,])){
      numnewcols<-length(yout[1,])-length(youtb[1,])
      youtbnew<-youtb[,1:numnewcols]*0
      youtb<-cbind(youtb,youtbnew)
    }
    updated <- rbind(youtb, yout[2:length(yout[,1]),])
    updated[,1]<-updated[,1]-Dburnin
    timeplot <-c(-0.1,tail(times,n=1))#+tail(youtb[,1],n=1)
    yout<-updated
    
    relevanttimes<-yout[,1]
    finaltime<-tail(relevanttimes,n=1)
    firsttime<-finaltime-my_p
    numbertimepoints<-length(relevanttimes[relevanttimes>=firsttime])
    lengthtimepoints<-length(relevanttimes)
    indices<-(lengthtimepoints-numbertimepoints+1):lengthtimepoints
    
    yout<-yout[indices,]
    
    ## define initial conditions for second MDA scenario
    tailyALL<-tail(yout,n=1)
    taily<-tailyALL[,2:9]
    lasttime<-tailyALL[,1]
    c<-my_coverage[2]
    if(my_zeta!=0){
      yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
                 c*tailyALL[,10:11] + (1-c)*tailyALL[,12:13],
                 c*tailyALL[,12:13] + (1-c)*0)
    }else{
      yinit <- tailyALL[,2:11]
    }
    times <- seq(0, MDAf, by = 0.01)
    
    ## Calculate lags
    # find time points for last two years
    t0<-yout[,1]
    lg0<-yout[,2]
    ls0<-yout[,6]
    
    ## run MDA
    MDAon<-1
    if(my_zeta!=0){
      system.time(
        yout2 <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, 
                      MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=2, 
                      postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI))
    }else{
      system.time(
        yout2 <- dede(y = yinit, times = times, func = derivs_MDA_noc, parms = parameters, 
                      MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=2, 
                      postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI))
    }
    
    ## combine time series to end of second MDA
    yout2[,1]<-yout2[,1]+tail(yout[,1],n=1)
    
    if(length(yout2[1,])>length(yout[1,])){
      numnewcols<-length(yout2[1,])-length(yout[1,])
      youtnew<-yout[,1:numnewcols]*0
      yout<-cbind(yout,youtnew)
    }
    
    updated <- rbind(yout, yout2[2:length(yout2[,1]),])
    timeplot <-updated[,1]
    yout<-updated
    
    for(i in 3:NumberMDAs){
      
      ## define initial conditions for current MDA
      tailyALL<-tail(yout,n=1)
      taily<-tailyALL[,2:9]
      lasttime<-tailyALL[,1]
      #yinit<-c(taily,0, 0, tailyALL[,c(10:(9+2*i))])
      
      c<-my_coverage[i]
      if(my_zeta!=0){
        yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
                   c*tailyALL[,10:(7+2*i)] + (1-c)*tailyALL[,12:(9+2*i)],
                   c*tailyALL[,(8+2*i):(9+2*i)] + (1-c)*0)
      }else{
        yinit <- tailyALL[,2:11]
      }
      times <- seq(0, MDAf, by = 0.01)
      
      ## Calculate lags
      # find time points for last two years
      lags <- calculateLags(yout,my_p)
      t0<-lags$timelag
      lg0<-lags$LGRlag
      ls0<-lags$LSORlag
      
      ## run MDA
      MDAon<-1
      if(my_zeta!=0){
        system.time(
          yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                        parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
        )}else{
          system.time(
            yout2 <- dede(y = yinit, times = times, func = derivs_MDA_noc, 
                          parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
          )
        }
      
      ##-----------------------------
      ## combine time series to end of this MDA
      ##-----------------------------
      yout<-combineDFs(yout,yout2)
    }
    
    RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
    RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
    
    RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
    RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
    
    RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
    RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
    RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
    
    RecFinalPrevalence<-tail(RecPrevalence,n=1)
    RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
    RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
    
    tempOutputM = c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                    RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence,
                    uniformsamples[j,1],uniformsamples[j,2],uniformsamples[j,3],uniformsamples[j,4])
    
    tempOutputM
    
  }
  ll[[n]]<-t(OutputM)
}

parallel::stopCluster(cl = my.cluster)

if(my_zeta==0.3){
  save(ll, file = filename)
} 
if(my_zeta==0){
  save(ll, file = filename)
}
