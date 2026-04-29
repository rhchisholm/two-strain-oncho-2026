# Calibrate GR-only from LHS samples of y to Asubende data

library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

ABR<-65000 # 80% prev

# Load LHS parameter samples
load("fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu.Rdata")
y<-ll[[1]][,9] # in equations.R, this parameter is labelled "r"
names(y)<-NULL
num_sims<-length(y)

#Intervention parameters
MDAi<-1 # I: MDAi=1, M: MDAi=0
MDAon<-0 # 0:no MDA, 1: MDA (starts off at zero for the burnin)
MDAf <-1 #Time between MDAs (years)
NumberMDAs <- 20
Dburnin<-300
DpostMDA<-10
# Coverage annual
my_coverage<-c(65,65,65,65,65,
               65,65,65,65,65,
               65,81.7,80.1,80,84.10,
               84.10,74.70,79.40,87.10,86.20)/100

# Load equations (define functions)
source("equations.R")

# Load fixed parameters
source("fixed-parameters.R")

# Specify inheritance parameter values for GR only model
my_CGR_GG<-1
my_CGR_SS<-0
my_CGR_SG<-0
# not relevant in GR-only model, set to anything:
my_dom_c_SS<-0
my_dom_c_SG<-0
my_psiSOR<-my_psiGR*0
my_r<-0 # not relevant during burn in, set to anything

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

yinit <- c(10, 10, 10, 10, 0, 0, 0, 0, 20, 0)
times <- seq(0, Dburnin, by = 0.1)

## Calculate lags
t0<-(0:100)*my_p/100-my_p
lg0<-0*t0+10
ls0<-0*lg0

##-----------------------------
## run burn in no MDA
##-----------------------------
system.time(
  youtb <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, MDAon=0, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=0, postMDA=0,fMDA=MDAf,icw=0)
)

SOR_W_Proportion <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
Final_SOR_proportion_worms<-tail(SOR_W_Proportion,n=1)
Final_SOR_proportion_worms_BI<-Final_SOR_proportion_worms

##----------------------------------------------------------
## For each value of the cumulative effect, run MDAi program
##----------------------------------------------------------
Cumulative_effect<-c(0.3,0)

for(ce in 1:2){
  PrevFinal<-c()
  OutputA<-vector("list",num_sims)
  my_zeta<-Cumulative_effect[ce]
  
  for(nn in 1:num_sims){
    
    parameters <- c(m=my_m, beta=my_beta, Omega=my_Omega,
                    muH=my_muH, p=my_p, deltaH0=my_deltaH0,
                    deltaHinfinity=my_deltaHinfinity, cH=my_cH,
                    r=y[nn], lambda0=my_lambda0, 
                    lambdaMAX=my_lambdaMAX, psiGR=my_psiGR,
                    nu=my_nu, omega=my_omega, psiSOR=my_psiSOR, 
                    sigmaW=my_sigmaW, deltaV0=my_deltaV0, cV=my_cV, 
                    CGR_GG=my_CGR_GG, CGR_SS=my_CGR_SS, CGR_SG=my_CGR_SG,
                    dom_c_SS=my_dom_c_SS, dom_c_SG=my_dom_c_SG,
                    kW=my_kW, zeta=my_zeta, sigmaM0=my_sigmaM0, 
                    alphaH=my_alphaH, g=my_g, sigmaL=my_sigmaL,
                    muL=my_muL, alphaV=my_alphaV, epsilon=my_epsilon, 
                    coverage=my_coverage)
    
    ##-----------------------------
    ## define initial conditions for MDA scenario from burn in
    ##-----------------------------
    tailyALL<-tail(youtb,n=1)
    taily<-tailyALL[,2:9]
    c<-my_coverage[1]
    if(my_zeta!=0){
      yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], c*tailyALL[,10:11])}else
      {
        yinit<-tailyALL[,2:11]
      }
    
    times <- seq(0, MDAf, by = 0.01)
    
    ## Calculate lags
    t0<-(0:100)*my_p/100-my_p
    lg0<-0*t0+taily[1]
    ls0<-0*t0+taily[5]
    
    ##-----------------------------
    ## run MDA
    ##-----------------------------
    MDAon<-1
    if(my_zeta!=0){
      system.time(
        yout <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, 
                     MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=1, postMDA=0, 
                     fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
      )}else{
        system.time(
          yout <- dede(y = yinit, times = times, func = derivs_MDA_noc, parms = parameters, 
                       MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=1, postMDA=0, 
                       fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
        )  
      }
    
    ##-----------------------------
    ## combine time series from burn in to end of first MDA
    ##-----------------------------
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
    
    ##-----------------------------
    ## define initial conditions for second MDA scenario
    ##-----------------------------
    tailyALL<-tail(yout,n=1)
    taily<-tailyALL[,2:9]
    lasttime<-tailyALL[,1]
    c<-my_coverage[2]

    if(my_zeta!=0){
      yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
                 c*tailyALL[,10:11] + (1-c)*tailyALL[,12:13],
                 c*tailyALL[,12:13] + (1-c)*0)
      times <- seq(0, MDAf, by = 0.01)}else{
        yinit<-tailyALL[,2:11]
      }
    
    ## Calculate lags
    # find time points for last two years
    t0<-yout[,1]
    lg0<-yout[,2]
    ls0<-yout[,6]
    
    ##-----------------------------
    ## run MDA
    ##-----------------------------
    MDAon<-1
    if(my_zeta!=0){
      system.time(
        yout2 <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, 
                      MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=2, postMDA=0, 
                      fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
      )}else{
        system.time(
          yout2 <- dede(y = yinit, times = times, func = derivs_MDA_noc, parms = parameters, 
                        MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=2, postMDA=0, 
                        fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
        )
      }
    
    ##-----------------------------
    ## combine time series to end of second MDA
    ##-----------------------------
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
      ##-----------------------------
      ## define initial conditions for current MDA
      ##-----------------------------
      tailyALL<-tail(yout,n=1)
      taily<-tailyALL[,2:9]
      lasttime<-tailyALL[,1]

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
      
      ##-----------------------------
      ## run MDA
      ##-----------------------------
      MDAon<-1
      if(my_zeta!=0){
        system.time(
          yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                        parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, 
                        expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
        )}else{
          system.time(
            yout2 <- dede(y = yinit, times = times, func = derivs_MDA_noc, 
                          parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, 
                          expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
          )
        }
      
      ##-----------------------------
      ## combine time series to end of this MDA
      ##-----------------------------
      yout<-combineDFs(yout,yout2)
    }
    
    RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
    PrevFinal[nn]<-tail(RecPrevalence,n=1)
    OutputA[[nn]] <-yout
    
  }
    
  l <- vector("list", 3)
  l[[1]]<-t(OutputA)
  l[[2]]<-PrevFinal
  l[[3]]<-y
  
  if(my_zeta!=0){
    save(l, file = "GR_only_model_LHS_simulations_annual_cumu.Rdata")
  }else{
    save(l, file = "GR_only_model_LHS_simulations_annual_no_cumu.Rdata")
  }
}

