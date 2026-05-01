# First re-run annual MDA from burnin for posterior samples and save

library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

# load equations (define functions)
source("equations.R")

ABR<-65000 # 80% prev

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

# load fixed parameters
source("fixed-parameters.R")

# Specify inheritance parameter values for GR only model
my_CGR_GG<-1
my_CGR_SS<-0
my_CGR_SG<-0
# not relevant in GR-only model, set to anything:
my_dom_c_SS<-0
my_dom_c_SG<-0
my_psiSOR<-my_psiGR*0


# Specify unknown ivermectin effect parameter
load("data_files/GR_only_model_LHS_simulations_annual_cumu.Rdata")
my_zeta<-0.3
PrevFinalc<-l[[2]]
y<-l[[3]]
posteriorsc<-matrix(c(y,PrevFinalc),ncol=2)
columns<-c('y','Final_mf_prevalence')
dfall <-as.data.frame(posteriorsc)
names(dfall) <- columns
data=subset(dfall, Final_mf_prevalence>0.1 & Final_mf_prevalence<0.12)
my_r<-median(data$y)

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

# Get variable values at final time of burn in
SOR_W_Proportion <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
Final_SOR_proportion_worms<-tail(SOR_W_Proportion,n=1)
Final_SOR_proportion_wormsBI<-Final_SOR_proportion_worms
tailyALLBI<-tail(youtb,n=1)
tailyBI<-tailyALLBI[,2:9]

## define initial conditions for MDA scenario from burn in
cBI<-my_coverage[1]
yinitBI <- c(tailyBI,0*cBI+(1-cBI)*tailyALLBI[,10], 0*cBI+(1-cBI)*tailyALLBI[,11], cBI*tailyALLBI[,10:11])
timesBI <- seq(0, MDAf, by = 0.01)

## Calculate lags
t0BI<-(0:100)*my_p/100-my_p
lg0BI<-0*t0BI+tailyBI[1]
ls0BI<-0*t0BI+tailyBI[5]

## run first MDA
MDAon<-1
system.time(
  yout <- dede(y = yinitBI, times = timesBI, func = derivs_MDA, parms = parameters, 
               MDAon=MDAon, t0=t0BI, lg0=lg0BI, ls0=ls0BI, expMIN=0, expMAX=1, 
               postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
)

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
yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
           c*tailyALL[,10:11] + (1-c)*tailyALL[,12:13],
           c*tailyALL[,12:13] + (1-c)*0)
times <- seq(0, MDAf, by = 0.01)

## Calculate lags
# find time points for last two years
t0<-yout[,1]
lg0<-yout[,2]
ls0<-yout[,6]

## run second MDA
MDAon<-1
system.time(
  yout2 <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, 
                MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=2, 
                postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
)

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

# run remaining annual MDAs
for(i in 3:NumberMDAs){
  
  ## define initial conditions for current MDA
  tailyALL<-tail(yout,n=1)
  taily<-tailyALL[,2:9]
  lasttime<-tailyALL[,1]
  #yinit<-c(taily,0, 0, tailyALL[,c(10:(9+2*i))])
  
  c<-my_coverage[i]
  yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
             c*tailyALL[,10:(7+2*i)] + (1-c)*tailyALL[,12:(9+2*i)],
             c*tailyALL[,(8+2*i):(9+2*i)] + (1-c)*0)
  times <- seq(0, MDAf, by = 0.01)
  
  ## Calculate lags
  # find time points for last two years
  lags <- calculateLags(yout,my_p)
  t0<-lags$timelag
  lg0<-lags$LGRlag
  ls0<-lags$LSORlag
  
  ## run MDA
  MDAon<-1
  system.time(
    yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                  parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
  )
  
  ##-----------------------------
  ## combine time series to end of this MDA
  ##-----------------------------
  yout<-combineDFs(yout,yout2)

}

yout_annual <- yout

################################################################
# Now run biannual MDA
################################################################

#Intervention parameters
MDAi<-1 # I: MDAi=1, M: MDAi=0
MDAon<-1 # 0:no MDA, 1: MDA (starts off at zero for the burnin)
MDAf <-0.5 #Time between MDAs (years)
NumberMDAsBiAnnual <- 34

# Coverage (yearly)
my_coverageALL<-c(92.00,92.00,87.30,88.00,0.00,74.20,
                  50.60,79.60,81.90,78.30,83.00,
                  79.50,82.10,78.70,73.10,77.50,
                  77.10,84.30)/100

# Coverage (annual followed by biannual MDAi)
my_coverage2<-rep(my_coverageALL[1:4],each=2)
my_coverage2<-c(my_coverage,my_coverage2)
my_coverage3<-rep(my_coverageALL[6:18],each=2)
my_coverage3<-c(my_coverage2,my_coverage3)

# Get variable values at final time of burn in
Final_SOR_proportion_worms_BI<-Final_SOR_proportion_wormsBI

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
                coverage=my_coverage2)

for(i in 21:28){
  ##-----------------------------
  ## define initial conditions for current MDA
  ##-----------------------------
  tailyALL<-tail(yout,n=1)
  lasttime<-tailyALL[,1]
  taily<-tailyALL[,2:9]
  
  c<-my_coverage2[i]
  
  yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
             c*tailyALL[,10:(7+2*i)] + (1-c)*tailyALL[,12:(9+2*i)],
             c*tailyALL[,(8+2*i):(9+2*i)] + (1-c)*0)
  
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
  system.time(
    yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                  parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
  )
  
  ##-----------------------------
  ## combine time series to end of this MDA
  ##-----------------------------
  yout<-combineDFs(yout,yout2)
}

# Simulate one year without MDA
## define new stable dominance parameters
Updated_SOR_W_Proportion <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
Updated_Final_SOR_proportion_worms<-tail(Updated_SOR_W_Proportion,n=1)
DifferenceSOR<-Final_SOR_proportion_worms_BI-Updated_Final_SOR_proportion_worms

##-----------------------------
## define initial conditions
##-----------------------------
tailyALL<-tail(yout,n=1)
numcomp<-length(tailyALL)
taily<-tailyALL[,2:numcomp]
lasttime<-tailyALL[,1]
yinit<-taily
times <- seq(0, 1, by = 0.01)

## Calculate lags
# find time points for last two years
lags <- calculateLags(yout,my_p)
t0<-lags$timelag
lg0<-lags$LGRlag
ls0<-lags$LSORlag

##-----------------------------
## run model no MDA
##-----------------------------
MDAon<-0
system.time(
  yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=1, fMDA=MDAf, icw=DifferenceSOR)
)

##-----------------------------
## combine time series to end of this year
##-----------------------------
yout<-combineDFs(yout,yout2)

# Simulate remainder biannual MDAi
# Coverage (biannually)

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
                coverage=my_coverage3)

for(i in 29:54){
  ##-----------------------------
  ## define initial conditions for current MDA
  ##-----------------------------
  tailyALL<-tail(yout,n=1)
  lasttime<-tailyALL[,1]
  taily<-tailyALL[,2:9]
  
  c<-my_coverage3[i]
  
  yinit <- c(taily,0*c+(1-c)*tailyALL[,10], 0*c+(1-c)*tailyALL[,11], 
             c*tailyALL[,10:(7+2*i)] + (1-c)*tailyALL[,12:(9+2*i)],
             c*tailyALL[,(8+2*i):(9+2*i)] + (1-c)*0)
  
  times <- seq(0, MDAf, by = 0.1)
  
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
  system.time(
    yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                  parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_worms_BI)
  )
  
  ##-----------------------------
  ## combine time series to end of this MDA
  ##-----------------------------
  yout<-combineDFs(yout,yout2)
}

yout_biannual <- yout

# Final values after burn in
tailybi<-tailyBI

# Final values after annual
tailyALL<-tail(yout_annual,n=1)
lasttime<-tailyALL[,1]
tailya<-tailyALL[,2:9]

# Final values after biannual
tailyALL<-tail(yout_biannual,n=1)
lasttime<-tailyALL[,1]
tailyba<-tailyALL[,2:9]

A<-c(tailybi,tailya,tailyba)
GR_only_model<-matrix(A,nrow=3,byrow = TRUE)

save(GR_only_model,file = "data_files/GR_only_model_final_stats.Rdata")
save(yout_biannual,file = "data_files/GR_only_model_full_trajectory.Rdata")
