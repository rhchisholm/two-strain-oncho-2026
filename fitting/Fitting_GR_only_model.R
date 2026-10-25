# Calibrate GR-only from LHS samples of y to Asubende data

library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

ABR<-65000 # 80% prev

LHSsamples<- 2000 # Number of LHS samples
y <- randomLHS(LHSsamples,1) # LHS samples from U[0,1], 
# in equations.R, this parameter is labelled "r"
names(y)<-NULL
num_sims<-LHSsamples

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

############################################################
## Set up parallel backend
############################################################

n.cores <- parallel::detectCores() - 1

my.cluster <- parallel::makeCluster(
  n.cores,
  type = "FORK"       # appropriate for Mac
)

doParallel::registerDoParallel(my.cluster)

cat(
  "Running simulations using",
  foreach::getDoParWorkers(),
  "cores\n"
)

for(ce in 2:2){

  my_zeta <- Cumulative_effect[ce]
  
  results <- foreach(
    nn = 1:num_sims,
    .packages = c("deSolve", "epiR")
  ) %dopar% {
    
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
    
    RecPrevalence <- prevalence(
      yout[,9] + yout[,5],
      my_k0,
      my_k1
    )
    
    ## Return both outputs from this simulation
    list(
      PrevFinal = tail(RecPrevalence, n = 1),
      trajectory = yout
    )
  }
  
  ## Extract results returned by the parallel workers
  PrevFinal <- vapply(
    results,
    function(x) x$PrevFinal,
    numeric(1)
  )
  
  OutputA <- lapply(
    results,
    function(x) x$trajectory
  )
    
  l <- vector("list", 3)
  l[[1]]<-t(OutputA)
  l[[2]]<-PrevFinal
  l[[3]]<-y

  if(my_zeta!=0){
    save(l, file = "data_files/GR_only_model_LHS_simulations_annual_cumu_2000.Rdata")
  }else{
    save(l, file = "data_files/GR_only_model_LHS_simulations_annual_no_cumu_2000.Rdata")
  }
}

parallel::stopCluster(my.cluster)

###############################################################################
## GR-only model: LHS sampling adequacy
##
## Assesses whether 1000 LHS samples of y are sufficient to characterise
## the GR-only calibration by examining:
##
##   A. cumulative acceptance proportion; and
##   B. stability of the median and 2.5th-97.5th percentiles of retained y.
##
## Results are shown for:
##   - Primary calibration: final mf prevalence = 0.10-0.12
##   - Wider prevalence sensitivity: final mf prevalence = 0.08-0.14
###############################################################################

## Load GR-only LHS results WITH cumulative ivermectin effect
load(
  "data_files/GR_only_model_LHS_simulations_annual_cumu_2000.Rdata"
)

## Saved object:
##   l[[2]] = final mf prevalence after annual MDAi
##   l[[3]] = sampled values of y

final_prev <- as.numeric(l[[2]])
y <- as.numeric(l[[3]])

stopifnot(length(final_prev) == length(y))

Ntotal <- length(y)

cat("Total LHS samples:", Ntotal, "\n")


###############################################################################
## Define calibration criteria
###############################################################################

accepted_primary <-
  final_prev >= 0.10 &
  final_prev <= 0.12

accepted_wider <-
  final_prev >= 0.08 &
  final_prev <= 0.14

cat(
  "Primary accepted:",
  sum(accepted_primary),
  "of", Ntotal,
  "=", mean(accepted_primary), "\n"
)

cat(
  "Wider prevalence accepted:",
  sum(accepted_wider),
  "of", Ntotal,
  "=", mean(accepted_wider), "\n"
)


###############################################################################
## Evaluate at increasing cumulative LHS sample sizes
###############################################################################

Ns <- seq(400, 2000, by=400)
Ns <- Ns[Ns <= Ntotal]

###############################################################################
## Function to calculate summaries at each N
###############################################################################

summarise_LHS <- function(accepted, criterion_name) {
  
  do.call(
    rbind,
    lapply(Ns, function(N) {
      
      retained_y <- y[1:N][accepted[1:N]]
      
      data.frame(
        criterion = criterion_name,
        N = N,
        n_accepted = length(retained_y),
        acceptance_proportion = mean(accepted[1:N]),
        y_025 = quantile(retained_y, 0.025, na.rm = TRUE),
        y_median = median(retained_y, na.rm = TRUE),
        y_975 = quantile(retained_y, 0.975, na.rm = TRUE)
      )
    })
  )
}


primary_summary <- summarise_LHS(
  accepted_primary,
  "Primary (10-12%)"
)

wider_summary <- summarise_LHS(
  accepted_wider,
  "Wider prevalence (8-14%)"
)

GR_LHS_summary <- rbind(
  primary_summary,
  wider_summary
)

print(GR_LHS_summary)


###############################################################################
## Save numerical results
###############################################################################

data_dir <- "stability_data"
figure_dir <- "stability_figures"

dir.create(
  data_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

dir.create(
  figure_dir,
  showWarnings = FALSE,
  recursive = TRUE
)

write.csv(
  GR_LHS_summary,
  file.path(
    data_dir,
    "GR_only_LHS_stability.csv"
  ),
  row.names = FALSE
)


###############################################################################
## Three-panel figure
###############################################################################

## Text sizes
axis_cex   <- 1.5
label_cex  <- 1.6
title_cex  <- 1.4
legend_cex <- 1.2
panel_cex  <- 1.5

png(
  file.path(
    figure_dir,
    "GR_only_LHS_stability.png"
  ),
  width = 4200,
  height = 1600,
  res = 300
)

## Three equal-width panels
layout(
  matrix(1:3, nrow = 1),
  widths = c(1, 1, 1)
)

## Slightly larger margins for larger text
par(mar = c(5.5, 6.2, 4, 1.5))


###############################################################################
## Panel A: cumulative acceptance proportion
###############################################################################

yrange_A <- range(
  c(
    0,
    wider_summary$acceptance_proportion
  )
)

plot(
  primary_summary$N,
  primary_summary$acceptance_proportion,
  type = "b",
  pch = 16,
  lty = 1,
  xlim = range(Ns),
  ylim = yrange_A,
  xlab = "Number of LHS samples",
  ylab = "Proportion of parameter sets accepted",
  cex.axis = axis_cex,
  cex.lab = label_cex
)

lines(
  wider_summary$N,
  wider_summary$acceptance_proportion,
  type = "b",
  pch = 1,
  lty = 2
)

legend(
  "topright",
  inset = c(0.02, 0.2),
  legend = c(
    "Primary (10-12%)",
    "Wider prevalence (8-14%)"
  ),
  pch = c(16, 1),
  lty = c(1, 2),
  bty = "n",
  horiz = FALSE,
  cex = legend_cex
)

mtext(
  "A",
  side = 3,
  adj = 0,
  line = 1,
  font = 2,
  cex = panel_cex
)


###############################################################################
## Panel B: retained y - primary calibration
###############################################################################

yrange_B <- range(
  primary_summary$y_025,
  primary_summary$y_975,
  na.rm = TRUE
)

plot(
  primary_summary$N,
  primary_summary$y_median,
  type = "b",
  pch = 16,
  lty = 1,
  xlim = range(Ns),
  ylim = yrange_B,
  xlab = "Number of LHS samples",
  ylab = expression(y),
  cex.axis = axis_cex,
  cex.lab = label_cex
)

## 2.5th percentile
lines(
  primary_summary$N,
  primary_summary$y_025,
  type = "b",
  pch = 1,
  lty = 2
)

## 97.5th percentile
lines(
  primary_summary$N,
  primary_summary$y_975,
  type = "b",
  pch = 1,
  lty = 2
)

mtext(
  "Primary (10-12%)",
  side = 3,
  line = 0.5,
  cex = title_cex
)

mtext(
  "B",
  side = 3,
  adj = 0,
  line = 1,
  font = 2,
  cex = panel_cex
)


###############################################################################
## Panel C: retained y - wider prevalence sensitivity analysis
###############################################################################

yrange_C <- range(
  wider_summary$y_025,
  wider_summary$y_975,
  na.rm = TRUE
)

plot(
  wider_summary$N,
  wider_summary$y_median,
  type = "b",
  pch = 16,
  lty = 1,
  xlim = range(Ns),
  ylim = yrange_C,
  xlab = "Number of LHS samples",
  ylab = expression(y),
  cex.axis = axis_cex,
  cex.lab = label_cex
)

## 2.5th percentile
lines(
  wider_summary$N,
  wider_summary$y_025,
  type = "b",
  pch = 1,
  lty = 2
)

## 97.5th percentile
lines(
  wider_summary$N,
  wider_summary$y_975,
  type = "b",
  pch = 1,
  lty = 2
)

mtext(
  "Wider prevalence (8-14%)",
  side = 3,
  line = 0.5,
  cex = title_cex
)

mtext(
  "C",
  side = 3,
  adj = 0,
  line = 1,
  font = 2,
  cex = panel_cex
)

## Common key for panels B and C
legend(
  "topright",
  inset = c(0.02, 0.05),
  legend = c(
    "Median",
    "2.5th/97.5th percentile"
  ),
  pch = c(16, 1),
  lty = c(1, 2),
  bty = "n",
  cex = legend_cex
)


###############################################################################
## Finish
###############################################################################

dev.off()

layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))