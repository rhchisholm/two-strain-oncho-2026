
single_VR <- function(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin,my_zetaX){
  
  source("../fitting/fixed-parameters.R", local = TRUE)
  source("../fitting/equations.R", local = TRUE)
  my_zeta<-my_zetaX
  
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
  my_coverage<-rep(0.8, NumberMDAs)
  
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
  
  # Run burnin for each sample
  OutputB <- foreach(
    j=1:num_samples#, 
    #.combine=cbind
  ) %dopar% {
    
    my_CGR_GG<-PostParms[j,1]
    my_CGR_SS<-PostParms[j,2]
    my_CGR_SG<-PostParms[j,3]
    
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
    youtb
    
  }
  
  # Run MDA for each sample
  MDAon<-1
  
  OutputM <- foreach(
    j=1:num_samples#, 
    #.combine=cbind
  ) %dopar% {
    
    my_CGR_GG<-PostParms[j,1]
    my_CGR_SS<-PostParms[j,2]
    my_CGR_SG<-PostParms[j,3]
    my_dom_c_SS<-PostParms[j,4]
    my_dom_c_SG<-PostParms[j,5]
    my_psiSOR<-my_psiGR*PostParms[j,6]
    my_r<-PostParms[j,7]
    
    youtb<-OutputB[[j]]
    SOR_W_ProportionBI <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
    Final_SOR_proportion_wormsBI<-tail(SOR_W_ProportionBI,n=1)
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
    
    ## run MDA
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
    
    for(i in 3:NumberMDAs){
      
      ## define initial conditions for current MDA
      tailyALL<-tail(yout,n=1)
      taily<-tailyALL[,2:9]
      lasttime<-tailyALL[,1]
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
      system.time(
        yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                      parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
      )
      
      ##-----------------------------
      ## combine time series to end of this MDA
      ##-----------------------------
      yout<-combineDFs(yout,yout2)
    }
  
    yout
    
  }
  save(OutputM, file = filename)
  parallel::stopCluster(cl = my.cluster)
}


single_GR <- function(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin){
  
  source("../fitting/fixed-parameters.R", local = TRUE)
  source("../fitting/equations.R", local = TRUE)
  
  # run burn in for each sample and store
  # Burn in details
  yinit <- c(10, 10, 10, 10, 0, 0, 0, 0, 20, 0)
  times <- seq(0, Dburnin, by = 0.1)
  t0<-(0:100)*my_p/100-my_p
  lg0<-0*t0+10
  ls0<-1*lg0
  # Not relevant for burn-in but set them anyway:
  my_dom_c_SS<-1
  my_dom_c_SG<-3
  my_psiSOR<-my_psiGR*(10^0.5) 
  my_r<-0
  my_coverage<-rep(0.8, NumberMDAs)
  
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
  
  # Run burnin for each sample
  OutputB <- foreach(
    j=1:num_samples#, 
    #.combine=cbind
  ) %dopar% {
    
    my_CGR_GG<-1
    my_CGR_SS<-PostParms[j,2]
    my_CGR_SG<-PostParms[j,3]
    
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
    youtb
    
  }
  
  # Run MDA for each sample
  MDAon<-1
  
  OutputM <- foreach(
    j=1:num_samples#, 
    #.combine=cbind
  ) %dopar% {
    
    my_CGR_GG<-1
    my_CGR_SS<-PostParms[j,2]
    my_CGR_SG<-PostParms[j,3]
    my_dom_c_SS<-PostParms[j,4]
    my_dom_c_SG<-PostParms[j,5]
    my_psiSOR<-my_psiGR*PostParms[j,6]
    my_r<-PostParms[j,7]
    
    youtb<-OutputB[[j]]
    SOR_W_ProportionBI <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
    Final_SOR_proportion_wormsBI<-tail(SOR_W_ProportionBI,n=1)
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
    
    ## run MDA
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
    
    for(i in 3:NumberMDAs){
      
      ## define initial conditions for current MDA
      tailyALL<-tail(yout,n=1)
      taily<-tailyALL[,2:9]
      lasttime<-tailyALL[,1]
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
      system.time(
        yout2 <- dede(y = yinit, times = times, func = derivs_MDA, 
                      parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
      )
      
      ##-----------------------------
      ## combine time series to end of this MDA
      ##-----------------------------
      yout<-combineDFs(yout,yout2)
    }
    # trim trajectories to save:
    #yout<-yout[abs(yout[,1] %% 0.05) < 1e-9, ]
    yout
    
  }
  
  save(OutputM, file = filename)
  parallel::stopCluster(cl = my.cluster)
}

single_VR_noc <- function(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin){
  
  source("../fitting/fixed-parameters.R", local = TRUE)
  source("../fitting/equations.R", local = TRUE)
  
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
  my_coverage<-rep(0.8, NumberMDAs)
  
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
  
  # Run burnin for each sample
  OutputB <- foreach(
    j=1:num_samples#, 
    #.combine=cbind
  ) %dopar% {
    
    my_CGR_GG<-PostParms[j,1]
    my_CGR_SS<-PostParms[j,2]
    my_CGR_SG<-PostParms[j,3]
    
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
                    kW=my_kW, zeta=0, sigmaM0=my_sigmaM0, 
                    alphaH=my_alphaH, g=my_g, sigmaL=my_sigmaL,
                    muL=my_muL, alphaV=my_alphaV, epsilon=my_epsilon, 
                    coverage=my_coverage)
    
    youtb <- dede(y = yinit, times = times, func = derivs_MDA, parms = parameters, 
                  MDAon=0, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=0, postMDA=0,
                  fMDA=MDAf,icw=0)
    youtb
    
  }
  
  # Run MDA for each sample
  MDAon<-1
  
  OutputM <- foreach(
    j=1:num_samples#, 
    #.combine=cbind
  ) %dopar% {
    
    my_CGR_GG<-PostParms[j,1]
    my_CGR_SS<-PostParms[j,2]
    my_CGR_SG<-PostParms[j,3]
    my_dom_c_SS<-PostParms[j,4]
    my_dom_c_SG<-PostParms[j,5]
    my_psiSOR<-my_psiGR*PostParms[j,6]
    my_r<-PostParms[j,7]
    
    youtb<-OutputB[[j]]
    SOR_W_ProportionBI <- 1 - (youtb[,3]+youtb[,4])/(youtb[,3]+youtb[,4]+youtb[,7]+youtb[,8])
    Final_SOR_proportion_wormsBI<-tail(SOR_W_ProportionBI,n=1)
    tailyALLBI<-tail(youtb,n=1)
    tailyBI<-tailyALLBI[,2:9]
    
    ## define initial conditions for MDA scenario from burn in
    cBI<-my_coverage[1]
    yinitBI <-tailyALLBI[,2:11]
    timesBI <- seq(0, MDAf, by = 0.01)
    
    ## Calculate lags
    t0BI<-(0:100)*my_p/100-my_p
    lg0BI<-0*t0BI+tailyBI[1]
    ls0BI<-0*t0BI+tailyBI[5]
    
    parameters <- c(m=my_m, beta=my_beta, Omega=my_Omega,
                    muH=my_muH, p=my_p, deltaH0=my_deltaH0,
                    deltaHinfinity=my_deltaHinfinity, cH=my_cH,
                    r=my_r, lambda0=my_lambda0, 
                    lambdaMAX=my_lambdaMAX, psiGR=my_psiGR,
                    nu=my_nu, omega=my_omega, psiSOR=my_psiSOR, 
                    sigmaW=my_sigmaW, deltaV0=my_deltaV0, cV=my_cV, 
                    CGR_GG=my_CGR_GG, CGR_SS=my_CGR_SS, CGR_SG=my_CGR_SG,
                    dom_c_SS=my_dom_c_SS, dom_c_SG=my_dom_c_SG,
                    kW=my_kW, zeta=0, sigmaM0=my_sigmaM0, 
                    alphaH=my_alphaH, g=my_g, sigmaL=my_sigmaL,
                    muL=my_muL, alphaV=my_alphaV, epsilon=my_epsilon, 
                    coverage=my_coverage)
    
    ## run MDA
    system.time(
      yout <- dede(y = yinitBI, times = timesBI, func = derivs_MDA_noc, parms = parameters, 
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
    yinit <- tailyALL[,2:11]
    times <- seq(0, MDAf, by = 0.01)
    
    ## Calculate lags
    # find time points for last two years
    t0<-yout[,1]
    lg0<-yout[,2]
    ls0<-yout[,6]
    
    ## run MDA
    system.time(
      yout2 <- dede(y = yinit, times = times, func = derivs_MDA_noc, parms = parameters, 
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
    
    for(i in 3:NumberMDAs){
      
      ## define initial conditions for current MDA
      tailyALL<-tail(yout,n=1)
      taily<-tailyALL[,2:9]
      lasttime<-tailyALL[,1]
      yinit <- tailyALL[,2:11]
      times <- seq(0, MDAf, by = 0.01)
      
      ## Calculate lags
      # find time points for last two years
      lags <- calculateLags(yout,my_p)
      t0<-lags$timelag
      lg0<-lags$LGRlag
      ls0<-lags$LSORlag
      
      ## run MDA
      system.time(
        yout2 <- dede(y = yinit, times = times, func = derivs_MDA_noc, 
                      parms = parameters, MDAon=MDAon, t0=t0, lg0=lg0, ls0=ls0, expMIN=0, expMAX=i, postMDA=0, fMDA=MDAf, icw=Final_SOR_proportion_wormsBI)
      )
      
      ##-----------------------------
      ## combine time series to end of this MDA
      ##-----------------------------
      yout<-combineDFs(yout,yout2)
    }
    
    yout
    
  }
  save(OutputM, file = filename)
  parallel::stopCluster(cl = my.cluster)
}

