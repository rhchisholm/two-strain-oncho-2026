library(deSolve)

# Functions for 2 strain 1 patch model

##-----------------------------
## the derivatives function (with cumulative effect)
##-----------------------------
derivs_MDA <- function(t,y,parameters,MDAon,t0,lg0,ls0,expMIN,expMAX,postMDA,fMDA,icw) {
  
  Historyg <- splinefun(t0,lg0,method = "natural")
  Historys <- splinefun(t0,ls0,method = "natural")
  
  m <-parameters["m"] 
  beta<-parameters["beta"] 
  Omega<-parameters["Omega"] 
  muH<-parameters["muH"]
  p<-parameters["p"]
  deltaH0<-parameters["deltaH0"]
  deltaHinfinity<-parameters["deltaHinfinity"]
  cH<-parameters["cH"] 
  lambda0<-parameters["lambda0"]
  lambdaMAX<-parameters["lambdaMAX"] 
  psiGR<-parameters["psiGR"] 
  nu<-parameters["nu"] 
  omega<-parameters["omega"] 
  psiSOR<-parameters["psiSOR"] 
  sigmaW<-parameters["sigmaW"] 
  deltaV0<-parameters["deltaV0"] 
  cV<-parameters["cV"] 
  CGR_GG<-parameters["CGR_GG"] 
  CGR_SS<-parameters["CGR_SS"] 
  CGR_SG<-parameters["CGR_SG"] 
  dom_c_SS<-parameters["dom_c_SS"] 
  dom_c_SG<-parameters["dom_c_SG"] 
  kW<-parameters["kW"] 
  zeta<-parameters["zeta"] 
  sigmaM0<-parameters["sigmaM0"]  
  alphaH<-parameters["alphaH"] 
  g<-parameters["g"] 
  sigmaL<-parameters["sigmaL"] 
  muL<-parameters["muL"] 
  alphaV<-parameters["alphaV"] 
  epsilon<-parameters["epsilon"] 
  
  cov<-0
  
  # Current coverage
  if (MDAon==1){
    cov<-parameters[paste("coverage",expMAX,sep="")]
  }else{
    if (postMDA==1){
      cov<-parameters[paste("coverage",expMAX,sep="")]
    }
  }
  
  r<-parameters["r"]
  omega0<-0.59
  omegabar<-omega0
  
  # Dominance parameters dependence on SOR worm proportion during MDA
  if(MDAon==1){
    propsor<-(y[6]+y[7])/(y[2]+y[3]+y[6]+y[7])
    CGR_SSf<-CGR_fn(propsor,icw,dom_c_SS,CGR_SS)
    CGR_SGf<-CGR_fn(propsor,icw,dom_c_SG,CGR_SG)
  }else{ #otherwise constant
    CGR_SSf<-CGR_SS
    CGR_SGf<-CGR_SG
  }
  
  
  # Set new dominance parameters post MDA based on new SOR proportion
  if(postMDA==1){
    CGR_SSf<-CGR_fn(0,icw,dom_c_SS,CGR_SS)
    CGR_SGf<-CGR_fn(0,icw,dom_c_SG,CGR_SG)
  }
  
  # Ensure probabilities are between 0.05, 0.5 and 1
  if(CGR_SSf<0.05){
    CGR_SSf<-0.05
  }
  if(CGR_SGf<0.5){
    CGR_SGf<-0.5
  }
  if(CGR_SSf>0.5){
    CGR_SSf<-0.5
  }
  if(CGR_SGf>1){
    CGR_SGf<-1
  }
  
  if (t < p){
    lag1 <- Historyg(t-p)
    lag2 <- Historys(t-p)
  }else{
    lag1 <- lagvalue(t - p,1)
    lag2 <- lagvalue(t - p,5)
  }
  
  if(expMAX>0){
    #dWvect = [dWG1 dWS1 dWG2 dWS2 ...]
    dWvec<-c()
  }
  
  # Specify worm equation for each exposure group and phenotype
  for(j in expMIN:expMAX){
    if(j==0){
      dWG <- 0.5 * m * beta * Omega * 
        (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
        (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag1 - sigmaW * y[9]
      dWS <- 0.5 * m * beta * Omega *  
        (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
        (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag2 - sigmaW * y[10]
    }else{
      dWvec[2*j-1] <- - sigmaW * y[9+2*j]
      dWvec[2*j] <- - sigmaW * y[10+2*j]
    }
  }
  
  # Total worm numbers for each phenotype
  SumWG <- 0
  SumWS <- 0
  for(j in expMIN:expMAX){
    SumWG <- SumWG + y[9+2*j]
    SumWS <- SumWS + y[10+2*j]
  }
  
  CumulativeEffectG <- 0
  CumulativeEffectS <- 0
  
  # Calculate average cumulative effects for each phenotype
  for(j in expMIN:expMAX){
    if(j==0){
      CE <- 1
    }else{
      CE <- (1-zeta)^(j-1)  
    }
    CumulativeEffectG <- CumulativeEffectG + CE * y[9+2*j] / SumWG
    CumulativeEffectS <- CumulativeEffectS + CE * y[10+2*j] / SumWS 
  }
  
  CumulativeEffectS <- 1 # Assume there is no cumulative effect for SOR worms
  
  dLG <- beta * Omega * deltaV0 / (1 + cV * (y[4] + y[8])) * y[4]  - 
    (alphaH/g + sigmaL + muL + alphaV * Omega * (y[4] + y[8])) * y[1]
  
  dNG <- 0.5 * m * beta * Omega *  
    (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
    (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag1 + 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiGR * (t+postMDA*fMDA))) * y[3] -
    (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiGR * (t+postMDA*fMDA)) + sigmaW) * y[2]
  
  dFG <-  (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiGR * (t+postMDA*fMDA))) * y[2] - 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiGR * (t+postMDA*fMDA)) + sigmaW) * y[3]
  
  totalworms<-y[2]+y[3]+y[6]+y[7]
  prob0GR<-(1+(y[2]+y[3])/kW)^(-kW-1)
  prob0SOR<-(1+(y[6]+y[7])/kW)^(-kW-1)
  prob0W<-(1+totalworms/kW)^(-kW-1)
  probanyW<-1-prob0W
  
  dMG <- epsilon * probanyW * CGR_GG * y[3] * CumulativeEffectG * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * CGR_SGf * y[7] * CumulativeEffectS * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * CGR_SGf * y[3] * CumulativeEffectG * (y[6]+y[7]) / totalworms +
    epsilon * probanyW * CGR_SSf * y[7] * CumulativeEffectS * (y[6]+y[7]) / totalworms -
    (sigmaM0 + (MDAon+postMDA) * cov * ((t+postMDA*fMDA) + nu)^(-omega)) * y[4]
  
  dLS <- beta * Omega * deltaV0 / (1 + cV * (y[4] + y[8])) * y[8] - 
    (alphaH/g + sigmaL + muL + alphaV * Omega * (y[4] + y[8])) * y[5]
  
  dNS <- 0.5 * m * beta * Omega *
    (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
    (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag2 + 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiSOR * (t+postMDA*fMDA))) * y[7] -
    (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiSOR * (t+postMDA*fMDA)) + sigmaW) * y[6]
  
  dFS <-  (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiSOR * (t+postMDA*fMDA))) * y[6] - 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiSOR * (t+postMDA*fMDA)) + sigmaW) * y[7]
  
  dMS <-  epsilon * probanyW * (1 - CGR_GG) * y[3] * CumulativeEffectG * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * (1 - CGR_SGf) * y[7] * CumulativeEffectS * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * (1 - CGR_SGf) * y[3] * CumulativeEffectG * (y[6]+y[7]) / totalworms+
    epsilon * probanyW * (1 - CGR_SSf) * y[7] * CumulativeEffectS * (y[6]+y[7]) / totalworms -
    (sigmaM0 + (MDAon+postMDA) * cov * ((t+postMDA*fMDA) + nu)^(-omega)) * y[8]
  
  
  if(expMAX==0){
    list(c(dLG, dNG, dFG, dMG, dLS, dNS, dFS, dMS, dWG, dWS))
  }else{
    list(c(dLG, dNG, dFG, dMG, dLS, dNS, dFS, dMS, dWG, dWS, dWvec))
  }
}

##-----------------------------
## the derivatives function (no cumulative effect)
##-----------------------------
derivs_MDA_noc <- function(t,y,parameters,MDAon,t0,lg0,ls0,expMIN,expMAX,postMDA,fMDA,icw) {
  
  Historyg <- splinefun(t0,lg0,method = "natural")
  Historys <- splinefun(t0,ls0,method = "natural")
  
  m <-parameters["m"] 
  beta<-parameters["beta"] 
  Omega<-parameters["Omega"] 
  muH<-parameters["muH"]
  p<-parameters["p"]
  deltaH0<-parameters["deltaH0"]
  deltaHinfinity<-parameters["deltaHinfinity"]
  cH<-parameters["cH"] 
  lambda0<-parameters["lambda0"]
  lambdaMAX<-parameters["lambdaMAX"] 
  psiGR<-parameters["psiGR"] 
  nu<-parameters["nu"] 
  omega<-parameters["omega"] 
  psiSOR<-parameters["psiSOR"] 
  sigmaW<-parameters["sigmaW"] 
  deltaV0<-parameters["deltaV0"] 
  cV<-parameters["cV"] 
  CGR_GG<-parameters["CGR_GG"] 
  CGR_SS<-parameters["CGR_SS"] 
  CGR_SG<-parameters["CGR_SG"] 
  dom_c_SS<-parameters["dom_c_SS"] 
  dom_c_SG<-parameters["dom_c_SG"] 
  kW<-parameters["kW"] 
  zeta<-parameters["zeta"] 
  sigmaM0<-parameters["sigmaM0"]  
  alphaH<-parameters["alphaH"] 
  g<-parameters["g"] 
  sigmaL<-parameters["sigmaL"] 
  muL<-parameters["muL"] 
  alphaV<-parameters["alphaV"] 
  epsilon<-parameters["epsilon"] 
  
  cov<-0
  
  # Current coverage
  if (MDAon==1){
    cov<-parameters[paste("coverage",expMAX,sep="")]
  }else{
    if (postMDA==1){
      cov<-parameters[paste("coverage",expMAX,sep="")]
    }
  }
  
  r<-parameters["r"]
  omega0<-0.59
  omegabar<-omega0
  
  # Dominance parameters dependence on SOR worm proportion during MDA
  if(MDAon==1){
    propsor<-(y[6]+y[7])/(y[2]+y[3]+y[6]+y[7])
    CGR_SSf<-CGR_fn(propsor,icw,dom_c_SS,CGR_SS)
    CGR_SGf<-CGR_fn(propsor,icw,dom_c_SG,CGR_SG)
  }else{ #otherwise constant
    CGR_SSf<-CGR_SS
    CGR_SGf<-CGR_SG
  }
  
  # Set new dominance parameters post MDA based on new SOR proportion
  if(postMDA==1){
    CGR_SSf<-CGR_fn(0,icw,dom_c_SS,CGR_SS)
    CGR_SGf<-CGR_fn(0,icw,dom_c_SG,CGR_SG)
  }
  
  # Ensure probabilities are between 0.05, 0.5 and 1
  if(CGR_SSf<0.05){
    CGR_SSf<-0.05
  }
  if(CGR_SGf<0.5){
    CGR_SGf<-0.5
  }
  if(CGR_SSf>0.5){
    CGR_SSf<-0.5
  }
  if(CGR_SGf>1){
    CGR_SGf<-1
  }
  
  if (t < p){
    lag1 <- Historyg(t-p)
    lag2 <- Historys(t-p)
  }else{
    lag1 <- lagvalue(t - p,1)
    lag2 <- lagvalue(t - p,5)
  }
  
  # Specify worm equations for each phenotype
  
  dWG <- 0.5 * m * beta * Omega * 
    (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
    (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag1 - sigmaW * y[9]
  dWS <- 0.5 * m * beta * Omega *  
    (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
    (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag2 - sigmaW * y[10]
  
  CumulativeEffectS <- 1 # Assume there is no cumulative effect for SOR worms
  CumulativeEffectG <- 1 # Assume there is no cumulative effect for GR worms
  
  dLG <- beta * Omega * deltaV0 / (1 + cV * (y[4] + y[8])) * y[4]  - 
    (alphaH/g + sigmaL + muL + alphaV * Omega * (y[4] + y[8])) * y[1]
  
  dNG <- 0.5 * m * beta * Omega *  
    (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
    (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag1 + 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiGR * (t+postMDA*fMDA))) * y[3] -
    (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiGR * (t+postMDA*fMDA)) + sigmaW) * y[2]
  
  dFG <-  (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiGR * (t+postMDA*fMDA))) * y[2] - 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiGR * (t+postMDA*fMDA)) + sigmaW) * y[3]
  
  totalworms<-y[2]+y[3]+y[6]+y[7]
  prob0GR<-(1+(y[2]+y[3])/kW)^(-kW-1)
  prob0SOR<-(1+(y[6]+y[7])/kW)^(-kW-1)
  prob0W<-(1+totalworms/kW)^(-kW-1)
  probanyW<-1-prob0W
  
  dMG <- epsilon * probanyW * CGR_GG * y[3] * CumulativeEffectG * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * CGR_SGf * y[7] * CumulativeEffectS * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * CGR_SGf * y[3] * CumulativeEffectG * (y[6]+y[7]) / totalworms +
    epsilon * probanyW * CGR_SSf * y[7] * CumulativeEffectS * (y[6]+y[7]) / totalworms -
    (sigmaM0 + (MDAon+postMDA) * cov * ((t+postMDA*fMDA) + nu)^(-omega)) * y[4]
  
  dLS <- beta * Omega * deltaV0 / (1 + cV * (y[4] + y[8])) * y[8] - 
    (alphaH/g + sigmaL + muL + alphaV * Omega * (y[4] + y[8])) * y[5]
  
  dNS <- 0.5 * m * beta * Omega *
    (deltaH0 + deltaHinfinity * cH * m * beta * (lag1 + lag2)) /
    (1 + cH * m * beta * (lag1 + lag2)) * exp(-muH * p) * lag2 + 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiSOR * (t+postMDA*fMDA))) * y[7] -
    (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiSOR * (t+postMDA*fMDA)) + sigmaW) * y[6]
  
  dFS <-  (omega0 - (MDAon+postMDA) * cov * omegabar * r * exp(-psiSOR * (t+postMDA*fMDA))) * y[6] - 
    (lambda0 + (MDAon+postMDA) * cov * r * lambdaMAX * exp(-psiSOR * (t+postMDA*fMDA)) + sigmaW) * y[7]
  
  dMS <-  epsilon * probanyW * (1 - CGR_GG) * y[3] * CumulativeEffectG * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * (1 - CGR_SGf) * y[7] * CumulativeEffectS * (y[2]+y[3]) / totalworms +
    epsilon * probanyW * (1 - CGR_SGf) * y[3] * CumulativeEffectG * (y[6]+y[7]) / totalworms+
    epsilon * probanyW * (1 - CGR_SSf) * y[7] * CumulativeEffectS * (y[6]+y[7]) / totalworms -
    (sigmaM0 + (MDAon+postMDA) * cov * ((t+postMDA*fMDA) + nu)^(-omega)) * y[8]
  
  list(c(dLG, dNG, dFG, dMG, dLS, dNS, dFS, dMS, dWG, dWS))
  
}

combineDFs<-function(df1,df2){
  df2[,1]<-df2[,1]+tail(df1[,1],n=1)
  
  if(length(df2[1,])>length(df1[1,])){
    numnewcols<-length(df2[1,])-length(df1[1,])
    df1new<-df1[,1:numnewcols]*0
    df1<-cbind(df1,df1new)
  }
  
  updated <- rbind(df1, df2[2:length(df2[,1]),])
  return(updated)
}

calculateLags<-function(df1,p){
  relevanttimes<-df1[,1]
  finaltime<-tail(relevanttimes,n=1)
  firsttime<-finaltime-p
  numbertimepoints<-length(relevanttimes[relevanttimes>=firsttime])
  lengthtimepoints<-length(relevanttimes)
  indices<-(lengthtimepoints-numbertimepoints+1):lengthtimepoints
  df1<-df1[indices,]
  t0<-df1[,1]
  lg0<-df1[,2]
  ls0<-df1[,6]
  # shift times to -p:0
  t0<-t0-finaltime
  lags<-list("timelag"=t0,"LGRlag"=lg0,"LSORlag"=ls0)
  return(lags)
}

prevalence<-function(M,k0,k1){
  kM<- k0*M/(1+k1*M)
  prev<-1-(1+M/(kM))^(-kM)
  return(prev)
}

# Function for dominance parameter dependence on SOR frequency during MDA programs
CGR_fn<-function(prop,ic,c1,c2){
  cgr<-c1*(ic-prop)+c2
  if (is.nan(cgr)){
    cgr<-0
  }
  return(cgr)
}