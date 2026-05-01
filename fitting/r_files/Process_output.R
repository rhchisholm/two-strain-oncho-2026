# Process burnin plus annual MDAi program simulations for 
# 
# Cumulative effect on: accepted parameter samples
# Creates file: posteriors_ABR_65000_UU_cumu.Rdata
#
# Cumulative effect off: all parameter samples (as none are accepted)
# Creates file: posteriors_ABR_65000_UU_no_cumu.Rdata

# First: with cumulative effect

# # #################################################
# # # Process output from burnin simulations ########
# # #################################################

# determine burnins that satisfy criteria 1
preMDA <- vector("list", length(OutputB))
indv2<-which(l[[2]][,1]-l[[2]][,3]>0 & l[[2]][,3]-l[[2]][,2]>0 & l[[1]][,1]< 0.15 & l[[1]][,1] > 0.05)
num_samples<-length(indv2)

# extract burnin data for those satisfying Criteria 1
for(n in 1:num_samples){
  
  yout<-OutputB[[n]]
  
  samp<-indv2[n]
  my_CGR_GG<-l[[2]][samp,1]
  my_CGR_SS<-l[[2]][samp,2]
  my_CGR_SG<-l[[2]][samp,3]
  
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
  
  preMDA[[n]] = c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                  RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence,
                  my_CGR_GG,my_CGR_SS,my_CGR_SG)
  
}

# Save data from fitting pre MDA in data frame
columns2<- c('c_GR,GR','c0_SOR,SOR','c0_GR,SOR', 
             'Initial_mf_Prevalence', 'Initial_SOR_worm_proportion')
PreMDAdata<- data.frame(matrix(nrow=0, ncol=length(columns2))) 
names(PreMDAdata) <- columns2

for(jj in 1:num_samples){
  
  youtb<-OutputB[[jj]]
  endstep<-length(OutputB[[1]][,1])
  Final_SOR_W_ProportionBI <- 1 - (youtb[endstep,3]+youtb[endstep,4])/
    (youtb[endstep,3]+youtb[endstep,4]+youtb[endstep,7]+youtb[endstep,8])
  Final_prevBI <- prevalence(youtb[endstep,9]+youtb[endstep,5],my_k0,my_k1)
  
  samp<-indv2[jj]
  C_GG<-l[[2]][samp,1]
  C_SS<-l[[2]][samp,2]
  C_SG<-l[[2]][samp,3]
  
  newrow<-c(C_GG,C_SS,C_SG, 
            Final_prevBI, Final_SOR_W_ProportionBI)
  PreMDAdata[nrow(PreMDAdata) + 1, ]<-newrow
  
}

#######################################################
## Process output from annual MDAi simulations ########
#######################################################

# LHS simulation data for annual MDAi saved in ll
postMDA <- ll

# Calculate summary statistics for each LHS sample
dSOR <- vector("list", num_samples)
fSOR <- vector("list", num_samples)
fPrev <- vector("list", num_samples)
params <- vector("list", num_samples)

for(n in 1:num_samples){
  FP <- postMDA[[n]][,5]
  ISOR <- preMDA[[n]][1]
  fSOR[[n]] <- postMDA[[n]][,1]
  DiffSOR <- fSOR[[n]] - ISOR 
  dSOR[[n]] <- DiffSOR
  fPrev[[n]] <- FP
  dom_css <- postMDA[[n]][,6]
  dom_csg <- postMDA[[n]][,7]
  ratef <- postMDA[[n]][,8]
  r <- postMDA[[n]][,9]
  params[[n]] <- matrix(c(dom_css,dom_csg,ratef,r), ncol=4)
}

# Now find posterior samples that satisfy Criteria 2 (final prevalence and worm proportion 
# constraints)
FinalParams <- c()
FinalPrev <- c()
FinaldSOR <- c()
FinalSOR <- c()
BIsamplen <- c()

if(SA==0){
  for(n in 1:num_samples){
    for(i in 1:50){
      if(dSOR[[n]][i]>0.1){ #original
        if(fPrev[[n]][i]>0.1){
          if(fPrev[[n]][i]<0.12){
            p <- c(preMDA[[n]][6:8],params[[n]][i,])
            FinalParams <-rbind(FinalParams, p)
            FinalPrev <- rbind(FinalPrev,fPrev[[n]][i])
            FinaldSOR <- rbind(FinaldSOR,dSOR[[n]][i])
            FinalSOR <- rbind(FinalSOR,fSOR[[n]][i])
            BIsamplen <- rbind(BIsamplen, n)
          }
        } 
      }
    }
  }
}
if(SA==1){
  for(n in 1:num_samples){
    for(i in 1:50){
      if(dSOR[[n]][i]>0.05){ 
        if(dSOR[[n]][i]<0.1){ 
          if(fPrev[[n]][i]>0.1){
            if(fPrev[[n]][i]<0.12){
              p <- c(preMDA[[n]][6:8],params[[n]][i,])
              FinalParams <-rbind(FinalParams, p)
              FinalPrev <- rbind(FinalPrev,fPrev[[n]][i])
              FinaldSOR <- rbind(FinaldSOR,dSOR[[n]][i])
              FinalSOR <- rbind(FinalSOR,fSOR[[n]][i])
              BIsamplen <- rbind(BIsamplen, n)
            }
          }
        } 
      }
    }
  }
}

# Save annual MDA data for simulations satisfying criteria 1 and 2
posteriors<-matrix(c(BIsamplen,FinalParams,FinalPrev,FinaldSOR,FinalSOR),ncol=11)
save(posteriors, file = filename) 
