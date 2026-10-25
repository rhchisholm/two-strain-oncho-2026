###############################################################################
## Stage 2 calibration: select retained parameter combinations
###############################################################################
##
## PURPOSE
##
## This script processes output from the Stage 1 burn-in simulations and
## Stage 2 annual-MDAi calibration simulations for the variable-response model.
##
## It reconstructs the pre-MDA quantities for the inheritance parameter sets
## retained following Stage 1, calculates the relevant outcomes following the
## 20-year annual MDAi program in Asubende, and identifies Stage 2 parameter
## combinations satisfying the specified calibration criteria.
##
## The script is sourced by Fitting_VR_model.R after the corresponding
## Stage 2 simulations have been generated.
##
## STAGE 1 CRITERIA
##
## Inheritance parameter sets are retained if:
##   - pre-MDA SOR worm proportion is between 0.05 and 0.15; and
##   - the inheritance parameters satisfy:
##         c_GR,GR > c0_GR,SOR > c0_SOR,SOR
##
## STAGE 2 CALIBRATION ANALYSES
##
## The calibration analysis is specified by the object SA:
##   SA = 0: Primary calibration
##           final mf prevalence = 0.10-0.12
##           increase in SOR worm proportion > 0.10
##   SA = 1: Moderate-increase sensitivity analysis
##           final mf prevalence = 0.10-0.12
##           increase in SOR worm proportion = 0.05-0.10
##   SA = 2: Weak-increase sensitivity analysis
##           final mf prevalence = 0.10-0.12
##           increase in SOR worm proportion < 0.05
##   SA = 3: Wider-prevalence sensitivity analysis
##           final mf prevalence = 0.08-0.14
##           increase in SOR worm proportion > 0.10
##   SA = 5: No offspring-phenotype selection sensitivity analysis
##           final mf prevalence = 0.10-0.12
##           no criterion is imposed on the change in SOR worm proportion
##
## For SA = 5, the offspring-phenotype selection parameters are fixed at zero
## in the corresponding Stage 2 simulations.
##
## INPUT OBJECTS
##
## The script assumes that the following objects have already been created or
## loaded by Fitting_VR_model.R:
##
##   OutputB
##       Stage 1 burn-in trajectories for retained inheritance parameter sets.
##
##   l
##       Stage 1 LHS results, including inheritance parameter values and
##       pre-MDA SOR worm proportions.
##
##   ll
##       Stage 2 annual-MDAi simulation summaries for the current calibration
##       analysis.
##
##   SA
##       Integer specifying the calibration analysis listed above.
##
##   filename
##       Filename to which the retained parameter combinations are saved.
##
## The prevalence() function and fixed model parameters must also be available
## in the current R environment.
##
## OUTPUT
##
## The script creates the object:
##
##   posteriors
##
## containing the parameter combinations satisfying the specified Stage 2
## calibration criteria, together with the corresponding final mf prevalence,
## change in SOR worm proportion, and final SOR worm proportion.
##
## The first column identifies the corresponding Stage 1 inheritance parameter
## set.
##
## The object is saved using:
##   save(posteriors, file = filename)
##
## TERMINOLOGY
##
## The historical object name "posteriors" is retained for compatibility with
## the existing analysis workflow. These are retained parameter combinations
## satisfying specified calibration criteria and are not Bayesian posterior
## samples.
##
###############################################################################

##################################################
## Process output from burnin simulations ########
##################################################

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
    for(i in 1:100){
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
    for(i in 1:100){
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
if(SA==2){
  for(n in 1:num_samples){
    for(i in 1:100){
      if(dSOR[[n]][i]<0.05){ 
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
if(SA==3){
  for(n in 1:num_samples){
    for(i in 1:100){
      if(dSOR[[n]][i]>0.1){ #original
        if(fPrev[[n]][i]>0.08){
          if(fPrev[[n]][i]<0.14){
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
if(SA==5){
  for(n in 1:num_samples){
    for(i in 1:100){
      #if(dSOR[[n]][i]>0.1){ #original
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
      #}
    }
  }
}

# Save annual MDA data for simulations satisfying criteria 1 and 2
posteriors<-matrix(c(BIsamplen,FinalParams,FinalPrev,FinaldSOR,FinalSOR),ncol=11)
save(posteriors, file = filename) 
