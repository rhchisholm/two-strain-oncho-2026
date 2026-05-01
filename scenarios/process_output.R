# Process output from MDAi and MDAm simulations for different endemic prev

library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)


if (file.exists("processed_outputs.RData")==FALSE){
  
  ABRall<-c(2942,5350,27000) # annual biting rates
  MDAfall<-c(0.5,1) # time between MDA rounds (years)
  MDAnall<-c(40,20) # total number of MDA rounds
  
  source("../Fitting/equations.R")
  
  calculate_final_data <- function(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1){
    
    RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
    RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
    RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
    RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
    RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
    
    if(MDAf==0.5){
      time10<-which(yout[,1]==10)
      if(length(time10)==0){
        a<-yout[,1]
        closest_value <- a[which.min(abs(a - 10))]
        time10<-which(a==closest_value)
      }
      RecFinal_SOR_proportion_worms <- RecSOR_proportion_worms[time10]
      RecFinal_SOR_proportion_worms20 <- tail(RecSOR_proportion_worms,n=1)
      RecFinal_SOR_proportion_mf <- RecSOR_proportion_mf[time10]
      RecFinalPrevalence<-RecPrevalence[time10]
      RecFinalPrevalenceGR<-RecPrevalenceGR[time10]
      RecFinalPrevalenceSOR<-RecPrevalenceSOR[time10]
    }else{
      RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
      RecFinal_SOR_proportion_worms20 <- RecFinal_SOR_proportion_worms
      RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
      RecFinalPrevalence<-tail(RecPrevalence,n=1)
      RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
      RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
    }
    
    RecInitial_SOR_proportion_worms <- RecSOR_proportion_worms[11]
    
    if(MDAf==0.5){
      obs_times<-seq(from = 0, to = 40*MDAf, by = MDAf)
    }else{
      obs_times<-seq(from = 0, to = 20*MDAf, by = MDAf)
    }
    time<-yout[,1]
    obs_ind<-0*(1:length(obs_times))
    for(ii in 1:length(obs_times)){
      oi<-which(time==obs_times[ii])
      if(length(oi)==0){
        closest_value <- time[which.min(abs(time - obs_times[ii]))]
        oi<-which(time==closest_value)
      }
      obs_ind[ii]<-oi
    }
    obs<-RecPrevalence[obs_ind]
    pOTTIS<-obs_times[min(which(obs<0.014))]
    # if(is.na(pOTTIS)==TRUE){
    #   pOTTIS<-0
    # }
    
    newrow <- c(RecFinal_SOR_proportion_worms, RecFinal_SOR_proportion_mf,
                      RecFinalPrevalenceGR, RecFinalPrevalenceSOR, RecFinalPrevalence,
                      ABR, MDAf, SORpresent,
                      RecInitial_SOR_proportion_worms, MDAi, RecFinal_SOR_proportion_worms20, pOTTIS)
    
    return(newrow)
  }
  
  columns<-c('Final_SOR_proportion_worms','Final_SOR_proportion_mf',
             'Final_Prevalence_GR','Final_Prevalence_SOR','Final_Prevalence',
             'Endemic_prevalence','MDA_frequency', 'SOR_present',
             'Initial_SOR_proportion_worms','MDAi','Final_SOR_proportion_worms20yr', 'pOTTIS')
  empty_matrix <- matrix(ncol = length(columns), nrow = 0)
  df <- as.data.frame(empty_matrix)
  
  ########################################################
  # Ivermectin MDAs, original parameterisation, VR model #
  ########################################################
  
  MDAi<-1 # I: MDAi=1, M: MDAi=0
  
  # Original model scenarios, with MDAi
  filenames<-c(list("output/MDAi_Low_B.Rdata", "output/MDAi_Med_B.Rdata","output/MDAi_Hig_B.Rdata"), 
               list("output/MDAi_Low_A.Rdata", "output/MDAi_Med_A.Rdata", "output/MDAi_Hig_A.Rdata"))
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-1
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  df0<-df
  
  ##################################################
  # Ivermectin MDAs, SA parameterisation, VR model #
  ##################################################
  df <- as.data.frame(empty_matrix)
  filenames<-c(list("output/MDAi_Low_BSA.Rdata", "output/MDAi_Med_BSA.Rdata","output/MDAi_Hig_BSA.Rdata"), 
               list("output/MDAi_Low_ASA.Rdata", "output/MDAi_Med_ASA.Rdata", "output/MDAi_Hig_ASA.Rdata"))
  
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-1
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  dfSA0<-df
  
  ########################################################
  # Moxidectin MDAs, original parameterisation, VR model #
  ########################################################
  df <- as.data.frame(empty_matrix)
  MDAi<-0 # I: MDAi=1, M: MDAi=0
  
  # Original model scenarios, with MDAm
  filenames<-c(list("output/MDAm_Low_B.Rdata", "output/MDAm_Med_B.Rdata","output/MDAm_Hig_B.Rdata"), 
               list("output/MDAm_Low_A.Rdata", "output/MDAm_Med_A.Rdata", "output/MDAm_Hig_A.Rdata"))
  
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-1
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  df1<-df
  
  ##################################################
  # Moxidectin MDAs, SA parameterisation, VR model #
  ##################################################
  df <- as.data.frame(empty_matrix)
  # Sensitity analysis model scenarios, with MDAm
  filenames<-c(list("output/MDAm_Low_BSA.Rdata", "output/MDAm_Med_BSA.Rdata","output/MDAm_Hig_BSA.Rdata"), 
               list("output/MDAm_Low_ASA.Rdata", "output/MDAm_Med_ASA.Rdata", "output/MDAm_Hig_ASA.Rdata"))
  
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-1
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  dfSA1<-df
  
  ########################################################
  # Ivermectin MDAs, original parameterisation, GR model #
  ########################################################
  MDAi<-1 # I: MDAi=1, M: MDAi=0
  df <- as.data.frame(empty_matrix)
  
  # Original model with GR only, with MDAi
  filenames<-c(list("output/MDAi_Low_B_GR.Rdata", "output/MDAi_Med_B_GR.Rdata","output/MDAi_Hig_B_GR.Rdata"), 
               list("output/MDAi_Low_A_GR.Rdata", "output/MDAi_Med_A_GR.Rdata", "output/MDAi_Hig_A_GR.Rdata"))
  
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-0
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  df2<-df
  
  ##################################################
  # Ivermectin MDAs, SA parameterisation, GR model #
  ##################################################
  df <- as.data.frame(empty_matrix)
  # SA model with GR only, with MDAi
  filenames<-c(list("output/MDAi_Low_B_GRSA.Rdata", "output/MDAi_Med_B_GRSA.Rdata","output/MDAi_Hig_B_GRSA.Rdata"), 
               list("output/MDAi_Low_A_GRSA.Rdata", "output/MDAi_Med_A_GRSA.Rdata", "output/MDAi_Hig_A_GRSA.Rdata"))
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-0
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  dfSA2<-df
  ########################################################
  # Moxidectin MDAs, original parameterisation, GR model #
  ########################################################
  MDAi<-0 # I: MDAi=1, M: MDAi=0
  df <- as.data.frame(empty_matrix)
  
  # Original model with GR only, with MDAm
  filenames<-c(list("output/MDAm_Low_B_GR.Rdata", "output/MDAm_Med_B_GR.Rdata","output/MDAm_Hig_B_GR.Rdata"), 
               list("output/MDAm_Low_A_GR.Rdata", "output/MDAm_Med_A_GR.Rdata", "output/MDAm_Hig_A_GR.Rdata"))
  
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-0
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  df3<-df
  
  ##################################################
  # Moxidectin MDAs, SA parameterisation, GR model #
  ##################################################
  df <- as.data.frame(empty_matrix)
  # SA model with GR only, with MDAm
  filenames<-c(list("output/MDAm_Low_B_GRSA.Rdata", "output/MDAm_Med_B_GRSA.Rdata","output/MDAm_Hig_B_GRSA.Rdata"), 
               list("output/MDAm_Low_A_GRSA.Rdata", "output/MDAm_Med_A_GRSA.Rdata", "output/MDAm_Hig_A_GRSA.Rdata"))
  
  load(filenames[[1]])
  num_samples<-length(OutputM)
  SORpresent<-0
  
  for(ii in 1:3){
    ABR<-ABRall[ii]
    source("../Fitting/fixed-parameters.R")
    for(k in 1:2){
      ind<-(k-1)*3+ii
      MDAf <-MDAfall[k]
      NumberMDAs <- MDAnall[k]
      filename<-filenames[[ind]]
      load(filename)
      RecFinal <- matrix(0,nrow=num_samples,ncol=length(columns))
      for(n in 1:num_samples){
        yout<-OutputM[[n]]
        RecFinal[n,] <- calculate_final_data(yout, ABR, MDAf, SORpresent, MDAi, my_k0, my_k1) 
      }
      dft <-as.data.frame(RecFinal)
      df <- rbind(df, dft)
    }
  }
  
  names(df)<-columns
  df$MDA_frequency<-as.factor(df$MDA_frequency)
  df$Endemic_prevalence<-as.factor(df$Endemic_prevalence)
  df$SOR_present<-as.factor(df$SOR_present)
  df$MDAi<-as.factor(df$MDAi)
  dfSA3<-df
  
  # Combine all data 
  df <- rbind(df0, df1, df2, df3)
  dfSA <- rbind(dfSA0, dfSA1, dfSA2, dfSA3)
  df$Final_Prevalence<-df$Final_Prevalence*100
  dfSA$Final_Prevalence<-dfSA$Final_Prevalence*100
  
  save(df, dfSA, file = "processed_outputs.RData")
  
}

