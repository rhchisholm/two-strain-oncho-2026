library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)


Dburnin<-300
ABRall<-c(2942,5350,27000)
MDAfall<-c(0.5,1)
MDAnall<-c(40,20)

# load equations (define functions) and fixed parameters
source("Single_scenario.R")

# SCENARIOS

# Original model scenarios, with MDAi
filenames<-c(list("MDAi_Low_B", "MDAi_Med_B","MDAi_Hig_B"), 
             list("MDAi_Low_A", "MDAi_Med_A", "MDAi_Hig_A"))
             
# Load posterior parameter sets
load("../fitting/accepted_trajectories.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-1 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}

# Sensitivity analysis model scenarios, with MDAi
filenames<-c(list("MDAi_Low_BSA", "MDAi_Med_BSA","MDAi_Hig_BSA"), 
             list("MDAi_Low_ASA", "MDAi_Med_ASA", "MDAi_Hig_ASA"))

# Load posterior parameter sets
load("../fitting/accepted_trajectoriesSA.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-1 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}


# Original model scenarios, with MDAm
filenames<-c(list("MDAm_Low_B", "MDAm_Med_B","MDAm_Hig_B"), 
             list("MDAm_Low_A", "MDAm_Med_A", "MDAm_Hig_A"))

# Load posterior parameter sets
load("../fitting/accepted_trajectories.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-0 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}

# Sensitity analysis model scenarios, with MDAm
filenames<-c(list("MDAm_Low_BSA", "MDAm_Med_BSA","MDAm_Hig_BSA"), 
             list("MDAm_Low_ASA", "MDAm_Med_ASA", "MDAm_Hig_ASA"))

# Load posterior parameter sets
load("../fitting/accepted_trajectoriesSA.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-0 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}

# Original model with GR only, with MDAi
filenames<-c(list("MDAi_Low_B_GR", "MDAi_Med_B_GR","MDAi_Hig_B_GR"), 
             list("MDAi_Low_A_GR", "MDAi_Med_A_GR", "MDAi_Hig_A_GR"))

# Load posterior parameter sets
load("../fitting/accepted_trajectories.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-1 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_GR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}

# SA model with GR only, with MDAi
filenames<-c(list("MDAi_Low_B_GRSA", "MDAi_Med_B_GRSA","MDAi_Hig_B_GRSA"), 
             list("MDAi_Low_A_GRSA", "MDAi_Med_A_GRSA", "MDAi_Hig_A_GRSA"))

# Load posterior parameter sets
load("../fitting/accepted_trajectoriesSA.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-1 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_GR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}

# Original model with GR only, with MDAm
filenames<-c(list("MDAm_Low_B_GR", "MDAm_Med_B_GR","MDAm_Hig_B_GR"), 
             list("MDAm_Low_A_GR", "MDAm_Med_A_GR", "MDAm_Hig_A_GR"))

# Load posterior parameter sets
load("../fitting/accepted_trajectories.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-0 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_GR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}

# SA model with GR only, with MDAm
filenames<-c(list("MDAm_Low_B_GRSA", "MDAm_Med_B_GRSA","MDAm_Hig_B_GRSA"), 
             list("MDAm_Low_A_GRSA", "MDAm_Med_A_GRSA", "MDAm_Hig_A_GRSA"))

# Load posterior parameter sets
load("../fitting/accepted_trajectoriesSA.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-0 # I: MDAi=1, M: MDAi=0

for(ii in 1:3){
  ABR<-ABRall[ii]
  for(k in 1:2){
    ind<-(k-1)*3+ii
    MDAf <-MDAfall[k]
    NumberMDAs <- MDAnall[k]
    filename<-filenames[[ind]]
    single_GR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
  }
}