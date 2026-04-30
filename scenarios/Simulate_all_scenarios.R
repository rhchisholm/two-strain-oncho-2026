library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

if (!dir.exists("output")) {
  dir.create("output")
}

Dburnin<-300
ABRall<-c(2942,5350,27000)
MDAfall<-c(0.5,1)
MDAnall<-c(40,20)

# load equations (define functions) and fixed parameters
source("Single_scenario.R")

# # SCENARIOS
# 
# # Original model scenarios, with MDAi
# filenames<-c(list("output/MDAi_Low_B.Rdata", "output/MDAi_Med_B.Rdata","output/MDAi_Hig_B.Rdata"), 
#              list("output/MDAi_Low_A.Rdata", "output/MDAi_Med_A.Rdata", "output/MDAi_Hig_A.Rdata"))
#              
# # Load posterior parameter sets
# load("../fitting/accepted_trajectories.Rdata")
# PostParms<-data[,2:8]
# num_samples<-length(PostParms[,1])
# MDAi<-1 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
#   }
# }
# 
# # Sensitivity analysis model scenarios, with MDAi
# filenames<-c(list("output/MDAi_Low_BSA.Rdata", "output/MDAi_Med_BSA.Rdata","output/MDAi_Hig_BSA.Rdata"), 
#              list("output/MDAi_Low_ASA.Rdata", "output/MDAi_Med_ASA.Rdata", "output/MDAi_Hig_ASA.Rdata"))
# 
# # Load posterior parameter sets
# load("../fitting/accepted_trajectoriesSA.Rdata")
# PostParms<-data[,2:8]
# num_samples<-length(PostParms[,1])
# MDAi<-1 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
#   }
# }
# 
# 
# # Original model scenarios, with MDAm
# filenames<-c(list("output/MDAm_Low_B.Rdata", "output/MDAm_Med_B.Rdata","output/MDAm_Hig_B.Rdata"), 
#              list("output/MDAm_Low_A.Rdata", "output/MDAm_Med_A.Rdata", "output/MDAm_Hig_A.Rdata"))
# 
# # Load posterior parameter sets
# load("../fitting/accepted_trajectories.Rdata")
# PostParms<-data[,2:8]
# num_samples<-length(PostParms[,1])
# MDAi<-0 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
#   }
# }
# 
# # Sensitity analysis model scenarios, with MDAm
# filenames<-c(list("output/MDAm_Low_BSA.Rdata", "output/MDAm_Med_BSA.Rdata","output/MDAm_Hig_BSA.Rdata"), 
#              list("output/MDAm_Low_ASA.Rdata", "output/MDAm_Med_ASA.Rdata", "output/MDAm_Hig_ASA.Rdata"))
# 
# # Load posterior parameter sets
# load("../fitting/accepted_trajectoriesSA.Rdata")
# PostParms<-data[,2:8]
# num_samples<-length(PostParms[,1])
# MDAi<-0 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
#   }
# }
# 
# # Original model with GR only, with MDAi
# filenames<-c(list("output/MDAi_Low_B_GR.Rdata", "output/MDAi_Med_B_GR.Rdata","output/MDAi_Hig_B_GR.Rdata"), 
#              list("output/MDAi_Low_A_GR.Rdata", "output/MDAi_Med_A_GR.Rdata", "output/MDAi_Hig_A_GR.Rdata"))
# 
# # Load posterior parameter sets
# load("../fitting/accepted_trajectories.Rdata")
# PostParms<-data[,2:8]
# num_samples<-length(PostParms[,1])
# MDAi<-1 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_GR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
#   }
# }
# 
# # SA model with GR only, with MDAi
# filenames<-c(list("output/MDAi_Low_B_GRSA.Rdata", "output/MDAi_Med_B_GRSA.Rdata","output/MDAi_Hig_B_GRSA.Rdata"), 
#              list("output/MDAi_Low_A_GRSA.Rdata", "output/MDAi_Med_A_GRSA.Rdata", "output/MDAi_Hig_A_GRSA.Rdata"))
# 
# # Load posterior parameter sets
# load("../fitting/accepted_trajectoriesSA.Rdata")
# PostParms<-data[,2:8]
# num_samples<-length(PostParms[,1])
# MDAi<-1 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_GR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
#   }
# }
# 


# Original model with GR only, with MDAm
filenames<-c(list("output/MDAm_Low_B_GR.Rdata", "output/MDAm_Med_B_GR.Rdata","output/MDAm_Hig_B_GR.Rdata"), 
             list("output/MDAm_Low_A_GR.Rdata", "output/MDAm_Med_A_GR.Rdata", "output/MDAm_Hig_A_GR.Rdata"))

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
filenames<-c(list("output/MDAm_Low_B_GRSA.Rdata", "output/MDAm_Med_B_GRSA.Rdata","output/MDAm_Hig_B_GRSA.Rdata"), 
             list("output/MDAm_Low_A_GRSA.Rdata", "output/MDAm_Med_A_GRSA.Rdata", "output/MDAm_Hig_A_GRSA.Rdata"))

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