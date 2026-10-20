library(foreach)
library(doParallel)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

if (!dir.exists("output")) {
  dir.create("output")
}

Dburnin<-300 # duration of burnin
ABRall<-c(2942,5350,27000) # annual biting rates
MDAfall<-c(0.5,1) # time between MDA rounds (years)
MDAnall<-c(40,20) # total number of MDA rounds
my_zetaX<-0.3

# load equations (define functions) and fixed parameters
source("Single_scenario.R")

# SCENARIOS

# Original model scenarios, with MDAi
filenames<-c(list("output/MDAi_Low_B.Rdata", "output/MDAi_Med_B.Rdata","output/MDAi_Hig_B.Rdata"),
             list("output/MDAi_Low_A.Rdata", "output/MDAi_Med_A.Rdata", "output/MDAi_Hig_A.Rdata"))

# Load posterior parameter sets
load("../fitting/data_files/primary_retained_parameter_outcomes_100.Rdata")
PostParms<-data[,2:8]
num_samples<-length(PostParms[,1])
MDAi<-1 # I: MDAi=1, M: MDAi=0

# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin,my_zetaX)
#   }
# }
# 
# # Original model scenarios, with MDAm
# filenames<-c(list("output/MDAm_Low_B.Rdata", "output/MDAm_Med_B.Rdata","output/MDAm_Hig_B.Rdata"),
#              list("output/MDAm_Low_A.Rdata", "output/MDAm_Med_A.Rdata", "output/MDAm_Hig_A.Rdata"))
# 
# MDAi<-0 # I: MDAi=1, M: MDAi=0
# 
# for(ii in 1:3){
#   ABR<-ABRall[ii]
#   for(k in 1:2){
#     ind<-(k-1)*3+ii
#     MDAf <-MDAfall[k]
#     NumberMDAs <- MDAnall[k]
#     filename<-filenames[[ind]]
#     single_VR(PostParms,num_samples,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin,my_zetaX)
#   }
# }
# 
# # Original model with GR only, with MDAi
# filenames<-c(list("output/MDAi_Low_B_GR.Rdata", "output/MDAi_Med_B_GR.Rdata","output/MDAi_Hig_B_GR.Rdata"),
#              list("output/MDAi_Low_A_GR.Rdata", "output/MDAi_Med_A_GR.Rdata", "output/MDAi_Hig_A_GR.Rdata"))
# 
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

# # Original model with GR only, with MDAm
# filenames<-c(list("output/MDAm_Low_B_GR.Rdata", "output/MDAm_Med_B_GR.Rdata","output/MDAm_Hig_B_GR.Rdata"),
#              list("output/MDAm_Low_A_GR.Rdata", "output/MDAm_Med_A_GR.Rdata", "output/MDAm_Hig_A_GR.Rdata"))
# 
# MDAi<-0 # I: MDAi=1, M: MDAi=0
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

###############################################################################
## SENSITIVITY ANALYSIS SCENARIOS
##
## Run the 70% baseline-prevalence scenario under ivermectin MDA for:
##
##   - annual MDAi for 20 years; and
##   - biannual MDAi for 20 years (40 rounds).
##
## Variable-response model:
##   1. Moderate SOR increase
##   2. Weak SOR increase
##   3. No offspring-phenotype selection
##   4. Wider prevalence calibration
##
## GR-only model:
##   - Primary prevalence calibration (10-12%)
##   - Wider prevalence calibration (8-14%)
##
## The SOR-change and no-selection sensitivity analyses do not define
## alternative GR-only models because these calibration assumptions concern
## SOR phenotype dynamics, which are absent from the GR-only model.
###############################################################################


## ---------------------------------------------------------------------------
## Common scenario settings
## ---------------------------------------------------------------------------

ABR <- 27000          # approximately 70% baseline mf prevalence
MDAi <- 1             # ivermectin
Dburnin <- 300
my_zetaX <- 0.3

## Treatment frequencies
MDAf_sens <- c(
  Biannual = 0.5,
  Annual = 1
)

## Number of rounds over 20 years
MDAn_sens <- c(
  Biannual = 40,
  Annual = 20
)


###############################################################################
## VARIABLE-RESPONSE MODEL
###############################################################################

## Define the four sensitivity analyses and corresponding retained
## parameter/outcome files.

VR_sensitivity_files <- c(
  
  Moderate =
    "../fitting/data_files/primary_retained_parameter_outcomesSA_100.Rdata",
  
  Weak =
    "../fitting/data_files/primary_retained_parameter_outcomesSA2_100.Rdata",
  
  No_selection =
    "../fitting/data_files/primary_retained_parameter_outcomesSA5_100.Rdata",
  
  Wider_prevalence =
    "../fitting/data_files/primary_retained_parameter_outcomesSA3_100.Rdata"
)

## Short labels used in output filenames
VR_labels <- c(
  Moderate = "Moderate",
  Weak = "Weak",
  No_selection = "NoSelection",
  Wider_prevalence = "WiderPrev"
)


for (analysis_name in names(VR_sensitivity_files)) {
  
  cat(
    "\nRunning VR sensitivity analysis:",
    analysis_name,
    "\n"
  )
  
  ## Load object 'data'
  load(VR_sensitivity_files[[analysis_name]])
  
  ## Parameter columns required by single_VR()
  PostParms <- data[, 2:8]
  
  num_samples <- nrow(PostParms)
  
  cat(
    "Number of retained parameter combinations:",
    num_samples,
    "\n"
  )
  
  
  ## Run annual and biannual MDAi
  for (frequency in names(MDAf_sens)) {
    
    MDAf <- MDAf_sens[[frequency]]
    NumberMDAs <- MDAn_sens[[frequency]]
    
    freq_label <- ifelse(
      frequency == "Annual",
      "A",
      "B"
    )
    
    filename <- paste0(
      "output/MDAi_Hig_",
      freq_label,
      "_",
      VR_labels[[analysis_name]],
      ".Rdata"
    )
    
    cat(
      "  ",
      frequency,
      "MDAi ->",
      filename,
      "\n"
    )
    
    single_VR(
      PostParms,
      num_samples,
      filename,
      MDAf,
      NumberMDAs,
      ABR,
      MDAi,
      Dburnin,
      my_zetaX
    )
  }
}

###############################################################################
## GR-ONLY MODEL
###############################################################################
##
## The moderate-SOR, weak-SOR and no-selection sensitivity analyses concern
## SOR dynamics and therefore do not generate different GR-only parameter
## sets.
##
## We therefore run the GR-only model under:
##
##   1. Primary prevalence calibration: 10-12%
##   2. Wider prevalence calibration: 8-14%
##
###############################################################################

GR_calibration_files <- c(
  
  Primary =
    "../fitting/data_files/primary_retained_parameter_outcomes_100.Rdata",
  
  Wider_prevalence =
    "../fitting/data_files/primary_retained_parameter_outcomesSA3_100.Rdata"
)

GR_labels <- c(
  Primary = "Primary",
  Wider_prevalence = "WiderPrev"
)

for (analysis_name in names(GR_calibration_files)) {
  
  cat(
    "\nRunning GR-only calibration:",
    analysis_name,
    "\n"
  )
  
  ## Load object 'data'
  load(GR_calibration_files[[analysis_name]])
  
  PostParms_GR <- data[, 2:8]
  
  num_samples_GR <- nrow(PostParms_GR)
  
  cat(
    "Number of retained GR-only parameter combinations:",
    num_samples_GR,
    "\n"
  )
  
  for (frequency in names(MDAf_sens)) {
    
    MDAf <- MDAf_sens[[frequency]]
    NumberMDAs <- MDAn_sens[[frequency]]
    
    freq_label <- ifelse(
      frequency == "Annual",
      "A",
      "B"
    )
    
    filename <- paste0(
      "output/MDAi_Hig_",
      freq_label,
      "_GR_",
      GR_labels[[analysis_name]],
      ".Rdata"
    )
    
    cat(
      "  ",
      frequency,
      "MDAi ->",
      filename,
      "\n"
    )
    
    single_GR(
      PostParms_GR,
      num_samples_GR,
      filename,
      MDAf,
      NumberMDAs,
      ABR,
      MDAi,
      Dburnin
    )
  }
}

