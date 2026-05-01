# Fit variable response model

ABR<-65000 # 80% prev

# First stage: generate 1000 LHS parameter sets for c_ij, 
# calculate and store summary stats and LHS samples of burnins in the file: 
# "fitting_burnin_65000_UU.Rdata"

filename<-"data_files/fitting_burnin_65000_UU.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_phenotype_parameters.R")
}else{
  load(filename)
}

# Run burnins for LHS samples that satisfy initial SOR proportion criteria
# of between 5-15% SOR worms, and c_gg > c0sg > c0_ss, and store burnins in the 
# file: "burnin_157_ABR_65000_UU.Rdata"

filename<-"data_files/burnin_157_ABR_65000_UU.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/create_burnins_from_posterior_samples.R")
} else {
  load(filename)
}

# Second stage: generate 50 LHS samples of ivermectin response parameters, and
# run 20 year MDAi program from burnin for each of these 50 LHS samples, for 
# each of the accepted LHS samples from step 1. 

# with cumulative effect
my_zetaX<-0.3
filename<-"data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_SOR_response_parameters_from_posterior_burnins.R")
} else {
  load(filename)
}

# process output, determine samples that satisfy criteria 2: change in SOR worm 
# proportion x,  so that x > 0.1
SA<-0
filename<-"data_files/posteriors_ABR_65000_UU_cumu.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that 0.05 < x < 0.1 
SA<-1
filename<-"data_files/posteriors_ABR_65000_UU_cumuSA.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# without cumulative effect
my_zetaX<-0
filename<-"data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_no_cumu.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_SOR_response_parameters_from_posterior_burnins.R")
} else {
  load(filename)
}

# process output, determine samples that satisfy criteria 2: change in SOR worm 
# proportion x,  so that x > 0.1
SA<-0
filename<-"data_files/posteriors_ABR_65000_UU_no_cumu.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that 0.05 < x < 0.1 
SA<-1
filename<-"data_files/posteriors_ABR_65000_UU_no_cumuSA.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}


