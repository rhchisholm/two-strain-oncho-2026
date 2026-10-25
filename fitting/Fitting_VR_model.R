###############################################################################
## Variable-response model calibration
###############################################################################
##
## PURPOSE
##
## This is the main script used to calibrate the variable-response (VR)
## onchocerciasis transmission model.
##
## Calibration is performed in two stages:
##
##   Stage 1:
##     Sample the three phenotype-inheritance parameters and retain parameter
##     sets that generate a pre-MDA SOR worm proportion between 0.05 and 0.15
##     and satisfy the specified inheritance-parameter ordering.
##
##   Stage 2:
##     Cross the inheritance parameter sets retained from Stage 1 with
##     Latin hypercube samples of the four MDA-effect parameters and simulate
##     the 20-year annual MDAi calibration period in Asubende.
##
## The script runs the primary calibration and four sensitivity analyses,
## first under the model with a cumulative ivermectin effect and then under
## the corresponding model without a cumulative ivermectin effect.
##
## Existing output files are loaded rather than regenerated. Therefore,
## deleting an output file causes the corresponding simulation or processing
## step to be rerun the next time this script is executed.
##
## STAGE 1: INHERITANCE PARAMETERS
##
## 1. Generate 1000 LHS samples of:
##
##      c_GR,GR
##      c0_SOR,SOR
##      c0_GR,SOR
##
## 2. Simulate each parameter set to pre-MDA endemic equilibrium.
##
## 3. Retain inheritance parameter sets satisfying:
##
##      Belief 1:
##        pre-MDA SOR worm proportion = 0.05-0.15
##
##      inheritance ordering:
##        c_GR,GR > c0_GR,SOR > c0_SOR,SOR
##
## The current primary design retains 157 Stage 1 parameter sets.
##
## STAGE 2: MDA-EFFECT PARAMETERS
##
## For each inheritance parameter set retained from Stage 1, sample 100
## combinations of the four MDA-effect parameters:
##
##      s_SOR,SOR
##      s_SOR,GR
##      f
##      y
##
## Each combination is used to simulate the 20-year annual MDAi programme
## using the Asubende treatment-coverage history.
##
## The primary analysis assumes a cumulative ivermectin effect on GR worm
## fertility:
##
##      my_zetaX = 0.3
##
## PRIMARY AND SENSITIVITY CALIBRATIONS
##
## The variable SA passed to Process_output.R determines the Stage 2
## calibration criteria:
##
##   SA = 0 : Primary calibration
##            final mf prevalence = 0.10-0.12
##            increase in SOR worm proportion > 0.10
##
##   SA = 1 : Moderate-increase sensitivity analysis
##            final mf prevalence = 0.10-0.12
##            increase in SOR worm proportion = 0.05-0.10
##
##   SA = 2 : Weak-increase sensitivity analysis
##            final mf prevalence = 0.10-0.12
##            increase in SOR worm proportion < 0.05
##
##   SA = 3 : Wider-prevalence sensitivity analysis
##            final mf prevalence = 0.08-0.14
##            increase in SOR worm proportion > 0.10
##
##   SA = 5 : No offspring-phenotype selection sensitivity analysis
##            s_SOR,SOR = s_SOR,GR = 0
##            final mf prevalence = 0.10-0.12
##            no SOR-frequency change criterion (Belief 3) is imposed
##
## Final retained combinations additionally satisfy:
##
##      c0_SOR,SOR < 0.5
##      c0_GR,SOR  > 0.5
##
## CUMULATIVE-EFFECT COMPARISON
##
## The complete Stage 2 analysis is also run without the cumulative
## ivermectin effect:
##
##      my_zetaX = 0
##
## The same primary and alternative calibration criteria are applied to these
## simulations to assess whether parameter combinations satisfying the
## calibration criteria can be identified without the cumulative effect.
##
## These analyses are used as a model-structure comparison; failure to retain
## parameter combinations under the no-cumulative-effect model should be
## interpreted conditional on the model structure, sampled parameter ranges
## and calibration criteria.
##
## MAIN OUTPUT FILES
##
## Stage 1:
##
##   data_files/fitting_burnin_65000_UU.Rdata
##     Results from the 1000-point inheritance-parameter LHS.
##
##   data_files/burnin_157_ABR_65000_UU.Rdata
##     Burn-in trajectories for the inheritance parameter sets retained
##     following Stage 1.
##
## Stage 2 simulations with cumulative effect:
##
##   fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100.Rdata
##     Primary 100-point Stage 2 LHS simulation output.
##
##   fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100_s0.Rdata
##     Corresponding Stage 2 simulations with the offspring-phenotype
##     selection coefficients fixed at zero.
##
## Retained Stage 2 parameter combinations with cumulative effect:
##
##   posteriors_ABR_65000_UU_cumu_100.Rdata
##     Primary calibration.
##
##   posteriors_ABR_65000_UU_cumuSA_100.Rdata
##     Moderate-increase sensitivity analysis.
##
##   posteriors_ABR_65000_UU_cumuSA2_100.Rdata
##     Weak-increase sensitivity analysis.
##
##   posteriors_ABR_65000_UU_cumuSA3_100.Rdata
##     Wider-prevalence sensitivity analysis.
##
##   posteriors_ABR_65000_UU_cumuSA5_100.Rdata
##     No offspring-phenotype selection sensitivity analysis.
##
## No-cumulative-effect simulations and retained combinations use the
## corresponding filenames containing "no_cumu".
##
## TERMINOLOGY
##
## The objects and filenames retain the historical name "posteriors".
## These are described as retained parameter combinations/distributions
## rather than Bayesian posterior samples or posterior distributions.
## 
###############################################################################

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

# Second stage: generate 50 or 100 LHS samples of ivermectin response parameters, and
# run 20 year MDAi program from burnin for each of these LHS samples, for 
# each of the accepted LHS samples from step 1. 
# set.seed(20261013) # seed for 100
# with cumulative effect
my_zetaX<-0.3
filename<-"data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_SOR_response_parameters_from_posterior_burnins.R")
} else {
  load(filename)
}

# process output, determine samples that satisfy criteria 2: change in SOR worm 
# proportion x,  so that x > 0.1
SA<-0
filename<-"data_files/posteriors_ABR_65000_UU_cumu_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that 0.05 < x < 0.1 
SA<-1
filename<-"data_files/posteriors_ABR_65000_UU_cumuSA_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that x < 0.05 
SA<-2
filename<-"data_files/posteriors_ABR_65000_UU_cumuSA2_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 3: final prev
# is between 0.08 and 0.14
SA<-3
filename<-"data_files/posteriors_ABR_65000_UU_cumuSA3_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# SA when selection coefficients are set to 0, with cumulative effect
my_zetaX<-0.3
filename<-"data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100_s0.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_2SOR_response_parameters_from_posterior_burnins.R")
} else {
  load(filename)
}

# Then process output, determine samples with not constraint on SOR  
# worm proportion x
SA<-5
filename<-"data_files/posteriors_ABR_65000_UU_cumuSA5_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# without cumulative effect
# set.seed(20261013) # seed for 100
my_zetaX<-0
filename<-"data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_no_cumu_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_SOR_response_parameters_from_posterior_burnins.R")
} else {
  load(filename)
}

# process output, determine samples that satisfy criteria 2: change in SOR worm 
# proportion x,  so that x > 0.1
SA<-0
filename<-"data_files/posteriors_ABR_65000_UU_no_cumu_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that 0.05 < x < 0.1 
SA<-1
filename<-"data_files/posteriors_ABR_65000_UU_no_cumuSA_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that  x < 0.05 
SA<-2
filename<-"data_files/posteriors_ABR_65000_UU_no_cumuSA2_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# process output for sensitivity analysis with different criteria 2: change in 
# SOR worm proportion x, so that  x < 0.05 
SA<-3
filename<-"data_files/posteriors_ABR_65000_UU_no_cumuSA3_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

# SA when selection coefficients are set to 0, with no cumulative effect
filename<-"data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_no_cumu_100_s0.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_2SOR_response_parameters_from_posterior_burnins.R")
} else {
  load(filename)
}

# Then process output, determine samples with not constraint on SOR  
# worm proportion x
SA<-5
filename<-"data_files/posteriors_ABR_65000_UU_no_cumuSA5_100.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Process_output.R")
} else {
  load(filename)
}

