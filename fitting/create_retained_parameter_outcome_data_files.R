###############################################################################
## Create retained parameter/outcome data files
## Primary calibration and four sensitivity analyses
###############################################################################
##
## PURPOSE
##
## This script creates processed data objects 
##
## These files contains a data frame ("data") giving the retained
## parameter combinations together with summary outcomes from:
##
##   - the 20-year annual MDAi calibration period; and
##   - the subsequent 18-year biannual MDAi projection.
##
## Files are generated for:
##
##   1. Primary calibration
##   2. Moderate-increase sensitivity analysis
##   3. Weak-increase sensitivity analysis
##   4. No offspring-phenotype selection sensitivity analysis
##   5. Wider-prevalence sensitivity analysis
##
###############################################################################

## ---------------------------------------------------------------------------
## Model functions and common data
## ---------------------------------------------------------------------------

source("equations.R")

MDAi <- 1
ABR <- 65000

source("fixed-parameters.R")

## Stage 1 burn-in trajectories
load("data_files/burnin_157_ABR_65000_UU.Rdata")


###############################################################################
## Define analyses
###############################################################################

analyses <- list(
  
  Primary = list(
    posterior_file =
      "data_files/posteriors_ABR_65000_UU_cumu_100.Rdata",
    
    trajectory_file =
      "data_files/Final_posterior_burnins_plusbiannualMDAs_cumu_100.Rdata",
    
    output_file =
      "data_files/primary_retained_parameter_outcomes_100.Rdata"
  ),
  
  
  Moderate = list(
    posterior_file =
      "data_files/posteriors_ABR_65000_UU_cumuSA_100.Rdata",
    
    trajectory_file =
      "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA_100.Rdata",
    
    output_file =
      "data_files/primary_retained_parameter_outcomesSA_100.Rdata"
  ),
  
  
  Weak = list(
    posterior_file =
      "data_files/posteriors_ABR_65000_UU_cumuSA2_100.Rdata",
    
    trajectory_file =
      "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA2_100.Rdata",
    
    output_file =
      "data_files/primary_retained_parameter_outcomesSA2_100.Rdata"
  ),
  
  
  No_selection = list(
    posterior_file =
      "data_files/posteriors_ABR_65000_UU_cumuSA5_100.Rdata",
    
    trajectory_file =
      "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA5_100.Rdata",
    
    output_file =
      "data_files/primary_retained_parameter_outcomesSA5_100.Rdata"
  ),
  
  
  Wider_prevalence = list(
    posterior_file =
      "data_files/posteriors_ABR_65000_UU_cumuSA3_100.Rdata",
    
    trajectory_file =
      "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA3_100.Rdata",
    
    output_file =
      "data_files/primary_retained_parameter_outcomesSA3_100.Rdata"
  )
)


###############################################################################
## Function to process one calibration analysis
###############################################################################

create_accepted_data <- function(
    posterior_file,
    trajectory_file,
    output_file,
    analysis_name) {
  
  cat("\n--------------------------------------------------\n")
  cat("Processing:", analysis_name, "\n")
  cat("--------------------------------------------------\n")
  
  
  ## -------------------------------------------------------------------------
  ## Load retained Stage 2 parameter combinations
  ## -------------------------------------------------------------------------
  
  load(posterior_file)
  ## loads object: posteriors
  
  ## Apply final inheritance-parameter constraints
  posteriors <- posteriors[
    posteriors[, 3] < 0.5 &
      posteriors[, 4] > 0.5,
    ,
    drop = FALSE
  ]
  
  num_sims <- nrow(posteriors)
  
  cat(
    "Retained parameter combinations:",
    num_sims,
    "\n"
  )
  
  
  ## -------------------------------------------------------------------------
  ## Load subsequent biannual MDAi trajectories
  ## -------------------------------------------------------------------------
  
  load(trajectory_file)
  ## loads object: OutputMM
  
  if (length(OutputMM) != num_sims) {
    
    stop(
      analysis_name,
      ": number of OutputMM trajectories (",
      length(OutputMM),
      ") does not match number of retained posterior rows (",
      num_sims,
      ")."
    )
  }
  
  ## -------------------------------------------------------------------------
  ## Extract outcomes at end of biannual projection
  ## -------------------------------------------------------------------------
  
  RecFinal <- matrix(
    NA_real_,
    nrow = num_sims,
    ncol = 6
  )
  
  for (n in seq_len(num_sims)) {
    
    yout <- OutputMM[[n]]
    
    ## SOR worm proportion
    RecSOR_proportion_worms <-
      1 -
      (yout[, 3] + yout[, 4]) /
      (yout[, 3] + yout[, 4] +
         yout[, 7] + yout[, 8])
    
    ## SOR mf proportion
    RecSOR_proportion_mf <-
      1 -
      yout[, 5] /
      (yout[, 5] + yout[, 9])
    
    ## Phenotype-specific and total mf prevalence
    RecPrevalenceGR <-
      prevalence(
        yout[, 5],
        my_k0,
        my_k1
      )
    
    RecPrevalenceSOR <-
      prevalence(
        yout[, 9],
        my_k0,
        my_k1
      )
    
    RecPrevalence <-
      prevalence(
        yout[, 9] + yout[, 5],
        my_k0,
        my_k1
      )
    
    
    ## Final outcomes
    RecFinal_SOR_proportion_worms <-
      tail(
        RecSOR_proportion_worms,
        1
      )
    
    RecFinal_dSOR_proportion_worms <-
      RecFinal_SOR_proportion_worms -
      RecSOR_proportion_worms[1]
    
    RecFinal_SOR_proportion_mf <-
      tail(
        RecSOR_proportion_mf,
        1
      )
    
    RecFinalPrevalence <-
      tail(
        RecPrevalence,
        1
      )
    
    RecFinalPrevalenceGR <-
      tail(
        RecPrevalenceGR,
        1
      )
    
    RecFinalPrevalenceSOR <-
      tail(
        RecPrevalenceSOR,
        1
      )
    
    
    RecFinal[n, ] <- c(
      RecFinal_SOR_proportion_worms,
      RecFinal_SOR_proportion_mf,
      RecFinalPrevalenceGR,
      RecFinalPrevalenceSOR,
      RecFinalPrevalence,
      RecFinal_dSOR_proportion_worms
    )
  }
  
  #############################################################################
  ## Construct data frame
  #############################################################################
  
  columns <- c(
    "SampleID",
    "c_GRGR",
    "c0_SORSOR",
    "c0_GRSOR",
    "s_SORSOR",
    "s_SORGR",
    "psi_SOR/psi_GR",
    "y",
    
    ## End of annual calibration
    "Final_mf_prevalence_annual",
    "Increase_in_SOR_worm_proportion",
    "Final_SOR_worm_proportion",
    
    ## End of subsequent biannual projection
    "Final_mf_prevalence_biannual",
    "Final_SOR_worm_proportion_biannual_MDAi",
    
    ## Pre-MDA endemic equilibrium
    "Initial_mf_prevalence",
    "Initial_SOR_worm_proportion",
    
    ## Change during biannual projection
    "Further_increase_in_SOR_worm_proportion"
  )
  
  ALLdata <- data.frame(
    matrix(
      NA_real_,
      nrow = num_sims,
      ncol = length(columns)
    )
  )
  
  names(ALLdata) <- columns
  
  #############################################################################
  ## Add pre-MDA and annual/biannual outcomes
  #############################################################################
  
  for (ii in seq_len(num_sims)) {
    
    ## Stage 1 inheritance sample associated with this retained combination
    sampleBI <- posteriors[ii, 1]
    
    youtb <- OutputB[[sampleBI]]
    
    endstep <- nrow(youtb)
    
    
    ## Pre-MDA SOR worm proportion
    Initial_SOR_W_Proportion <-
      1 -
      (youtb[endstep, 3] + youtb[endstep, 4]) /
      (
        youtb[endstep, 3] +
          youtb[endstep, 4] +
          youtb[endstep, 7] +
          youtb[endstep, 8]
      )
    
    
    ## Pre-MDA mf prevalence
    Initial_prev <-
      prevalence(
        youtb[endstep, 9] +
          youtb[endstep, 5],
        my_k0,
        my_k1
      )
    
    ## Increase in SOR during subsequent biannual projection
    Further_increase_SOR <-
      RecFinal[ii, 1] -
      posteriors[ii, 11]
    
    ## posteriors already contains:
    ##
    ##  1  SampleID
    ##  2  c_GRGR
    ##  3  c0_SORSOR
    ##  4  c0_GRSOR
    ##  5  s_SORSOR
    ##  6  s_SORGR
    ##  7  psi_SOR/psi_GR
    ##  8  y
    ##  9  final annual mf prevalence
    ## 10  increase in SOR during annual calibration
    ## 11  final annual SOR worm proportion
    
    ALLdata[ii, ] <- c(
      posteriors[ii, ],
      
      ## biannual outcomes
      RecFinal[ii, 5],
      RecFinal[ii, 1],
      
      ## pre-MDA outcomes
      Initial_prev,
      Initial_SOR_W_Proportion,
      
      ## change during biannual projection
      Further_increase_SOR
    )
  }
  
  #############################################################################
  ## Historical object name
  #############################################################################
  
  ## Retain the object name "data" for compatibility with existing scripts
  
  data <- ALLdata
  
  #############################################################################
  ## Save
  #############################################################################
  
  save(
    data,
    file = output_file
  )
  
  cat(
    "Saved:",
    output_file,
    "\n"
  )
  
  cat(
    "Rows saved:",
    nrow(data),
    "\n"
  )
  
  ## Return invisibly in case the function is run interactively
  invisible(data)
}


###############################################################################
## Run all five analyses
###############################################################################

for (analysis_name in names(analyses)) {
  
  x <- analyses[[analysis_name]]
  
  create_accepted_data(
    posterior_file = x$posterior_file,
    trajectory_file = x$trajectory_file,
    output_file = x$output_file,
    analysis_name = analysis_name
  )
}


###############################################################################
## Finished
###############################################################################

cat(
  "\nFinished creating retained parameter/outcome files for all analyses.\n"
)