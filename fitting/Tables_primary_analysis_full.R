###############################################################################
## Supplementary tables: primary-calibration outcomes
###############################################################################
##
## PURPOSE
##
## This script reconstructs the retained parameter combinations and
## Asubende outcomes required for the supplementary summary tables from the
## primary 100-point Stage 2 calibration.
##
## ANALYSES IN THIS SCRIPT
##
## 1. Reconstruct the 157 inheritance parameter sets retained after Stage 1.
## 2. Reconstruct the primary Stage 2 annual-MDAi simulation quantities.
## 3. Load the primary retained parameter combinations and subsequent biannual
##    MDAi simulations for Asubende.
## 4. Extract the final retained primary parameter combinations after applying
##    the inheritance-parameter constraints c0_SORSOR < 0.5 and c0_GRSOR > 0.5.
## 5. Calculate 2.5th, 25th, 50th, 75th and 97.5th percentiles for the seven
##    retained model parameters and seven Asubende summary statistics.
##
## OUTPUT
##
## Summary-statistic objects are saved in:
##   tables/Tables_Asubende_100.Rdata
##
## The saved objects are:
##   ParameterStatsO  - quantiles of retained primary parameter values.
##   SummaryStatsO    - quantiles of primary Asubende calibration/projection
##                      summary statistics.
##
## INTERPRETATION
##
## These are distributions among parameter combinations retained under the
## specified calibration criteria and parameter constraints. They are not
## Bayesian posterior distributions, and the 2.5th-97.5th percentile intervals
## should not be interpreted as Bayesian credible intervals.
##
###############################################################################

if (!dir.exists("tables")) {
  dir.create("tables")
}

## -----------------------------------------------------------------------------
## 1. Load common model objects and Stage 1 retained inheritance sets
## -----------------------------------------------------------------------------

load("data_files/burnin_157_ABR_65000_UU.Rdata")
load("data_files/fitting_burnin_65000_UU.Rdata")

indv2 <- which(
  l[[2]][,1] - l[[2]][,3] > 0 &
  l[[2]][,3] - l[[2]][,2] > 0 &
  l[[1]][,1] > 0.05 &
  l[[1]][,1] < 0.15
)
num_samples <- length(indv2)

source("equations.R")
MDAi <- 1
ABR <- 65000
source("fixed-parameters.R")

## OutputB contains the Stage 1 burn-in simulations corresponding to the
## retained inheritance parameter sets and is used below to reconstruct
## pre-MDA prevalence and SOR worm proportion.

## -----------------------------------------------------------------------------
## 2. Reconstruct retained primary Asubende parameter/output data
## -----------------------------------------------------------------------------
columns<-c('SampleID','c_GRGR','c0_SORSOR','c0_GRSOR', 
           's_SORSOR', 's_SORGR','psi_SOR/psi_GR','y',
           'Final_mf_prevalence_annual','Increase_in_SOR_worm_proportion','Final_SOR_worm_proportion', # annual
           'Final_mf_prevalence_biannual','Final_SOR_worm_proportion_biannual_MDAi', # biannual
           'Initial_mf_prevalence', 'Initial_SOR_worm_proportion', # pre annual MDA
           'Further_increase_in_SOR_worm_proportion')  # end biannual - end annual

ALLdataAccepted<- data.frame(matrix(nrow=0, ncol=length(columns))) 
names(ALLdataAccepted) <- columns

## Load primary cumulative-effect calibration and subsequent biannual MDAi output
load("data_files/Final_posterior_burnins_plusbiannualMDAs_cumu_100.Rdata")
load("data_files/posteriors_ABR_65000_UU_cumu_100.Rdata")
posteriors <- posteriors[posteriors[,3]<0.5 & posteriors[,4]>0.5,]
num_sims<-length(posteriors[,1])
RecFinal <- matrix(0,nrow=num_sims,ncol=6)

for(n in 1:num_sims){
  
  yout<-OutputMM[[n]]
  
  RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
  RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
  RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
  RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
  RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
  
  RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
  RecFinal_dSOR_proportion_worms <- RecFinal_SOR_proportion_worms - RecSOR_proportion_worms[1] #whole program (a+b)

  RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
  RecFinalPrevalence<-tail(RecPrevalence,n=1)
  RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
  RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
  
  RecFinal[n,] <- c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                    RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence, 
                    RecFinal_dSOR_proportion_worms)
}

## Reconstruct pre-MDA quantities and combine annual-calibration and biannual-projection outcomes
for(ii in 1:num_sims){
  
  sampleBI<-posteriors[ii,1]
  youtb<-OutputB[[sampleBI]]
  endstep<-length(OutputB[[1]][,1])
  Final_SOR_W_ProportionBI <- 1 - (youtb[endstep,3]+youtb[endstep,4])/
    (youtb[endstep,3]+youtb[endstep,4]+youtb[endstep,7]+youtb[endstep,8])
  Final_prevBI <- prevalence(youtb[endstep,9]+youtb[endstep,5],my_k0,my_k1)
  IncreaseSORAB <- RecFinal[ii,1] - posteriors[ii,11]
  newrow<-c(posteriors[ii,],RecFinal[ii,5],RecFinal[ii,1], 
            Final_prevBI, Final_SOR_W_ProportionBI, IncreaseSORAB)
  ALLdataAccepted[nrow(ALLdataAccepted) + 1, ]<-newrow
  
}

ALLdataAccepted<-subset(ALLdataAccepted, c0_SORSOR<0.5 & c0_GRSOR>0.5)
ALLdataAcceptedO<-ALLdataAccepted

## -----------------------------------------------------------------------------
## 3. Create tables of summary statistics for primary Asubende simulations
## -----------------------------------------------------------------------------

param_names<-c('c_GRGR',
               'c0_SORSOR',
               'c0_GRSOR',
               's_SOR,SOR',
               's_SOR,GR',
               'psi_SOR/psi_GR',
               'y')

stats_names<-c('Final_mf_prevalence_annual',
               'Increase_in_SOR_worm_proportion_annual',
               'Final_SOR_worm_proportion_annual',
               'Final_mf_prevalence_biannual',
               'Further_increase_in_SOR_worm_proportion_biannual',
               'Final_SOR_worm_proportion_biannual',
               'Initial_SOR_worm_proportion')

ParameterStats<-matrix(ncol = 5, nrow = 7, 
                       dimnames = list(param_names,
                                       c("0.025", "0.25", "0.5", "0.75","0.975")))
SummaryStats<-matrix(ncol = 5, nrow = 7, 
                     dimnames = list(stats_names,
                                     c("0.025", "0.25", "0.5", "0.75","0.975")))

## Quantiles of retained primary parameter values
ParameterStatsO<-ParameterStats
for(i in 1:7){
  ParameterStatsO[i,]<-c(quantile(ALLdataAcceptedO[,i+1], 
                                 probs=c(0.025,0.25,0.5,0.75,0.975)))
}

# Create matrix of distributions of retained summary statistics
SummaryStatsO<-SummaryStats
inds<-c(9,10,11,12,16,13)
for(i in 1:7){
  if(i==7){
    dataSORi<-ALLdataAcceptedO[,11]-ALLdataAcceptedO[,10]
    SummaryStatsO[i,]<-c(quantile(dataSORi, 
                                  probs=c(0.025,0.25,0.5,0.75,0.975)))
  }else{
    SummaryStatsO[i,]<-c(quantile(ALLdataAcceptedO[,inds[i]], 
                                 probs=c(0.025,0.25,0.5,0.75,0.975)))
  }
}

save(ParameterStatsO,SummaryStatsO,file = "tables/Tables_Asubende_100.Rdata")


###############################################################################
## 4. GR-only model: primary calibration
###############################################################################
##
## The GR-only model has a single calibrated MDA-effect parameter, y.
## The primary calibration retains y values for which final mf prevalence
## following the 20-year annual MDAi calibration period is 10-12%.
###############################################################################

## Load annual-MDAi calibration results
load("data_files/GR_only_model_LHS_simulations_annual_cumu_2000.Rdata")

## Expected contents:
##   l[[2]] = final mf prevalence after annual MDAi
##   l[[3]] = sampled y values
GR_final_prev_annual <- as.numeric(l[[2]])
GR_y <- as.numeric(l[[3]])

stopifnot(length(GR_final_prev_annual) == length(GR_y))

## Identify wider and primary retained sets.
## The primary set is nested within the wider 8-14% set.
GR_wider_ind <- which(
  GR_final_prev_annual >= 0.08 &
    GR_final_prev_annual <= 0.14
)

GR_primary_ind <- which(
  GR_final_prev_annual >= 0.10 &
    GR_final_prev_annual <= 0.12
)

primary_within_wider <- which(
  GR_wider_ind %in% GR_primary_ind
)

GR_y_primary <- GR_y[GR_primary_ind]
GR_final_prev_annual_primary <-
  GR_final_prev_annual[GR_primary_ind]

cat(
  "\nGR-only primary calibration:\n",
  "Total LHS samples: ", length(GR_y), "\n",
  "Retained samples: ", length(GR_primary_ind), "\n",
  "Acceptance proportion: ",
  length(GR_primary_ind) / length(GR_y), "\n"
)

## Retained y distribution
GR_ParameterStatsO <- matrix(
  NA_real_,
  nrow = 1,
  ncol = 5,
  dimnames = list(
    "y",
    c("0.025", "0.25", "0.5", "0.75", "0.975")
  )
)

GR_ParameterStatsO["y", ] <- quantile(
  GR_y_primary,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

## Annual and biannual prevalence summaries
GR_SummaryStatsO <- matrix(
  NA_real_,
  nrow = 2,
  ncol = 5,
  dimnames = list(
    c(
      "Final_mf_prevalence_annual",
      "Final_mf_prevalence_biannual"
    ),
    c("0.025", "0.25", "0.5", "0.75", "0.975")
  )
)

GR_SummaryStatsO["Final_mf_prevalence_annual", ] <- quantile(
  GR_final_prev_annual_primary,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

## Load subsequent biannual-MDAi trajectories.
## This file contains the trajectories for the wider 8-14% retained set.
load(
  "data_files/Final_posterior_burnins_plusbiannualMDAs_GR_cumu_2000.Rdata"
)

## Select the nested primary 10-12% subset.
OutputMM_primary <- OutputMM[primary_within_wider]

stopifnot(
  length(OutputMM_primary) == length(GR_y_primary)
)

GR_final_prev_biannual <- numeric(length(OutputMM_primary))

for (n in seq_along(OutputMM_primary)) {
  yout <- OutputMM_primary[[n]]

  RecPrevalence <- prevalence(
    yout[,5] + yout[,9],
    my_k0,
    my_k1
  )

  GR_final_prev_biannual[n] <- tail(RecPrevalence, 1)
}

GR_SummaryStatsO["Final_mf_prevalence_biannual", ] <- quantile(
  GR_final_prev_biannual,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

save(
  GR_ParameterStatsO,
  GR_SummaryStatsO,
  GR_y_primary,
  GR_final_prev_annual_primary,
  GR_final_prev_biannual,
  file = "tables/Tables_Asubende_GR_primary.Rdata"
)

cat("\nGR-only retained parameter distribution:\n")
print(GR_ParameterStatsO)

cat("\nGR-only Asubende summary statistics:\n")
print(GR_SummaryStatsO)


###############################################################################
## 5. GR-only model: wider-prevalence sensitivity analysis
###############################################################################
##
## The wider-prevalence analysis retains y values for which final mf
## prevalence following the 20-year annual MDAi calibration period is 8-14%.
## All trajectories in the loaded GR-only biannual file correspond to this
## wider retained set.
###############################################################################

GR_y_wider <- GR_y[GR_wider_ind]

GR_final_prev_annual_wider <-
  GR_final_prev_annual[GR_wider_ind]

cat(
  "\nGR-only wider-prevalence calibration:\n",
  "Total LHS samples: ", length(GR_y), "\n",
  "Retained samples: ", length(GR_wider_ind), "\n",
  "Acceptance proportion: ",
  length(GR_wider_ind) / length(GR_y), "\n"
)

stopifnot(length(OutputMM) == length(GR_y_wider))

## Retained y distribution
GR_ParameterStats_wider <- matrix(
  NA_real_,
  nrow = 1,
  ncol = 5,
  dimnames = list(
    "y",
    c("0.025", "0.25", "0.5", "0.75", "0.975")
  )
)

GR_ParameterStats_wider["y", ] <- quantile(
  GR_y_wider,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

## Annual and biannual prevalence summaries
GR_SummaryStats_wider <- matrix(
  NA_real_,
  nrow = 2,
  ncol = 5,
  dimnames = list(
    c(
      "Final_mf_prevalence_annual",
      "Final_mf_prevalence_biannual"
    ),
    c("0.025", "0.25", "0.5", "0.75", "0.975")
  )
)

GR_SummaryStats_wider["Final_mf_prevalence_annual", ] <- quantile(
  GR_final_prev_annual_wider,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

GR_final_prev_biannual_wider <- numeric(length(OutputMM))

for (n in seq_along(OutputMM)) {
  yout <- OutputMM[[n]]

  RecPrevalence <- prevalence(
    yout[,5] + yout[,9],
    my_k0,
    my_k1
  )

  GR_final_prev_biannual_wider[n] <- tail(RecPrevalence, 1)
}

GR_SummaryStats_wider["Final_mf_prevalence_biannual", ] <- quantile(
  GR_final_prev_biannual_wider,
  probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  na.rm = TRUE
)

save(
  GR_ParameterStats_wider,
  GR_SummaryStats_wider,
  GR_y_wider,
  GR_final_prev_annual_wider,
  GR_final_prev_biannual_wider,
  file = "tables/Tables_Asubende_GR_wider_prevalence.Rdata"
)

cat("\nGR-only wider-prevalence retained parameter distribution:\n")
print(GR_ParameterStats_wider)

cat("\nGR-only wider-prevalence Asubende summary statistics:\n")
print(GR_SummaryStats_wider)


###############################################################################
## 6. Export compact data for repository table reproduction
###############################################################################
##
## This file contains only the summary matrices and compact retained values
## required to reproduce the supplementary calibration table objects.
## It does not contain raw model trajectories.
###############################################################################

dir.create(
  "data_processed",
  showWarnings = FALSE,
  recursive = TRUE
)

save(
  ## Variable-response primary calibration
  ParameterStatsO,
  SummaryStatsO,

  ## GR-only primary calibration
  GR_ParameterStatsO,
  GR_SummaryStatsO,
  GR_y_primary,
  GR_final_prev_annual_primary,
  GR_final_prev_biannual,

  ## GR-only wider-prevalence calibration
  GR_ParameterStats_wider,
  GR_SummaryStats_wider,
  GR_y_wider,
  GR_final_prev_annual_wider,
  GR_final_prev_biannual_wider,

  file = "data_processed/primary_tables_repo_data.Rdata"
)

message(
  "Saved compact repository table data: ",
  "data_processed/primary_tables_repo_data.Rdata"
)
