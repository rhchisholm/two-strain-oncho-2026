###############################################################################
## Create compact plotting data for calibration figures
###############################################################################
## Run this script once in the full analysis environment, where the large raw
## simulation files are available. The resulting compact RData files can be
## committed to the repository; the large raw trajectory files are not needed
## thereafter for reproducing Figures 2 and S2-S4.
###############################################################################

if (!dir.exists("data_processed")) dir.create("data_processed", recursive = TRUE)

source("equations.R")
MDAi <- 1
ABR <- 65000
source("fixed-parameters.R")

###############################################################################
## A. Compact data for Figures S2-S4
###############################################################################

load("data_files/burnin_157_ABR_65000_UU.Rdata")
load("data_files/fitting_burnin_65000_UU.Rdata")

indv2 <- which(
  l[[2]][,1] - l[[2]][,3] > 0 &
  l[[2]][,3] - l[[2]][,2] > 0 &
  l[[1]][,1] > 0.05 & l[[1]][,1] < 0.15
)
num_samples <- length(indv2)

preMDA <- vector("list", num_samples)
for (n in seq_len(num_samples)) {
  yout <- OutputB[[n]]
  samp <- indv2[n]
  sor_worm <- 1 - (yout[,3] + yout[,4]) /
    (yout[,3] + yout[,4] + yout[,7] + yout[,8])
  preMDA[[n]] <- c(
    tail(sor_worm, 1),
    l[[2]][samp,1], l[[2]][samp,2], l[[2]][samp,3]
  )
}

stage2_dataframe <- function(postMDA, cumulative_effect) {
  z <- vector("list", num_samples * nrow(postMDA[[1]]))
  k <- 1
  for (n in seq_len(num_samples)) {
    initial_sor <- preMDA[[n]][1]
    for (i in seq_len(nrow(postMDA[[n]]))) {
      z[[k]] <- data.frame(
        SampleID = n,
        c_GRGR = preMDA[[n]][2],
        c0_SORSOR = preMDA[[n]][3],
        c0_GRSOR = preMDA[[n]][4],
        s_SORSOR = postMDA[[n]][i,6],
        s_SORGR = postMDA[[n]][i,7],
        f = postMDA[[n]][i,8],
        y = postMDA[[n]][i,9],
        Final_mf_prevalence = postMDA[[n]][i,5],
        Increase_in_SOR_worm_proportion = postMDA[[n]][i,1] - initial_sor,
        Cumulative_Effect = factor(cumulative_effect, levels = c(0,1)),
        Final_SOR_worm_proportion = postMDA[[n]][i,1]
      )
      k <- k + 1
    }
  }
  do.call(rbind, z)
}

load("data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100.Rdata")
dfc <- stage2_dataframe(ll, 1)
load("data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_no_cumu_100.Rdata")
dfnoc <- stage2_dataframe(ll, 0)
S2_VR_data <- subset(rbind(dfc, dfnoc), c0_SORSOR < 0.5 & c0_GRSOR > 0.5)

## GR-only S2 data. Update filenames here if the final archived files use a
## different LHS size/name.
load("data_files/GR_only_model_LHS_simulations_annual_cumu_2000.Rdata")
S2_GR_cumu <- data.frame(
  y = as.numeric(l[[3]]),
  Final_mf_prevalence = as.numeric(l[[2]]) * 100,
  Cumulative_effect = factor(1, levels = c(0,1))
)
load("data_files/GR_only_model_LHS_simulations_annual_no_cumu_2000.Rdata")
S2_GR_nocumu <- data.frame(
  y = as.numeric(l[[3]]),
  Final_mf_prevalence = as.numeric(l[[2]]) * 100,
  Cumulative_effect = factor(0, levels = c(0,1))
)
S2_GR_data <- rbind(S2_GR_cumu, S2_GR_nocumu)

## S3/S4: use the already processed primary retained parameter/outcome file.
load("data_files/primary_retained_parameter_outcomes_100.Rdata")
S34_primary_data <- data

save(
  S2_VR_data, S2_GR_data, S34_primary_data,
  file = "data_processed/primary_calibration_figure_data.Rdata",
  compress = "xz"
)

###############################################################################
## B. Compact data for Figure 2
###############################################################################

## VR trajectories and primary retained outcomes
load("data_files/Final_posterior_burnins_plusbiannualMDAs_cumu_100.Rdata")
OutputMM_VR <- OutputMM
load("data_files/primary_retained_parameter_outcomes_100.Rdata")
VR_primary <- data

## GR trajectories saved under wider 8-14% criterion
load("data_files/Final_posterior_burnins_GR_2000.Rdata")
OutputM_GR_all <- OutputM
load("data_files/Final_posterior_burnins_plusbiannualMDAs_GR_cumu_2000.Rdata")
OutputMM_GR_all <- OutputMM

GR_prev_annual_all <- vapply(
  OutputM_GR_all,
  function(yout) tail(prevalence(yout[,9] + yout[,5], my_k0, my_k1), 1),
  numeric(1)
)
keep_GR <- which(GR_prev_annual_all >= 0.10 & GR_prev_annual_all <= 0.12)
OutputM_GR <- OutputM_GR_all[keep_GR]
OutputMM_GR <- OutputMM_GR_all[keep_GR]
stopifnot(length(OutputM_GR) == length(OutputMM_GR))

## Panel A: save only time, trajectory identifier and total mf prevalence.
VR_long <- lapply(seq_along(OutputMM_VR), function(i) {
  yout <- OutputMM_VR[[i]]
  data.frame(
    time = yout[,1],
    trajectory = paste0("VR_", i),
    value = 100 * prevalence(yout[,9] + yout[,5], my_k0, my_k1)
  )
})
Figure2_VR_trajectories <- do.call(rbind, VR_long)

GR_long <- lapply(seq_along(OutputMM_GR), function(i) {
  yout <- OutputMM_GR[[i]]
  data.frame(
    time = yout[,1],
    trajectory = paste0("GR_", i),
    value = 100 * prevalence(yout[,9] + yout[,5], my_k0, my_k1)
  )
})
Figure2_GR_trajectories <- do.call(rbind, GR_long)

## Panel B distributions
GR_prev_annual <- vapply(
  OutputM_GR,
  function(yout) tail(prevalence(yout[,9] + yout[,5], my_k0, my_k1), 1),
  numeric(1)
)
GR_prev_biannual <- vapply(
  OutputMM_GR,
  function(yout) tail(prevalence(yout[,9] + yout[,5], my_k0, my_k1), 1),
  numeric(1)
)

Figure2_panelB_data <- rbind(
  data.frame(period="Annual calibration", model="Variable response",
             prevalence=VR_primary$Final_mf_prevalence_annual * 100),
  data.frame(period="Annual calibration", model="GR-only",
             prevalence=GR_prev_annual * 100),
  data.frame(period="Biannual projection", model="Variable response",
             prevalence=VR_primary$Final_mf_prevalence_biannual * 100),
  data.frame(period="Biannual projection", model="GR-only",
             prevalence=GR_prev_biannual * 100)
)
Figure2_panelB_data$period <- factor(
  Figure2_panelB_data$period,
  levels=c("Annual calibration", "Biannual projection")
)
Figure2_panelB_data$model <- factor(
  Figure2_panelB_data$model,
  levels=c("Variable response", "GR-only")
)

save(
  Figure2_VR_trajectories,
  Figure2_GR_trajectories,
  Figure2_panelB_data,
  file="data_processed/Figure2_plot_data.Rdata",
  compress="xz"
)

message("Compact repository plotting files created in data_processed/")
