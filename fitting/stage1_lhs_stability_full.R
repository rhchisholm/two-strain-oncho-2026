###############################################################################
## Stage 1 calibration: LHS sampling stability and independent-design checks
###############################################################################
##
## PURPOSE
##
## This script assesses the adequacy and reproducibility of the Latin
## hypercube sampling (LHS) used for Stage 1 of the variable-response model
## calibration. Stage 1 samples the three inheritance parameters and retains
## parameter sets satisfying Belief 1 and the specified inheritance ordering.
##
## Stage 1 acceptance requires:
##
##   Belief 1:
##     pre-MDA SOR worm proportion = 0.05-0.15
##
##   Inheritance ordering:
##     c_GR,GR > c0_GR,SOR > c0_SOR,SOR
##
## Internal parameter names used in the first part of this script are:
##     p1 = c_GR,GR
##     p2 = c0_SOR,SOR
##     p3 = c0_GR,SOR
##
## ANALYSES
##
## 1. Primary Stage 1 LHS sampling adequacy
##    - Uses the primary 1000-point inheritance-parameter LHS design.
##    - Examines the cumulative acceptance proportion at increasing sample
##      sizes (N = 200, 400, 600, 800 and 1000).
##    - Examines stability of the median and 2.5th-97.5th percentiles of
##      each inheritance parameter among accepted samples.
##    - Compares the joint accepted parameter space identified among the
##      first 800 LHS samples with that identified among the final 200
##      samples, to assess whether additional sampling reveals distinct
##      plausible regions.
##
## 2. Independent LHS design / random-seed analysis
##    - Compares the primary 1000-point LHS design with two independently
##      generated 250-point LHS designs using different random seeds.
##    - Applies the same Stage 1 acceptance criteria to all three designs.
##    - Compares:
##        * acceptance proportions;
##        * medians and 2.5th-97.5th percentiles of accepted parameter values;
##        * pairwise distributions of accepted inheritance parameters.
##
##    This analysis assesses whether the Stage 1 results are sensitive to the
##    particular LHS design/random seed used for the primary calibration.
##
## OUTPUT
##
## Figures are saved in:
##   stability_figures/
##
## Main figure outputs:
##
##   Stage1_acceptance_stability.png
##     Cumulative Stage 1 acceptance proportion at selected LHS sample sizes.
##
##   Stage1_parameter_convergence.png
##     Stability of the median and 2.5th-97.5th percentiles of the three
##     accepted inheritance parameters as the LHS sample size increases.
##
##   Stage1_pairwise_parameter_space.png
##     Pairwise accepted parameter distributions comparing samples identified
##     among the first 800 and final 200 points of the primary LHS design.
##
##   Stage1_independent_LHS_parameter_stability.png
##     Comparison of accepted inheritance-parameter distributions across the
##     primary 1000-point design and two independent 250-point designs.
##
##   Stage1_independent_LHS_pairwise.png
##     Pairwise comparison of the accepted parameter regions identified by
##     the three independent LHS designs.
##
## Data summaries are saved in:
##   stability_data/
##
## Main data outputs:
##
##   Stage1_independent_LHS_acceptance.csv
##     Sample size, number accepted and acceptance proportion for each
##     independent LHS design.
##
##   Stage1_independent_LHS_parameter_summary.csv
##     Median and 2.5th-97.5th percentiles of each inheritance parameter
##     among accepted samples for each independent LHS design.
##
## INTERPRETATION
##
## Stability of the acceptance proportion and accepted parameter distributions
## as the primary LHS sample size increases provides evidence that the sampled
## Stage 1 parameter space has been adequately characterised.
##
## Similar acceptance proportions and accepted parameter regions across the
## independently generated LHS designs provide evidence that the calibration
## results are not strongly dependent on the particular LHS design/random seed.
##
## Reported intervals are the 2.5th-97.5th percentiles of accepted parameter
## values and should not be interpreted as Bayesian credible intervals.
##
###############################################################################

## Output directories
figure_dir <- "stability_figures"
data_dir   <- "stability_data"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

filename <- "data_files/fitting_burnin_65000_UU.Rdata"
load(filename)

figure1name <- file.path(figure_dir, "Stage1_acceptance_stability.png")
figure2name <- file.path(figure_dir, "Stage1_parameter_convergence.png")
figure3name <- file.path(figure_dir, "Stage1_pairwise_parameter_space.png")

############################################################
## 1. Extract model outputs and parameter values
############################################################

## Model output used as part of the plausibility criterion:
## final proportion from the model for each LHS parameter set.
final_prop <- l[[1]][, 1]

## Store the three sampled parameters in a data frame.
## Simple internal names are retained for convenient coding;
## mathematical notation is used in the figures below.
pars <- data.frame(
  p1 = l[[2]][, 1],  # c_gg
  p2 = l[[2]][, 2],  # c^0_ss
  p3 = l[[2]][, 3]   # c^0_sg
)

############################################################
## 2. Identify plausible (accepted) parameter sets
############################################################

## A parameter set is accepted if:
##   1. the final model proportion is between 0.05 and 0.15;
##   2. c_gg > c^0_sg; and
##   3. c^0_sg > c^0_ss.
##
## Hence the required parameter ordering is:
##              c_gg > c^0_sg > c^0_ss

accepted <- final_prop >= 0.05 &
            final_prop <= 0.15 &
            pars$p1 > pars$p3 &
            pars$p3 > pars$p2

## Report the total number and overall proportion accepted.
sum(accepted)
mean(accepted)

############################################################
## 3. Cumulative acceptance rate
############################################################

## LHS sample number.
N <- seq_along(accepted)

## Calculate the proportion of all samples accepted up to
## each sample number.
cum_accept_rate <- cumsum(accepted) / N

############################################################
## Plot at selected LHS sample sizes
############################################################

Ns <- c(200, 400, 600, 800, 1000)
acceptance_rate_N <- cum_accept_rate[Ns]

############################################################
## Save figure
############################################################

png(
  figure1name,
  width = 1800,
  height = 1500,
  res = 300
)

plot(
  Ns,
  acceptance_rate_N,
  type = "b",
  xlim = c(200, 1000),
  ylim = c(0, 0.2),
  xlab = "Number of inheritance LHS samples",
  ylab = "Proportion of parameter sets accepted",
  cex.axis = 1.3,
  cex.lab = 1.4
)

dev.off()

############################################################
## 4. Stability of marginal accepted parameter distributions
############################################################

## Sample sizes at which the accepted parameter distributions
## will be compared.
Ns <- c(200, 400, 600, 800, 1000)

## Retain only sample sizes that are available in the data.
## This makes the code usable if fewer than 1000 samples are
## analysed in future.
Ns <- Ns[Ns <= nrow(pars)]

## At each sample size, calculate:
##   - the number of accepted parameter sets;
##   - the median of each accepted parameter; and
##   - the 2.5th and 97.5th percentiles of each accepted parameter.
##
## Stability of these summaries as N increases indicates that
## the marginal plausible parameter distributions are being
## characterised consistently.
parameter_summary <- do.call(
  rbind,
  lapply(Ns, function(N) {

    ## Retain the first N parameter samples.
    x <- pars[1:N, , drop = FALSE]

    ## Retain only samples satisfying the plausibility criteria.
    x <- x[accepted[1:N], , drop = FALSE]

    data.frame(
      N = N,
      n_accepted = nrow(x),
      p1_025 = quantile(x$p1, 0.025),
      p1_med = median(x$p1),
      p1_975 = quantile(x$p1, 0.975),
      
      p2_025 = quantile(x$p2, 0.025),
      p2_med = median(x$p2),
      p2_975 = quantile(x$p2, 0.975),
      
      p3_025 = quantile(x$p3, 0.025),
      p3_med = median(x$p3),
      p3_975 = quantile(x$p3, 0.975)
    )
  })
)

## Display the numerical summaries.
parameter_summary

############################################################
## 5. Plot convergence of marginal parameter distributions
############################################################

## Function to plot the median and 2.5th/97.5th percentiles
## of one parameter against the number of LHS samples.
##
## 'parameter' gives the internal parameter name (p1, p2, p3).
## 'label' gives the mathematical notation shown on the figure.
plot_convergence <- function(summary, parameter, label,
                             xlim = NULL, ylim = NULL) {

  med <- summary[[paste0(parameter, "_med")]]
  lo  <- summary[[paste0(parameter, "_025")]]
  hi  <- summary[[paste0(parameter, "_975")]]
  
  ## If limits are not supplied, use the data range
  if (is.null(xlim))
    xlim <- range(summary$N)
  
  if (is.null(ylim))
    ylim <- range(lo, hi)

  ## Plot the median.
  plot(summary$N, med,
       type = "b",
       ylim = ylim,
       xlab = "Number of LHS samples",
       ylab = label,
       cex.axis = 1.3,
       cex.lab = 1.4)

  ## Add the 2.5th percentile.
  lines(summary$N, lo,
        type = "b",
        lty = 2)

  ## Add the 97.5th percentile.
  lines(summary$N, hi,
        type = "b",
        lty = 2)
}

## Create four layout regions:
##   panels 1-3 = convergence plots for the three parameters;
##   panel 4   = common legend.
##

png(figure2name,
    width = 3600,
    height = 1000,
    res = 300)

## The fourth panel is narrower so that the legend sits to
## the right of the three main plotting panels.
layout(matrix(1:4, nrow = 1),
       widths = c(1, 1, 1, 0.7))

## Margins for the three parameter plots.
par(mar = c(4.5, 6.5, 2, 1))

plot_convergence(parameter_summary, "p1",
                 expression(c[gg]),
                 ylim = c(0.9, 1))

plot_convergence(parameter_summary, "p2",
                 expression(c[ss]^0),
                 ylim = c(0, 1))

plot_convergence(parameter_summary, "p3",
                 expression(c[sg]^0),
                 ylim = c(0, 1))

## Use the fourth panel for the common legend.
par(mar = c(0, 0, 0, 0))
plot.new()

legend("center",
       legend = c("Median",
                  "2.5th/97.5th percentile"),
       lty = c(1, 2),
       pch = c(1, 1),
       bty = "n")

dev.off()
## Reset to the standard single-panel layout and margins.
layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))

############################################################
## 6. Stability of the joint plausible parameter region
############################################################

## Extract all accepted parameter sets.
accpars <- pars[accepted, , drop = FALSE]

## Recover the original LHS sample number associated with
## each accepted parameter set.
accN <- which(accepted)

## Compare accepted parameter sets obtained among the first
## 800 LHS samples with those obtained among the final 200.
##
## If 1000 samples provide adequate coverage, the accepted
## samples from the final 200 should generally occupy the same
## regions of parameter space as those from the first 800,
## rather than revealing distinct new plausible regions.
early <- accN <= 800
late  <- accN > 800

############################################################
## 7. Pairwise plots of accepted parameter values
############################################################

## Again use three plotting panels followed by a narrower
## fourth panel containing the common legend.
png(figure3name,
    width = 3600,
    height = 1000,
    res = 300)

layout(matrix(1:4, nrow = 1),
       widths = c(1, 1, 1, 0.7))

par(mar = c(4.5, 6.5, 2, 1))

## ---- c_gg versus c^0_ss ----

plot(accpars$p1[early],
     accpars$p2[early],
     pch = 1,
     xlim = c(0.9, 1),
     ylim = c(0, 1),
     xlab = expression(c[gg]),
     ylab = expression(c[ss]^0),
     cex.axis = 1.3,
     cex.lab = 1.4)

## Overlay accepted samples arising from the final 200
## LHS samples using filled points.
points(accpars$p1[late],
       accpars$p2[late],
       pch = 16)

## ---- c_gg versus c^0_sg ----

plot(accpars$p1[early],
     accpars$p3[early],
     pch = 1,
     xlim = c(0.9, 1),
     ylim = c(0, 1),
     xlab = expression(c[gg]),
     ylab = expression(c[sg]^0),
     cex.axis = 1.3,
     cex.lab = 1.4)

points(accpars$p1[late],
       accpars$p3[late],
       pch = 16)

## ---- c^0_ss versus c^0_sg ----

plot(accpars$p2[early],
     accpars$p3[early],
     pch = 1,
     xlim = c(0, 1),
     ylim = c(0, 1),
     xlab = expression(c[ss]^0),
     ylab = expression(c[sg]^0),
     cex.axis = 1.3,
     cex.lab = 1.4)

points(accpars$p2[late],
       accpars$p3[late],
       pch = 16)

## Use the fourth panel for the common legend.
par(mar = c(0, 0, 0, 0))
plot.new()

legend("center",
       legend = c("First 800",
                  "Final 200"),
       pch = c(1, 16),
       bty = "n")

dev.off()

## Reset graphical layout and margins.
layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))

############################################################
## 8. Run two additional 250-point LHS designs
############################################################
# Before running this, change the number of LHSsamples <- 250 
# in Fitting_phenotype_parameters.R

seednumber<-123
ABR<-65000 # 80% prev

# First stage: generate 250 LHS parameter sets for c_ij, 
# calculate and store summary stats

filename<-"data_files/fitting_burnin_65000_additional_LHS1.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_phenotype_parameters.R")
}else{
  load(filename)
}

# Run burnins for LHS samples that satisfy initial SOR proportion criteria
# of between 5-15% SOR worms, and c_gg > c0sg > c0_ss
# 
seednumber<-123456
filename<-"data_files/burnin_157_ABR_65000_additional_LHS1.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/create_burnins_from_posterior_samples.R")
} else {
  load(filename)
}

# And repeat for another random seed

seednumber<-456
ABR<-65000 # 80% prev

filename<-"data_files/fitting_burnin_65000_additional_LHS2.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/Fitting_phenotype_parameters.R")
}else{
  load(filename)
}

seednumber<-456123
filename<-"data_files/burnin_157_ABR_65000_additional_LHS2.Rdata"
if (file.exists(filename)==FALSE){
  source("r_files/create_burnins_from_posterior_samples.R")
} else {
  load(filename)
}

############################################################
## 9. Random-seed / independent LHS design analysis
############################################################

load("data_files/fitting_burnin_65000_UU.Rdata")
l_primary <- l

############################################################
## 9.1 Function to process a Stage 1 LHS design
############################################################

process_stage1_LHS <- function(l, design_name) {
  
  ## Model output used for Belief 1
  final_prop <- l[[1]][, 1]
  
  ## Inheritance parameter values
  pars <- data.frame(
    c_gg  = l[[2]][, 1],
    c0_ss = l[[2]][, 2],
    c0_sg = l[[2]][, 3]
  )
  
  ## Stage 1 acceptance criterion:
  ##   Belief 1: 0.05 < SOR proportion < 0.15
  ##   inheritance ordering: c_gg > c0_sg > c0_ss
  accepted <- final_prop > 0.05 &
    final_prop < 0.15 &
    pars$c_gg > pars$c0_sg &
    pars$c0_sg > pars$c0_ss
  
  ## Extract accepted parameter sets
  accpars <- pars[accepted, , drop = FALSE]
  
  ## Calculate parameter summaries
  parameter_summary <- data.frame(
    parameter = c("c_gg", "c0_ss", "c0_sg"),
    
    lower = c(
      quantile(accpars$c_gg,  0.025),
      quantile(accpars$c0_ss, 0.025),
      quantile(accpars$c0_sg, 0.025)
    ),
    
    median = c(
      median(accpars$c_gg),
      median(accpars$c0_ss),
      median(accpars$c0_sg)
    ),
    
    upper = c(
      quantile(accpars$c_gg,  0.975),
      quantile(accpars$c0_ss, 0.975),
      quantile(accpars$c0_sg, 0.975)
    )
  )
  
  ## Return all information needed for comparisons
  list(
    design = design_name,
    N = length(accepted),
    n_accepted = sum(accepted),
    acceptance_rate = mean(accepted),
    accepted_parameters = accpars,
    parameter_summary = parameter_summary
  )
}

############################################################
## 9.2 Process the original 1000-point LHS design
############################################################

primary <- process_stage1_LHS(
  l_primary,
  "Primary (N=1000)"
)

############################################################
## 9.3 Load and process additional 250-point design 1
############################################################

load("data_files/fitting_burnin_65000_additional_LHS1.Rdata")

## The load() command replaces l with the data from this file
l_LHS1 <- l

LHS1 <- process_stage1_LHS(
  l_LHS1,
  "Independent 1 (N=250)"
)

############################################################
## 9.4 Load and process additional 250-point design 2
############################################################

load("data_files/fitting_burnin_65000_additional_LHS2.Rdata")

l_LHS2 <- l

LHS2 <- process_stage1_LHS(
  l_LHS2,
  "Independent 2 (N=250)"
)

############################################################
## 9.5 Compare acceptance proportions
############################################################

stage1_acceptance_summary <- data.frame(
  
  design = c(
    primary$design,
    LHS1$design,
    LHS2$design
  ),
  
  N = c(
    primary$N,
    LHS1$N,
    LHS2$N
  ),
  
  n_accepted = c(
    primary$n_accepted,
    LHS1$n_accepted,
    LHS2$n_accepted
  ),
  
  acceptance_rate = c(
    primary$acceptance_rate,
    LHS1$acceptance_rate,
    LHS2$acceptance_rate
  )
)

print(stage1_acceptance_summary)

############################################################
## 9.6 Save acceptance summary
############################################################

write.csv(
  stage1_acceptance_summary,
  file.path(data_dir, "Stage1_independent_LHS_acceptance.csv"),
  row.names = FALSE
)

############################################################
## 9.7 Compare accepted parameter distributions
############################################################

## Add design labels to each parameter summary
primary$parameter_summary$design <- primary$design
LHS1$parameter_summary$design <- LHS1$design
LHS2$parameter_summary$design <- LHS2$design

stage1_parameter_summary <- rbind(
  primary$parameter_summary,
  LHS1$parameter_summary,
  LHS2$parameter_summary
)

print(stage1_parameter_summary)

write.csv(
  stage1_parameter_summary,
  file.path(data_dir, "Stage1_independent_LHS_parameter_summary.csv"),
  row.names = FALSE
)

############################################################
## 9.8 Plot parameter distributions across LHS designs
############################################################

## Function to plot median and 95% uncertainty interval
## for one parameter across the three independent designs
plot_design_comparison <- function(parameter, label, ylim) {
  
  x <- stage1_parameter_summary[
    stage1_parameter_summary$parameter == parameter,
    ,
    drop = FALSE
  ]
  
  ## Positions shifted well away from panel boundaries
  xpos <- c(1.4, 2.0, 2.6)
  
  plot(
    xpos,
    x$median,
    pch = 16,
    ylim = ylim,
    xlim = c(1, 3),
    xaxs = "i",
    xaxt = "n",
    xlab = "LHS design",
    ylab = label,
    cex.axis = 1.3,
    cex.lab = 1.4,
    mgp = c(2.7, 0.8, 0)
  )
  
  ## 2.5th-97.5th percentile intervals
  arrows(
    xpos,
    x$lower,
    xpos,
    x$upper,
    angle = 90,
    code = 3,
    length = 0.05
  )
  
  axis(
    1,
    at = xpos,
    labels = c(
      "Primary\nN=1000",
      "Indep. 1\nN=250",
      "Indep. 2\nN=250"
    ),
    cex.axis = 1.1
  )
}

############################################################
## 9.9 Three-panel comparison figure
############################################################

png(
  file.path(figure_dir, "Stage1_independent_LHS_parameter_stability.png"),
  width = 3200,
  height = 1200,
  res = 300
)

par(
  mfrow = c(1, 3),
  
  ## bottom, left, top, right
  ## Extra left margin prevents y labels being clipped
  mar = c(6, 6.5, 2, 1.5),
  
  ## Outer margin gives the whole figure some breathing room
  oma = c(0, 1, 0, 1)
)

plot_design_comparison(
  "c_gg",
  expression(c[gg]),
  ylim = c(0.9, 1)
)

plot_design_comparison(
  "c0_ss",
  expression(c[ss]^0),
  ylim = c(0.01, 1)
)

plot_design_comparison(
  "c0_sg",
  expression(c[sg]^0),
  ylim = c(0.01, 1)
)

dev.off()

par(
  mfrow = c(1, 1),
  mar = c(5.1, 4.1, 4.1, 2.1),
  oma = c(0, 0, 0, 0)
)

############################################################
## 9.10 Pairwise distributions across independent designs
############################################################

## Add design identifier to accepted parameter sets
acc_primary <- primary$accepted_parameters
acc_primary$design <- "Primary"

acc_LHS1 <- LHS1$accepted_parameters
acc_LHS1$design <- "Independent 1"

acc_LHS2 <- LHS2$accepted_parameters
acc_LHS2$design <- "Independent 2"

############################################################
## Plot all three independent designs using different
## plotting symbols
############################################################

png(
  file.path(figure_dir, "Stage1_independent_LHS_pairwise.png"),
  width = 3600,
  height = 1200,
  res = 300
)

## Three plots plus a narrow fourth panel for the legend
layout(
  matrix(1:4, nrow = 1),
  widths = c(1, 1, 1, 0.8)
)

par(
  mar = c(4.5, 4.5, 2, 1)
)


## ---- c_gg versus c0_ss ----

plot(
  acc_primary$c_gg,
  acc_primary$c0_ss,
  pch = 1,
  xlim = range(
    c(acc_primary$c_gg,
      acc_LHS1$c_gg,
      acc_LHS2$c_gg)
  ),
  ylim = range(
    c(acc_primary$c0_ss,
      acc_LHS1$c0_ss,
      acc_LHS2$c0_ss)
  ),
  xlab = expression(c[gg]),
  ylab = expression(c[ss]^0),
  cex.axis = 1.3,
  cex.lab = 1.4
)

points(
  acc_LHS1$c_gg,
  acc_LHS1$c0_ss,
  pch = 16
)

points(
  acc_LHS2$c_gg,
  acc_LHS2$c0_ss,
  pch = 2
)

## ---- c_gg versus c0_sg ----

plot(
  acc_primary$c_gg,
  acc_primary$c0_sg,
  pch = 1,
  xlim = range(
    c(acc_primary$c_gg,
      acc_LHS1$c_gg,
      acc_LHS2$c_gg)
  ),
  ylim = range(
    c(acc_primary$c0_sg,
      acc_LHS1$c0_sg,
      acc_LHS2$c0_sg)
  ),
  xlab = expression(c[gg]),
  ylab = expression(c[sg]^0),
  cex.axis = 1.3,
  cex.lab = 1.4
)

points(
  acc_LHS1$c_gg,
  acc_LHS1$c0_sg,
  pch = 16
)

points(
  acc_LHS2$c_gg,
  acc_LHS2$c0_sg,
  pch = 2
)

## ---- c0_ss versus c0_sg ----

plot(
  acc_primary$c0_ss,
  acc_primary$c0_sg,
  pch = 1,
  xlim = range(
    c(acc_primary$c0_ss,
      acc_LHS1$c0_ss,
      acc_LHS2$c0_ss)
  ),
  ylim = range(
    c(acc_primary$c0_sg,
      acc_LHS1$c0_sg,
      acc_LHS2$c0_sg)
  ),
  xlab = expression(c[ss]^0),
  ylab = expression(c[sg]^0),
  cex.axis = 1.3,
  cex.lab = 1.4
)

points(
  acc_LHS1$c0_ss,
  acc_LHS1$c0_sg,
  pch = 16
)

points(
  acc_LHS2$c0_ss,
  acc_LHS2$c0_sg,
  pch = 2
)

############################################################
## Common legend
############################################################

par(mar = c(0, 0, 0, 0))

plot.new()

legend(
  "center",
  legend = c(
    "Primary (N=1000)",
    "Independent 1 (N=250)",
    "Independent 2 (N=250)"
  ),
  pch = c(1, 16, 2),
  bty = "n",
  cex = 1.2
)

dev.off()

## Reset plotting layout
layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))


###############################################################################
## REPOSITORY DATA EXPORT
## Run after all Stage 1 analyses above.
###############################################################################
dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)

stage1_primary_pairwise <- data.frame(
  c_gg = accpars$p1,
  c0_ss = accpars$p2,
  c0_sg = accpars$p3,
  sampling_block = ifelse(accN <= 800, "First 800", "Final 200")
)

stage1_independent_pairwise <- rbind(
  transform(acc_primary, design = "Primary"),
  transform(acc_LHS1, design = "Independent 1"),
  transform(acc_LHS2, design = "Independent 2")
)

stage1_acceptance_stability <- data.frame(
  N = Ns,
  acceptance_rate = acceptance_rate_N
)

save(
  stage1_acceptance_stability,
  parameter_summary,
  stage1_primary_pairwise,
  stage1_acceptance_summary,
  stage1_parameter_summary,
  stage1_independent_pairwise,
  file = "data_processed/stage1_stability_repo_data.Rdata"
)

message("Saved data_processed/stage1_stability_repo_data.Rdata")
