###############################################################################
## Stage 2 calibration diagnostics and sensitivity analyses
###############################################################################
##
## PURPOSE
##
## This script performs diagnostic and sensitivity analyses for Stage 2 of the
## variable-response model calibration. Stage 2 combines the inheritance
## parameter sets retained following Stage 1 with Latin hypercube samples (LHS)
## of the four MDA-effect parameters.
##
## The primary Stage 2 calibration requires:
##   Belief 2: final mf prevalence = 0.10-0.12
##   Belief 3: increase in SOR worm proportion > 0.10
## together with the final inheritance-parameter constraints:
##   c0_SORSOR < 0.5
##   c0_GRSOR  > 0.5
##
## ANALYSES
##
## 1. Primary Stage 2 LHS sampling adequacy
##    - Evaluates cumulative LHS sample sizes N = 20, 40, 60, 80 and 100.
##    - Examines stability of:
##        * calibration acceptance proportions;
##        * final retention proportions;
##        * medians and 2.5th-97.5th percentiles of retained parameters.
##    - Compares parameter space represented by the first 80 and final 20
##      MDA-effect LHS points.
##
## 2. Parameter dependence and practical identifiability
##    - Calculates pairwise Spearman correlations among the seven calibrated
##      parameters in the final retained primary parameter combinations.
##    - Compares the central 95% interval of retained values for each parameter
##      with its prior range as a diagnostic of practical identifiability.
##
## 3. Independent Stage 2 LHS design comparison
##    - Compares the primary 100-point MDA-effect LHS design with the
##      independently generated 50-point LHS design.
##    - Compares acceptance proportions, retained parameter distributions,
##      and pairwise parameter space.
##
## 4. Sensitivity analyses
##    The primary calibration is compared with four alternative calibrations:
##
##      Moderate increase:
##        final mf prevalence = 0.10-0.12
##        increase in SOR worm proportion = 0.05-0.10
##
##      Weak increase:
##        final mf prevalence = 0.10-0.12
##        increase in SOR worm proportion < 0.05
##
##      No offspring-phenotype selection:
##        s_SORSOR = s_SORGR = 0
##        final mf prevalence = 0.10-0.12
##        no Belief 3 criterion is imposed
##
##      Wider prevalence:
##        final mf prevalence = 0.08-0.14
##        increase in SOR worm proportion > 0.10
##
##    For these analyses the script:
##      - assesses stability of acceptance proportions with increasing LHS size;
##      - calculates final acceptance/retention proportions;
##      - compares retained parameter distributions across the primary and
##        sensitivity calibrations.
##
## OUTPUT
##
## Figures are saved in:
##   stability_figures/
##
## Main figure outputs:
##   Stage2_acceptance_stability_primary.png
##     Primary Stage 2 acceptance-rate stability.
##
##   Stage2_parameter_stability_100_all_parameters.png
##     Stability of medians and 2.5th-97.5th percentiles of all seven
##     retained parameters as the number of MDA-effect LHS samples increases.
##
##   Stage2_pairwise_stability_100.png
##     Pairwise comparison of retained MDA-effect parameters arising from
##     the first 80 versus final 20 LHS points.
##
##   Stage2_independent_LHS_parameter_stability.png
##     Comparison of retained parameter distributions between the independent
##     50-point and primary 100-point Stage 2 LHS designs.
##
##   Stage2_independent_LHS_pairwise.png
##     Pairwise parameter-space comparison of the independent 50- and
##     100-point Stage 2 LHS designs.
##
##   Stage2_acceptance_stability_sensitivity_analyses.png
##     Acceptance-rate stability for the four sensitivity calibrations.
##
##   Stage2_sensitivity_parameter_distributions.png
##     Comparison of retained parameter distributions across the primary
##     calibration and four sensitivity analyses.
##
## Data summaries are saved in:
##   stability_data/
##
## Main data outputs:
##   Stage2_acceptance_stability_100.csv
##   Stage2_parameter_stability_100.csv
##   Stage2_parameter_Spearman_correlations.csv
##   Stage2_parameter_identifiability_summary.csv
##   Stage2_independent_LHS_parameter_summary.csv
##   Stage2_sensitivity_parameter_summary.csv
##   Stage2_sensitivity_acceptance_summary.csv
##
## INTERPRETATION
##
## These analyses assess numerical sampling stability, parameter dependence,
## practical identifiability and sensitivity to the calibration assumptions.
## The identifiability analysis is conditional on the specified model,
## parameter prior ranges and calibration criteria; it is not a structural
## identifiability analysis.
##
## All reported parameter intervals are the 2.5th-97.5th percentiles of
## retained parameter values and should not be interpreted as Bayesian
## credible intervals.
##
###############################################################################

# ---- Output directories ----
figure_dir <- "stability_figures"
data_dir   <- "stability_data"
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Load Stage 1 data ----
load("data_files/burnin_157_ABR_65000_UU.Rdata")
load("data_files/fitting_burnin_65000_UU.Rdata")

indv2 <- which(
  l[[2]][,1] - l[[2]][,3] > 0 &
  l[[2]][,3] - l[[2]][,2] > 0 &
  l[[1]][,1] > 0.05 &
  l[[1]][,1] < 0.15
)
num_samples <- length(indv2)
message("Stage 1 inheritance sets: ", num_samples)

source("equations.R")
MDAi <- 1
ABR <- 65000
source("fixed-parameters.R")

# Reconstruct pre-MDA SOR proportion and inheritance parameters
preMDA <- vector("list", num_samples)
for (n in 1:num_samples) {
  yout <- OutputB[[n]]
  samp <- indv2[n]
  sor <- 1 - (yout[,3] + yout[,4]) /
    (yout[,3] + yout[,4] + yout[,7] + yout[,8])

  preMDA[[n]] <- c(
    tail(sor, 1), NA, NA, NA, NA,
    l[[2]][samp,1], l[[2]][samp,2], l[[2]][samp,3]
  )
}

# ---- Load Stage 2 cumulative-effect simulations ----
load("data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100.Rdata")
postMDA <- ll

# Number of unique Stage 2 MDA-effect LHS points examined
Ns <- c(20, 40, 60, 80, 100)

stage2_summary <- function(N) {
  z <- vector("list", num_samples * N)
  k <- 1

  for (n in 1:num_samples) {
    initial_SOR <- preMDA[[n]][1]
    c_GRGR <- preMDA[[n]][6]
    c0_SORSOR <- preMDA[[n]][7]
    c0_GRSOR <- preMDA[[n]][8]

    for (i in 1:N) {
      final_SOR <- postMDA[[n]][i,1]
      final_prev <- postMDA[[n]][i,5]
      delta_SOR <- final_SOR - initial_SOR

      belief23 <- final_prev > 0.10 &
        final_prev < 0.12 &
        delta_SOR > 0.10
        #delta_SOR > 0.05 & delta_SOR < 0.10
        #delta_SOR < 0.05

      constraints <- c0_SORSOR < 0.5 & c0_GRSOR > 0.5

      z[[k]] <- data.frame(
        lhs_sample = i,
        inheritance_sample = n,
        c_GRGR = c_GRGR,
        c0_SORSOR = c0_SORSOR,
        c0_GRSOR = c0_GRSOR,
        s_SORSOR = postMDA[[n]][i,6],
        s_SORGR = postMDA[[n]][i,7],
        f = postMDA[[n]][i,8],
        y = postMDA[[n]][i,9],
        final_prev = final_prev,
        delta_SOR = delta_SOR,
        belief23 = belief23,
        retained = belief23 & constraints
      )
      k <- k + 1
    }
  }

  d <- do.call(rbind, z)
  acc <- d[d$retained, , drop = FALSE]
  pars <- c("c_GRGR","c0_SORSOR","c0_GRSOR",
            "s_SORSOR","s_SORGR","f","y")

  q <- if (nrow(acc) > 0) {
    sapply(acc[,pars,drop=FALSE], quantile,
           probs=c(0.025,0.5,0.975), na.rm=TRUE)
  } else {
    matrix(NA, 3, length(pars),
           dimnames=list(c("2.5%","50%","97.5%"), pars))
  }

  list(
    N=N,
    n_evaluated=nrow(d),
    n_belief23=sum(d$belief23),
    n_retained=sum(d$retained),
    acceptance_belief23=mean(d$belief23),
    acceptance_final=mean(d$retained),
    quantiles=q,
    accepted=acc
  )
}

results <- lapply(Ns, stage2_summary)

# ---- Acceptance summary ----
acceptance_summary <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    N_MDA_LHS=x$N,
    combinations_evaluated=x$n_evaluated,
    n_Beliefs_2_3=x$n_belief23,
    acceptance_Beliefs_2_3=x$acceptance_belief23,
    n_final_retained=x$n_retained,
    final_retention_proportion=x$acceptance_final
  )
}))
print(acceptance_summary)

# At N = 100 this should reproduce:
# 15700 combinations evaluated
# 962 satisfying Beliefs 2 and 3
# 724 finally retained

write.csv(acceptance_summary, file.path(data_dir, "Stage2_acceptance_stability_100.csv"),
          row.names=FALSE)

# ---- Acceptance stability plot ----
png(file.path(figure_dir, "Stage2_acceptance_stability_100.png"),
    width=1800, height=1500, res=300)

plot(acceptance_summary$N_MDA_LHS,
     acceptance_summary$acceptance_Beliefs_2_3,
     type="b",
     ylim=range(c(acceptance_summary$acceptance_Beliefs_2_3,
                  acceptance_summary$final_retention_proportion)),
     xlab="Number of MDA-effect LHS samples",
     ylab="Proportion of parameter combinations",
     cex.axis=1.3, cex.lab=1.4)

lines(acceptance_summary$N_MDA_LHS,
      acceptance_summary$final_retention_proportion,
      type="b", lty=2)

legend("topright",
       legend=c("Satisfy Beliefs 2 and 3","Final retained"),
       lty=c(1,2), pch=c(1,1), bty="n", cex=1.2)
dev.off()

# ---- Parameter quantile stability ----
parameters <- c("c_GRGR","c0_SORSOR","c0_GRSOR",
                "s_SORSOR","s_SORGR","f","y")

parameter_stability <- NULL
for (j in seq_along(results)) {
  q <- results[[j]]$quantiles
  for (p in parameters) {
    parameter_stability <- rbind(
      parameter_stability,
      data.frame(
        N=results[[j]]$N,
        parameter=p,
        lower=q["2.5%",p],
        median=q["50%",p],
        upper=q["97.5%",p]
      )
    )
  }
}
print(parameter_stability)
write.csv(parameter_stability, file.path(data_dir, "Stage2_parameter_stability_100.csv"),
          row.names=FALSE)



plot_stage2_convergence <- function(parameter, label, ylim) {
  
  x <- parameter_stability[
    parameter_stability$parameter == parameter,
  ]
  
  plot(
    x$N,
    x$median,
    type = "b",
    ylim = ylim,
    xlab = "Number of MDA-effect LHS samples",
    ylab = label,
    cex.axis = 1.3,
    cex.lab = 1.4
  )
  
  lines(
    x$N,
    x$lower,
    type = "b",
    lty = 2
  )
  
  lines(
    x$N,
    x$upper,
    type = "b",
    lty = 2
  )
}

############################################################
## Stability plots for all seven retained parameters
############################################################

png(
  file.path(figure_dir, "Stage2_parameter_stability_100_all_parameters.png"),
  width = 3600,
  height = 2200,
  res = 300
)

## Top row:
## c_GRGR | c0_SORSOR | c0_GRSOR | legend
##
## Bottom row:
## s_SORSOR | s_SORGR | f | y

layout(
  matrix(
    c(
      1, 2, 3, 4,
      5, 6, 7, 8
    ),
    nrow = 2,
    byrow = TRUE
  ),
  widths = c(1, 1, 1, 1)
)

par(
  mar = c(5, 6, 2, 1)
)

############################################################
## Top row: inheritance parameters
############################################################

plot_stage2_convergence(
  "c_GRGR",
  expression(c[GR*",GR"]),
  c(0.9, 1)
)

plot_stage2_convergence(
  "c0_SORSOR",
  expression(c[SOR*",SOR"]^0),
  c(0.01, 1)
)

plot_stage2_convergence(
  "c0_GRSOR",
  expression(c[GR*",SOR"]^0),
  c(0.01, 1)
)

############################################################
## Legend in top-right panel
############################################################

par(mar = c(0, 0, 0, 0))

plot.new()

legend(
  "center",
  legend = c(
    "Median",
    "2.5th/97.5th percentile"
  ),
  lty = c(1, 2),
  pch = c(1, 1),
  bty = "n",
  cex = 1.2
)

############################################################
## Bottom row: MDA-effect parameters
############################################################

par(
  mar = c(5, 6, 2, 1)
)

plot_stage2_convergence(
  "s_SORSOR",
  expression(s[SOR*",SOR"]),
  c(0, 20)
)

plot_stage2_convergence(
  "s_SORGR",
  expression(s[SOR*",GR"]),
  c(0, 20)
)

plot_stage2_convergence(
  "f",
  expression(f),
  c(1, 10)
)

plot_stage2_convergence(
  "y",
  expression(y),
  c(0, 1)
)

dev.off()

layout(1)

par(
  mar = c(5.1, 4.1, 4.1, 2.1)
)

# ---- Pairwise stability: first 80 vs final 20 unique LHS points ----
acc100_stability <- results[[which(Ns==100)]]$accepted
early <- acc100_stability$lhs_sample <= 80
late <- acc100_stability$lhs_sample > 80

message("Retained combinations from first 80 LHS points: ", sum(early))
message("Retained combinations from final 20 LHS points: ", sum(late))

pairs_to_plot <- list(
  c("s_SORSOR","s_SORGR"),
  c("s_SORSOR","f"),
  c("s_SORSOR","y"),
  c("s_SORGR","f"),
  c("s_SORGR","y"),
  c("f","y")
)

labs <- list(
  s_SORSOR=expression(s[SOR*",SOR"]),
  s_SORGR=expression(s[SOR*",GR"]),
  f=expression(f),
  y=expression(y)
)

png(file.path(figure_dir, "Stage2_pairwise_stability_100.png"),
    width=3600, height=2000, res=300)

layout(matrix(c(1, 2, 3, 7,
                4, 5, 6, 8),nrow=2,byrow=TRUE),
       widths=c(1,1,1,0.7))
par(mar=c(4.5,4.5,2,1))

for (z in pairs_to_plot) {
  xv <- z[1]
  yv <- z[2]

  plot(acc100_stability[[xv]][early], acc100_stability[[yv]][early],
       pch=1,
       xlim=range(acc100_stability[[xv]],na.rm=TRUE),
       ylim=range(acc100_stability[[yv]],na.rm=TRUE),
       xlab=labs[[xv]], ylab=labs[[yv]],
       cex.axis=1.2, cex.lab=1.3)

  points(acc100_stability[[xv]][late], acc100_stability[[yv]][late], pch=16)
}

par(mar=c(0,0,0,0))
plot.new()
legend("center", legend=c("First 80","Final 20"),
       pch=c(1,16), bty="n", cex=1.2)
plot.new()
dev.off()

layout(1)
par(mar=c(5.1,4.1,4.1,2.1))

# End of Stage 2 stability analysis

############################################################
## Parameter correlations and practical identifiability
############################################################

## Final retained parameter combinations from the complete
## 100-point Stage 2 LHS design
acc100 <- results[[which(Ns == 100)]]$accepted

## Extract the seven calibrated parameters
param_data <- acc100[, c(
  "c_GRGR",
  "c0_SORSOR",
  "c0_GRSOR",
  "s_SORSOR",
  "s_SORGR",
  "f",
  "y"
)]

############################################################
## Pairwise Spearman correlations
############################################################

spearman_cor <- cor(
  param_data,
  method = "spearman",
  use = "complete.obs"
)

print(round(spearman_cor, 2))

write.csv(
  round(spearman_cor, 3),
  file.path(data_dir, "Stage2_parameter_Spearman_correlations.csv")
)

############################################################
## Practical identifiability: retained vs prior ranges
############################################################

## These names must exactly match the columns in param_data
prior_ranges <- data.frame(
  parameter = c(
    "c_GRGR",
    "c0_SORSOR",
    "c0_GRSOR",
    "s_SORSOR",
    "s_SORGR",
    "f",
    "y"
  ),
  
  prior_min = c(
    0.9,   # c_GRGR
    0.01,  # c0_SORSOR
    0.01,  # c0_GRSOR
    0,     # s_SORSOR
    0,     # s_SORGR
    1,     # f
    0      # y
  ),
  
  prior_max = c(
    1,     # c_GRGR
    1,     # c0_SORSOR
    1,     # c0_GRSOR
    20,    # s_SORSOR
    20,    # s_SORGR
    10,    # f
    1      # y
  )
)


############################################################
## Calculate retained parameter distributions
############################################################

identifiability_summary <- data.frame(
  parameter = names(param_data),
  
  q025 = sapply(
    param_data,
    quantile,
    probs = 0.025,
    na.rm = TRUE
  ),
  
  median = sapply(
    param_data,
    median,
    na.rm = TRUE
  ),
  
  q975 = sapply(
    param_data,
    quantile,
    probs = 0.975,
    na.rm = TRUE
  ),
  
  row.names = NULL
)


############################################################
## Add prior ranges
############################################################

## match() is preferable here because it preserves the
## parameter order and makes mismatches easier to detect
ii <- match(
  identifiability_summary$parameter,
  prior_ranges$parameter
)

## Check that every parameter has been matched
if (any(is.na(ii))) {
  stop(
    "Prior ranges missing for: ",
    paste(
      identifiability_summary$parameter[is.na(ii)],
      collapse = ", "
    )
  )
}

identifiability_summary$prior_min <-
  prior_ranges$prior_min[ii]

identifiability_summary$prior_max <-
  prior_ranges$prior_max[ii]


############################################################
## Compare retained and prior widths
############################################################

## Width of prior range
identifiability_summary$prior_width <-
  identifiability_summary$prior_max -
  identifiability_summary$prior_min

## Width of central 95% of retained parameter values
identifiability_summary$retained_95_width <-
  identifiability_summary$q975 -
  identifiability_summary$q025

## Width of retained 95% interval relative to prior width
identifiability_summary$relative_width <-
  identifiability_summary$retained_95_width /
  identifiability_summary$prior_width

## Express as a percentage for easier interpretation
identifiability_summary$relative_width_percent <-
  100 * identifiability_summary$relative_width


############################################################
## Display and save
############################################################

print(identifiability_summary)

write.csv(
  identifiability_summary,
  file.path(data_dir, "Stage2_parameter_identifiability_summary.csv"),
  row.names = FALSE
)

############################################################
## Compare retained parameter distributions between the
## independent 50-point and 100-point Stage 2 LHS designs
############################################################


############################################################
## Primary 100-point retained set is already stored in acc100
############################################################

############################################################
## 2. Load original independent 50-point Stage 2 design
############################################################

load(
  "data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu.Rdata"
)

postMDA50 <- ll


############################################################
## 3. Process the 50-point design using the same criteria
############################################################

results50 <- NULL

for (n in 1:num_samples) {
  
  ## Pre-MDA SOR worm proportion
  initial_SOR <- preMDA[[n]][1]
  
  ## Inheritance parameters
  c_GRGR    <- preMDA[[n]][6]
  c0_SORSOR <- preMDA[[n]][7]
  c0_GRSOR  <- preMDA[[n]][8]
  
  for (i in 1:50) {
    
    ## Model outputs
    final_SOR  <- postMDA50[[n]][i, 1]
    final_prev <- postMDA50[[n]][i, 5]
    
    ## Change in SOR worm proportion
    delta_SOR <- final_SOR - initial_SOR
    
    ## Beliefs 2 and 3
    belief23 <-
      final_prev > 0.10 &
      final_prev < 0.12 &
      delta_SOR > 0.10
    
    ## Additional parameter constraints
    constraints <-
      c0_SORSOR < 0.5 &
      c0_GRSOR > 0.5
    
    ## Store final retained combinations only
    if (belief23 & constraints) {
      
      results50 <- rbind(
        results50,
        data.frame(
          c_GRGR = c_GRGR,
          c0_SORSOR = c0_SORSOR,
          c0_GRSOR = c0_GRSOR,
          s_SORSOR = postMDA50[[n]][i, 6],
          s_SORGR  = postMDA50[[n]][i, 7],
          f         = postMDA50[[n]][i, 8],
          y         = postMDA50[[n]][i, 9]
        )
      )
    }
  }
}

acc50_independent <- results50

## Check that this reproduces the expected number retained
nrow(acc50_independent)

## Should be 399


############################################################
## 4. Function to calculate parameter summaries
############################################################

summarise_design <- function(x, design) {
  
  parameters <- c(
    "c_GRGR",
    "c0_SORSOR",
    "c0_GRSOR",
    "s_SORSOR",
    "s_SORGR",
    "f",
    "y"
  )
  
  output <- NULL
  
  for (p in parameters) {
    
    output <- rbind(
      output,
      data.frame(
        design = design,
        parameter = p,
        
        lower = quantile(
          x[[p]],
          0.025,
          na.rm = TRUE
        ),
        
        median = median(
          x[[p]],
          na.rm = TRUE
        ),
        
        upper = quantile(
          x[[p]],
          0.975,
          na.rm = TRUE
        )
      )
    )
  }
  
  output
}


############################################################
## 5. Calculate summaries for both independent designs
############################################################

summary50 <- summarise_design(
  acc50_independent,
  "Independent\nN=50"
)

summary100 <- summarise_design(
  acc100,
  "Primary\nN=100"
)

design_comparison <- rbind(
  summary50,
  summary100
)

print(design_comparison)

write.csv(
  design_comparison,
  file.path(data_dir, "Stage2_independent_LHS_parameter_summary.csv"),
  row.names = FALSE
)


############################################################
## Plot comparison of independent Stage 2 LHS designs
############################################################

## Function to plot the two independent designs
plot_design_comparison <- function(parameter, label,
                                   ylim) {
  
  x <- design_comparison[
    design_comparison$parameter == parameter,
    ,
    drop = FALSE
  ]
  
  ## Put primary design first
  x <- x[
    match(
      c("Primary\nN=100",
        "Independent\nN=50"),
      x$design
    ),
  ]
  
  ## Move points away from panel borders
  xpos <- c(1.25, 1.75)
  
  ## Plot medians
  plot(
    xpos,
    x$median,
    pch = 16,
    ylim = ylim,
    xlim = c(1, 2),
    xaxt = "n",
    xlab = "LHS design",
    ylab = label,
    cex.axis = 1.3,
    cex.lab = 1.4
  )
  
  ## Add 2.5th-97.5th percentile intervals
  arrows(
    xpos,
    x$lower,
    xpos,
    x$upper,
    angle = 90,
    code = 3,
    length = 0.05
  )
  
  ## Design labels
  axis(
    1,
    at = xpos,
    labels = c(
      "Primary\nN=100",
      "Independent\nN=50"
    ),
    cex.axis = 1.1
  )
}


############################################################
## Two-row figure
############################################################

png(
  file.path(figure_dir, "Stage2_independent_LHS_parameter_stability.png"),
  width = 3200,
  height = 2200,
  res = 300
)

## 4 columns x 2 rows
##
## Top:
## c_GRGR | c0_SORSOR | c0_GRSOR | blank
##
## Bottom:
## s_SORSOR | s_SORGR | f | y

layout(
  matrix(
    c(
      1, 2, 3, 8,
      4, 5, 6, 7
    ),
    nrow = 2,
    byrow = TRUE
  ),
  widths = c(1, 1, 1, 1)
)

par(
  mar = c(6, 5, 2, 1)
)


############################################################
## Top row: inheritance parameters
############################################################

plot_design_comparison(
  "c_GRGR",
  expression(c[GR*",GR"]),
  ylim = c(0.9, 1)
)

plot_design_comparison(
  "c0_SORSOR",
  expression(c[SOR*",SOR"]^0),
  ylim = c(0.01, 1)
)

plot_design_comparison(
  "c0_GRSOR",
  expression(c[GR*",SOR"]^0),
  ylim = c(0.01, 1)
)


############################################################
## Bottom row: MDA-effect parameters
############################################################

plot_design_comparison(
  "s_SORSOR",
  expression(s[SOR*",SOR"]),
  ylim = c(0, 20)
)

plot_design_comparison(
  "s_SORGR",
  expression(s[SOR*",GR"]),
  ylim = c(0, 20)
)

plot_design_comparison(
  "f",
  expression(f),
  ylim = c(1, 10)
)

plot_design_comparison(
  "y",
  expression(y),
  ylim = c(0, 1)
)


############################################################
## Blank top-right panel
############################################################

par(mar = c(0, 0, 0, 0))
plot.new()


############################################################
## Finish figure
############################################################

dev.off()

layout(1)
par(mar = c(5.1, 4.1, 4.1, 2.1))

############################################################
## Acceptance-rate stability under primary and sensitivity
## calibration criteria
############################################################

## Number of MDA-effect LHS samples
Ns <- c(20, 40, 60, 80, 100)

############################################################
## Load no-selection Stage 2 simulations
############################################################

load(
  "data_files/fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu_100_s0.Rdata"
)

postMDA_s0 <- ll

############################################################
## Function to calculate acceptance rates for a threshold-
## based scenario using the primary Stage 2 simulations
############################################################

stage2_acceptance_scenario <- function(N,
                                       prev_lo,
                                       prev_hi,
                                       sor_lo = NULL,
                                       sor_hi = NULL) {

  n_evaluated <- 0
  n_calibration <- 0
  n_retained <- 0

  for (n in 1:num_samples) {

    initial_SOR <- preMDA[[n]][1]
    c0_SORSOR <- preMDA[[n]][7]
    c0_GRSOR  <- preMDA[[n]][8]

    for (i in 1:N) {

      final_SOR  <- postMDA[[n]][i, 1]
      final_prev <- postMDA[[n]][i, 5]
      delta_SOR  <- final_SOR - initial_SOR

      ## Belief 2
      belief2 <- final_prev > prev_lo & final_prev < prev_hi

      ## Belief 3 / sensitivity criterion
      if (!is.null(sor_lo) && !is.null(sor_hi)) {
        belief3 <- delta_SOR > sor_lo & delta_SOR < sor_hi
      } else if (!is.null(sor_lo)) {
        belief3 <- delta_SOR > sor_lo
      } else if (!is.null(sor_hi)) {
        belief3 <- delta_SOR < sor_hi
      } else {
        belief3 <- TRUE
      }

      calibration <- belief2 & belief3
      constraints <- c0_SORSOR < 0.5 & c0_GRSOR > 0.5

      n_evaluated  <- n_evaluated + 1
      n_calibration <- n_calibration + calibration
      n_retained   <- n_retained + (calibration & constraints)
    }
  }

  data.frame(
    N = N,
    n_evaluated = n_evaluated,
    acceptance = n_calibration / n_evaluated,
    final_retention = n_retained / n_evaluated
  )
}


############################################################
## Function for no offspring-phenotype selection scenario
## Uses separate s = 0 simulations and does NOT impose
## Belief 3; acceptance is based on Belief 2 only.
############################################################

stage2_acceptance_no_selection <- function(N,
                                           prev_lo = 0.10,
                                           prev_hi = 0.12) {

  n_evaluated <- 0
  n_calibration <- 0
  n_retained <- 0

  for (n in 1:num_samples) {

    c0_SORSOR <- preMDA[[n]][7]
    c0_GRSOR  <- preMDA[[n]][8]

    for (i in 1:N) {

      final_prev <- postMDA_s0[[n]][i, 5]
      belief2 <- final_prev > prev_lo & final_prev < prev_hi
      constraints <- c0_SORSOR < 0.5 & c0_GRSOR > 0.5

      n_evaluated  <- n_evaluated + 1
      n_calibration <- n_calibration + belief2
      n_retained   <- n_retained + (belief2 & constraints)
    }
  }

  data.frame(
    N = N,
    n_evaluated = n_evaluated,
    acceptance = n_calibration / n_evaluated,
    final_retention = n_retained / n_evaluated
  )
}


############################################################
## Run primary and four sensitivity scenarios
############################################################

## Primary: prevalence 0.10-0.12; SOR increase >0.10
primary <- do.call(
  rbind,
  lapply(
    Ns,
    stage2_acceptance_scenario,
    prev_lo = 0.10,
    prev_hi = 0.12,
    sor_lo = 0.10
  )
)

## Moderate increase: prevalence 0.10-0.12; SOR increase 0.05-0.10
moderate <- do.call(
  rbind,
  lapply(
    Ns,
    stage2_acceptance_scenario,
    prev_lo = 0.10,
    prev_hi = 0.12,
    sor_lo = 0.05,
    sor_hi = 0.10
  )
)

## Weak increase: prevalence 0.10-0.12; SOR increase <0.05
low <- do.call(
  rbind,
  lapply(
    Ns,
    stage2_acceptance_scenario,
    prev_lo = 0.10,
    prev_hi = 0.12,
    sor_hi = 0.05
  )
)

## Wider prevalence: prevalence 0.08-0.14; SOR increase >0.10
wide <- do.call(
  rbind,
  lapply(
    Ns,
    stage2_acceptance_scenario,
    prev_lo = 0.08,
    prev_hi = 0.14,
    sor_lo = 0.10
  )
)

## No offspring-phenotype selection: s = 0; Belief 2 only
no_selection <- do.call(
  rbind,
  lapply(
    Ns,
    stage2_acceptance_no_selection,
    prev_lo = 0.10,
    prev_hi = 0.12
  )
)


############################################################
## Primary acceptance stability as a standalone figure
############################################################

yr_primary <- range(
  c(primary$acceptance, primary$final_retention),
  na.rm = TRUE
)
ypad_primary <- 0.05 * diff(yr_primary)

png(
  file.path(figure_dir, "Stage2_acceptance_stability_primary.png"),
  width = 1800,
  height = 1500,
  res = 300
)

par(mar = c(5, 5, 2, 1))

plot(
  primary$N,
  primary$acceptance,
  type = "b",
  ylim = c(0,0.08
    #max(0, yr_primary[1] - ypad_primary),
    #yr_primary[2] + ypad_primary
  ),
  xlim = range(Ns),
  xlab = "Number of MDA-effect LHS samples",
  ylab = "Proportion of parameter combinations",
  cex.axis = 1.3,
  cex.lab = 1.4
)

lines(
  primary$N,
  primary$final_retention,
  type = "b",
  lty = 2
)

legend(
  "topright",
  legend = c("Satisfy calibration criteria", "Final retained"),
  lty = c(1, 2),
  pch = c(1, 1),
  bty = "n",
  cex = 1.1
)

dev.off()


############################################################
## Four-panel acceptance stability figure for sensitivity
## analyses only (primary omitted)
############################################################

## Common y-axis across all four sensitivity panels
yr_sensitivity <- range(
  c(
    moderate$acceptance,
    moderate$final_retention,
    low$acceptance,
    low$final_retention,
    no_selection$acceptance,
    no_selection$final_retention,
    wide$acceptance,
    wide$final_retention
  ),
  na.rm = TRUE
)

ypad_sensitivity <- 0.05 * diff(yr_sensitivity)

ylim_sensitivity <- c(0,
  #max(0, yr_sensitivity[1] - ypad_sensitivity),
  yr_sensitivity[2] + ypad_sensitivity
)

plot_acceptance_sensitivity <- function(x, title) {
  plot(
    x$N,
    x$acceptance,
    type = "b",
    ylim = ylim_sensitivity,
    xlim = range(Ns),
    xlab = "Number of MDA-effect LHS samples",
    ylab = "Proportion of parameter combinations",
    main = title,
    cex.axis = 1.2,
    cex.lab = 1.3,
    cex.main = 1.2
  )

  lines(
    x$N,
    x$final_retention,
    type = "b",
    lty = 2
  )
}

png(
  file.path(figure_dir, "Stage2_acceptance_stability_sensitivity_analyses.png"),
  width = 2600,
  height = 2200,
  res = 300
)

par(
  mfrow = c(2, 2),
  mar = c(5, 5, 3, 1)
)

plot_acceptance_sensitivity(
  moderate,
  "A.  Moderate increase"
)

legend(
  "topright",
  legend = c("Satisfy calibration criteria", "Final retained"),
  lty = c(1, 2),
  pch = c(1, 1),
  bty = "n",
  cex = 0.9
)

plot_acceptance_sensitivity(
  low,
  "B.  Weak increase"
)

plot_acceptance_sensitivity(
  no_selection,
  "C.  No offspring-phenotype selection"
)

plot_acceptance_sensitivity(
  wide,
  "D.  Wider prevalence"
)

dev.off()

par(mfrow = c(1, 1))
par(mar = c(5.1, 4.1, 4.1, 2.1))

############################################################
## Pairwise comparison of independent 100-point and
## 50-point Stage 2 LHS designs
############################################################

## acc100:
##   final retained combinations from primary N = 100 design
##
## acc50_independent:
##   final retained combinations from independent N = 50 design


############################################################
## Parameter pairs to plot
############################################################

pairs_to_plot <- list(
  c("s_SORSOR", "s_SORGR"),
  c("s_SORSOR", "f"),
  c("s_SORSOR", "y"),
  c("s_SORGR", "f"),
  c("s_SORGR", "y"),
  c("f", "y")
)


############################################################
## Mathematical labels
############################################################

labs <- list(
  s_SORSOR = expression(s[SOR*",SOR"]),
  s_SORGR  = expression(s[SOR*",GR"]),
  f        = expression(f),
  y        = expression(y)
)


############################################################
## Prior ranges
##
## Using fixed prior ranges makes the two designs directly
## comparable and gives all panels consistent axis limits.
############################################################

prior_limits <- list(
  s_SORSOR = c(0, 20),
  s_SORGR  = c(0, 20),
  f        = c(1, 10),
  y        = c(0, 1)
)


############################################################
## Save figure
############################################################

png(
  file.path(figure_dir, "Stage2_independent_LHS_pairwise.png"),
  width = 3600,
  height = 2000,
  res = 300
)

## Six pairwise plots + legend + blank panel
layout(
  matrix(
    c(
      1, 2, 3, 7,
      4, 5, 6, 8
    ),
    nrow = 2,
    byrow = TRUE
  ),
  widths = c(1, 1, 1, 0.8)
)

par(
  mar = c(4.5, 5.5, 2, 1)
)


############################################################
## Draw pairwise plots
############################################################

for (z in pairs_to_plot) {
  
  xv <- z[1]
  yv <- z[2]
  
  ## Primary 100-point design
  plot(
    acc100[[xv]],
    acc100[[yv]],
    pch = 1,
    xlim = prior_limits[[xv]],
    ylim = prior_limits[[yv]],
    xlab = labs[[xv]],
    ylab = labs[[yv]],
    cex.axis = 1.2,
    cex.lab = 1.3,
    mgp = c(2.8, 0.8, 0)
  )
  
  ## Independent 50-point design
  points(
    acc50_independent[[xv]],
    acc50_independent[[yv]],
    pch = 16
  )
}


############################################################
## Legend
############################################################

par(mar = c(0, 0, 0, 0))

plot.new()

legend(
  "center",
  legend = c(
    "Primary (N=100)",
    "Independent (N=50)"
  ),
  pch = c(1, 16),
  bty = "n",
  cex = 1.2
)


############################################################
## Blank final panel
############################################################

plot.new()


############################################################
## Finish
############################################################

dev.off()

layout(1)

par(
  mar = c(5.1, 4.1, 4.1, 2.1)
)


############################################################
## Compare retained parameter distributions across
## primary and sensitivity analyses
############################################################

############################################################
## 1. Construct the complete Stage 2 parameter/output dataset
############################################################

all_stage2 <- vector("list", num_samples * 100)
k <- 1

for (n in 1:num_samples) {
  
  ## Pre-MDA SOR proportion and Stage 1 inheritance parameters
  initial_SOR <- preMDA[[n]][1]
  
  c_GRGR    <- preMDA[[n]][6]
  c0_SORSOR <- preMDA[[n]][7]
  c0_GRSOR  <- preMDA[[n]][8]
  
  for (i in 1:100) {
    
    ## Stage 2 model outputs
    final_SOR  <- postMDA[[n]][i, 1]
    final_prev <- postMDA[[n]][i, 5]
    
    ## Change in SOR worm proportion
    delta_SOR <- final_SOR - initial_SOR
    
    all_stage2[[k]] <- data.frame(
      
      ## Stage 1 parameters
      c_GRGR    = c_GRGR,
      c0_SORSOR = c0_SORSOR,
      c0_GRSOR  = c0_GRSOR,
      
      ## Stage 2 parameters
      s_SORSOR = postMDA[[n]][i, 6],
      s_SORGR  = postMDA[[n]][i, 7],
      f        = postMDA[[n]][i, 8],
      y        = postMDA[[n]][i, 9],
      
      ## Calibration outputs
      final_prev = final_prev,
      delta_SOR  = delta_SOR
    )
    
    k <- k + 1
  }
}

all_stage2 <- do.call(rbind, all_stage2)


############################################################
## 2. Common final parameter constraints
############################################################

constraints <-
  all_stage2$c0_SORSOR < 0.5 &
  all_stage2$c0_GRSOR > 0.5


############################################################
## 3. Load no-selection sensitivity analysis
############################################################

## In this scenario the two offspring-phenotype selection coefficients
## are fixed at zero. The saved Stage 2 output retains the same
## nine-column structure as the primary simulations:
## columns 6-7 = s_SORSOR, s_SORGR (= 0); columns 8-9 = f, y.

all_stage2_s0 <- vector("list", num_samples * 100)
k <- 1

for (n in 1:num_samples) {
  initial_SOR <- preMDA[[n]][1]
  c_GRGR      <- preMDA[[n]][6]
  c0_SORSOR   <- preMDA[[n]][7]
  c0_GRSOR    <- preMDA[[n]][8]

  for (i in 1:100) {
    final_SOR  <- postMDA_s0[[n]][i, 1]
    final_prev <- postMDA_s0[[n]][i, 5]

    all_stage2_s0[[k]] <- data.frame(
      c_GRGR = c_GRGR,
      c0_SORSOR = c0_SORSOR,
      c0_GRSOR = c0_GRSOR,
      s_SORSOR = 0,
      s_SORGR = 0,
      f = postMDA_s0[[n]][i, 8],
      y = postMDA_s0[[n]][i, 9],
      final_prev = final_prev,
      delta_SOR = final_SOR - initial_SOR
    )
    k <- k + 1
  }
}

all_stage2_s0 <- do.call(rbind, all_stage2_s0)
constraints_s0 <-
  all_stage2_s0$c0_SORSOR < 0.5 &
  all_stage2_s0$c0_GRSOR > 0.5

## No Belief 3 threshold is imposed here: selection through the
## offspring-phenotype functions is removed mechanistically.
keep_no_selection <-
  all_stage2_s0$final_prev > 0.10 &
  all_stage2_s0$final_prev < 0.12 &
  #all_stage2_s0$delta_SOR >0.1 &
  constraints_s0

acc_no_selection <- all_stage2_s0[keep_no_selection, , drop = FALSE]
message("No-selection combinations retained: ", nrow(acc_no_selection))

############################################################
## 4. Define the four threshold-based calibration scenarios
############################################################

## Primary:
## prevalence 10-12%, increase in SOR > 0.10
keep_primary <-
  all_stage2$final_prev > 0.10 &
  all_stage2$final_prev < 0.12 &
  all_stage2$delta_SOR > 0.10 &
  constraints


## Moderate increase:
## prevalence 10-12%, increase in SOR 0.05-0.10
keep_moderate <-
  all_stage2$final_prev > 0.10 &
  all_stage2$final_prev < 0.12 &
  all_stage2$delta_SOR > 0.05 &
  all_stage2$delta_SOR < 0.10 &
  constraints


## Low increase:
## prevalence 10-12%, increase in SOR < 0.05
keep_low <-
  all_stage2$final_prev > 0.10 &
  all_stage2$final_prev < 0.12 &
  all_stage2$delta_SOR < 0.05 &
  constraints


## Wider prevalence:
## prevalence 8-14%, increase in SOR > 0.10
keep_wide <-
  all_stage2$final_prev > 0.08 &
  all_stage2$final_prev < 0.14 &
  all_stage2$delta_SOR > 0.10 &
  constraints


############################################################
## 4. Extract retained parameter combinations
############################################################

acc_primary  <- all_stage2[keep_primary, , drop = FALSE]
acc_moderate <- all_stage2[keep_moderate, , drop = FALSE]
acc_low      <- all_stage2[keep_low, , drop = FALSE]
acc_wide     <- all_stage2[keep_wide, , drop = FALSE]

## Check numbers retained
c(
  Primary = nrow(acc_primary),
  Moderate = nrow(acc_moderate),
  Low = nrow(acc_low),
  No_selection = nrow(acc_no_selection),
  Wider_prevalence = nrow(acc_wide)
)

############################################################
## 5. Summarise retained parameter distributions
############################################################

parameters <- c(
  "c_GRGR",
  "c0_SORSOR",
  "c0_GRSOR",
  "s_SORSOR",
  "s_SORGR",
  "f",
  "y"
)


summarise_scenario <- function(x, scenario) {
  
  do.call(
    rbind,
    lapply(parameters, function(p) {
      
      data.frame(
        scenario = scenario,
        parameter = p,
        
        lower = quantile(
          x[[p]],
          0.025,
          na.rm = TRUE
        ),
        
        median = median(
          x[[p]],
          na.rm = TRUE
        ),
        
        upper = quantile(
          x[[p]],
          0.975,
          na.rm = TRUE
        )
      )
    })
  )
}


sensitivity_parameter_summary <- rbind(
  
  summarise_scenario(
    acc_primary,
    "Primary"
  ),
  
  summarise_scenario(
    acc_moderate,
    "Moderate"
  ),
  
  summarise_scenario(
    acc_low,
    "Low"
  ),

  summarise_scenario(
    acc_no_selection,
    "No selection"
  ),
  
  summarise_scenario(
    acc_wide,
    "Wider prevalence"
  )
)


print(sensitivity_parameter_summary)

write.csv(
  sensitivity_parameter_summary,
  file.path(data_dir, "Stage2_sensitivity_parameter_summary.csv"),
  row.names = FALSE
)

############################################################
## 6. Plotting function
############################################################

plot_sensitivity_comparison <- function(parameter,
                                        label,
                                        ylim,
                                        show_xlab = FALSE) {
  
  x <- sensitivity_parameter_summary[
    sensitivity_parameter_summary$parameter == parameter,
    ,
    drop = FALSE
  ]
  
  scenario_order <- c(
    "Primary",
    "Moderate",
    "Low",
    "No selection",
    "Wider prevalence"
  )
  
  x <- x[
    match(scenario_order, x$scenario),
    ,
    drop = FALSE
  ]
  
  xpos <- 1:5
  
  plot(
    xpos,
    x$median,
    pch = 16,
    ylim = ylim,
    xlim = c(0.5, 5.5),
    xaxt = "n",
    xlab = "",
    ylab = label,
    cex.axis = 1.2,
    cex.lab = 1.3,
    mgp = c(2.8, 0.8, 0)
  )
  
  arrows(
    xpos,
    x$lower,
    xpos,
    x$upper,
    angle = 90,
    code = 3,
    length = 0.04
  )
  
  ## Tick marks
  axis(
    1,
    at = xpos,
    labels = FALSE,
    tck = -0.02
  )
  
  ## Move labels below the tick marks
  label_y <- par("usr")[3] -
    0.045 * diff(par("usr")[3:4])
  
  text(
    x = xpos,
    y = label_y,
    labels = c(
      "Primary",
      "Moderate",
      "Weak",
      "No selection",
      "Wider prevalence"
    ),
    srt = 45,
    adj = 1,
    xpd = TRUE,
    cex = 0.85
  )
  
  if (show_xlab) {
    mtext(
      "Calibration analysis",
      side = 1,
      line = 5.3,
      cex = 1.3
    )
  }
}

############################################################
## 7. Seven-panel sensitivity comparison
############################################################

png(
  file.path(figure_dir, "Stage2_sensitivity_parameter_distributions.png"),
  width = 3800,
  height = 1800,
  res = 300
)

## Top:
## c_GRGR | c0_SORSOR | c0_GRSOR | key
##
## Bottom:
## s_SORSOR | s_SORGR | f | y

layout(
  matrix(
    1:8,
    nrow = 2,
    byrow = TRUE
  ),
  widths = c(1, 1, 1, 1)
)


############################################################
## Top row: inheritance parameters
############################################################

par(mar = c(5.5, 6, 1, 1))

plot_sensitivity_comparison(
  "c_GRGR",
  expression(c[GR*",GR"]),
  c(0.9, 1)
)

plot_sensitivity_comparison(
  "c0_SORSOR",
  expression(c[SOR*",SOR"]^0),
  c(0.01, 1)
)

plot_sensitivity_comparison(
  "c0_GRSOR",
  expression(c[GR*",SOR"]^0),
  c(0.01, 1)
)


############################################################
## Top-right: explanation of intervals
############################################################

par(mar = c(0, 0, 0, 0))
plot.new()

legend(
  "center",
  legend = c(
    "Median",
    "2.5th/97.5th percentile"
  ),
  lty = c(NA, 1),
  pch = c(16, NA),
  bty = "n",
  cex = 1.1
)


############################################################
## Bottom row: MDA-effect parameters
############################################################

par(mar = c(7.5, 6, 1, 1))

plot_sensitivity_comparison(
  "s_SORSOR",
  expression(s[SOR*",SOR"]),
  c(0, 20),
  show_xlab = TRUE
)

plot_sensitivity_comparison(
  "s_SORGR",
  expression(s[SOR*",GR"]),
  c(0, 20),
  show_xlab = TRUE
)

plot_sensitivity_comparison(
  "f",
  expression(f),
  c(1, 10),
  show_xlab = TRUE
)

plot_sensitivity_comparison(
  "y",
  expression(y),
  c(0, 1),
  show_xlab = TRUE
)


############################################################
## Finish
############################################################

dev.off()

layout(1)

par(
  mar = c(5.1, 4.1, 4.1, 2.1)
)

############################################################
## Acceptance proportions for sensitivity analyses
############################################################

## For the four threshold-based analyses, counts come from the
## common 100-point Stage 2 design. For the no-selection analysis,
## Belief 3 is not imposed; acceptance is based on Belief 2 and the
## common inheritance-parameter constraints.

get_n100 <- function(x, field) x[[field]][x$N == 100]

no_selection_n_eval <- nrow(all_stage2_s0)
no_selection_n_belief2 <- sum(
  all_stage2_s0$final_prev > 0.10 &
  all_stage2_s0$final_prev < 0.12
)
no_selection_n_retained <- nrow(acc_no_selection)

sensitivity_acceptance_summary <- data.frame(
  analysis = c(
    "Primary",
    "Moderate increase",
    "Low increase",
    "No selection",
    "Wider prevalence"
  ),
  combinations_evaluated = c(
    get_n100(primary, "n_evaluated"),
    get_n100(moderate, "n_evaluated"),
    get_n100(low, "n_evaluated"),
    no_selection_n_eval,
    get_n100(wide, "n_evaluated")
  ),
  n_satisfying_calibration = c(
    round(get_n100(primary, "acceptance") * get_n100(primary, "n_evaluated")),
    round(get_n100(moderate, "acceptance") * get_n100(moderate, "n_evaluated")),
    round(get_n100(low, "acceptance") * get_n100(low, "n_evaluated")),
    no_selection_n_belief2,
    round(get_n100(wide, "acceptance") * get_n100(wide, "n_evaluated"))
  ),
  n_final_retained = c(
    round(get_n100(primary, "final_retention") * get_n100(primary, "n_evaluated")),
    round(get_n100(moderate, "final_retention") * get_n100(moderate, "n_evaluated")),
    round(get_n100(low, "final_retention") * get_n100(low, "n_evaluated")),
    no_selection_n_retained,
    round(get_n100(wide, "final_retention") * get_n100(wide, "n_evaluated"))
  )
)

sensitivity_acceptance_summary$acceptance_proportion <-
  sensitivity_acceptance_summary$n_satisfying_calibration /
  sensitivity_acceptance_summary$combinations_evaluated

sensitivity_acceptance_summary$final_retention_proportion <-
  sensitivity_acceptance_summary$n_final_retained /
  sensitivity_acceptance_summary$combinations_evaluated

sensitivity_acceptance_summary$acceptance_percent <-
  100 * sensitivity_acceptance_summary$acceptance_proportion

sensitivity_acceptance_summary$final_retention_percent <-
  100 * sensitivity_acceptance_summary$final_retention_proportion

print(sensitivity_acceptance_summary)

write.csv(
  sensitivity_acceptance_summary,
  file.path(data_dir, "Stage2_sensitivity_acceptance_summary.csv"),
  row.names = FALSE
)


###############################################################################
## REPOSITORY DATA EXPORT
## Run after all Stage 2 analyses above.
###############################################################################
dir.create("data_processed", showWarnings = FALSE, recursive = TRUE)

## Accepted primary combinations labelled by the part of the 100-point
## MDA-effect LHS design from which they arose.
stage2_primary_pairwise <- acc100_stability
stage2_primary_pairwise$sampling_block <- ifelse(
  stage2_primary_pairwise$lhs_sample <= 80,
  "First 80",
  "Final 20"
)

## Accepted combinations for the independent-design comparison.
## Retain only the parameter columns common to both designs.

pairwise_parameters <- c(
  "c_GRGR",
  "c0_SORSOR",
  "c0_GRSOR",
  "s_SORSOR",
  "s_SORGR",
  "f",
  "y"
)

stage2_independent_pairwise <- rbind(
  
  data.frame(
    acc100[
      ,
      pairwise_parameters,
      drop = FALSE
    ],
    design = "Primary N=100"
  ),
  
  data.frame(
    acc50_independent[
      ,
      pairwise_parameters,
      drop = FALSE
    ],
    design = "Independent N=50"
  )
)

## Acceptance-stability objects for the four sensitivity analyses.
stage2_sensitivity_stability <- list(
  Moderate = moderate,
  Weak = low,
  No_selection = no_selection,
  Wider_prevalence = wide
)

save(
  acceptance_summary,
  parameter_stability,
  stage2_primary_pairwise,
  spearman_cor,
  identifiability_summary,
  design_comparison,
  stage2_independent_pairwise,
  sensitivity_parameter_summary,
  sensitivity_acceptance_summary,
  stage2_sensitivity_stability,
  file = "data_processed/stage2_stability_repo_data.Rdata"
)

message("Saved data_processed/stage2_stability_repo_data.Rdata")
