###############################################################################
## Asubende projections: comparison across calibration sensitivity analyses
###############################################################################
##
## PURPOSE
##
## This script compares the Asubende calibration outcomes and subsequent
## biannual MDAi projections across the primary calibration and four
## sensitivity analyses.
##
## The analyses are:
##
##   1. Primary calibration
##        final mf prevalence = 0.10-0.12
##        increase in SOR worm proportion > 0.10
##
##   2. Moderate-increase sensitivity analysis
##        final mf prevalence = 0.10-0.12
##        increase in SOR worm proportion = 0.05-0.10
##
##   3. Weak-increase sensitivity analysis
##        final mf prevalence = 0.10-0.12
##        increase in SOR worm proportion < 0.05
##
##   4. No offspring-phenotype selection sensitivity analysis
##        s_SORSOR = s_SORGR = 0
##        final mf prevalence = 0.10-0.12
##        no Belief 3 criterion is imposed
##
##   5. Wider-prevalence sensitivity analysis
##        final mf prevalence = 0.08-0.14
##        increase in SOR worm proportion > 0.10
##
## All analyses use the common final inheritance-parameter constraints:
##        c0_SORSOR < 0.5
##        c0_GRSOR  > 0.5
##
## ANALYSES
##
## For each calibration analysis, the script loads:
##
##   - the final retained parameter combinations following the 20-year
##     annual MDAi calibration period; and
##   - the corresponding simulations of the subsequent 18-year biannual
##     MDAi period in Asubende.
##
## For every retained parameter combination, the following quantities are
## extracted:
##
##   Annual MDAi calibration period:
##     - final mf prevalence;
##     - increase in SOR worm proportion;
##     - final SOR worm proportion.
##
##   Subsequent biannual MDAi projection:
##     - final mf prevalence;
##     - further increase in SOR worm proportion;
##     - final SOR worm proportion.
##
##   The pre-MDA SOR worm proportion is also reconstructed for each
##   parameter combination.
##
## For each outcome and calibration analysis, the script calculates:
##     - 2.5th percentile;
##     - 25th percentile;
##     - median;
##     - 75th percentile;
##     - 97.5th percentile;
##     - number of retained parameter combinations.
##
## FIGURE
##
## The script produces a six-panel comparison across the five calibration
## analyses:
##
##   A. Final mf prevalence at the end of annual MDAi calibration
##   B. Increase in SOR worm proportion during annual MDAi calibration
##   C. Final SOR worm proportion at the end of annual MDAi calibration
##   D. Final mf prevalence at the end of the subsequent biannual
##      MDAi projection
##   E. Further increase in SOR worm proportion during the subsequent
##      biannual MDAi projection
##   F. Final SOR worm proportion at the end of the subsequent biannual
##      MDAi projection
##
## Points show medians and vertical intervals show the 2.5th-97.5th
## percentiles across retained parameter combinations.
##
## OUTPUT
##
## Figure output is saved in:
##   stability_figures/
##
##   Asubende_sensitivity_outcomes.png
##     Six-panel comparison of Asubende outcomes across the primary
##     calibration and four sensitivity analyses.
##
## Data outputs are saved in:
##   stability_data/
##
##   Asubende_sensitivity_all_outcomes.csv
##     Individual outcomes for every retained parameter combination under
##     each calibration analysis.
##
##   Asubende_sensitivity_summary.csv
##     Quantile summaries of each Asubende outcome under each calibration
##     analysis.
##
## INTERPRETATION
##
## The annual-MDAi outcomes are partly constrained by the corresponding
## calibration criteria and should not be interpreted as independent model
## predictions.
##
## Outcomes during the subsequent biannual MDAi period were not used as
## calibration criteria and are model projections conditional on the
## corresponding retained parameter combinations.
##
## Comparison across the five analyses is used to assess which Asubende
## projections are quantitatively sensitive to the calibration assumptions
## and which qualitative conclusions are maintained when the assumed
## magnitude or mechanism of SOR selection is changed.
##
## Reported intervals are the 2.5th-97.5th percentiles of outcomes across
## retained parameter combinations and should not be interpreted as Bayesian
## credible intervals.
##
###############################################################################

## Assumes working directory is the project root.
figure_dir <- "stability_figures"
data_dir   <- "stability_data"

dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(data_dir, showWarnings = FALSE, recursive = TRUE)

############################################################
## 1. Load common model objects
############################################################

load("data_files/burnin_157_ABR_65000_UU.Rdata")

source("equations.R")

MDAi <- 1
ABR  <- 65000

source("fixed-parameters.R")

############################################################
## 2. Files for the five calibration analyses
############################################################

scenario_files <- data.frame(
  
  scenario = c(
    "Primary",
    "Moderate",
    "Weak",
    "No selection",
    "Wider prevalence"
  ),
  
  biannual_file = c(
    "data_files/Final_posterior_burnins_plusbiannualMDAs_cumu_100.Rdata",
    "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA_100.Rdata",
    "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA2_100.Rdata",
    "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA5_100.Rdata",
    "data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA3_100.Rdata"
  ),
  
  posterior_file = c(
    "data_files/posteriors_ABR_65000_UU_cumu_100.Rdata",
    "data_files/posteriors_ABR_65000_UU_cumuSA_100.Rdata",
    "data_files/posteriors_ABR_65000_UU_cumuSA2_100.Rdata",
    "data_files/posteriors_ABR_65000_UU_cumuSA5_100.Rdata",
    "data_files/posteriors_ABR_65000_UU_cumuSA3_100.Rdata"
  ),
  
  stringsAsFactors = FALSE
)

############################################################
## 3. Function to extract Asubende results
############################################################

extract_asubende_results <- function(biannual_file,
                                     posterior_file,
                                     scenario) {
  
  ##########################################################
  ## Load biannual simulation output
  ##########################################################
  
  e1 <- new.env()
  load(biannual_file, envir = e1)
  
  if (!exists("OutputMM", envir = e1)) {
    stop(
      "OutputMM not found in: ",
      biannual_file
    )
  }
  
  OutputMM_local <- e1$OutputMM
  
  ##########################################################
  ## Load corresponding accepted parameter combinations
  ##########################################################
  
  e2 <- new.env()
  load(posterior_file, envir = e2)
  
  if (!exists("posteriors", envir = e2)) {
    stop(
      "posteriors not found in: ",
      posterior_file
    )
  }
  
  post <- e2$posteriors
  
  ##########################################################
  ## Apply final inheritance-parameter constraints
  ##
  ## columns:
  ## 3 = c0_SORSOR
  ## 4 = c0_GRSOR
  ##########################################################
  
  keep <-
    post[, 3] < 0.5 &
    post[, 4] > 0.5
  
  post <- post[keep, , drop = FALSE]
  
  ##########################################################
  ## Check alignment between posteriors and OutputMM
  ##########################################################
  
  if (length(OutputMM_local) != nrow(post)) {
    
    stop(
      scenario,
      ": number of OutputMM simulations (",
      length(OutputMM_local),
      ") does not match number of retained posterior rows (",
      nrow(post),
      ")."
    )
  }
  
  ##########################################################
  ## Extract final biannual-MDA outcomes
  ##########################################################
  
  num_sims <- nrow(post)
  
  out <- vector("list", num_sims)
  
  for (n in seq_len(num_sims)) {
    
    yout <- OutputMM_local[[n]]
    
    ## SOR worm proportion over time
    sor_worm <-
      1 -
      (yout[, 3] + yout[, 4]) /
      (
        yout[, 3] + yout[, 4] +
          yout[, 7] + yout[, 8]
      )
    
    ## Total mf prevalence over time
    prev <-
      prevalence(
        yout[, 5] + yout[, 9],
        my_k0,
        my_k1
      )
    
    ## End of biannual period
    final_sor_biannual <-
      tail(sor_worm, 1)
    
    final_prev_biannual <-
      tail(prev, 1)
    
    ########################################################
    ## Values at end of annual-MDA calibration period
    ##
    ## Existing code:
    ## column 9  = final annual mf prevalence
    ## column 10 = increase in SOR during annual MDA
    ## column 11 = final SOR proportion after annual MDA
    ########################################################
    
    final_prev_annual <- post[n, 9]
    
    increase_sor_annual <- post[n, 10]
    
    final_sor_annual <- post[n, 11]
    
    ########################################################
    ## Initial pre-MDA SOR proportion
    ########################################################
    
    initial_sor <-
      final_sor_annual -
      increase_sor_annual
    
    ########################################################
    ## Further change during biannual period
    ########################################################
    
    further_increase_sor <-
      final_sor_biannual -
      final_sor_annual
    
    ########################################################
    ## Store
    ########################################################
    
    out[[n]] <- data.frame(
      
      scenario = scenario,
      
      SampleID = post[n, 1],
      
      final_prev_annual =
        final_prev_annual,
      
      increase_sor_annual =
        increase_sor_annual,
      
      final_sor_annual =
        final_sor_annual,
      
      final_prev_biannual =
        final_prev_biannual,
      
      further_increase_sor_biannual =
        further_increase_sor,
      
      final_sor_biannual =
        final_sor_biannual,
      
      initial_sor =
        initial_sor
    )
  }
  
  do.call(rbind, out)
}

############################################################
## 4. Extract all five calibration analyses
############################################################

asubende_list <- vector(
  "list",
  nrow(scenario_files)
)

for (i in seq_len(nrow(scenario_files))) {
  
  message(
    "Processing ",
    scenario_files$scenario[i]
  )
  
  asubende_list[[i]] <-
    extract_asubende_results(
      
      biannual_file =
        scenario_files$biannual_file[i],
      
      posterior_file =
        scenario_files$posterior_file[i],
      
      scenario =
        scenario_files$scenario[i]
    )
}

asubende_all <-
  do.call(
    rbind,
    asubende_list
  )

## Check retained sample numbers
print(
  table(asubende_all$scenario)
)

## Save raw extracted outcomes
write.csv(
  asubende_all,
  file.path(
    data_dir,
    "Asubende_sensitivity_all_outcomes.csv"
  ),
  row.names = FALSE
)

############################################################
## 5. Calculate summary statistics
############################################################

outcomes <- c(
  
  "final_prev_annual",
  
  "increase_sor_annual",
  
  "final_sor_annual",
  
  "final_prev_biannual",
  
  "further_increase_sor_biannual",
  
  "final_sor_biannual",
  
  "initial_sor"
)

summarise_outcomes <- function(d) {
  
  output <- NULL
  
  for (s in unique(d$scenario)) {
    
    ds <- d[
      d$scenario == s,
      ,
      drop = FALSE
    ]
    
    for (v in outcomes) {
      
      q <- quantile(
        ds[[v]],
        probs = c(
          0.025,
          0.25,
          0.5,
          0.75,
          0.975
        ),
        na.rm = TRUE
      )
      
      output <- rbind(
        output,
        data.frame(
          scenario = s,
          outcome = v,
          
          q025 = q[1],
          q25  = q[2],
          median = q[3],
          q75  = q[4],
          q975 = q[5],
          
          n = sum(
            !is.na(ds[[v]])
          ),
          
          row.names = NULL
        )
      )
    }
  }
  
  output
}

asubende_summary <-
  summarise_outcomes(
    asubende_all
  )

write.csv(
  asubende_summary,
  file.path(
    data_dir,
    "Asubende_sensitivity_summary.csv"
  ),
  row.names = FALSE
)

print(asubende_summary)

############################################################
## 6. Comparative Asubende figure
############################################################

scenario_order <- c(
  "Primary",
  "Moderate",
  "Weak",
  "No selection",
  "Wider prevalence"
)

plot_asubende_outcome <- function(outcome,
                                  label,
                                  ylim = NULL,
                                  show_xlab = FALSE) {
  
  x <- asubende_summary[
    asubende_summary$outcome == outcome,
    ,
    drop = FALSE
  ]
  
  x <- x[
    match(
      scenario_order,
      x$scenario
    ),
    ,
    drop = FALSE
  ]
  
  xpos <- 1:5
  
  
  ## Automatic y range if not supplied
  if (is.null(ylim)) {
    
    yr <- range(
      c(
        x$q025,
        x$q975
      ),
      na.rm = TRUE
    )
    
    pad <- 0.05 * diff(yr)
    
    ylim <- c(
      max(0, yr[1] - pad),
      yr[2] + pad
    )
  }
  
  
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
    cex.lab = 1.3
  )
  
  arrows(
    xpos,
    x$q025,
    xpos,
    x$q975,
    angle = 90,
    code = 3,
    length = 0.04
  )
  
    axis(
    1,
    at = xpos,
    labels = FALSE,
    tck = -0.02
  )
  
  ## Move labels below ticks
  label_y <-
    par("usr")[3] -
    0.045 *
    diff(par("usr")[3:4])
  
  text(
    xpos,
    label_y,
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
      cex = 1.2
    )
  }
}

############################################################
## 7. Six-panel figure
############################################################

png(
  file.path(
    figure_dir,
    "Asubende_sensitivity_outcomes.png"
  ),
  width = 3300,
  height = 2200,
  res = 300
)

layout(
  matrix(
    1:6,
    nrow = 2,
    byrow = TRUE
  )
)

############################################################
## Top row: end of annual MDAi calibration period
############################################################

par(
  mar = c(7.5, 6, 2.5, 1)
)

plot_asubende_outcome(
  "final_prev_annual",
  "Final mf prevalence",
  ylim = c(0, 0.15)
)

mtext(
  "A.  End of annual MDAi calibration",
  side = 3,
  line = 0.5,
  cex = 1.1
)

plot_asubende_outcome(
  "increase_sor_annual",
  "Increase in SOR worm proportion",
  ylim = NULL
)

mtext(
  "B.  End of annual MDAi calibration",
  side = 3,
  line = 0.5,
  cex = 1.1
)

plot_asubende_outcome(
  "final_sor_annual",
  "Final SOR worm proportion",
  ylim = c(0, 0.6)
)

mtext(
  "C.  End of annual MDAi calibration",
  side = 3,
  line = 0.5,
  cex = 1.1
)

############################################################
## Bottom row: subsequent biannual MDAi projection period
############################################################

par(
  mar = c(8.5, 6, 2.5, 1)
)

plot_asubende_outcome(
  "final_prev_biannual",
  "Final mf prevalence",
  ylim = NULL,
  show_xlab = TRUE
)

mtext(
  "D.  End of biannual MDAi projection",
  side = 3,
  line = 0.5,
  cex = 1.1
)

plot_asubende_outcome(
  "further_increase_sor_biannual",
  "Further increase in SOR worm proportion",
  ylim = NULL,
  show_xlab = TRUE
)

mtext(
  "E.  End of biannual MDAi projection",
  side = 3,
  line = 0.5,
  cex = 1.1
)

plot_asubende_outcome(
  "final_sor_biannual",
  "Final SOR worm proportion",
  ylim = c(0, 0.6),
  show_xlab = TRUE
)

mtext(
  "F.  End of biannual MDAi projection",
  side = 3,
  line = 0.5,
  cex = 1.1
)

dev.off()

layout(1)

par(
  mar = c(5.1, 4.1, 4.1, 2.1)
)


###############################################################################
## REPOSITORY DATA EXPORT
## The CSV files written above are the compact repository inputs.
###############################################################################
message(
  "Repository sensitivity inputs written to stability_data/: ",
  "Asubende_sensitivity_all_outcomes.csv and ",
  "Asubende_sensitivity_summary.csv"
)
