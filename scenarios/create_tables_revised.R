###############################################################################
## Create scenario summary tables
###############################################################################
##
## Uses processed_all_scenarios.RData created by process_all_scenario_outputs.R.
##
## Outputs:
##   1. Primary scenario tables across 30%, 50%, and 70% baseline prevalence.
##   2. Comparative sensitivity tables at ~70% baseline prevalence for:
##        Primary, Moderate, Weak, No selection, Wider prevalence.
##
## Summaries are median and 2.5th-97.5th percentiles.
###############################################################################

library(dplyr)
library(tidyr)

load("processed_all_scenarios.RData")

if (!dir.exists("tables")) dir.create("tables")

qfun <- function(x) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(c(q025=NA, median=NA, q975=NA))
  c(q025=unname(quantile(x,0.025)),
    median=median(x),
    q975=unname(quantile(x,0.975)))
}

fmt_interval <- function(med, lo, hi, digits=2) {
  ifelse(
    is.na(med),
    NA_character_,
    paste0(
      formatC(med, digits=digits, format="f"),
      " (",
      formatC(lo, digits=digits, format="f"),
      "-",
      formatC(hi, digits=digits, format="f"),
      ")"
    )
  )
}


summarise_metric <- function(dat, metric, groups) {
  dat %>%
    group_by(across(all_of(groups))) %>%
    summarise(
      n = sum(is.finite(.data[[metric]])),
      q025 = qfun(.data[[metric]])["q025"],
      median = qfun(.data[[metric]])["median"],
      q975 = qfun(.data[[metric]])["q975"],
      .groups="drop"
    )
}

## Labels for original baseline settings
prev_labels <- c(`2942`="30%", `5350`="50%", `27000`="70%")

###############################################################################
## PRIMARY SCENARIO TABLES
###############################################################################

primary <- processed_all %>%
  filter(Calibration == "Primary") %>%
  mutate(
    Baseline_prevalence = unname(prev_labels[as.character(ABR)])
  )

## Final mf prevalence (%)
primary_final_prev <- summarise_metric(
  primary, "Final_Prevalence_percent",
  c("Baseline_prevalence","Drug","Model","Treatment_frequency")
)

###############################################################################
## Years to pOTTIS - primary scenarios
###############################################################################

## Distribution of time to pOTTIS among simulations that reached pOTTIS
primary_pOTTIS <- summarise_metric(
  primary,
  "pOTTIS",
  c(
    "Baseline_prevalence",
    "Drug",
    "Model",
    "Treatment_frequency"
  )
)

## Total number of simulations in each scenario
primary_pOTTIS_total <- primary %>%
  group_by(
    Baseline_prevalence,
    Drug,
    Model,
    Treatment_frequency
  ) %>%
  summarise(
    n_total = n(),
    .groups = "drop"
  )

## Add number and proportion reaching pOTTIS
primary_pOTTIS <- primary_pOTTIS %>%
  left_join(
    primary_pOTTIS_total,
    by = c(
      "Baseline_prevalence",
      "Drug",
      "Model",
      "Treatment_frequency"
    )
  ) %>%
  rename(
    n_reached = n
  ) %>%
  mutate(
    proportion_reached = n_reached / n_total,
    percent_reached = 100 * proportion_reached
  )

primary_pOTTIS_formatted <- primary_pOTTIS %>%
  mutate(
    `Median years (2.5th-97.5th percentile)` =
      fmt_interval(
        median,
        q025,
        q975,
        1
      ),
    
    `Reached pOTTIS, n/N (%)` =
      paste0(
        n_reached,
        "/",
        n_total,
        " (",
        round(percent_reached, 1),
        "%)"
      )
  ) %>%
  
  select(
    Baseline_prevalence,
    Drug,
    Model,
    Treatment_frequency,
    `Median years (2.5th-97.5th percentile)`,
    `Reached pOTTIS, n/N (%)`
  )

## Final SOR worm proportion after 20 years: VR only
primary_SOR <- summarise_metric(
  primary %>% filter(Model == "VR"),
  "Final_SOR_proportion_worms20yr",
  c("Baseline_prevalence","Drug","Treatment_frequency")
)

## Initial SOR worm proportion: VR only
primary_initial_SOR <- summarise_metric(
  primary %>% filter(Model == "VR", Drug == "Ivermectin"),
  "Initial_SOR_proportion_worms",
  c("Baseline_prevalence")
)

###############################################################################
## SENSITIVITY TABLES - 70% BASELINE, IVERMECTIN
###############################################################################

sens <- processed_all %>%
  filter(
    ABR == 27000,
    Drug == "Ivermectin",
    Calibration %in% c(
      "Primary","Moderate","Weak","No selection","Wider prevalence"
    )
  )

## Final mf prevalence (%): all calibrations, VR and GR-only
sensitivity_final_prev <- summarise_metric(
  sens, "Final_Prevalence_percent",
  c("Calibration","Model","Treatment_frequency")
)

## Years to pOTTIS: all calibrations, VR and GR-only
sensitivity_pOTTIS <- summarise_metric(
  sens, "pOTTIS",
  c("Calibration","Model","Treatment_frequency")
)

## Also record the total number of trajectories, so the proportion reaching
## pOTTIS can be reported rather than conditioning silently on those that do.
pottis_total <- sens %>%
  group_by(Calibration, Model, Treatment_frequency) %>%
  summarise(n_total=n(), .groups="drop")

sensitivity_pOTTIS <- sensitivity_pOTTIS %>%
  left_join(pottis_total,
            by=c("Calibration","Model","Treatment_frequency")) %>%
  rename(n_reached=n) %>%
  mutate(proportion_reached=n_reached/n_total)

## Final SOR worm proportion after 20 years: VR only
sensitivity_SOR <- summarise_metric(
  sens %>% filter(Model == "VR"),
  "Final_SOR_proportion_worms20yr",
  c("Calibration","Treatment_frequency")
)

## Initial SOR worm proportion: VR only
sensitivity_initial_SOR <- summarise_metric(
  sens %>% filter(Model == "VR"),
  "Initial_SOR_proportion_worms",
  c("Calibration")
)

###############################################################################
## Formatted manuscript-friendly versions
###############################################################################

sensitivity_final_prev_formatted <- sensitivity_final_prev %>%
  mutate(`Median (2.5th-97.5th percentile)` =
           fmt_interval(median,q025,q975,2)) %>%
  select(Calibration, Model, Treatment_frequency,
         `Median (2.5th-97.5th percentile)`)

sensitivity_pOTTIS_formatted <- sensitivity_pOTTIS %>%
  mutate(
    `Median years (2.5th-97.5th percentile)` =
      fmt_interval(median,q025,q975,1),
    `Reached pOTTIS` = paste0(n_reached,"/",n_total)
  ) %>%
  select(Calibration, Model, Treatment_frequency,
         `Median years (2.5th-97.5th percentile)`,
         `Reached pOTTIS`, proportion_reached)

sensitivity_SOR_formatted <- sensitivity_SOR %>%
  mutate(`Median (2.5th-97.5th percentile)` =
           fmt_interval(median,q025,q975,3)) %>%
  select(Calibration, Treatment_frequency,
         `Median (2.5th-97.5th percentile)`)

###############################################################################
## Save
###############################################################################

save(
  primary_final_prev,
  primary_pOTTIS,
  primary_SOR,
  primary_initial_SOR,
  sensitivity_final_prev,
  sensitivity_pOTTIS,
  sensitivity_SOR,
  sensitivity_initial_SOR,
  sensitivity_final_prev_formatted,
  sensitivity_pOTTIS_formatted,
  sensitivity_SOR_formatted,
  file="tables/Tables_scenarios_updated.RData"
)

write.csv(primary_final_prev,
          "tables/primary_final_prevalence.csv", row.names=FALSE)
write.csv(primary_pOTTIS,
          "tables/primary_pOTTIS.csv", row.names=FALSE)
write.csv(primary_pOTTIS_formatted,"tables/primary_pOTTIS_formatted.csv",
  row.names = FALSE
)
write.csv(primary_SOR,
          "tables/primary_final_SOR.csv", row.names=FALSE)

write.csv(sensitivity_final_prev_formatted,
          "tables/sensitivity_final_prevalence.csv", row.names=FALSE)
write.csv(sensitivity_pOTTIS_formatted,
          "tables/sensitivity_pOTTIS.csv", row.names=FALSE)
write.csv(sensitivity_SOR_formatted,
          "tables/sensitivity_final_SOR.csv", row.names=FALSE)

message("Updated scenario tables created in tables/")
