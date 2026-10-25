###############################################################################
## Create scenario figures
###############################################################################
##
## Uses processed_all_scenarios.RData created by process_all_scenario_outputs.R.
##
## Primary figures reproduce the main scenario comparisons across baseline
## prevalence settings. Sensitivity figures compare all five calibration
## analyses at ~70% baseline prevalence under ivermectin.
###############################################################################

library(ggplot2)
library(dplyr)
library(patchwork)

load("processed_all_scenarios.RData")

processed_all$Treatment_frequency <- factor(
  processed_all$Treatment_frequency,
  levels = c("Biannual", "Annual")
)

model_labels <- c(
  "VR" = "SOR present",
  "GR-only" = "No SOR present"
)

prepare_box_positions <- function(dat) {
  
  dat %>%
    mutate(
      x_base = as.numeric(Baseline_prevalence),
      
      x_plot = x_base +
        ifelse(
          Treatment_frequency == "Biannual",
          -0.18,
          0.18
        )
    )
}

if (!dir.exists("figures")) dir.create("figures")

prev_labels <- c(`2942`="30%", `5350`="50%", `27000`="70%")

primary <- processed_all %>%
  filter(Calibration == "Primary") %>%
  mutate(
    Baseline_prevalence=factor(
      unname(prev_labels[as.character(ABR)]),
      levels=c("30%","50%","70%")
    )
  )

primary <- prepare_box_positions(primary)
###############################################################################
## F3: Final mf prevalence - primary scenarios
###############################################################################

make_prev_panel <- function(dat, drug_name, show_legend = TRUE) {
  
  ggplot(
    dat %>% filter(Drug == drug_name),
    aes(
      x = Baseline_prevalence,
      y = Final_Prevalence_percent,
      colour = Treatment_frequency
    )
  ) +
    
    geom_boxplot(
      aes(
        group = interaction(
          Baseline_prevalence,
          Treatment_frequency
        )
      ),
      position = position_dodge(width = 0.75),
      width = 0.6,
      outlier.size = 0.7
    ) +
    
    scale_y_log10() +
    
    facet_wrap(
      ~Model,
      nrow = 2,
      labeller = labeller(Model = model_labels)
    ) +
    
    geom_hline(
      yintercept = 1.4,
      linetype = "dashed",
      colour = "blue"
    ) +
    
    annotate(
      "text",
      x = 0.75,
      y = 2.5,
      label = "pOTTIS",
      colour = "blue"
    ) +
    
    coord_cartesian(
      ylim = c(0.0001, 10)
    ) +
    
    scale_colour_manual(
      values = c(
        "Biannual" = "#9ebcda",
        "Annual" = "#8856a7"
      ),
      breaks = c(
        "Biannual",
        "Annual"
      )
    ) +
    
    labs(
      x = "Endemic mf prevalence",
      y = "Final mf prevalence (%)",
      colour = "MDA frequency",
      title = drug_name
    ) +
    
    theme_bw() +
    
    theme(
      text = element_text(size = 14),
      strip.text = element_text(face = "bold"),
      legend.position =
        if (show_legend) "right" else "none"
    )
}

p3a <- make_prev_panel(primary,"Ivermectin",FALSE)
p3b <- make_prev_panel(primary,"Moxidectin",TRUE) + ylab(NULL)

ggsave(
  "figures/F3_updated.png",
  p3a+p3b+plot_annotation(tag_levels="A"),
  width=10, height=6, dpi=300
)

###############################################################################
## F4: Years to pOTTIS - primary scenarios
###############################################################################

make_pottis_panel <- function(dat, drug_name, show_legend = TRUE) {
  
  plot_data <- dat %>%
    filter(Drug == drug_name)
  
  ## Identify groups where no trajectories reached pOTTIS
  missing_groups <- plot_data %>%
    group_by(
      Model,
      Baseline_prevalence,
      Treatment_frequency,
      x_base,
      x_plot
    ) %>%
    summarise(
      all_missing = all(is.na(pOTTIS)),
      .groups = "drop"
    ) %>%
    filter(all_missing)
  
  ggplot(
    plot_data,
    aes(
      x = x_plot,
      y = pOTTIS,
      colour = Treatment_frequency,
      group = interaction(
        x_base,
        Treatment_frequency
      )
    )
  ) +
    
    geom_boxplot(
      width = 0.28,
      outlier.size = 0.7,
      na.rm = TRUE
    ) +
    
    ## Star where no trajectories reached pOTTIS
    geom_point(
      data = missing_groups,
      aes(
        x = x_plot,
        y = 5,
        colour = Treatment_frequency
      ),
      inherit.aes = FALSE,
      shape = 8,
      size = 4,
      show.legend = FALSE
    ) +
    
    facet_wrap(
      ~Model,
      nrow = 2,
      labeller = labeller(
        Model = model_labels
      )
    ) +
    
    scale_x_continuous(
      breaks = 1:3,
      labels = c(
        "30%",
        "50%",
        "70%"
      )
    ) +
    
    scale_colour_manual(
      values = c(
        "Biannual" = "#9ebcda",
        "Annual" = "#8856a7"
      ),
      breaks = c(
        "Biannual",
        "Annual"
      )
    ) +
    
    coord_cartesian(
      xlim = c(0.5, 3.5),
      ylim = c(0, 20)
    ) +
    
    labs(
      x = "Endemic mf prevalence",
      y = "Years to pOTTIS",
      colour = "MDA frequency",
      title = drug_name
    ) +
    
    theme_bw() +
    
    theme(
      text = element_text(size = 14),
      strip.text = element_text(face = "bold"),
      legend.position =
        if (show_legend) "right" else "none"
    )
}

p4a <- make_pottis_panel(primary,"Ivermectin",FALSE)
p4b <- make_pottis_panel(primary,"Moxidectin",TRUE) + ylab(NULL)

ggsave(
  "figures/F4_updated.png",
  p4a+p4b+plot_annotation(tag_levels="A"),
  width=10, height=6, dpi=300
)

###############################################################################
## F5: Final SOR worm proportion - primary VR scenarios
##
## Pink shading shows the distribution of the initial SOR worm proportion
## before MDA:
##
##   darker band = interquartile range (Q1-Q3)
##   lighter band = values within the boxplot whiskers (excluding outliers)
##
###############################################################################

primary_vr <- primary %>%
  filter(Model == "VR") %>%
  mutate(
    Treatment_frequency = factor(
      Treatment_frequency,
      levels = c("Biannual", "Annual")
    )
  )

###############################################################################
## Calculate initial SOR distribution for each baseline prevalence
###############################################################################

initial_SOR_stats <- primary_vr %>%
  
  ## Initial SOR distribution is independent of treatment frequency.
  ## Use one treatment scenario to avoid duplicating each parameter set.
  filter(
    Drug == "Ivermectin",
    Treatment_frequency == "Annual"
  ) %>%
  
  group_by(Baseline_prevalence) %>%
  
  summarise(
    
    Q1 = quantile(
      Initial_SOR_proportion_worms,
      0.25,
      na.rm = TRUE
    ),
    
    Q3 = quantile(
      Initial_SOR_proportion_worms,
      0.75,
      na.rm = TRUE
    ),
    
    IQR_value = IQR(
      Initial_SOR_proportion_worms,
      na.rm = TRUE
    ),
    
    ## Standard boxplot whiskers:
    ## most extreme values within 1.5 x IQR of Q1/Q3
    
    lower_whisker = min(
      Initial_SOR_proportion_worms[
        Initial_SOR_proportion_worms >=
          Q1 - 1.5 * IQR_value
      ],
      na.rm = TRUE
    ),
    
    upper_whisker = max(
      Initial_SOR_proportion_worms[
        Initial_SOR_proportion_worms <=
          Q3 + 1.5 * IQR_value
      ],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

###############################################################################
## Convert prevalence categories to plotting positions
###############################################################################

initial_SOR_stats$x <-
  as.numeric(
    initial_SOR_stats$Baseline_prevalence
  )

p5a <- ggplot(
  primary_vr %>%
    filter(Drug == "Ivermectin"),
  
  aes(
    x = Baseline_prevalence,
    y = Final_SOR_proportion_worms20yr,
    colour = Treatment_frequency
  )
) +
  
  ## Light band: values within whiskers
  geom_rect(
    data = initial_SOR_stats,
    aes(
      xmin = x - 0.5,
      xmax = x + 0.5,
      ymin = lower_whisker,
      ymax = upper_whisker
    ),
    inherit.aes = FALSE,
    fill = "pink",
    alpha = 0.20
  ) +
  
  ## Dark band: interquartile range
  geom_rect(
    data = initial_SOR_stats,
    aes(
      xmin = x - 0.5,
      xmax = x + 0.5,
      ymin = Q1,
      ymax = Q3
    ),
    inherit.aes = FALSE,
    fill = "pink",
    alpha = 0.40
  ) +
  
  geom_boxplot(
    aes(
      group = interaction(
        Baseline_prevalence,
        Treatment_frequency
      )
    ),
    position = position_dodge(width = 0.75),
    width = 0.6
  ) +
  
  scale_colour_manual(
    values = c(
      "Biannual" = "#9ebcda",
      "Annual" = "#8856a7"
    ),
    breaks = c(
      "Biannual",
      "Annual"
    )
  ) +
  
  coord_cartesian(
    ylim = c(0, 0.3)
  ) +
  
  labs(
    x = "Endemic mf prevalence",
    y = "Final SOR worm proportion",
    colour = "MDA frequency",
    title = "Ivermectin"
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(size = 14),
    legend.position = "none"
  )


p5b <- ggplot(
  primary_vr %>%
    filter(Drug == "Moxidectin"),
  
  aes(
    x = Baseline_prevalence,
    y = Final_SOR_proportion_worms20yr,
    colour = Treatment_frequency
  )
) +
  
  ## Light band: values within whiskers
  geom_rect(
    data = initial_SOR_stats,
    aes(
      xmin = x - 0.5,
      xmax = x + 0.5,
      ymin = lower_whisker,
      ymax = upper_whisker
    ),
    inherit.aes = FALSE,
    fill = "pink",
    alpha = 0.20
  ) +
  
  ## Dark band: interquartile range
  geom_rect(
    data = initial_SOR_stats,
    aes(
      xmin = x - 0.5,
      xmax = x + 0.5,
      ymin = Q1,
      ymax = Q3
    ),
    inherit.aes = FALSE,
    fill = "pink",
    alpha = 0.40
  ) +
  
  geom_boxplot(
    aes(
      group = interaction(
        Baseline_prevalence,
        Treatment_frequency
      )
    ),
    position = position_dodge(width = 0.75),
    width = 0.6
  ) +
  
  scale_colour_manual(
    values = c(
      "Biannual" = "#9ebcda",
      "Annual" = "#8856a7"
    ),
    breaks = c(
      "Biannual",
      "Annual"
    )
  ) +
  
  coord_cartesian(
    ylim = c(0, 0.3)
  ) +
  
  labs(
    x = "Endemic mf prevalence",
    y = NULL,
    colour = "MDA frequency",
    title = "Moxidectin"
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(size = 14)
  )

p5 <-p5a + p5b + plot_annotation(tag_levels = "A")

ggsave(
  "figures/F5_updated.png",
  plot = p5,
  width = 10,
  height = 4.5,
  dpi = 300
)

###############################################################################
## Sensitivity data: 70% baseline prevalence, ivermectin
###############################################################################

sens <- processed_all %>%
  filter(
    ABR==27000,
    Drug=="Ivermectin",
    Calibration %in% c(
      "Primary","Moderate","Weak","No selection","Wider prevalence"
    )
  ) %>%
  mutate(
    Calibration=factor(
      Calibration,
      levels=c(
        "Primary","Moderate","Weak","No selection","Wider prevalence"
      )
    )
  )

###############################################################################
## Sensitivity Figure 1: Final mf prevalence
###############################################################################

s_prev <- ggplot(
  sens,
  aes(
    x = Calibration,
    y = Final_Prevalence_percent,
    colour = Treatment_frequency
  )
) +
  geom_boxplot(
    aes(
      group = interaction(
        Calibration,
        Treatment_frequency
      )
    ),
    position = position_dodge(width = 0.75),
    outlier.size = 0.6
  ) +
  scale_colour_manual(
    values = c(
      "Biannual" = "#9ebcda",
      "Annual" = "#8856a7"
    ),
    breaks = c(
      "Biannual",
      "Annual"
    )
  ) +
  scale_y_log10() +
  facet_wrap(
    ~Model,
    nrow = 1,
    labeller = labeller(Model = model_labels)
  ) +
  geom_hline(
    yintercept = 1.4,
    linetype = "dashed",
    colour = "blue"
  ) +
  annotate(
    "text",
    x = 3.25,
    y = 1.1,
    label = "pOTTIS",
    colour = "blue"
  ) +
  coord_cartesian(
    ylim = c(0.01, 10)
  ) +
  labs(
    x = NULL,
    y = "Final mf prevalence (%)",
    colour = "MDA frequency"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 14),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

ggsave(
  "figures/Sensitivity_final_prevalence.png",
  s_prev,
  width=11, height=5.5, dpi=300
)

###############################################################################
## Sensitivity Figure 2: Years to pOTTIS
###############################################################################
sens_pottis <- sens %>%
  mutate(
    x_base = as.numeric(Calibration),
    x_plot = x_base +
      ifelse(
        Treatment_frequency == "Biannual",
        -0.18,
        0.18
      )
  )

## Identify groups where no simulations reach pOTTIS
missing_pottis_sens <- sens_pottis %>%
  group_by(
    Calibration,
    Model,
    Treatment_frequency,
    x_base,
    x_plot
  ) %>%
  summarise(
    all_missing = all(is.na(pOTTIS)),
    .groups = "drop"
  ) %>%
  filter(all_missing)

s_pottis <- ggplot(
  sens_pottis,
  aes(
    x = x_plot,
    y = pOTTIS,
    colour = Treatment_frequency,
    group = interaction(
      Calibration,
      Treatment_frequency
    )
  )
) +
  
  geom_boxplot(
    width = 0.28,
    outlier.size = 0.6,
    na.rm = TRUE
  ) +
  
  ## Star where no simulations reach pOTTIS
  geom_point(
    data = missing_pottis_sens,
    aes(
      x = x_plot,
      y = 5,
      colour = Treatment_frequency
    ),
    inherit.aes = FALSE,
    shape = 8,
    size = 4,
    show.legend = FALSE
  ) +
  
  facet_wrap(
    ~Model,
    nrow = 1,
    labeller = labeller(Model = model_labels)
  ) +
  
  scale_x_continuous(
    breaks = 1:5,
    labels = c(
      "Primary",
      "Moderate",
      "Weak",
      "No selection",
      "Wider prevalence"
    )
  ) +
  
  scale_colour_manual(
    values = c(
      "Biannual" = "#9ebcda",
      "Annual" = "#8856a7"
    ),
    breaks = c(
      "Biannual",
      "Annual"
    )
  ) +
  
  coord_cartesian(
    xlim = c(0.5, 5.5),
    ylim = c(0, 20)
  ) +
  
  labs(
    x = NULL,
    y = "Years to pOTTIS",
    colour = "MDA frequency"
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(size = 14),
    strip.text = element_text(face = "bold"),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

ggsave(
  "figures/Sensitivity_pOTTIS.png",
  s_pottis,
  width=11, height=5.5, dpi=300
)

###############################################################################
## Sensitivity Figure 3: Final SOR worm proportion - VR only
###############################################################################

###############################################################################
## Initial SOR distributions for sensitivity analyses
###############################################################################

initial_SOR_stats_sens <- sens %>%
  filter(
    Model == "VR",
    Treatment_frequency == "Annual"
  ) %>%
  group_by(Calibration) %>%
  summarise(
    Q1 = quantile(
      Initial_SOR_proportion_worms,
      0.25,
      na.rm = TRUE
    ),
    
    Q3 = quantile(
      Initial_SOR_proportion_worms,
      0.75,
      na.rm = TRUE
    ),
    
    IQR_value = IQR(
      Initial_SOR_proportion_worms,
      na.rm = TRUE
    ),
    
    lower_whisker = min(
      Initial_SOR_proportion_worms[
        Initial_SOR_proportion_worms >=
          Q1 - 1.5 * IQR_value
      ],
      na.rm = TRUE
    ),
    
    upper_whisker = max(
      Initial_SOR_proportion_worms[
        Initial_SOR_proportion_worms <=
          Q3 + 1.5 * IQR_value
      ],
      na.rm = TRUE
    ),
    
    .groups = "drop"
  )

## Numeric x positions of the five calibration categories
initial_SOR_stats_sens$x <-
  as.numeric(initial_SOR_stats_sens$Calibration)

###############################################################################
## Sensitivity Figure: Final SOR worm proportion
###############################################################################

s_sor <- ggplot(
  sens %>%
    filter(Model == "VR"),
  
  aes(
    x = Calibration,
    y = Final_SOR_proportion_worms20yr,
    colour = Treatment_frequency
  )
) +
  
  ## Light pink:
  ## initial SOR values within boxplot whiskers
  geom_rect(
    data = initial_SOR_stats_sens,
    aes(
      xmin = x - 0.5,
      xmax = x + 0.5,
      ymin = lower_whisker,
      ymax = upper_whisker
    ),
    inherit.aes = FALSE,
    fill = "pink",
    alpha = 0.20
  ) +
  
  ## Dark pink:
  ## initial SOR interquartile range
  geom_rect(
    data = initial_SOR_stats_sens,
    aes(
      xmin = x - 0.5,
      xmax = x + 0.5,
      ymin = Q1,
      ymax = Q3
    ),
    inherit.aes = FALSE,
    fill = "pink",
    alpha = 0.40
  ) +
  
  geom_boxplot(
    aes(
      group = interaction(
        Calibration,
        Treatment_frequency
      )
    ),
    position = position_dodge(
      width = 0.75
    ),
    width = 0.6,
    outlier.size = 0.6
  ) +
  
  scale_colour_manual(
    values = c(
      "Biannual" = "#9ebcda",
      "Annual" = "#8856a7"
    ),
    breaks = c(
      "Biannual",
      "Annual"
    )
  ) +
  
  labs(
    x = NULL,
    y = "Final SOR worm proportion",
    colour = "MDA frequency"
  ) +
  
  theme_bw() +
  
  theme(
    text = element_text(size = 14),
    axis.text.x = element_text(
      angle = 30,
      hjust = 1
    )
  )

ggsave(
  "figures/Sensitivity_final_SOR.png",
  s_sor,
  width=8, height=5.5, dpi=300
)

message("Updated scenario figures created in figures/")
