###############################################################################
## Figure S1: MDAi coverage in Asubende, Ghana, 1987-2024
###############################################################################
##
## PURPOSE
##
## Plot MDAi coverage through time
##
## Coverage values unavailable for 1987-1997 are set to 65% and the unavailable
## 1998 coverage value is set to 80%. No MDAi was implemented in 2011.
##
## This script is self-contained and does not require model simulation outputs.
##
## OUTPUT
##
##   figures/SF1.png
##
###############################################################################

library(ggplot2)

dir.create("figures", showWarnings = FALSE, recursive = TRUE)

coverage <- c(
  65, 65, 65, 65, 65,
  65, 65, 65, 65, 65,
  65, 81.7, 80.1, 80, 84.10,
  84.10, 74.70, 79.40, 87.10, 86.20,
  92.00, 92, 87.30, 88.00, 0.00, 74.20,
  50.60, 79.60, 81.90, 78.30, 83.00,
  79.50, 82.10, 78.70, 73.10, 77.50,
  77.10, 84.30
) / 100

coverage_data <- data.frame(
  Year = 1987:2024,
  Coverage = coverage
)

coverage_data$Data_source <- ifelse(
  coverage_data$Year %in% 1987:1997 |
    coverage_data$Year == 1998,
  "Estimated",
  "Measured"
)

coverage_data$Data_source <- factor(
  coverage_data$Data_source,
  levels = c("Measured", "Estimated")
)

p <- ggplot(
  coverage_data,
  aes(
    x = Year,
    y = Coverage,
    colour = Data_source
  )
) +
  geom_point() +
  scale_colour_manual(
    values = c(
      "Measured" = "#009988",
      "Estimated" = "#EE3377"
    )
  ) +
  coord_cartesian(
    xlim = c(1987, 2024),
    ylim = c(0, 1)
  ) +
  labs(
    x = "Year",
    y = "MDAi coverage (proportion)",
    colour = "Data source"
  ) +
  theme_bw() +
  theme(
    text = element_text(size = 18)
  )

ggsave(
  "figures/SF1.png",
  plot = p,
  width = 8,
  height = 6,
  dpi = 300
)
