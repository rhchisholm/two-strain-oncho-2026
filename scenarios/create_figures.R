# create figures summarising different endemic prevalence scenarios
library(ggplot2)
library(dplyr)
library(patchwork)
library(tidyr)

load("processed_outputs.RData")

if (!dir.exists("figures")) {
  dir.create("figures")
}

#####################################
# Figures original parameterisation #
#####################################

MDAi.labs <- c("Ivermectin", "Moxidectin")
names(MDAi.labs) <- c("1", "0")
df$MDAi_f <- factor(df$MDAi, levels = c('1', '0'))
# Reorder facets
df_reordered <- df 

# Final mf prevalence 
df_reorderedI<-df_reordered %>% filter(MDAi == "1")
df_reorderedM<-df_reordered %>% filter(MDAi == "0")
SOR.labs <- c("SOR present", "No SOR")
names(SOR.labs) <- c("1", "0")

p3 <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                y = Final_Prevalence, color=factor(MDA_frequency))) +
  geom_boxplot() + scale_y_log10() + facet_wrap(~as.factor(SOR_present), nrow=2) + 
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + coord_cartesian(ylim = c(0.0001, 10)) + 
  #ylim(c(0,6)) +
  facet_wrap(~SOR_present, labeller = labeller(SOR_present=SOR.labs), nrow=2) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  annotate("segment", x = 0.33, xend = 3.66, y = 0.014*100, yend = 0.014*100,
           colour = "blue", linetype = "dashed") +
  annotate("text", label = "pOTTIS", x = 0.75, y = 2.5, color = "blue") +
  labs(x="Endemic mf prevalence", y="Final mf prevalence (%)") + ggtitle('Ivermectin') +
  theme_bw() +  
  theme(legend.position = "none", text = element_text(size = 14))

p4 <- ggplot(df_reorderedM, aes(x = Endemic_prevalence, 
                                y = Final_Prevalence, color=factor(MDA_frequency))) +
  geom_boxplot() + scale_y_log10() + facet_wrap(~as.factor(SOR_present), nrow=2) + 
  theme(strip.text = element_text(face = "bold", size = 14)) +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + coord_cartesian(ylim = c(0.0001, 10)) + #
  #ylim(c(0,6)) +
  facet_wrap(~SOR_present, labeller = labeller(SOR_present=SOR.labs), nrow=2) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  annotate("segment", x = 0.33, xend = 3.66, y = 0.014*100, yend = 0.014*100,
           colour = "blue", linetype = "dashed") +
  annotate("text", label = "pOTTIS", x = 0.75, y = 2.5, color = "blue") +
  labs(x="Endemic mf prevalence", y="Final mf prevalence (%)") + 
  ggtitle('Moxidectin') + ylab(NULL) + theme_bw() + 
  theme(text = element_text(size = 14))

p<- p3+p4+ plot_annotation(tag_levels = 'A') 

ggsave("figures/F3.png", plot = p, width = 10, height = 6, dpi=300)

# SOR figures

# Reorder facets
df_reordered <- df 
df_reorderedI<-df_reordered %>% filter(MDAi == "1" & SOR_present == "1")
df_reorderedM<-df_reordered %>% filter(MDAi == "0" & SOR_present == "1")

p5 <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                y = Final_SOR_proportion_worms20yr, color=factor(MDA_frequency))) +
  geom_boxplot() +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0.0,0.3)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Final SOR worm proportion") + ggtitle('Ivermectin') +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14))

p6 <- ggplot(df_reorderedM, aes(x = Endemic_prevalence, 
                                y = Final_SOR_proportion_worms20yr, color=factor(MDA_frequency))) +
  geom_boxplot() +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0.0,0.3)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Final SOR worm proportion") + 
  ggtitle('Moxidectin') + ylab(NULL) + theme_bw() + theme(text = element_text(size = 14))

bp <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                y = Initial_SOR_proportion_worms, color=factor(MDAi))) +
  geom_boxplot() 

bp

# Build the plot to get the calculated statistics
ggplot_data <- ggplot_build(bp)$data[[1]]

# Access the values (these match the output of boxplot.stats for hinges)
# ymin_interactive is lower whisker, lower is Q1, middle is median, upper is Q3, ymax_interactive is upper whisker
q1_ggplot <- ggplot_data$lower
q3_ggplot <- ggplot_data$upper
q0_ggplot <- ggplot_data$ymin
q4_ggplot <- ggplot_data$ymax

patchwork2 <- p5 +
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q1_ggplot, ymax=q3_ggplot, fill="pink", alpha=0.4) + 
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q0_ggplot, ymax=q1_ggplot, fill="pink", alpha=0.2) + 
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q3_ggplot, ymax=q4_ggplot, fill="pink", alpha=0.2) + 
  p6 +
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q1_ggplot, ymax=q3_ggplot, fill="pink", alpha=0.4) +
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q0_ggplot, ymax=q1_ggplot, fill="pink", alpha=0.2) + 
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q3_ggplot, ymax=q4_ggplot, fill="pink", alpha=0.2) 

p<-patchwork2 + plot_annotation(tag_levels = 'A')

ggsave("figures/F5.png", plot = p, width = 10, height = 4, dpi=300)

# Final pottis figures
# Reorder facets
df_reordered <- df 
df_reorderedI<-df_reordered %>% filter(MDAi == "1")
df_reorderedM<-df_reordered %>% filter(MDAi == "0")

p1 <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                y = pOTTIS, color=factor(MDA_frequency)))+ #, na.rm = TRUE) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  facet_wrap(~as.factor(SOR_present), nrow=2) +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0,20)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Years to pOTTIS") + ggtitle('Ivermectin') +
  facet_wrap(~factor(SOR_present, levels=c('1','0'), labels = SOR.labs), nrow=2) +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14))+
  # Position star at missing mox plot location
  geom_point(
    data = data.frame(
      group = 3,
      value = 5,
      subgroup = "Y",
      SOR_present = 1
    ),
    aes(x = group, y = value),
    shape = 8,
    size = 4,
    color = "#8856a7",
    inherit.aes = FALSE,
    position = position_dodge(width = 0.8)
  )

p2 <- ggplot(df_reorderedM, aes(x = Endemic_prevalence, 
                                y = pOTTIS, color=factor(MDA_frequency))) +
  geom_boxplot() + facet_wrap(~as.factor(SOR_present), nrow=2) + 
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0,20)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Years to pOTTIS") + 
  facet_wrap(~factor(SOR_present, levels=c('1','0'), labels = SOR.labs), nrow=2) +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  ggtitle('Moxidectin') + ylab(NULL) + theme_bw() + 
  theme(text = element_text(size = 14))

p <- p1 + p2 + plot_annotation(tag_levels = 'A')

ggsave("figures/F4.png", plot = p, width = 10, height = 6, dpi=300)

###############################
# Figures SA parameterisation #
###############################

dfSA$MDAi_f <- factor(dfSA$MDAi, levels = c('1', '0'))
# Reorder facets
df_reordered <- dfSA 

# Final mf prevalence 
df_reorderedI<-df_reordered %>% filter(MDAi == "1")
df_reorderedM<-df_reordered %>% filter(MDAi == "0")
SOR.labs <- c("SOR present", "No SOR")
names(SOR.labs) <- c("1", "0")

p3s <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                 y = Final_Prevalence, color=factor(MDA_frequency))) +
  geom_boxplot() + scale_y_log10() + facet_wrap(~as.factor(SOR_present), nrow=2) + 
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + coord_cartesian(ylim = c(0.0001, 10)) + 
  #ylim(c(0,6)) +
  facet_wrap(~SOR_present, labeller = labeller(SOR_present=SOR.labs), nrow=2) +
  theme(strip.text = element_text(face = "bold", size = 14)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  annotate("segment", x = 0.33, xend = 3.66, y = 0.014*100, yend = 0.014*100,
           colour = "blue", linetype = "dashed") +
  annotate("text", label = "pOTTIS", x = 0.75, y = 2.5, color = "blue") +
  labs(x="Endemic mf prevalence", y="Final mf prevalence (%)") + ggtitle('Ivermectin') +
  theme_bw() +  
  theme(legend.position = "none", text = element_text(size = 14))

p4s <- ggplot(df_reorderedM, aes(x = Endemic_prevalence, 
                                 y = Final_Prevalence, color=factor(MDA_frequency))) +
  geom_boxplot() + scale_y_log10() + facet_wrap(~as.factor(SOR_present), nrow=2) + 
  theme(strip.text = element_text(face = "bold", size = 14)) +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + coord_cartesian(ylim = c(0.0001, 10)) + #
  #ylim(c(0,6)) +
  facet_wrap(~SOR_present, labeller = labeller(SOR_present=SOR.labs), nrow=2) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  annotate("segment", x = 0.33, xend = 3.66, y = 0.014*100, yend = 0.014*100,
           colour = "blue", linetype = "dashed") +
  annotate("text", label = "pOTTIS", x = 0.75, y = 2.5, color = "blue") +
  labs(x="Endemic mf prevalence", y="Final mf prevalence (%)") + 
  ggtitle('Moxidectin') + ylab(NULL) + theme_bw() + 
  theme(text = element_text(size = 14))

p<-p3s+p4s+ plot_annotation(tag_levels = 'A') 

ggsave("figures/F3_SA.png", plot = p, width = 10, height = 6, dpi=300)

# SOR figures

# Reorder facets
df_reordered <- dfSA 
df_reorderedI<-df_reordered %>% filter(MDAi == "1" & SOR_present == "1")
df_reorderedM<-df_reordered %>% filter(MDAi == "0" & SOR_present == "1")

p5s <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                 y = Final_SOR_proportion_worms20yr, color=factor(MDA_frequency))) +
  geom_boxplot() +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0.0,0.3)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Final SOR worm proportion") + ggtitle('Ivermectin') +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14))

p6s <- ggplot(df_reorderedM, aes(x = Endemic_prevalence, 
                                 y = Final_SOR_proportion_worms20yr, color=factor(MDA_frequency))) +
  geom_boxplot() +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0.0,0.3)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Final SOR worm proportion") + 
  ggtitle('Moxidectin') + ylab(NULL) + theme_bw() + theme(text = element_text(size = 14))

bps <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                 y = Initial_SOR_proportion_worms, color=factor(MDAi))) +
  geom_boxplot() 

bps

# Build the plot to get the calculated statistics
ggplot_data <- ggplot_build(bps)$data[[1]]

# Access the values (these match the output of boxplot.stats for hinges)
# ymin_interactive is lower whisker, lower is Q1, middle is median, upper is Q3, ymax_interactive is upper whisker
q1_ggplot <- ggplot_data$lower
q3_ggplot <- ggplot_data$upper
q0_ggplot <- ggplot_data$ymin
q4_ggplot <- ggplot_data$ymax

patchwork2 <- p5s +
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q1_ggplot, ymax=q3_ggplot, fill="pink", alpha=0.4) + 
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q0_ggplot, ymax=q1_ggplot, fill="pink", alpha=0.2) + 
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q3_ggplot, ymax=q4_ggplot, fill="pink", alpha=0.2) + 
  p6s +
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q1_ggplot, ymax=q3_ggplot, fill="pink", alpha=0.4) +
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q0_ggplot, ymax=q1_ggplot, fill="pink", alpha=0.2) + 
  annotate("rect", xmin=c(0.5,1.5,2.5), xmax=c(1.5,2.5,3.5), ymin=q3_ggplot, ymax=q4_ggplot, fill="pink", alpha=0.2) 

p<- patchwork2 + plot_annotation(tag_levels = 'A')

ggsave("figures/F5_SA.png", plot = p, width = 10, height = 4, dpi=300)

# Final pottis figures
# Reorder facets
df_reordered <- dfSA 
df_reorderedI<-df_reordered %>% filter(MDAi == "1")
df_reorderedM<-df_reordered %>% filter(MDAi == "0")


p1s <- ggplot(df_reorderedI, aes(x = Endemic_prevalence, 
                                 y = pOTTIS, color=factor(MDA_frequency)))+ #, na.rm = TRUE) +
  geom_boxplot(position = position_dodge2(preserve = "single")) +
  facet_wrap(~as.factor(SOR_present), nrow=2) +
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0,20)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Years to pOTTIS") + ggtitle('Ivermectin') +
  facet_wrap(~factor(SOR_present, levels=c('1','0'), labels = SOR.labs), nrow=2) +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14))+
  # Position star at missing mox plot location
  geom_point(
    data = data.frame(
      group = 3,
      value = 5,
      subgroup = "Y",
      SOR_present = 1
    ),
    aes(x = group, y = value),
    shape = 8,
    size = 4,
    color = "#8856a7",
    inherit.aes = FALSE,
    position = position_dodge(width = 0.8)
  )

p2s <- ggplot(df_reorderedM, aes(x = Endemic_prevalence, 
                                 y = pOTTIS, color=factor(MDA_frequency))) +
  geom_boxplot() + facet_wrap(~as.factor(SOR_present), nrow=2) + 
  scale_color_manual(
    values = c("1" = "#8856a7", "0.5" = "#9ebcda"), # Assign specific colors
    labels = c("1" = "Annual", "0.5" = "Biannual") # Assign custom labels
  ) + ylim(c(0,20)) +
  scale_x_discrete(labels = c("30%", "50%", "70%")) +
  labs(color = "MDA frequency") +
  labs(x="Endemic mf prevalence", y="Years to pOTTIS") + 
  facet_wrap(~factor(SOR_present, levels=c('1','0'), labels = SOR.labs), nrow=2) +
  theme(strip.text = element_text(face = "bold", size = 12)) +
  ggtitle('Moxidectin') + ylab(NULL) + theme_bw() + 
  theme(text = element_text(size = 14))

p<- p1s + p2s + plot_annotation(tag_levels = 'A')

ggsave("figures/F4_SA.png", plot = p, width = 10, height = 6, dpi=300)




