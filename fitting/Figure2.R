###############################################################################
## Repository version: Figure 2
###############################################################################
## Requires only data_processed/Figure2_plot_data.Rdata.
## No raw simulation trajectories or model equations are required.
###############################################################################

library(ggplot2)
library(patchwork)
library(grid)

if (!dir.exists("figures")) dir.create("figures")
load("data_processed/Figure2_plot_data.Rdata")

f2a <- ggplot() +
  geom_line(data=Figure2_VR_trajectories,
            aes(x=time+1987,y=value,group=trajectory),alpha=0.025) +
  geom_line(data=Figure2_GR_trajectories,
            aes(x=time+1987,y=value,group=trajectory),colour="red",alpha=0.08) +
  annotate("segment",x=2007,xend=2007,y=0,yend=80,colour="darkcyan") +
  annotate("text",x=2006.8,y=40,colour="darkcyan",label="Calibration   Projection",size=6) +
  annotate("segment",x=2014,y=40,xend=2017,yend=40,colour="darkcyan",
           arrow=arrow(length=unit(0.2,"cm"),type="closed")) +
  annotate("text",label="Variable response model",x=2019.3,y=75,colour="black",size=6) +
  annotate("text",label="GR-only model",x=2016.6,y=68,colour="black",size=6) +
  annotate("segment",x=2010,xend=2012.2,y=75,yend=75,colour="black") +
  annotate("segment",x=2010,xend=2012.2,y=68,yend=68,colour="red") +
  labs(x="Time (years)",y="mf prevalence (%)") + theme_bw() +
  theme(legend.position="none",text=element_text(size=18))

f2b <- ggplot(Figure2_panelB_data,
              aes(x=period,y=prevalence,colour=model)) +
  geom_boxplot(width=0.55,outlier.size=0.7) +
  geom_hline(yintercept=1.4,colour="blue",linetype="dashed") +
  annotate("text",x=0.6,y=1.8,label="pOTTIS",colour="blue") +
  scale_colour_manual(values=c("Variable response"="black","GR-only"="red"),
                      labels=c("Variable response"="Variable-response","GR-only"="GR-only")) +
  scale_x_discrete(labels=c("Annual calibration"="Annual\ncalibration",
                            "Biannual projection"="Biannual\nprojection")) +
  labs(x=NULL,y="Final mf prevalence (%)",colour=NULL) + theme_bw() +
  theme(text=element_text(size=16),axis.title=element_text(size=17),
        axis.text=element_text(size=15),axis.text.x=element_text(size=14),
        legend.title=element_text(size=15),legend.text=element_text(size=14),
        legend.position="inside",legend.position.inside=c(0.78,0.92),
        legend.background=element_rect(fill="white",colour=NA))

p <- f2a+f2b+plot_annotation(tag_levels="A")+plot_layout(widths=c(2,1.25))
ggsave("figures/F2_100.png",plot=p,width=15,height=6,dpi=300)
