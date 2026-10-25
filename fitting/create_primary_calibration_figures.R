###############################################################################
## Repository version: create supplementary Figures S2-S4
###############################################################################
## Requires only:
##   data_processed/primary_calibration_figure_data.Rdata
## No raw simulation trajectories are required.
###############################################################################

library(ggplot2)
library(latex2exp)
library(patchwork)

figure_dir <- "figures"
dir.create(figure_dir, showWarnings=FALSE, recursive=TRUE)
load("data_processed/primary_calibration_figure_data.Rdata")

median_95 <- function(x) data.frame(
  y=median(x,na.rm=TRUE),
  ymin=quantile(x,0.025,na.rm=TRUE),
  ymax=quantile(x,0.975,na.rm=TRUE)
)

## Figure S2
b5 <- ggplot(S2_VR_data,
             aes(x=Final_mf_prevalence*100,
                 y=Increase_in_SOR_worm_proportion,
                 color=Cumulative_Effect)) +
  geom_point(alpha=0.4) +
  scale_color_manual(values=c("0"="#EE3377","1"="#009988"),
                     labels=c("0"="None","1"=TeX(r"(Present with $\zeta=0.3$)"))) +
  labs(x="Final mf prevalence (%)", y="Final increase in SOR worm proportion",
       title="Variable response model", color="Cumulative effect") +
  annotate("rect",xmin=10,xmax=12,ymin=0.10,ymax=0.20,fill="blue",alpha=0.20,color="blue") +
  annotate("rect",xmin=10,xmax=12,ymin=0.05,ymax=0.10,fill="#EE7733",alpha=0.20,color="#EE7733") +
  annotate("rect",xmin=10,xmax=12,ymin=0.00,ymax=0.05,fill="#EE0000",alpha=0.20,color="#EE0000") +
  annotate("rect",xmin=8,xmax=14,ymin=0.10,ymax=0.20,fill=NA,color="black",linetype="dashed") +
  theme_bw() + theme(legend.position="none",text=element_text(size=18))

b6 <- ggplot(S2_GR_data,
             aes(x=y,y=Final_mf_prevalence,color=Cumulative_effect)) +
  geom_point(alpha=0.4,size=1.5) +
  scale_color_manual(values=c("0"="#EE3377","1"="#009988"),
                     labels=c("0"="None","1"=TeX(r"(Present with $\zeta=0.3$)"))) +
  annotate("rect",xmin=0,xmax=1,ymin=10,ymax=12,fill="blue",alpha=0.20,color="blue") +
  annotate("rect",xmin=0,xmax=1,ymin=8,ymax=14,fill=NA,color="black",linetype="dashed") +
  labs(x="y",y="Final mf prevalence (%)",title="GR-only model",color="Cumulative effect") +
  theme_bw() + theme(text=element_text(size=18))

p <- b5+b6+plot_annotation(tag_levels="A")
ggsave(file.path(figure_dir,"SF2_100.png"),p,width=15,height=6,dpi=300)

## Figure S3
D <- S34_primary_data
S3_A <- data.frame(parameter=rep(c("c_GRGR","c0_SORSOR","c0_GRSOR","y"),each=nrow(D)),
                   value=c(D$c_GRGR,D$c0_SORSOR,D$c0_GRSOR,D$y))
S3_B <- data.frame(parameter=rep(c("s_SORSOR","s_SORGR","f"),each=nrow(D)),
                   value=c(D$s_SORSOR,D$s_SORGR,D[["psi_SOR/psi_GR"]]))
S3_A$parameter <- factor(S3_A$parameter,levels=c("c_GRGR","c0_SORSOR","c0_GRSOR","y"))
S3_B$parameter <- factor(S3_B$parameter,levels=c("s_SORSOR","s_SORGR","f"))

S3a <- ggplot(S3_A,aes(parameter,value)) +
  stat_summary(fun.data=median_95,geom="errorbar",width=0.12,linewidth=0.8) +
  stat_summary(fun=median,geom="point",size=2.8) +
  scale_x_discrete(labels=c("c_GRGR"=expression(c[GR*",GR"]),
                            "c0_SORSOR"=expression(c[SOR*",SOR"]^0),
                            "c0_GRSOR"=expression(c[GR*",SOR"]^0),"y"=expression(y))) +
  scale_y_continuous(limits=c(0,1),expand=expansion(mult=c(0,0.03))) +
  labs(x=NULL,y="Parameter value") + theme_bw() +
  theme(text=element_text(size=16),axis.text.x=element_text(size=15),panel.grid.minor=element_blank())

S3b <- ggplot(S3_B,aes(parameter,value)) +
  stat_summary(fun.data=median_95,geom="errorbar",width=0.12,linewidth=0.8) +
  stat_summary(fun=median,geom="point",size=2.8) +
  scale_x_discrete(labels=c("s_SORSOR"=expression(s[SOR*",SOR"]),
                            "s_SORGR"=expression(s[SOR*",GR"]),"f"=expression(f))) +
  scale_y_continuous(limits=c(0,20),expand=expansion(mult=c(0,0.03))) +
  labs(x=NULL,y="Parameter value") + theme_bw() +
  theme(text=element_text(size=16),axis.text.x=element_text(size=15),panel.grid.minor=element_blank())

S3 <- S3a+S3b+plot_annotation(tag_levels="A")
ggsave(file.path(figure_dir,"SF3_100.png"),S3,width=13,height=5.5,dpi=300)

## Figure S4
S4_A <- data.frame(outcome=rep(c("Final mf prevalence","Increase in SOR"),each=nrow(D)),
                   value=c(D$Final_mf_prevalence_annual,D$Increase_in_SOR_worm_proportion))
S4_B <- data.frame(outcome=rep(c("Final mf prevalence","Further increase in SOR"),each=nrow(D)),
                   value=c(D$Final_mf_prevalence_biannual,D$Further_increase_in_SOR_worm_proportion))
S4_A$outcome <- factor(S4_A$outcome,levels=c("Final mf prevalence","Increase in SOR"))
S4_B$outcome <- factor(S4_B$outcome,levels=c("Final mf prevalence","Further increase in SOR"))

S4a <- ggplot(S4_A,aes(outcome,value)) +
  stat_summary(fun.data=median_95,geom="errorbar",width=0.10,linewidth=0.8) +
  stat_summary(fun=median,geom="point",size=2.5) + coord_cartesian(ylim=c(0,0.30)) +
  scale_x_discrete(labels=c("Final mf prevalence"="mf prevalence","Increase in SOR"="Increase in\nSOR worm proportion")) +
  labs(x=NULL,y="Proportion",title="End of annual MDAi calibration") + theme_bw() +
  theme(text=element_text(size=16),axis.text.x=element_text(size=14),panel.grid.minor=element_blank())

S4b <- ggplot(S4_B,aes(outcome,value)) +
  stat_summary(fun.data=median_95,geom="errorbar",width=0.10,linewidth=0.8) +
  stat_summary(fun=median,geom="point",size=2.5) + coord_cartesian(ylim=c(0,0.30)) +
  scale_x_discrete(labels=c("Final mf prevalence"="mf prevalence","Further increase in SOR"="Further increase in\nSOR worm proportion")) +
  labs(x=NULL,y="Proportion",title="End of subsequent biannual MDAi projection") + theme_bw() +
  theme(text=element_text(size=16),axis.text.x=element_text(size=14),panel.grid.minor=element_blank())

S4 <- S4a+S4b+plot_annotation(tag_levels="A")
ggsave(file.path(figure_dir,"SF4_100.png"),S4,width=13,height=5.5,dpi=300)
