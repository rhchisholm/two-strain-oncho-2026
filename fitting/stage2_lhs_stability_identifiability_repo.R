###############################################################################
## Stage 2 stability/sensitivity figures: repository version
## Requires only data_processed/stage2_stability_repo_data.Rdata
###############################################################################
figure_dir="stability_figures"; dir.create(figure_dir,showWarnings=FALSE,recursive=TRUE)
load("data_processed/stage2_stability_repo_data.Rdata")

labs <- list(c_GRGR=expression(c[GR*",GR"]),c0_SORSOR=expression(c[SOR*",SOR"]^0),
 c0_GRSOR=expression(c[GR*",SOR"]^0),s_SORSOR=expression(s[SOR*",SOR"]),
 s_SORGR=expression(s[SOR*",GR"]),f=expression(f),y=expression(y))
lims <- list(c_GRGR=c(.9,1),c0_SORSOR=c(.01,1),c0_GRSOR=c(.01,1),
 s_SORSOR=c(0,20),s_SORGR=c(0,20),f=c(1,10),y=c(0,1))
pars=c("c_GRGR","c0_SORSOR","c0_GRSOR","s_SORSOR","s_SORGR","f","y")

## Primary acceptance stability
png(file.path(figure_dir,"Stage2_acceptance_stability_100.png"),1800,1500,res=300)
plot(acceptance_summary$N_MDA_LHS,acceptance_summary$acceptance_Beliefs_2_3,type="b",
 ylim=range(c(acceptance_summary$acceptance_Beliefs_2_3,acceptance_summary$final_retention_proportion)),
 xlab="Number of MDA-effect LHS samples",ylab="Proportion of parameter combinations",cex.axis=1.3,cex.lab=1.4)
lines(acceptance_summary$N_MDA_LHS,acceptance_summary$final_retention_proportion,type="b",lty=2)
legend("topright",c("Satisfy Beliefs 2 and 3","Final retained"),lty=c(1,2),pch=1,bty="n",cex=1.2)
dev.off()

## Seven-parameter convergence
plot_conv=function(p){x=parameter_stability[parameter_stability$parameter==p,]
 plot(x$N,x$median,type="b",ylim=lims[[p]],xlab="Number of MDA-effect LHS samples",
      ylab=labs[[p]],cex.axis=1.3,cex.lab=1.4);lines(x$N,x$lower,type="b",lty=2);lines(x$N,x$upper,type="b",lty=2)}
png(file.path(figure_dir,"Stage2_parameter_stability_100_all_parameters.png"),3600,2200,res=300)
layout(matrix(1:8,nrow=2,byrow=TRUE));par(mar=c(5,6,2,1))
for(p in pars[1:3]) plot_conv(p)
par(mar=c(0,0,0,0));plot.new();legend("center",c("Median","2.5th/97.5th percentile"),lty=c(1,2),pch=1,bty="n",cex=1.2)
par(mar=c(5,6,2,1));for(p in pars[4:7]) plot_conv(p)
dev.off();layout(1)

## First 80 versus final 20 pairwise
d=stage2_primary_pairwise; e=d$sampling_block=="First 80"; l=!e
pairs=list(c("s_SORSOR","s_SORGR"),c("s_SORSOR","f"),c("s_SORSOR","y"),c("s_SORGR","f"),c("s_SORGR","y"),c("f","y"))
png(file.path(figure_dir,"Stage2_pairwise_stability_100.png"),3600,2000,res=300)
layout(matrix(c(1,2,3,7,4,5,6,8),2,byrow=TRUE),widths=c(1,1,1,.7));par(mar=c(4.5,4.5,2,1))
for(z in pairs){x=z[1];y=z[2];plot(d[[x]][e],d[[y]][e],pch=1,xlim=range(d[[x]],na.rm=T),ylim=range(d[[y]],na.rm=T),xlab=labs[[x]],ylab=labs[[y]],cex.axis=1.2,cex.lab=1.3);points(d[[x]][l],d[[y]][l],pch=16)}
par(mar=c(0,0,0,0));plot.new();legend("center",c("First 80","Final 20"),pch=c(1,16),bty="n",cex=1.2);plot.new()
dev.off();layout(1)

## Independent 100 vs 50 parameter summaries
plot_design=function(p){x=design_comparison[design_comparison$parameter==p,]
 x=x[match(c("Primary\nN=100","Independent\nN=50"),x$design),];xp=c(1.25,1.75)
 plot(xp,x$median,pch=16,ylim=lims[[p]],xlim=c(1,2),xaxt="n",xlab="LHS design",ylab=labs[[p]],cex.axis=1.3,cex.lab=1.4)
 arrows(xp,x$lower,xp,x$upper,angle=90,code=3,length=.05);axis(1,at=xp,labels=c("Primary\nN=100","Independent\nN=50"),cex.axis=1.1)}
png(file.path(figure_dir,"Stage2_independent_LHS_parameter_stability.png"),3200,2200,res=300)
layout(matrix(c(1,2,3,8,4,5,6,7),2,byrow=T));par(mar=c(6,5,2,1))
for(p in pars[1:3])plot_design(p);for(p in pars[4:7])plot_design(p);par(mar=c(0,0,0,0));plot.new()
dev.off();layout(1)

## Independent-design pairwise
d=stage2_independent_pairwise;a=d[d$design=="Primary N=100",];b=d[d$design=="Independent N=50",]
png(file.path(figure_dir,"Stage2_independent_LHS_pairwise.png"),3600,2000,res=300)
layout(matrix(c(1,2,3,7,4,5,6,8),2,byrow=T),widths=c(1,1,1,.8));par(mar=c(4.5,5.5,2,1))
for(z in pairs){x=z[1];y=z[2];plot(a[[x]],a[[y]],pch=1,xlim=lims[[x]],ylim=lims[[y]],xlab=labs[[x]],ylab=labs[[y]],cex.axis=1.2,cex.lab=1.3);points(b[[x]],b[[y]],pch=16)}
par(mar=c(0,0,0,0));plot.new();legend("center",c("Primary (N=100)","Independent (N=50)"),pch=c(1,16),bty="n",cex=1.2);plot.new()
dev.off();layout(1)

## Sensitivity acceptance stability
ss=stage2_sensitivity_stability
yr=range(unlist(lapply(ss,function(x)c(x$acceptance,x$final_retention))),na.rm=T);ylim=c(0,yr[2]+.05*diff(yr))
plot_acc=function(x,title){plot(x$N,x$acceptance,type="b",ylim=ylim,xlim=range(x$N),xlab="Number of MDA-effect LHS samples",ylab="Proportion of parameter combinations",main=title,cex.axis=1.2,cex.lab=1.3,cex.main=1.2);lines(x$N,x$final_retention,type="b",lty=2)}
png(file.path(figure_dir,"Stage2_acceptance_stability_sensitivity_analyses.png"),2600,2200,res=300)
par(mfrow=c(2,2),mar=c(5,5,3,1));plot_acc(ss$Moderate,"A.  Moderate increase");legend("topright",c("Satisfy calibration criteria","Final retained"),lty=c(1,2),pch=1,bty="n",cex=.9);plot_acc(ss$Weak,"B.  Weak increase");plot_acc(ss$No_selection,"C.  No offspring-phenotype selection");plot_acc(ss$Wider_prevalence,"D.  Wider prevalence")
dev.off();par(mfrow=c(1,1))

## Retained parameter distributions across calibration analyses
plot_sens=function(p,show_xlab=FALSE){x=sensitivity_parameter_summary[sensitivity_parameter_summary$parameter==p,]
 ord=c("Primary","Moderate","Low","No selection","Wider prevalence");x=x[match(ord,x$scenario),];xp=1:5
 plot(xp,x$median,pch=16,ylim=lims[[p]],xlim=c(.5,5.5),xaxt="n",xlab="",ylab=labs[[p]],cex.axis=1.2,cex.lab=1.3,mgp=c(2.8,.8,0));arrows(xp,x$lower,xp,x$upper,angle=90,code=3,length=.04);axis(1,at=xp,labels=F,tck=-.02);yy=par("usr")[3]-.045*diff(par("usr")[3:4]);text(xp,yy,c("Primary","Moderate","Weak","No selection","Wider prevalence"),srt=45,adj=1,xpd=T,cex=.85);if(show_xlab)mtext("Calibration analysis",1,5.3,cex=1.3)}
png(file.path(figure_dir,"Stage2_sensitivity_parameter_distributions.png"),3800,1800,res=300)
layout(matrix(1:8,2,byrow=T));par(mar=c(5.5,6,1,1));for(p in pars[1:3])plot_sens(p)
par(mar=c(0,0,0,0));plot.new();legend("center",c("Median","2.5th/97.5th percentile"),lty=c(NA,1),pch=c(16,NA),bty="n",cex=1.1)
par(mar=c(7.5,6,1,1));for(p in pars[4:7])plot_sens(p,TRUE)
dev.off()
