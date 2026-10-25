###############################################################################
## Stage 1 LHS stability figures: repository version
## Requires only data_processed/stage1_stability_repo_data.Rdata
###############################################################################
figure_dir <- "stability_figures"
dir.create(figure_dir, showWarnings=FALSE, recursive=TRUE)
load("data_processed/stage1_stability_repo_data.Rdata")

## 1. Acceptance stability
png(file.path(figure_dir,"Stage1_acceptance_stability.png"),1800,1500,res=300)
plot(stage1_acceptance_stability$N,stage1_acceptance_stability$acceptance_rate,
     type="b",xlim=c(200,1000),ylim=c(0,0.2),
     xlab="Number of inheritance LHS samples",
     ylab="Proportion of parameter sets accepted",cex.axis=1.3,cex.lab=1.4)
dev.off()

## 2. Parameter convergence
plot_convergence <- function(summary, parameter, label, ylim) {
  plot(summary$N,summary[[paste0(parameter,"_med")]],type="b",ylim=ylim,
       xlab="Number of LHS samples",ylab=label,cex.axis=1.3,cex.lab=1.4)
  lines(summary$N,summary[[paste0(parameter,"_025")]],type="b",lty=2)
  lines(summary$N,summary[[paste0(parameter,"_975")]],type="b",lty=2)
}
png(file.path(figure_dir,"Stage1_parameter_convergence.png"),3600,1000,res=300)
layout(matrix(1:4,nrow=1),widths=c(1,1,1,0.7)); par(mar=c(4.5,6.5,2,1))
plot_convergence(parameter_summary,"p1",expression(c[gg]),c(0.9,1))
plot_convergence(parameter_summary,"p2",expression(c[ss]^0),c(0,1))
plot_convergence(parameter_summary,"p3",expression(c[sg]^0),c(0,1))
par(mar=c(0,0,0,0)); plot.new()
legend("center",legend=c("Median","2.5th/97.5th percentile"),lty=c(1,2),pch=c(1,1),bty="n")
dev.off(); layout(1)

## 3. First 800 versus final 200
d <- stage1_primary_pairwise
early <- d$sampling_block=="First 800"; late <- !early
png(file.path(figure_dir,"Stage1_pairwise_parameter_space.png"),3600,1000,res=300)
layout(matrix(1:4,nrow=1),widths=c(1,1,1,0.7)); par(mar=c(4.5,6.5,2,1))
plot(d$c_gg[early],d$c0_ss[early],pch=1,xlim=c(.9,1),ylim=c(0,1),
     xlab=expression(c[gg]),ylab=expression(c[ss]^0),cex.axis=1.3,cex.lab=1.4)
points(d$c_gg[late],d$c0_ss[late],pch=16)
plot(d$c_gg[early],d$c0_sg[early],pch=1,xlim=c(.9,1),ylim=c(0,1),
     xlab=expression(c[gg]),ylab=expression(c[sg]^0),cex.axis=1.3,cex.lab=1.4)
points(d$c_gg[late],d$c0_sg[late],pch=16)
plot(d$c0_ss[early],d$c0_sg[early],pch=1,xlim=c(0,1),ylim=c(0,1),
     xlab=expression(c[ss]^0),ylab=expression(c[sg]^0),cex.axis=1.3,cex.lab=1.4)
points(d$c0_ss[late],d$c0_sg[late],pch=16)
par(mar=c(0,0,0,0));plot.new();legend("center",c("First 800","Final 200"),pch=c(1,16),bty="n")
dev.off();layout(1)

## 4. Independent-design parameter summaries
plot_design <- function(parameter,label,ylim){
 x=stage1_parameter_summary[stage1_parameter_summary$parameter==parameter,]
 x=x[match(c("Primary (N=1000)","Independent 1 (N=250)","Independent 2 (N=250)"),x$design),]
 xp=c(1.4,2,2.6)
 plot(xp,x$median,pch=16,ylim=ylim,xlim=c(1,3),xaxs="i",xaxt="n",
      xlab="LHS design",ylab=label,cex.axis=1.3,cex.lab=1.4,mgp=c(2.7,.8,0))
 arrows(xp,x$lower,xp,x$upper,angle=90,code=3,length=.05)
 axis(1,at=xp,labels=c("Primary\nN=1000","Indep. 1\nN=250","Indep. 2\nN=250"),cex.axis=1.1)
}
png(file.path(figure_dir,"Stage1_independent_LHS_parameter_stability.png"),3200,1200,res=300)
par(mfrow=c(1,3),mar=c(6,6.5,2,1.5),oma=c(0,1,0,1))
plot_design("c_gg",expression(c[gg]),c(.9,1))
plot_design("c0_ss",expression(c[ss]^0),c(.01,1))
plot_design("c0_sg",expression(c[sg]^0),c(.01,1))
dev.off();par(mfrow=c(1,1))

## 5. Independent-design pairwise plots
d=stage1_independent_pairwise
a=d[d$design=="Primary",]; b=d[d$design=="Independent 1",]; c=d[d$design=="Independent 2",]
rng <- function(v) range(c(a[[v]],b[[v]],c[[v]]),na.rm=TRUE)
png(file.path(figure_dir,"Stage1_independent_LHS_pairwise.png"),3600,1200,res=300)
layout(matrix(1:4,nrow=1),widths=c(1,1,1,.8));par(mar=c(4.5,4.5,2,1))
plot(a$c_gg,a$c0_ss,pch=1,xlim=rng("c_gg"),ylim=rng("c0_ss"),xlab=expression(c[gg]),ylab=expression(c[ss]^0),cex.axis=1.3,cex.lab=1.4);points(b$c_gg,b$c0_ss,pch=16);points(c$c_gg,c$c0_ss,pch=2)
plot(a$c_gg,a$c0_sg,pch=1,xlim=rng("c_gg"),ylim=rng("c0_sg"),xlab=expression(c[gg]),ylab=expression(c[sg]^0),cex.axis=1.3,cex.lab=1.4);points(b$c_gg,b$c0_sg,pch=16);points(c$c_gg,c$c0_sg,pch=2)
plot(a$c0_ss,a$c0_sg,pch=1,xlim=rng("c0_ss"),ylim=rng("c0_sg"),xlab=expression(c[ss]^0),ylab=expression(c[sg]^0),cex.axis=1.3,cex.lab=1.4);points(b$c0_ss,b$c0_sg,pch=16);points(c$c0_ss,c$c0_sg,pch=2)
par(mar=c(0,0,0,0));plot.new();legend("center",c("Primary (N=1000)","Independent 1 (N=250)","Independent 2 (N=250)"),pch=c(1,16,2),bty="n",cex=1.2)
dev.off()
