###############################################################################
## Asubende calibration sensitivity figure: repository version
## Requires stability_data/Asubende_sensitivity_summary.csv only.
###############################################################################
figure_dir <- "stability_figures"
dir.create(figure_dir,showWarnings=FALSE,recursive=TRUE)
asubende_summary <- read.csv("stability_data/Asubende_sensitivity_summary.csv",
                             stringsAsFactors=FALSE)
scenario_order <- c("Primary","Moderate","Weak","No selection","Wider prevalence")

plot_asubende_outcome <- function(outcome,label,ylim=NULL,show_xlab=FALSE){
 x=asubende_summary[asubende_summary$outcome==outcome,,drop=FALSE]
 x=x[match(scenario_order,x$scenario),,drop=FALSE]; xpos=1:5
 if(is.null(ylim)){
   yr=range(c(x$q025,x$q975),na.rm=TRUE); pad=.05*diff(yr)
   ylim=c(max(0,yr[1]-pad),yr[2]+pad)
 }
 plot(xpos,x$median,pch=16,ylim=ylim,xlim=c(.5,5.5),xaxt="n",xlab="",ylab=label,cex.axis=1.2,cex.lab=1.3)
 arrows(xpos,x$q025,xpos,x$q975,angle=90,code=3,length=.04)
 axis(1,at=xpos,labels=FALSE,tck=-.02)
 label_y=par("usr")[3]-.045*diff(par("usr")[3:4])
 text(xpos,label_y,labels=scenario_order,srt=45,adj=1,xpd=TRUE,cex=.85)
 if(show_xlab)mtext("Calibration analysis",side=1,line=5.3,cex=1.2)
}
png(file.path(figure_dir,"Asubende_sensitivity_outcomes.png"),3300,2200,res=300)
layout(matrix(1:6,nrow=2,byrow=TRUE));par(mar=c(7.5,6,2.5,1))
plot_asubende_outcome("final_prev_annual","Final mf prevalence",c(0,.15));mtext("A.  End of annual MDAi calibration",3,.5,cex=1.1)
plot_asubende_outcome("increase_sor_annual","Increase in SOR worm proportion");mtext("B.  End of annual MDAi calibration",3,.5,cex=1.1)
plot_asubende_outcome("final_sor_annual","Final SOR worm proportion",c(0,.6));mtext("C.  End of annual MDAi calibration",3,.5,cex=1.1)
par(mar=c(8.5,6,2.5,1))
plot_asubende_outcome("final_prev_biannual","Final mf prevalence",show_xlab=TRUE);mtext("D.  End of biannual MDAi projection",3,.5,cex=1.1)
plot_asubende_outcome("further_increase_sor_biannual","Further increase in SOR worm proportion",show_xlab=TRUE);mtext("E.  End of biannual MDAi projection",3,.5,cex=1.1)
plot_asubende_outcome("final_sor_biannual","Final SOR worm proportion",c(0,.6),TRUE);mtext("F.  End of biannual MDAi projection",3,.5,cex=1.1)
dev.off()
