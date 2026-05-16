# Creates Figure 2 (and corresponding figure for the sensitivity analysis)
library(ggplot2)
library(patchwork)
library(reshape2)

if (!dir.exists("figures")) {
  dir.create("figures")
}

# load equations and parameters
source("equations.R")
MDAi<-1
ABR<-65000
source("fixed-parameters.R")

# load stats and params from 1000 LHS samples
load("data_files/fitting_burnin_65000_UU.Rdata")

# Determine simulations that match criteria 1
indv2<-which(l[[2]][,1]-l[[2]][,3]>0 & l[[2]][,3]-l[[2]][,2]>0 & l[[1]][,1]< 0.15 & l[[1]][,1] > 0.05)
num_samples<-length(indv2)

# loads OutputB (posterior burnin data) 
load("data_files/burnin_157_ABR_65000_UU.Rdata")

for(sa in 1:2){
  
  if(sa==1){
    filename1<-"data_files/Final_posterior_burnins_plusbiannualMDAs_cumu.Rdata"
    filename2<-"data_files/posteriors_ABR_65000_UU_cumu.Rdata"
    filename3<-"data_files/1patch2phenotypeALLdata.Rdata"
    filename4<-"data_files/accepted_trajectories.Rdata"
  }else{
    filename1<-"data_files/Final_posterior_burnins_plusbiannualMDAs_cumuSA.Rdata"
    filename2<-"data_files/posteriors_ABR_65000_UU_cumuSA.Rdata"
    filename3<-"data_files/1patch2phenotypeALLdataSA.Rdata"
    filename4<-"data_files/accepted_trajectoriesSA.Rdata"
  }
  
  # loads OutputMM (yout for biannual period)
  load(filename1)
  
  # loads posteriors (stats and params after annual MDA)
  load(filename2)
  posteriors <- posteriors[posteriors[,3]<0.5 & posteriors[,4]>0.5,]
  num_sims<-length(posteriors[,1])
  
  #########################################################
  ## Process output from biannual MDAi simulations ########
  #########################################################
  
  RecFinal <- matrix(0,nrow=num_sims,ncol=6)
  
  for(n in 1:num_sims){
    
    yout<-OutputMM[[n]]
    
    RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
    RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
    RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
    RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
    RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
    
    RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
    RecFinal_dSOR_proportion_worms <- RecFinal_SOR_proportion_worms - 
      RecSOR_proportion_worms[1]
    RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
    RecFinalPrevalence<-tail(RecPrevalence,n=1)
    RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
    RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
    
    RecFinal[n,] <- c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                      RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence, 
                      RecFinal_dSOR_proportion_worms)
  }
  
  # Turn matrix into data frame so we can plot using ggplot
  columns<-c('SampleID','c_GRGR','c0_SORSOR','c0_GRSOR', 
             's_SORSOR', 's_SORGR','psi_SOR/psi_GR','y',
             'Final_mf_prevalence_annual','Increase_in_SOR_worm_proportion', 'Final_SOR_worm_proportion', # annual
             'Final_mf_prevalence_biannual','Final_SOR_worm_proportion_biannual_MDAi', # biannual
             'Initial mf prevalence', 'Initial SOR worm proportion', # pre annual MDA
             'Further_increase_in_SOR_worm proportion')  # end biannual - end annual
  
  ALLdata<- data.frame(matrix(nrow=0, ncol=length(columns))) 
  names(ALLdata) <- columns
  
  for(ii in 1:num_sims){
    
    sampleBI<-posteriors[ii,1]
    youtb<-OutputB[[sampleBI]]
    endstep<-length(OutputB[[1]][,1])
    Final_SOR_W_ProportionBI <- 1 - (youtb[endstep,3]+youtb[endstep,4])/
      (youtb[endstep,3]+youtb[endstep,4]+youtb[endstep,7]+youtb[endstep,8])
    Final_prevBI <- prevalence(youtb[endstep,9]+youtb[endstep,5],my_k0,my_k1)
    
    newrow<-c(posteriors[ii,],RecFinal[ii,5],RecFinal[ii,1], 
              Final_prevBI, Final_SOR_W_ProportionBI, RecFinal[ii,6])
    ALLdata[nrow(ALLdata) + 1, ]<-newrow
      
  }
   
  # Data from fitting pre MDA
  columns2<- c('c_GR,GR','c0_SOR,SOR','c0_GR,SOR', 
            'Initial_mf_Prevalence', 'Initial_SOR_worm_proportion')
  PreMDAdata<- data.frame(matrix(nrow=0, ncol=length(columns2))) 
  names(PreMDAdata) <- columns2
  
  for(jj in 1:num_samples){
   
   youtb<-OutputB[[jj]]
   endstep<-length(OutputB[[1]][,1])
   Final_SOR_W_ProportionBI <- 1 - (youtb[endstep,3]+youtb[endstep,4])/
     (youtb[endstep,3]+youtb[endstep,4]+youtb[endstep,7]+youtb[endstep,8])
   Final_prevBI <- prevalence(youtb[endstep,9]+youtb[endstep,5],my_k0,my_k1)
   
   samp<-indv2[jj]
   C_GG<-l[[2]][samp,1]
   C_SS<-l[[2]][samp,2]
   C_SG<-l[[2]][samp,3]
   
   newrow<-c(C_GG,C_SS,C_SG, 
             Final_prevBI, Final_SOR_W_ProportionBI)
   PreMDAdata[nrow(PreMDAdata) + 1, ]<-newrow
   
  }
  
  if (file.exists(filename3)==FALSE){
    save(ALLdata, PreMDAdata, file = filename3) 
  }
  
  data=subset(ALLdata, c0_SORSOR<0.5 & c0_GRSOR>0.5)
  
  if (file.exists(filename4)==FALSE){
    save(data, file = filename4)
  }
  
  # Load GR-only model data
  load("data_files/GR_only_model_final_stats.Rdata")
  
  Final_GRprevBI <- prevalence(GR_only_model[1,4],my_k0,my_k1)
  Final_GRprevA <- prevalence(GR_only_model[2,4],my_k0,my_k1)
  Final_GRprevBA <- prevalence(GR_only_model[3,4],my_k0,my_k1)
  
  datan <- data
  names(datan)[9] <- "After annual MDAi"
  names(datan)[12] <- "After biannual MDAi"
  
  # Box plots (panel B)
  f2b <- ggplot(stack(datan[c(9,12)]*100), aes(x = ind, y = values)) + 
    geom_boxplot() +
    labs(x=NULL,
         y="Final mf prevalence (%)") +
    ylim(0,0.15*100) +
    annotate("segment", x = 0.5, xend = 1.5, y = Final_GRprevA*100, yend = Final_GRprevA*100,
             colour = "red", linetype = "solid") +
    annotate("segment", x = 1.5, xend = 2.5, y = Final_GRprevBA*100, yend = Final_GRprevBA*100,
             colour = "red", linetype = "dotdash") +
    annotate("segment", x = 0.5, xend = 2.5, y = 0.014*100, yend = 0.014*100,
             colour = "blue", linetype = "dashed") + theme_bw() + 
    annotate("text", label = "pOTTIS", x = 0.65, y = 0.02*100, color = "blue") +
    theme(text = element_text(size = 18))
  
  if(sa==1){
    f2bo<-f2b
  }else{
    f2bsa<-f2b
  }
  
  # Time series (panel A)
  row_indices <- which(ALLdata$`c0_SOR,SOR`<0.5 & ALLdata$`c0_GR,SOR`>0.5)
  nt<-nrow(posteriors)
  ntt<-nrow(OutputMM[[1]])
  
  StoreTrajectories<-matrix(1:(ntt*nt), nrow=ntt)
  
  for(i in 1:nt){
    index<-row_indices[i]
    if(nt==length(OutputMM)){
      yout<-OutputMM[[i]]}
    else{
      yout<-OutputMM[[index]]
    }
    Prevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
    StoreTrajectories[,i]<-Prevalence
    
  }
  dataf<- as.data.frame(StoreTrajectories[,]*100)
  dataf$time<-yout[,1]
  final_data <- melt(dataf, id='time')
  names(final_data) <- c('id','func', 'value')
  
  load("data_files/GR_only_model_full_trajectory.Rdata")
  youtgr<-yout_biannual
  Prevalencegr<-prevalence(youtgr[,9]+youtgr[,5],my_k0,my_k1)*100
  datafg<-as.data.frame(Prevalencegr)
  datafg$time<-youtgr[,1]
  final_datag <- melt(datafg, id='time')
  names(final_datag) <- c('id','func', 'value')
  
  f2a<- ggplot() + geom_line(data = final_datag, aes(x = id+1987, y = value, color = "red")) + 
    geom_line(data = final_data, aes(x = id+1987, y = value, group = func), alpha=0.0075) + xlab("Time (years)") +
    ylab("mf prevalence (%)") + theme_bw() + 
    annotate("text", label = "Variable response model", x = 32.3+1987, y = 0.75*100, color = "black", size=6) +
    annotate("text", label = "GR-only model", x = 29.6+1987, y = 0.68*100, color = "black", size=6) +
    annotate("segment", x = 23+1987, xend = 25.2+1987, y = 0.68*100, yend = 0.68*100,colour = "red") +
    annotate("segment", x = 23+1987, xend = 25.2+1987, y = 0.75*100, yend = 0.75*100,colour = "black") +
    annotate("segment", x = 2007, xend = 2007, y = 0, yend = 80,colour = "darkcyan") +
    annotate("text", x = 2008.1, y = 40, color = "darkcyan", label = "Fitting  Prediction", size=6) +
    annotate("segment", x = 2014, y = 40, xend=2017, yend=40, color = "darkcyan", arrow = arrow(length = unit(0.2, "cm"), type = "closed")) +
    theme(legend.position = "none", text = element_text(size = 18))
  
  if(sa==1){
    f2ao<-f2a
    patchworko <- f2ao+f2bo
  }else{
    f2asa<-f2a
    patchworksa <- f2asa+f2bsa
  }
  
}

# Figure 2
p <- patchworko + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2, 1.25))
ggsave("figures/F2.png", plot = p, width = 15, height = 5, dpi=300)

# Analogous supplementary figure from sensitivity analysis
p <- patchworksa + plot_annotation(tag_levels = 'A') + plot_layout(widths = c(2, 1.25))
ggsave("figures/F2_SA.png", plot = p, width = 15, height = 5, dpi=300)