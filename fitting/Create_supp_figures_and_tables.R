# Create Figures S1, S2, S3, S4, S6, S8 and matrices storing statistics
library(ggplot2)
library(latex2exp)
library(patchwork)

##-----------------------------------------------------
## Figure S1 ## Prevalence and MDA Coverage in Asubende
##-----------------------------------------------------

# Real data plots
columns1 <- c('Year','Coverage','Missing')
columns2 <- c('Year','Prevalence','Missing')

Realdata1<- data.frame(matrix(nrow=0, ncol=length(columns1))) 
names(Realdata1) <- columns1
Realdata2<- data.frame(matrix(nrow=0, ncol=length(columns2))) 
names(Realdata2) <- columns2

my_coverageA<-c(65,65,65,65,65,
                65,65,65,65,65,
                65,81.7,80.1,80,84.10,
                84.10,74.70,79.40,87.10,86.20)/100
my_coverageB<-c(92.00,92, 87.30,88.00,0.00,74.20,
                50.60,79.60,81.90,78.30,83.00,
                79.50,82.10,78.70,73.10,77.50,
                77.10,84.30)/100
my_coverageALL<-c(my_coverageA,my_coverageB)

cov.time<-seq(1987,2024,by=1)
prev.time<-c(1980,1987,1987,1987,1988,1988,1989,1989,1989,1990,1990,1990,1991,
             1992,1992,1992,1993,1994,1994,1995,1996,1997,1998,1999,2000,2001,
             2002,2003,2004,2004,2005,2006,2007,2007,2008,2009,2010,2012,2013,
             2014,2015,2016,2017,2018,2019,2020,2021,2022,2023,2024)

prevalence<-c(71.61,88.00,82.61,79.41,69.18,77.35,79.50,65.21,72.75,77.20,64.24,
              71.34,67.67,67.10,52.94,62.11,NA,NA,38.48,NA,NA,25.56,NA,NA,4.43,NA,
              NA,NA,NA,7.69,NA,NA,NA,11.36,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,NA,
              NA,NA,NA,NA)

for(i in 1:length(cov.time)){
  cov<-my_coverageALL[i]
  missing<-0
  if(cov==0.65){
    missing<-1}
  if(cov==0.8){
    missing<-1
  }
  time<-cov.time[i]
  newrow<-c(time, cov, missing)
  Realdata1[nrow(Realdata1) + 1, ]<-newrow
}

Realdata1$Missing<-as.factor(Realdata1$Missing)

for(i in 1:length(prev.time)){
  missing<-0
  time<-prev.time[i]
  prev<-prevalence[i]/100
  newrow<-c(time, prev,missing)
  Realdata2[nrow(Realdata2) + 1, ]<-newrow
  
}

Realdata2$Missing<-as.factor(Realdata2$Missing)

# Figure S1a
c1 <- ggplot(Realdata1, aes(x = Year, y = Coverage, color=Missing)) + 
  geom_point() +
  scale_color_manual(
    values = c("1" = "#EE3377", "0" = "#009988"), # Assign specific colors
    labels = c("0" = "Measured", "1" = "Estimated")) + # Assign custom labels
  labs(x="Year",
       y="MDAi coverage (proportion)", title=NULL) + xlim(1987,2024) +
  labs(color = "Data source") +
  ylim(0,1) + theme_bw() + 
  theme(text = element_text(size = 14))

# Figure S1b
c2 <- ggplot(Realdata2, aes(x = Year, y = Prevalence, color=Missing)) + 
  geom_point() +
  scale_color_manual(
    values = c("1" = "#EE3377", "0" = "#009988"), # Assign specific colors
    labels = c("0" = "Measured", "1" = "Estimated")) + # Assign custom labels
  labs(x="Year",
       y="mf prevalence", title=NULL) + xlim(1987,2024) +
  labs(color = "Data source") +
  ylim(0,1) + theme_bw() + 
  theme(legend.position = "none",text = element_text(size = 14))

dev.new()
c2+c1+plot_annotation(tag_levels = 'A') 

##-----------------------------------------------
## Load annual MDAi data from fitting simulations
##-----------------------------------------------

# Load burn ins and LHS simulation output
load("burnin_157_ABR_65000_UU.Rdata")
preMDA <- vector("list", 157)

load("fitting_burnin_65000_UU.Rdata")
indv2<-which(l[[2]][,1]-l[[2]][,3]>0 & l[[2]][,3]-l[[2]][,2]>0 & l[[1]][,1]< 0.15 & l[[1]][,1] > 0.05)
num_samples<-length(indv2)

# load equations and parameters
source("equations.R")
MDAi<-1
ABR<-65000
source("fixed-parameters.R")

# Load burnins
for(n in 1:num_samples){
  
  yout<-OutputB[[n]]
  
  samp<-indv2[n]
  my_CGR_GG<-l[[2]][samp,1]
  my_CGR_SS<-l[[2]][samp,2]
  my_CGR_SG<-l[[2]][samp,3]
  
  RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
  RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
  
  RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
  RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
  
  RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
  RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
  RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
  
  RecFinalPrevalence<-tail(RecPrevalence,n=1)
  RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
  RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
  
  preMDA[[n]] = c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                  RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence,
                  my_CGR_GG,my_CGR_SS,my_CGR_SG)
}

# Cumulative effect on
# Load LHS simulation data
load("fitting_MDA_ABR_65000_UU_from_posterior_samples_cumu.Rdata")
postMDA <- ll

# Calculate summary statistics for each LHS sample
dSOR <- vector("list", num_samples)
fPrev <- vector("list", num_samples)
fSOR <- vector("list", num_samples)
params <- vector("list", num_samples)

for(n in 1:num_samples){
  FP <- postMDA[[n]][,5]
  ISOR <- preMDA[[n]][1]
  fSOR[[n]] <- postMDA[[n]][,1]
  DiffSOR <- fSOR[[n]] - ISOR 
  dSOR[[n]] <- DiffSOR
  fPrev[[n]] <- FP
  dom_css <- postMDA[[n]][,6]
  dom_csg <- postMDA[[n]][,7]
  ratef <- postMDA[[n]][,8]
  r <- postMDA[[n]][,9]
  params[[n]] <- matrix(c(dom_css,dom_csg,ratef,r), ncol=4)
}

# Now save info for plots in dataframe 
FinalParams <- c()
FinalPrev <- c()
FinaldSOR <- c()
FinalSOR <- c()
BIsamplen <- c()
Cumu <- c()

for(n in 1:num_samples){
  for(i in 1:50){
      p <- c(preMDA[[n]][6:8],params[[n]][i,])
      FinalParams <-rbind(FinalParams, p)
      FinalPrev <- rbind(FinalPrev,fPrev[[n]][i])
      FinaldSOR <- rbind(FinaldSOR,dSOR[[n]][i])
      FinalSOR <- rbind(FinalSOR,fSOR[[n]][i])
      BIsamplen <- rbind(BIsamplen, n)
      Cumu <- rbind(Cumu, 1)
  }
}
columns<-c('SampleID','c_GRGR','c0_SORSOR','c0_GRSOR', 
           's_SOR,SOR', 's_SOR,GR','psi_SOR/psi_GR','y',
           'Final_mf_prevalence','Increase_in_SOR_worm_proportion','Cumulative_Effect','Final_SOR_worm_proportion')

posteriors<-matrix(c(BIsamplen,FinalParams,FinalPrev,FinaldSOR,Cumu,FinalSOR),ncol=12)
dfc <-as.data.frame(posteriors)
names(dfc) <- columns

# Cumulative effect off
# Load LHS simulation data
load("fitting_MDA_ABR_65000_UU_from_posterior_samples_no_cumu.Rdata")
postMDA <- ll

# Calculate summary statistics for each LHS sample
dSOR <- vector("list", num_samples)
fSOR <- vector("list", num_samples)
fPrev <- vector("list", num_samples)
params <- vector("list", num_samples)

for(n in 1:num_samples){
  FP <- postMDA[[n]][,5]
  ISOR <- preMDA[[n]][1]
  fSOR[[n]] <- postMDA[[n]][,1]
  DiffSOR <- fSOR[[n]] - ISOR 
  dSOR[[n]] <- DiffSOR
  fPrev[[n]] <- FP
  dom_css <- postMDA[[n]][,6]
  dom_csg <- postMDA[[n]][,7]
  ratef <- postMDA[[n]][,8]
  r <- postMDA[[n]][,9]
  params[[n]] <- matrix(c(dom_css,dom_csg,ratef,r), ncol=4)
}

# Now save info for plots in dataframe 
FinalParams <- c()
FinalPrev <- c()
FinaldSOR <- c()
FinalSOR <- c()
BIsamplen <- c()
Cumu <- c()

for(n in 1:num_samples){
  for(i in 1:50){
    p <- c(preMDA[[n]][6:8],params[[n]][i,])
    FinalParams <-rbind(FinalParams, p)
    FinalPrev <- rbind(FinalPrev,fPrev[[n]][i])
    FinaldSOR <- rbind(FinaldSOR,dSOR[[n]][i])
    FinalSOR <- rbind(FinalSOR,fSOR[[n]][i])
    BIsamplen <- rbind(BIsamplen, n)
    Cumu <- rbind(Cumu, 0)
  }
}
columns<-c('SampleID','c_GRGR','c0_SORSOR','c0_GRSOR', 
           's_SOR,SOR', 's_SOR,GR','psi_SOR/psi_GR','y',
           'Final_mf_prevalence','Increase_in_SOR_worm_proportion','Cumulative_Effect','Final_SOR_worm_proportion')

posteriorsn<-matrix(c(BIsamplen,FinalParams,FinalPrev,FinaldSOR,Cumu,FinalSOR),ncol=12)
dfnoc <-as.data.frame(posteriorsn)
names(dfnoc) <- columns

# Combine all data
dfall <- rbind(dfc, dfnoc)
dfall$Cumulative_Effect<-as.factor(dfall$Cumulative_Effect)
names(dfall) <- columns

##-----------------------------------------------------
## Figure S2 : VR and GR-only model calibration results
##-----------------------------------------------------

# Figure S1a
# Add a rectangle annotation
# Define the coordinates for the subset
x_min_subset <- 0.1
x_max_subset <- 0.12
y_min_subset <- 0.1
y_max_subset <- 0.2

y_min_subset2 <- 0.05
y_max_subset2 <- 0.1

data=subset(dfall, c0_SORSOR<0.5 & c0_GRSOR>0.5)

# VR model calibration, panel A
b5 <- ggplot(data, aes(x = Final_mf_prevalence, 
                       y = Increase_in_SOR_worm_proportion, color=factor(Cumulative_Effect))) +
  geom_point(alpha = 0.4) +
  scale_color_manual(
    values = c("0" = "#EE3377", "1" = "#009988"), # Assign specific colors
    labels = c("0" = "None", "1" = TeX(r"(Present with $\zeta=0.3$)")) # Assign custom labels
  ) + ggtitle('Variable response model') +
  labs(color = "Cumulative effect") +
  annotate("rect",
           xmin = x_min_subset, xmax = x_max_subset,
           ymin = y_min_subset, ymax = y_max_subset,
           fill = "blue", alpha = 0.2, color = "blue", linetype = "solid") +
  annotate("rect",
           xmin = x_min_subset, xmax = x_max_subset,
           ymin = y_min_subset2, ymax = y_max_subset2,
           fill = "#EE7733", alpha = 0.2, color = "#EE7733", linetype = "solid") +
  #ylim(0, 0.2) +
  labs(x="Final mf prevalence", y="Final increase in SOR worm proportion")  + theme_bw() + 
  theme(legend.position = "none", text = element_text(size = 14)) 

# GR-only model, panel B
load("GR_only_model_LHS_simulations_annual_no_cumu.Rdata")
PrevFinalnc<-l[[2]]
y<-l[[3]]
load("GR_only_model_LHS_simulations_annual_cumu.Rdata")
PrevFinalc<-l[[2]]

cumu<-(1:50)*0+1
posteriorsc<-matrix(c(y,PrevFinalc,cumu),ncol=3)
cumu<-(1:50)*0
posteriorsnc<-matrix(c(y,PrevFinalnc,cumu),ncol=3)
columns<-c('y','Final_mf_prevalence','Cumulative_effect')
df <-as.data.frame(posteriorsc)
dfn <-as.data.frame(posteriorsnc)
names(df) <- columns
names(dfn) <- columns
dfallgr <- rbind(df, dfn)
dfallgr$Cumulative_effect<-as.factor(dfallgr$Cumulative_effect)

b6 <- ggplot(dfallgr, aes(x = y, 
                          y = Final_mf_prevalence, color=factor(Cumulative_effect))) +
  geom_point(alpha = 0.8) +
  scale_color_manual(
    values = c("0" = "#EE3377", "1" = "#009988"), # Assign specific colors
    labels = c("0" = "None", "1" = TeX(r"(Present with $\zeta=0.3$)")) # Assign custom labels
  ) + ggtitle('GR-only model') +
  labs(color = "Cumulative effect") +
  labs(x="y", y="Final mf prevalence") +
  annotate("rect",
           xmin = 0, xmax = 1,
           ymin = 0.1, ymax = 0.12,
           fill = "blue", alpha = 0.2, color = "blue", linetype = "solid") + theme_bw() + 
  theme(text = element_text(size = 14))

dev.new()
b5+b6+plot_annotation(tag_levels = 'A') 

##------------------------------------------------------------------
## Figure S3 and S6 : VR model posterior distributions of parameters
##------------------------------------------------------------------

# Extract posterior distributions of parameters from accepted samples
# Without cumulative effect
datan=subset(dfall, Cumulative_Effect==0 & c0_SORSOR<0.5 & c0_GRSOR>0.5
            & Final_mf_prevalence>0.1 & Final_mf_prevalence<0.12 
            & Increase_in_SOR_worm_proportion>0.1)
# No accepted samples so don't plot results

# With cumulative effect
data=subset(dfall, Cumulative_Effect==1 & c0_SORSOR<0.5 & c0_GRSOR>0.5
             & Final_mf_prevalence>0.1 & Final_mf_prevalence<0.12 
             & Increase_in_SOR_worm_proportion>0.1)

# With cumulative effect sensitivity analysis
datas=subset(dfall, Cumulative_Effect==1 & c0_SORSOR<0.5 & c0_GRSOR>0.5
            & Final_mf_prevalence>0.1 & Final_mf_prevalence<0.12 
            & Increase_in_SOR_worm_proportion>0.05 & Increase_in_SOR_worm_proportion <0.1 )

## Figure S3, panel A
b1 <- ggplot(stack(data[c(2,3,4,8)]), aes(x = ind, y = values)) +
  geom_boxplot() + #geom_violin() + #
  labs(x=NULL, y=NULL, title="After annual MDAi program") +
  ylim(0, 1) + theme_bw() + 
  theme(text = element_text(size = 14))

## Figure S3, panel B
b2 <- ggplot(stack(data[c(5,6,7)]), aes(x = ind, y = values)) +
  geom_boxplot() + #geom_violin() + #
  labs(x=NULL, y=NULL, title="After annual MDAi program") +
  ylim(0, 20) + theme_bw() + 
  theme(text = element_text(size = 14))

dev.new()
b1+b2+plot_annotation(tag_levels = 'A') 

## Sensitivity analysis
## Figure S6, panel A
b1s <- ggplot(stack(datas[c(2,3,4,8)]), aes(x = ind, y = values)) +
  geom_boxplot() + #geom_violin() + #
  labs(x=NULL, y=NULL, title="After annual MDAi program") +
  ylim(0, 1) + theme_bw() + 
  theme(text = element_text(size = 14))

## Figure S6, panel B
b2s <- ggplot(stack(datas[c(5,6,7)]), aes(x = ind, y = values)) +
  geom_boxplot() + #geom_violin() + #
  labs(x=NULL, y=NULL, title="After annual MDAi program") +
  ylim(0, 20) + theme_bw() + 
  theme(text = element_text(size = 14))

dev.new()
b1s+b2s+plot_annotation(tag_levels = 'A') 

##------------------------------------------------------------------------------
## Figure S4 and S8 : Calibrated VR model, distributions of final mf prev & 
## increases in SOR worm proportion after annual and biannual MDA, Asubende
##------------------------------------------------------------------------------

# Create data frames summarizing simulations for accepted samples
columns<-c('SampleID','c_GRGR','c0_SORSOR','c0_GRSOR', 
           's_SORSOR', 's_SORGR','psi_SOR/psi_GR','y',
           'Final_mf_prevalence_annual','Increase_in_SOR_worm_proportion','Final_SOR_worm_proportion', # annual
           'Final_mf_prevalence_biannual','Final_SOR_worm_proportion_biannual_MDAi', # biannual
           'Initial_mf_prevalence', 'Initial_SOR_worm_proportion', # pre annual MDA
           'Further_increase_in_SOR_worm_proportion')  # end biannual - end annual

ALLdataAccepted<- data.frame(matrix(nrow=0, ncol=length(columns))) 
names(ALLdataAccepted) <- columns

# With cumulative effect, original fitting results
# Extract data from biannual MDAi
load("Final_posterior_burnins_plusbiannualMDAs_cumu.Rdata")
load("posteriors_ABR_65000_UU_cumu.Rdata")
posteriors <- posteriors[posteriors[,3]<0.5 & posteriors[,4]>0.5,]
num_sims<-length(posteriors[,1])
RecFinal <- matrix(0,nrow=num_sims,ncol=6)

for(n in 1:num_sims){
  
  yout<-OutputMM[[n]]
  
  RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
  RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
  RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
  RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
  RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
  
  RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
  RecFinal_dSOR_proportion_worms <- RecFinal_SOR_proportion_worms - RecSOR_proportion_worms[1] #whole program (a+b)

  RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
  RecFinalPrevalence<-tail(RecPrevalence,n=1)
  RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
  RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
  
  RecFinal[n,] <- c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                    RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence, 
                    RecFinal_dSOR_proportion_worms)
}

# Extract data from burnin and annual MDAi
for(ii in 1:num_sims){
  
  sampleBI<-posteriors[ii,1]
  youtb<-OutputB[[sampleBI]]
  endstep<-length(OutputB[[1]][,1])
  Final_SOR_W_ProportionBI <- 1 - (youtb[endstep,3]+youtb[endstep,4])/
    (youtb[endstep,3]+youtb[endstep,4]+youtb[endstep,7]+youtb[endstep,8])
  Final_prevBI <- prevalence(youtb[endstep,9]+youtb[endstep,5],my_k0,my_k1)
  IncreaseSORAB <- RecFinal[ii,1] - posteriors[ii,11]
  newrow<-c(posteriors[ii,],RecFinal[ii,5],RecFinal[ii,1], 
            Final_prevBI, Final_SOR_W_ProportionBI, IncreaseSORAB)
  ALLdataAccepted[nrow(ALLdataAccepted) + 1, ]<-newrow
  
}

ALLdataAccepted<-subset(ALLdataAccepted, c0_SORSOR<0.5 & c0_GRSOR>0.5)
ALLdataAcceptedO<-ALLdataAccepted

ALLdataAccepted<- data.frame(matrix(nrow=0, ncol=length(columns))) 
names(ALLdataAccepted) <- columns

# With cumulative effect, sensitivity analysis fitting results
# Extract data from biannual MDAi
load("Final_posterior_burnins_plusbiannualMDAs_cumuSA.Rdata")
load("posteriors_ABR_65000_UU_cumuSA.Rdata")
posteriors <- posteriors[posteriors[,3]<0.5 & posteriors[,4]>0.5,]
num_sims<-length(posteriors[,1])
RecFinal <- matrix(0,nrow=num_sims,ncol=6)

for(n in 1:num_sims){
  
  yout<-OutputMM[[n]]
  
  RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
  RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
  RecPrevalenceGR<-prevalence(yout[,5],my_k0,my_k1)
  RecPrevalenceSOR<-prevalence(yout[,9],my_k0,my_k1)
  RecPrevalence<-prevalence(yout[,9]+yout[,5],my_k0,my_k1)
  
  RecFinal_SOR_proportion_worms <- tail(RecSOR_proportion_worms,n=1)
  RecFinal_dSOR_proportion_worms <- RecFinal_SOR_proportion_worms - RecSOR_proportion_worms[1] #whole program (a+b)
  
  RecFinal_SOR_proportion_mf <- tail(RecSOR_proportion_mf,n=1)
  RecFinalPrevalence<-tail(RecPrevalence,n=1)
  RecFinalPrevalenceGR<-tail(RecPrevalenceGR,n=1)
  RecFinalPrevalenceSOR<-tail(RecPrevalenceSOR,n=1)
  
  RecFinal[n,] <- c(RecFinal_SOR_proportion_worms,RecFinal_SOR_proportion_mf,
                    RecFinalPrevalenceGR,RecFinalPrevalenceSOR,RecFinalPrevalence, 
                    RecFinal_dSOR_proportion_worms)
}

#  Extract data from burnin and annual MDAi
for(ii in 1:num_sims){
  
  sampleBI<-posteriors[ii,1]
  youtb<-OutputB[[sampleBI]]
  endstep<-length(OutputB[[1]][,1])
  Final_SOR_W_ProportionBI <- 1 - (youtb[endstep,3]+youtb[endstep,4])/
    (youtb[endstep,3]+youtb[endstep,4]+youtb[endstep,7]+youtb[endstep,8])
  Final_prevBI <- prevalence(youtb[endstep,9]+youtb[endstep,5],my_k0,my_k1)
  IncreaseSORAB <- RecFinal[ii,1] - posteriors[ii,11]
  newrow<-c(posteriors[ii,],RecFinal[ii,5],RecFinal[ii,1], 
            Final_prevBI, Final_SOR_W_ProportionBI, IncreaseSORAB)
  ALLdataAccepted[nrow(ALLdataAccepted) + 1, ]<-newrow
  
}

ALLdataAccepted<-subset(ALLdataAccepted, c0_SORSOR<0.5 & c0_GRSOR>0.5)
ALLdataAcceptedSA<-ALLdataAccepted

# Figure S4, panel A
b3 <- ggplot(stack(ALLdataAcceptedO[c(9,10)]), aes(x = ind, y = values)) + 
  geom_boxplot() +
  labs(x=NULL,
       y=NULL, title="After annual MDAi program") +
  ylim(0,0.3) + theme_bw() + 
  scale_x_discrete(labels = c("mf prevalence", "increase in\n proportion of SOR worms")) +
  theme(text = element_text(size = 14))

# Figure S4, panel B
b4 <- ggplot(stack(ALLdataAcceptedO[c(12,16)]), aes(x = ind, y = values)) + 
  geom_boxplot() +
  labs(x=NULL,
       y=NULL, title="After biannual MDAi program") +
  ylim(0,0.3) + theme_bw() + 
  scale_x_discrete(labels = c("mf prevalence", "further increase in\n proportion of SOR worms")) +
  theme(text = element_text(size = 14))

dev.new()
b3+b4+plot_annotation(tag_levels = 'A') 

# Figure S8, panel A
b3s <- ggplot(stack(ALLdataAcceptedSA[c(9,10)]), aes(x = ind, y = values)) + 
  geom_boxplot() +
  labs(x=NULL,
       y=NULL, title="After annual MDAi program") +
  ylim(0,0.3) + theme_bw() + 
  scale_x_discrete(labels = c("mf prevalence", "increase in\n proportion of SOR worms")) +
  theme(text = element_text(size = 14))

# Figure S8, panel B
b4s <- ggplot(stack(ALLdataAcceptedSA[c(12,16)]), aes(x = ind, y = values)) + 
  geom_boxplot() +
  labs(x=NULL,
       y=NULL, title="After biannual MDAi program") +
  ylim(0,0.3) + theme_bw() + 
  scale_x_discrete(labels = c("mf prevalence", "further increase in\n proportion of SOR worms")) +
  theme(text = element_text(size = 14))

dev.new()
b3s+b4s+plot_annotation(tag_levels = 'A')

##--------------------------------------------------------
## Create tables of summary stats for Asubende simulations
##--------------------------------------------------------

param_names<-c('c_GRGR',
               'c0_SORSOR',
               'c0_GRSOR',
               's_SOR,SOR',
               's_SOR,GR',
               'psi_SOR/psi_GR',
               'y')

stats_names<-c('Final_mf_prevalence_annual',
               'Increase_in_SOR_worm_proportion_annual',
               'Final_SOR_worm_proportion_annual',
               'Final_mf_prevalence_biannual',
               'Further_increase_in_SOR_worm_proportion_biannual',
               'Final_SOR_worm_proportion_biannual')

ParameterStats<-matrix(ncol = 5, nrow = 7, 
                       dimnames = list(param_names,
                                       c("0.05", "0.25", "0.5", "0.75","0.95")))
SummaryStats<-matrix(ncol = 5, nrow = 6, 
                     dimnames = list(stats_names,
                                     c("0.05", "0.25", "0.5", "0.75","0.95")))

# Original model, with cumulative effect
# Create matrix of posterior distributions of parameters
ParameterStatsO<-ParameterStats
for(i in 1:7){
  ParameterStatsO[i,]<-c(quantile(ALLdataAcceptedO[,i+1], 
                                 probs=c(0.05,0.25,0.5,0.75,0.95)))
}

# Create matrix of posterior distributions of summary stats
SummaryStatsO<-SummaryStats
inds<-c(9,10,11,12,16,13)
for(i in 1:6){
  SummaryStatsO[i,]<-c(quantile(ALLdataAcceptedO[,inds[i]], 
                                 probs=c(0.05,0.25,0.5,0.75,0.95)))
}

save(ParameterStatsO,SummaryStatsO,file = "Tables_Asubende.Rdata")

# Sensitivity analysis model, with cumulative effect
# Create matrix of posterior distributions of parameters
ParameterStatsSA<-ParameterStats
for(i in 1:7){
  ParameterStatsSA[i,]<-c(quantile(ALLdataAcceptedSA[,i+1], 
                                  probs=c(0.05,0.25,0.5,0.75,0.95)))
}

# Create matrix of posterior distributions of summary stats
SummaryStatsSA<-SummaryStats
inds<-c(9,10,11,12,16,13)
for(i in 1:6){
  SummaryStatsSA[i,]<-c(quantile(ALLdataAcceptedSA[,inds[i]], 
                                probs=c(0.05,0.25,0.5,0.75,0.95)))
}

save(ParameterStatsSA,SummaryStatsSA,file = "Tables_AsubendeSA.Rdata")