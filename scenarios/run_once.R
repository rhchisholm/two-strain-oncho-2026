# Run the VR model once for a given set of parameters, and plot some output

library(foreach)
library(grid)
library(epiR)
library(lhs)
library(deSolve)

# load equations (define functions) and fixed parameters
source("Single_scenario.R")

# SCENARIO parameters
MDAi<-1 # I: MDAi=1, M: MDAi=0
Dburnin<-300 # duration of burnin
ABRall<-c(27000) # annual biting rates
MDAfall<-c(1) # time between MDA rounds (years)
MDAnall<-c(20) # total number of MDA rounds
ABR<-ABRall[1]
MDAf <-MDAfall[1]
NumberMDAs <- MDAnall[1]

my_CGR_GG<-0.95
my_CGR_SS<-0.2
my_CGR_SG<-0.75
my_dom_c_SS<-9
my_dom_c_SG<-9
my_ratio<-5 # relative rate of loss of embryostatic effect ivermectin in SOR compared to GR worms
my_r<-0.55 
my_zeta<-0.3 # set cumulative effect of ivermectin on GR worm fertility here

PostParms<-matrix(c(my_CGR_GG, my_CGR_SS, my_CGR_SG,
                    my_dom_c_SS, my_dom_c_SG, my_ratio, my_r), 
                  nrow=1, ncol=7)

filename<-"test.Rdata"
if(my_zeta!=0){
  single_VR(PostParms,1,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin,my_zeta)
}else{
  single_VR_noc(PostParms,1,filename,MDAf,NumberMDAs,ABR,MDAi,Dburnin)
}

load("test.Rdata")
yout<-OutputM[[1]]

source("../fitting/fixed-parameters.R")
source("../fitting/equations.R")


tt <- yout[,1]
RecSOR_proportion_worms <- 1 - (yout[,3]+yout[,4])/(yout[,3]+yout[,4]+yout[,7]+yout[,8])
RecSOR_proportion_mf <- 1 - (yout[,5])/(yout[,5]+yout[,9])
RecPrevalenceGR <- prevalence(yout[,5],my_k0,my_k1)
RecPrevalenceSOR <- prevalence(yout[,9],my_k0,my_k1)
RecPrevalence <- prevalence(yout[,9]+yout[,5],my_k0,my_k1)

data <- c(tt,RecSOR_proportion_worms,RecSOR_proportion_mf,RecPrevalenceGR,RecPrevalenceSOR, RecPrevalence)
data <- matrix(data, nrow=length(RecSOR_proportion_worms),length(data))
dft <- as.data.frame(data)
names(dft) <- c("time","SOR_worm_proportion","SOR_mf_proportion","Prevalence_GR","Prevalence_SOR","Prevalence")

f <- ggplot(dft, aes(x = time)) + 
  geom_line(aes(y=Prevalence, color = "Total")) +
  geom_line(aes(y = Prevalence_GR, color = "GR")) +
  geom_line(aes(y = Prevalence_SOR, color = "SOR")) +
  labs(color = "") +
  xlab("Time (years)") +
  ylab("mf prevalence (%)") + theme_bw() + theme(legend.position = c(0.9,0.85))

f

g <- ggplot(dft, aes(x = time)) + 
  geom_line(aes(y=SOR_worm_proportion, color = "worms")) +
  geom_line(aes(y = SOR_mf_proportion, color = "mf")) +
  labs(color = "") +
  xlab("Time (years)") +
  ylab("Proportion SOR") + theme_bw() + theme(legend.position = c(0.9,0.85))

g