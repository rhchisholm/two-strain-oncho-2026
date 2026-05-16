# create tables summmarising different endemic prevalence scenarios

load("processed_outputs.RData")
ABRall<-c(2942,5350,27000)

if (!dir.exists("tables")) {
  dir.create("tables")
}

################################################################################
################################################################################
# Final prevalence tables 
################################################################################
#Table S3
S3Annual<-matrix(ncol = 5, nrow = 3, 
                         dimnames = list(c("30%", "50%","70%"),
                                         c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==1)
  S3Annual[i,]<-c(quantile(data$Final_Prevalence, 
                                     probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S3Biannual<-matrix(ncol = 5, nrow = 3, 
                         dimnames = list(c("30%", "50%","70%"),
                                         c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==1)
  S3Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                                   probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S4
S4Annual<-matrix(ncol = 5, nrow = 3, 
                         dimnames = list(c("30%", "50%","70%"),
                                         c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==0)
  S4Annual[i,]<-c(quantile(data$Final_Prevalence, 
                                   probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S4Biannual<-matrix(ncol = 5, nrow = 3, 
                           dimnames = list(c("30%", "50%","70%"),
                                           c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==0)
  S4Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                                     probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S5

S5Annual<-matrix(ncol = 5, nrow = 3, 
                         dimnames = list(c("30%", "50%","70%"),
                                         c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==1)
  S5Annual[i,]<-c(quantile(data$Final_Prevalence, 
                                     probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S5Biannual<-matrix(ncol = 5, nrow = 3, 
                           dimnames = list(c("30%", "50%","70%"),
                                           c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==1)
  S5Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                                       probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
#Table S6

S6Annual<-matrix(ncol = 5, nrow = 3, 
                         dimnames = list(c("30%", "50%","70%"),
                                         c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==0)
  S6Annual[i,]<-c(quantile(data$Final_Prevalence, 
                                   probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S6Biannual<-matrix(ncol = 5, nrow = 3, 
                           dimnames = list(c("30%", "50%","70%"),
                                           c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==0)
  S6Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                                     probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
################################################################################
# Time to pOTTIS
################################################################################
#Table S7
S7Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==1 & is.na(pOTTIS)==0)
  S7Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S7Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==1 & is.na(pOTTIS)==0)
  S7Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S8
S8Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==0 & is.na(pOTTIS)==0)
  S8Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S8Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==0 & is.na(pOTTIS)==0)
  S8Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S9

S9Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==1 & is.na(pOTTIS)==0)
  S9Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S9Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==1 & is.na(pOTTIS)==0)
  S9Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
#Table S10

S10Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==0 & is.na(pOTTIS)==0)
  S10Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S10Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==0 & is.na(pOTTIS)==0)
  S10Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
################################################################################
# Actual SOR worm proportion tables
################################################################################
#Table S11
S11Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==1 )
  S11Annual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S11Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==1 )
  S11Biannual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}



################################################################################
#Table S12

S12Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==1 )
  S12Annual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S12Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==1 )
  S12Biannual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}



################################################################################
################################################################################
# Final prevalence tables Sensitivity Analysis
################################################################################
df<-dfSA

#Table S14
S14Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==1)
  S14Annual[i,]<-c(quantile(data$Final_Prevalence, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S14Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==1)
  S14Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S15
S15Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==0)
  S15Annual[i,]<-c(quantile(data$Final_Prevalence, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S15Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==0)
  S15Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S16

S16Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==1)
  S16Annual[i,]<-c(quantile(data$Final_Prevalence, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S16Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==1)
  S16Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
#Table S17

S17Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==0)
  S17Annual[i,]<-c(quantile(data$Final_Prevalence, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S17Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==0)
  S17Biannual[i,]<-c(quantile(data$Final_Prevalence, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
################################################################################
# Time to pOTTIS Sensitivity Analysis
################################################################################
#Table S18
S18Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==1 & is.na(pOTTIS)==0)
  S18Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S18Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==1 & is.na(pOTTIS)==0)
  S18Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S19
S19Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==0 & is.na(pOTTIS)==0)
  S19Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S19Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==0 & is.na(pOTTIS)==0)
  S19Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}


################################################################################
#Table S20

S20Annual<-matrix(ncol = 5, nrow = 3, 
                 dimnames = list(c("30%", "50%","70%"),
                                 c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==1 & is.na(pOTTIS)==0)
  S20Annual[i,]<-c(quantile(data$pOTTIS, 
                           probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S20Biannual<-matrix(ncol = 5, nrow = 3, 
                   dimnames = list(c("30%", "50%","70%"),
                                   c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==1 & is.na(pOTTIS)==0)
  S20Biannual[i,]<-c(quantile(data$pOTTIS, 
                             probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
#Table S21

S21Annual<-matrix(ncol = 5, nrow = 3, 
                  dimnames = list(c("30%", "50%","70%"),
                                  c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==0 & is.na(pOTTIS)==0)
  S21Annual[i,]<-c(quantile(data$pOTTIS, 
                            probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S21Biannual<-matrix(ncol = 5, nrow = 3, 
                    dimnames = list(c("30%", "50%","70%"),
                                    c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==0 & is.na(pOTTIS)==0)
  S21Biannual[i,]<-c(quantile(data$pOTTIS, 
                              probs=c(0.05,0.25,0.5,0.75,0.95)))
}

################################################################################
################################################################################
# Actual SOR worm proportion tables Sensitivity Analysis
################################################################################
#Table S22
S22Annual<-matrix(ncol = 5, nrow = 3, 
                  dimnames = list(c("30%", "50%","70%"),
                                  c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==1 & SOR_present==1 )
  S22Annual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                            probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S22Biannual<-matrix(ncol = 5, nrow = 3, 
                    dimnames = list(c("30%", "50%","70%"),
                                    c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==1 & SOR_present==1 )
  S22Biannual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                              probs=c(0.05,0.25,0.5,0.75,0.95)))
}



################################################################################
#Table S23

S23Annual<-matrix(ncol = 5, nrow = 3, 
                  dimnames = list(c("30%", "50%","70%"),
                                  c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==1 & MDAi==0 & SOR_present==1 )
  S23Annual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                            probs=c(0.05,0.25,0.5,0.75,0.95)))
}

S23Biannual<-matrix(ncol = 5, nrow = 3, 
                    dimnames = list(c("30%", "50%","70%"),
                                    c("0.05", "0.25", "0.5", "0.75","0.95")))

for(i in 1:3){
  data=subset(df, Endemic_prevalence==ABRall[i] & MDA_frequency==0.5 & MDAi==0 & SOR_present==1 )
  S23Biannual[i,]<-c(quantile(data$Final_SOR_proportion_worms20yr, 
                              probs=c(0.05,0.25,0.5,0.75,0.95)))
}


s3<-rbind(S3Annual, S3Biannual)
s4<-rbind(S4Annual, S4Biannual)
s5<-rbind(S5Annual, S5Biannual)
s6<-rbind(S6Annual, S6Biannual)
s7<-rbind(S7Annual, S7Biannual)
s8<-rbind(S8Annual, S8Biannual)
s9<-rbind(S9Annual, S9Biannual)
s10<-rbind(S10Annual, S10Biannual)
s11<-rbind(S11Annual, S11Biannual)
s12<-rbind(S12Annual, S12Biannual)
s14<-rbind(S14Annual, S14Biannual)
s15<-rbind(S15Annual, S15Biannual)
s16<-rbind(S16Annual, S16Biannual)
s17<-rbind(S17Annual, S17Biannual)
s18<-rbind(S18Annual, S18Biannual) 
s19<-rbind(S19Annual, S19Biannual) 
s20<-rbind(S20Annual, S20Biannual) 
s21<-rbind(S21Annual, S21Biannual) 
s22<-rbind(S22Annual, S22Biannual)
s23<-rbind(S23Annual, S23Biannual)


save(s3,s4,s5,s6,
     s7,s8,s9,s10,
     s11,s12,
     s14,s15,s16,s17,
     s18,s19,s20,s21,
     s22,s23,
     file = "tables/Tables_scenarios_supp.Rdata")
