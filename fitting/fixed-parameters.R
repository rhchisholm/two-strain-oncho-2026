# Define fixed parameters

#Transmission to humans
#ABR<-65000 #80%prev
#ABR<-27000 #70%prev
#ABR<-5350 #50%prev
#ABR<-2942 #30%prev

my_beta<-31.2  # ABR = beta*m 
my_m<-ABR/my_beta #
my_Omega<-1
my_muH <-0.04 #
my_p<-2 #
my_deltaH0<- 0.0854 #
my_deltaHinfinity<-0.00299  #
my_cH<-5.86*(10^(-3)) #

# Flow between worms and ivermectin effects
#F->N: lambda0+lambda(tau), lambda(tau)=r*lambdaMAX*exp(-psi_i*tau)
#N->F: omega0+omega(tau), omega(tau)=-omegabar*r*exp(-psi_i*tau) 
#my_r<-0.4
my_lambda0 <-0.33
if(MDAi==1){
  my_lambdaMAX<-32.4  #
  my_psiGR<-19.6 #
  my_nu<-0.0096 #
  my_omega<- 1.25 #
}else{
  my_lambdaMAX<-462 #
  my_psiGR<-4.83 #
  my_nu<-0.04 #
  my_omega<- 1.82 #
}

#Worm loss
my_sigmaW<-0.1 #

#Mf production
my_deltaV0<-0.0207 #or 0.005 #
my_cV<-0.0148 #or 0 #
my_kW<-0.35 #
my_zeta<-0.3 #
my_epsilon<- 1.1538/2 #0.67/2 #

#Mf loss
my_sigmaM0<-0.8 #
my_alphaH<-0.5 #
my_g<-0.0096 #

#L loss
my_sigmaL<-104 #
my_muL<-52 #
my_alphaV<-0.6 #0.597 #

# Prevalence
my_k0<-0.013 #Karen, Turner 0.0553 (different form kM)
my_k1<-0.025 #Karen, Turner 0.4910