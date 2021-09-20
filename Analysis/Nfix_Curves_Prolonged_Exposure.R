###############################################################################################################
###############################################################################################################
#This script simulates N-fixation temperature response curves for different heating rates
###############################################################################################################
###############################################################################################################

####
#Set Working Directory as the "SNF_temporal_decline" folder
####

####
#Read in data that were generated in Temporal_Decline_Nfix_Minus_Temporal_Effect.R
####

#Morella
MOCE21a <- read.csv("MOCE21_022018_inst.csv")
MOCE21b <- read.csv("MOCE21_072219_inst.csv")
MOCE21c <- read.csv("MOCE21_100319_inst.csv")
MOCE26a <- read.csv("MOCE26_121219_inst.csv")
MOCE26b <- read.csv("MOCE26_121719_inst.csv")
MOCE26c <- read.csv("MOCE26_020820_inst.csv")
MOCE31a <- read.csv("MOCE31_022118_inst.csv")
MOCE31b <- read.csv("MOCE31_030619_inst.csv")
MOCE31c <- read.csv("MOCE31_031519_inst.csv")
M.inst<-read.csv("MOCE_inst_beta_fit.csv")

#Alnus
ALRU21a <- read.csv("ALRU21_070118_inst.csv")
ALRU21b <- read.csv("ALRU21_072218_inst.csv")
ALRU21c <- read.csv("ALRU21_073118_inst.csv")
ALRU26a <- read.csv("ALRU26_102519_inst.csv")
ALRU26b <- read.csv("ALRU26_111519_inst.csv")
ALRU26c <- read.csv("ALRU26_122019_inst.csv")
ALRU31a <- read.csv("ALRU31_021618_inst.csv")
ALRU31b <- read.csv("ALRU31_070818_inst.csv")
ALRU31c <- read.csv("ALRU31_073018_inst.csv")
A.inst<-read.csv("ALRU_inst_beta_fit.csv")

#Gliricidia
GLSE21a <- read.csv("GLSE21_021518_inst.csv")
GLSE21b <- read.csv("GLSE21_070218_inst.csv")
GLSE21c <- read.csv("GLSE21_090219_inst.csv")
GLSE26a <- read.csv("GLSE26_121419_inst.csv")
GLSE26b <- read.csv("GLSE26_121819_inst.csv")
GLSE26c <- read.csv("GLSE26_012120_inst.csv")
GLSE31a <- read.csv("GLSE31_071118_inst.csv")
GLSE31b <- read.csv("GLSE31_071618_inst.csv")
GLSE31c <- read.csv("GLSE31_072418_inst.csv")
G.inst<-read.csv("GLSE_inst_beta_fit.csv")

#Robinia
ROPS21a <- read.csv("ROPS21_021318_inst.csv")
ROPS21b <- read.csv("ROPS21_072318_inst.csv")
ROPS21c <- read.csv("ROPS21_080218_inst.csv")
ROPS26a <- read.csv("ROPS26_120319_inst.csv")
ROPS26b <- read.csv("ROPS26_121619_inst.csv")
ROPS26c <- read.csv("ROPS26_020620_inst.csv")
ROPS31a <- read.csv("ROPS31_020918_inst.csv")
ROPS31b <- read.csv("ROPS31_071418_inst.csv")
ROPS31c <- read.csv("ROPS31_072618_inst.csv")
R.inst<-read.csv("ROPS_inst_beta_fit.csv")

####
#Define modified beta function (equation 5)
####

beta <- function(ymax,Tmin,Topt,Tmax,T){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

####
#Parameter values for equation S5 and S6, which are part of equation S7
#See "Temporal_Decline_EqS7_Parameter_Calculations" script for derivation of these values
####

alpha<-0.0304076
beta2<-0.2923246
p<-0.6549813 
q<-0.2526679
c<-0.9572269

####
#Simulation of N-fixation temperature response curves at 1x, 2x slower, and 3x slower than
#the heating rate at which data were collected
####

###
#Morella
###

#Modified beta function (equation 5) parameters without decline due to the time spent at each temperature
Tmin21.i<-18.5*M.inst[12,2]+M.inst[11,2]
Tmin26.i<-23.5*M.inst[12,2]+M.inst[11,2]
Tmin31.i<-28.5*M.inst[12,2]+M.inst[11,2]
Topt21.i<-18.5*M.inst[14,2]+M.inst[13,2]
Topt26.i<-23.5*M.inst[14,2]+M.inst[13,2]
Topt31.i<-28.5*M.inst[14,2]+M.inst[13,2]
Tmax21.i<-18.5*M.inst[16,2]+M.inst[15,2]
Tmax26.i<-23.5*M.inst[16,2]+M.inst[15,2]
Tmax31.i<-28.5*M.inst[16,2]+M.inst[15,2]

#Measured modified beta function parameters
Topt21.o<-29.03
Topt26.o<-32.94
Topt31.o<-36.86

##
#21:15 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin21.i,Topt21.i,Tmax21.i,25)
ins.fit<-beta(1/ins.25,Tmin21.i,Topt21.i,Tmax21.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

M21.Topt.x1<-temp.sim[which.max(x1.fit)]
M21.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
M21.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

M21.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
M21.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
M21.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#26:20 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin26.i,Topt26.i,Tmax26.i,25)
ins.fit<-beta(1/ins.25,Tmin26.i,Topt26.i,Tmax26.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

M26.Topt.x1<-temp.sim[which.max(x1.fit)]
M26.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
M26.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

M26.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
M26.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
M26.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#31:25 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin31.i,Topt31.i,Tmax31.i,25)
ins.fit<-beta(1/ins.25,Tmin31.i,Topt31.i,Tmax31.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

M31.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
M31.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
M31.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

M31.Topt.x1<-temp.sim[which.max(x1.fit)]
M31.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
M31.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

###
#Alnus
###

#Modified beta function (equation 5) parameters without decline due to the time spent at each temperature
Tmin21.i<-18.5*A.inst[12,2]+A.inst[11,2]
Tmin26.i<-23.5*A.inst[12,2]+A.inst[11,2]
Tmin31.i<-28.5*A.inst[12,2]+A.inst[11,2]
Topt21.i<-18.5*A.inst[14,2]+A.inst[13,2]
Topt26.i<-23.5*A.inst[14,2]+A.inst[13,2]
Topt31.i<-28.5*A.inst[14,2]+A.inst[13,2]
Tmax21.i<-18.5*A.inst[16,2]+A.inst[15,2]
Tmax26.i<-23.5*A.inst[16,2]+A.inst[15,2]
Tmax31.i<-28.5*A.inst[16,2]+A.inst[15,2]

#Measured modified beta function parameters
Topt21.o<-32.42
Topt26.o<-32.74
Topt31.o<-33.07

##
#21:15 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin21.i,Topt21.i,Tmax21.i,25)
ins.fit<-beta(1/ins.25,Tmin21.i,Topt21.i,Tmax21.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

A21.Topt.x1<-temp.sim[which.max(x1.fit)]
A21.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
A21.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

A21.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
A21.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
A21.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#26:20 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin26.i,Topt26.i,Tmax26.i,25)
ins.fit<-beta(1/ins.25,Tmin26.i,Topt26.i,Tmax26.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

A26.Topt.x1<-temp.sim[which.max(x1.fit)]
A26.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
A26.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

A26.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
A26.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
A26.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#31:25 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin31.i,Topt31.i,Tmax31.i,25)
ins.fit<-beta(1/ins.25,Tmin31.i,Topt31.i,Tmax31.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

A31.Topt.x1<-temp.sim[which.max(x1.fit)]
A31.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
A31.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

A31.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
A31.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
A31.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

###
#Gliricidia
###

#Modified beta function (equation 5) parameters without decline due to the time spent at each temperature
Tmin21.i<-18.5*G.inst[13,2]+G.inst[12,2]
Tmin26.i<-23.5*G.inst[13,2]+G.inst[12,2]
Tmin31.i<-28.5*G.inst[13,2]+G.inst[12,2]
Topt21.i<-18.5*G.inst[15,2]+G.inst[14,2]
Topt26.i<-23.5*G.inst[15,2]+G.inst[14,2]
Topt31.i<-28.5*G.inst[15,2]+G.inst[14,2]
Tmax21.i<-G.inst[11,2]
Tmax26.i<-G.inst[11,2]
Tmax31.i<-G.inst[11,2]

#Measured modified beta function parameters
Topt21.o<-31.58
Topt26.o<-33.66
Topt31.o<-35.74

##
#21:15 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin21.i,Topt21.i,Tmax21.i,25)
ins.fit<-beta(1/ins.25,Tmin21.i,Topt21.i,Tmax21.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

G21.Topt.x1<-temp.sim[which.max(x1.fit)]
G21.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
G21.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

G21.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
G21.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
G21.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#26:20 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin26.i,Topt26.i,Tmax26.i,25)
ins.fit<-beta(1/ins.25,Tmin26.i,Topt26.i,Tmax26.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

G26.Topt.x1<-temp.sim[which.max(x1.fit)]
G26.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
G26.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

G26.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
G26.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
G26.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#31:25 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin31.i,Topt31.i,Tmax31.i,25)
ins.fit<-beta(1/ins.25,Tmin31.i,Topt31.i,Tmax31.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

G31.Topt.x1<-temp.sim[which.max(x1.fit)]
G31.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
G31.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

G31.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
G31.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
G31.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

###
#Robinia
###

#Modified beta function (equation 5) parameters without decline due to the time spent at each temperature
Tmin21.i<-18.5*R.inst[12,2]+R.inst[11,2]
Tmin26.i<-23.5*R.inst[12,2]+R.inst[11,2]
Tmin31.i<-28.5*R.inst[12,2]+R.inst[11,2]
Topt21.i<-18.5*R.inst[14,2]+R.inst[13,2]
Topt26.i<-23.5*R.inst[14,2]+R.inst[13,2]
Topt31.i<-28.5*R.inst[14,2]+R.inst[13,2]
Tmax21.i<-18.5*R.inst[16,2]+R.inst[15,2]
Tmax26.i<-23.5*R.inst[16,2]+R.inst[15,2]
Tmax31.i<-28.5*R.inst[16,2]+R.inst[15,2]

#Measured modified beta function parameters
Topt21.o<-31.89
Topt26.o<-32.02
Topt31.o<-32.14

##
#21:15 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin21.i,Topt21.i,Tmax21.i,25)
ins.fit<-beta(1/ins.25,Tmin21.i,Topt21.i,Tmax21.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

R21.Topt.x1<-temp.sim[which.max(x1.fit)]
R21.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
R21.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

R21.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
R21.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
R21.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#26:20 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours) 
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin26.i,Topt26.i,Tmax26.i,25)
ins.fit<-beta(1/ins.25,Tmin26.i,Topt26.i,Tmax26.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

R26.Topt.x1<-temp.sim[which.max(x1.fit)]
R26.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
R26.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

R26.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
R26.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
R26.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

##
#31:25 deg. C growing temperature
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.o
tau.sim <- seq(0,50,0.01) - Topt 

#
#1x heating rate
#

#Set run time (in hours)
runtime <- 8.333

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.a <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.a[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.a[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.a[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.a[i+1] <- SNF_realized.a[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#2x slower heating rate
#

#Set run time (in hours)
runtime <- 16.667

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.b <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.b[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.b[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.b[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.b[i+1] <- SNF_realized.b[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#
#3x slower heating rate
#

#Set run time (in hours)
runtime <- 25

#Create time vector for given run time
time.sim <- seq(0,runtime,(runtime-0)/(length(tau.sim)-1))

SNF_realized.c <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time.sim)-1)
SNF_realized.c[1] <- 1

jvec <- kvec <- rep(NA,length(time.sim)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time.sim)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau.sim[i]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  b <- exp(alpha*exp(beta2*tau.sim[i+1]))-1 
  d <- 1/(1+p*exp(q*tau.sim[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time.sim))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_realized.c[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_realized.c[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_realized.c[i+1] <- SNF_realized.c[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

#Calculate Topt and relative N-fixation rate at 15 and 40 deg. C for 1x and 3x slower heating rate
temp.sim<-tau.sim+Topt

ins.25<-beta(1,Tmin31.i,Topt31.i,Tmax31.i,25)
ins.fit<-beta(1/ins.25,Tmin31.i,Topt31.i,Tmax31.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x3.slow.fit<-ins.fit*SNF_realized.c

R31.Topt.x1<-temp.sim[which.max(x1.fit)]
R31.15.x1<-x1.fit[which(temp.sim==15)]/max(x1.fit)
R31.40.x1<-x1.fit[which(temp.sim==40)]/max(x1.fit)

R31.Topt.x3.slow<-temp.sim[which.max(x3.slow.fit)]
R31.15.x3.slow<-x3.slow.fit[which(temp.sim==15)]/max(x3.slow.fit)
R31.40.x3.slow<-x3.slow.fit[which(temp.sim==40)]/max(x3.slow.fit)

####
#Exports parameter estimates (available in "SNF_temporal_decline" folder)
####

#For 1x heating rate
Tr<-c("M21","M26","M31","A21","A26","A31","G21","G26","G31","R21","R26","R31")
Topt<-c(M21.Topt.x1,M26.Topt.x1,M31.Topt.x1,A21.Topt.x1,A26.Topt.x1,A31.Topt.x1,
        G21.Topt.x1,G26.Topt.x1,G31.Topt.x1,R21.Topt.x1,R26.Topt.x1,R31.Topt.x1)
r15<-c(M21.15.x1,M26.15.x1,M31.15.x1,A21.15.x1,A26.15.x1,A31.15.x1,
       G21.15.x1,G26.15.x1,G31.15.x1,R21.15.x1,R26.15.x1,R31.15.x1)
r40<-c(M21.40.x1,M26.40.x1,M31.40.x1,A21.40.x1,A26.40.x1,A31.40.x1,
       G21.40.x1,G26.40.x1,G31.40.x1,R21.40.x1,R26.40.x1,R31.40.x1)
write.csv(data.frame(Treat=Tr,Topt=Topt,r15=r15,r40=r40),"x1.Topt.r15.r40.csv")

#For 3x slower heating rate
Tr<-c("M21","M26","M31","A21","A26","A31","G21","G26","G31","R21","R26","R31")
Topt<-c(M21.Topt.x3.slow,M26.Topt.x3.slow,M31.Topt.x3.slow,A21.Topt.x3.slow,A26.Topt.x3.slow,A31.Topt.x3.slow,
        G21.Topt.x3.slow,G26.Topt.x3.slow,G31.Topt.x3.slow,R21.Topt.x3.slow,R26.Topt.x3.slow,R31.Topt.x3.slow)
r15<-c(M21.15.x3.slow,M26.15.x3.slow,M31.15.x3.slow,A21.15.x3.slow,A26.15.x3.slow,A31.15.x3.slow,
       G21.15.x3.slow,G26.15.x3.slow,G31.15.x3.slow,R21.15.x3.slow,R26.15.x3.slow,R31.15.x3.slow)
r40<-c(M21.40.x3.slow,M26.40.x3.slow,M31.40.x3.slow,A21.40.x3.slow,A26.40.x3.slow,A31.40.x3.slow,
       G21.40.x3.slow,G26.40.x3.slow,G31.40.x3.slow,R21.40.x3.slow,R26.40.x3.slow,R31.40.x3.slow)
write.csv(data.frame(Treat=Tr,Topt=Topt,r15=r15,r40=r40),"slow.3x.Topt.r15.r40.csv")

####
#Calculate simulated change in Topt and rates at 15 and 40 deg. C due to prolonged exposure (3x slower heating rate)
####

#Read in data (generated above) from "SNF_temporal_decline" folder
x1.dat<-read.csv("x1.Topt.r15.r40.csv")
slow.3x.dat<-read.csv("slow.3x.Topt.r15.r40.csv")

#Topt and rates at 15 and 40 deg. C due to prolonged exposure
perc.dif.time<-(slow.3x.dat/x1.dat)*100-100

###
#Calculate mean, min, and max
###

#Topt (deg. C difference in Topt, not %)
mean(x1.dat[,3]-slow.3x.dat[,3]) #1.8
min(x1.dat[,3]-slow.3x.dat[,3]) #1.3
max(x1.dat[,3]-slow.3x.dat[,3]) #2.2

#Rate at 15 deg. C
mean(perc.dif.time[,4]) #+5.1
min(perc.dif.time[,4]) #+4.8
max(perc.dif.time[,4]) #+5.8

#Rate at 40 deg. C
mean(perc.dif.time[,5]) #-48.0
min(perc.dif.time[,5]) #-65.5
max(perc.dif.time[,5]) #-19.0
