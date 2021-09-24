###############################################################################################################
###############################################################################################################
#This script simulates N-fixation data had they not been affected by the time spent at each temperature and
#refits the modified beta function (Equation 5) to these simulated data
###############################################################################################################
###############################################################################################################

####
#Load Necessary Package
####

library(bbmle)

####
#Read in data
####

###
#Morella
###

##
#21:15 deg. C growing temperature
##

M21a.dat<-read.csv("MOCE21_022018_Vmax_temp.csv")
M21b.dat<-read.csv("MOCE21_072219_Vmax_temp.csv")
M21c.dat<-read.csv("MOCE21_100319_Vmax_temp.csv")

##
#26:20 deg. C growing temperature
##

M26a.dat<-read.csv("MOCE26_121219_Vmax_temp.csv")
M26b.dat<-read.csv("MOCE26_121719_Vmax_temp.csv")
M26c.dat<-read.csv("MOCE26_020820_Vmax_temp.csv")

##
#31:25 deg. C growing temperature
##

M31a.dat<-read.csv("MOCE31_022118_Vmax_temp.csv")
M31b.dat<-read.csv("MOCE31_030619_Vmax_temp.csv")
M31c.dat<-read.csv("MOCE31_031519_Vmax_temp.csv")

###
#Alnus
###

##
#21:15 deg. C growing temperature
##

A21a.dat<-read.csv("ALRU21_070118_Vmax_temp.csv")
A21b.dat<-read.csv("ALRU21_072218_Vmax_temp.csv")
A21c.dat<-read.csv("ALRU21_073118_Vmax_temp.csv")

##
#26:20 deg. C growing temperature
##

A26a.dat<-read.csv("ALRU26_102519_Vmax_temp.csv")
A26b.dat<-read.csv("ALRU26_111519_Vmax_temp.csv")
A26c.dat<-read.csv("ALRU26_122019_Vmax_temp.csv")

##
#31:25 deg. C growing temperature
##

A31a.dat<-read.csv("ALRU31_070818_Vmax_temp.csv")
A31b.dat<-read.csv("ALRU31_073018_Vmax_temp.csv")
A31c.dat<-read.csv("ALRU31_021618_Vmax_temp.csv")

###
#Gliricidia
###

##
#21:15 deg. C growing temperature
##

G21a.dat<-read.csv("GLSE21_021518_Vmax_temp.csv")
G21b.dat<-read.csv("GLSE21_070218_Vmax_temp.csv")
G21c.dat<-read.csv("GLSE21_090219_Vmax_temp.csv")

##
#26:20 deg. C growing temperature
##

G26a.dat<-read.csv("GLSE26_121419_Vmax_temp.csv")
G26b.dat<-read.csv("GLSE26_121819_Vmax_temp.csv")
G26c.dat<-read.csv("GLSE26_012120_Vmax_temp.csv")

##
#31:25 deg. C growing temperature
##

G31a.dat<-read.csv("GLSE31_071118_Vmax_temp.csv")
G31b.dat<-read.csv("GLSE31_071618_Vmax_temp.csv")
G31c.dat<-read.csv("GLSE31_072418_Vmax_temp.csv")

###
#Robinia
###

##
#21:15 deg. C growing temperature
##

R21a.dat<-read.csv("ROPS21_021318_Vmax_temp.csv")
R21b.dat<-read.csv("ROPS21_072318_Vmax_temp.csv")
R21c.dat<-read.csv("ROPS21_080218_Vmax_temp.csv")

##
#26:20 deg. C growing temperature
##

R26a.dat<-read.csv("ROPS26_120319_Vmax_temp.csv")
R26b.dat<-read.csv("ROPS26_121619_Vmax_temp.csv")
R26c.dat<-read.csv("ROPS26_020620_Vmax_temp.csv")

##
#31:25 deg. C growing temperature
##

R31a.dat<-read.csv("ROPS31_020918_Vmax_temp.csv")
R31b.dat<-read.csv("ROPS31_071418_Vmax_temp.csv")
R31c.dat<-read.csv("ROPS31_072618_Vmax_temp.csv")

####
#Simulate N-fixation data had they not been affected by the time spent at each temperature
####

#Parameter values for equation S5 and S6
#See "Temporal_Decline_EqS7_Parameter_Calculations" script for derivation of these values
alpha<-0.0304076
beta2<-0.2923246
p<-0.6549813 
q<-0.2526679
c<-0.9572269

###
#Morella
###

##
#21:15 deg. C growing temperature
##

#Replicate 1
time<-M21a.dat$Time
time<-time-time[1]

Topt <- 29.03
tau<-M21a.dat$Temperature-Topt

Nfix<-M21a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the N-fixation rate at the 
  # next tau value that corresponds to the same N-fixation rate.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE21_022018_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-M21b.dat$Time
time<-time-time[1]

tau<-M21b.dat$Temperature-Topt

Nfix<-M21b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE21_072219_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-M21c.dat$Time
time<-time-time[1]

tau<-M21c.dat$Temperature-Topt

Nfix<-M21c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE21_100319_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#26:20 deg. C growing temperature
##

#Replicate 1
time<-M26a.dat$Time
time<-time-time[1]

Topt <- 32.94
tau<-M26a.dat$Temperature-Topt

Nfix<-M26a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE26_121219_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-M26b.dat$Time
time<-time-time[1]

tau<-M26b.dat$Temperature-Topt

Nfix<-M26b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE26_121719_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-M26c.dat$Time
time<-time-time[1]

tau<-M26c.dat$Temperature-Topt

Nfix<-M26c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE26_020820_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#31:25 deg. C growing temperature
##

#Replicate 1
time<-M31a.dat$Time
time<-time-time[1]

Topt <- 36.86
tau<-M31a.dat$Temperature-Topt

Nfix<-M31a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE31_022118_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-M31b.dat$Time
time<-time-time[1]

tau<-M31b.dat$Temperature-Topt

Nfix<-M31b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE31_030619_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-M31c.dat$Time
time<-time-time[1]

tau<-M31c.dat$Temperature-Topt

Nfix<-M31c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

MOCE31_031519_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

###
#Alnus
###

##
#21:15 deg. C growing temperature
##

#Replicate 1
time<-A21a.dat$Time
time<-time-time[1]

Topt <- 32.42
tau<-A21a.dat$Temperature-Topt

Nfix<-A21a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU21_070118_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-A21b.dat$Time
time<-time-time[1]

tau<-A21b.dat$Temperature-Topt

Nfix<-A21b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU21_072218_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-A21c.dat$Time
time<-time-time[1]

tau<-A21c.dat$Temperature-Topt

Nfix<-A21c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU21_073118_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#26:20 deg. C growing temperature
##

#Replicate 1
time<-A26a.dat$Time
time<-time-time[1]

Topt <- 32.74
tau<-A26a.dat$Temperature-Topt

Nfix<-A26a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU26_102519_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-A26b.dat$Time
time<-time-time[1]

tau<-A26b.dat$Temperature-Topt

Nfix<-A26b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU26_111519_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-A26c.dat$Time
time<-time-time[1]

tau<-A26c.dat$Temperature-Topt

Nfix<-A26c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU26_122019_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#31:25 deg. C growing temperature
##

#Replicate 1
time<-A31a.dat$Time
time<-time-time[1]

Topt <- 33.07
tau<-A31a.dat$Temperature-Topt

Nfix<-A31a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU31_070818_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-A31b.dat$Time
time<-time-time[1]

tau<-A31b.dat$Temperature-Topt

Nfix<-A31b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU31_073018_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-A31c.dat$Time
time<-time-time[1]

tau<-A31c.dat$Temperature-Topt

Nfix<-A31c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ALRU31_021618_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

###
#Gliricidia
###

##
#21:15 deg. C growing temperature
##

#Replicate 1
time<-G21a.dat$Time
time<-time-time[1]

Topt <- 31.58
tau<-G21a.dat$Temperature-Topt

Nfix<-G21a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE21_021518_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-G21b.dat$Time
time<-time-time[1]

tau<-G21b.dat$Temperature-Topt

Nfix<-G21b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE21_070218_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-G21c.dat$Time
time<-time-time[1]

tau<-G21c.dat$Temperature-Topt

Nfix<-G21c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE21_090219_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#26:20 deg. C growing temperature
##

#Replicate 1
time<-G26a.dat$Time
time<-time-time[1]

Topt <- 33.66
tau<-G26a.dat$Temperature-Topt

Nfix<-G26a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE26_121419_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-G26b.dat$Time
time<-time-time[1]

tau<-G26b.dat$Temperature-Topt

Nfix<-G26b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE26_121819_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-G26c.dat$Time
time<-time-time[1]

tau<-G26c.dat$Temperature-Topt

Nfix<-G26c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE26_012120_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#31:25 deg. C growing temperature
##

#Replicate 1
time<-G31a.dat$Time
time<-time-time[1]

Topt <- 35.74
tau<-G31a.dat$Temperature-Topt

Nfix<-G31a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE31_071118_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-G31b.dat$Time
time<-time-time[1]

tau<-G31b.dat$Temperature-Topt

Nfix<-G31b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE31_071618_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-G31c.dat$Time
time<-time-time[1]

tau<-G31c.dat$Temperature-Topt

Nfix<-G31c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

GLSE31_072418_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

###
#Robinia
###

##
#21:15 deg. C growing temperature
##

#Replicate 1
time<-R21a.dat$Time
time<-time-time[1]

Topt <- 31.89
tau<-R21a.dat$Temperature-Topt

Nfix<-R21a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS21_021318_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-R21b.dat$Time
time<-time-time[1]

tau<-R21b.dat$Temperature-Topt

Nfix<-R21b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS21_072318_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-R21c.dat$Time
time<-time-time[1]

tau<-R21c.dat$Temperature-Topt

Nfix<-R21c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS21_080218_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#26:20 deg. C growing temperature
##

#Replicate 1
time<-R26a.dat$Time
time<-time-time[1]

Topt <- 32.02
tau<-R26a.dat$Temperature-Topt

Nfix<-R26a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS26_120319_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-R26b.dat$Time
time<-time-time[1]

tau<-R26b.dat$Temperature-Topt

Nfix<-R26b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS26_121619_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-R26c.dat$Time
time<-time-time[1]

tau<-R26c.dat$Temperature-Topt

Nfix<-R26c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS26_020620_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

##
#31:25 deg. C growing temperature
##

#Replicate 1
time<-R31a.dat$Time
time<-time-time[1]

Topt <- 32.14
tau<-R31a.dat$Temperature-Topt

Nfix<-R31a.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS31_020918_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 2
time<-R31b.dat$Time
time<-time-time[1]

tau<-R31b.dat$Temperature-Topt

Nfix<-R31b.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS31_071418_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

#Replicate 3
time<-R31c.dat$Time
time<-time-time[1]

tau<-R31c.dat$Temperature-Topt

Nfix<-R31c.dat$Vmax

SNF_inst <- SNF_currenttime <- SNF_nexttime <- rep(NA,length(time)-1)
SNF_inst[1] <- 1

jvec <- kvec <- rep(NA,length(time)-1)
jvec[1] <- kvec[1] <- 1
for (i in 1:(length(time)-1)){
  # Initialize the temporal curve for the current and next tau values
  b <- exp(alpha*exp(beta2*tau[i]))-1
  d <- 1/(1+p*exp(q*tau[i]))
  SNF_currenttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  b <- exp(alpha*exp(beta2*tau[i+1]))-1
  d <- 1/(1+p*exp(q*tau[i+1]))
  SNF_nexttemp <- ((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*time))+d)
  
  # Now find the change in tau corresponding to the average
  # of the value at the current tau value and the value at the 
  # next tau value that corresponds to the same SNF value.
  
  # Index of the value that most closely matches the current relative N-fixation rate (SNF)
  j <- which.min(abs(SNF_inst[i] - SNF_currenttemp[1:i]))
  jvec[i+1] <- j
  dSNF_currenttemp <- SNF_currenttemp[j+1] - SNF_currenttemp[j]
  # Same but for the next tau value
  k <- which.min(abs(SNF_inst[i] - SNF_nexttemp[1:i]))
  kvec[i+1] <- k
  dSNF_nexttemp <- SNF_nexttemp[k+1] - SNF_nexttemp[k]
  SNF_inst[i+1] <- SNF_inst[i] + mean(c(dSNF_currenttemp,dSNF_nexttemp))
}

ROPS31_072618_inst<-data.frame(Time=time,Temperature=tau+Topt,Vmax.og=Nfix,Vmax.inst=Nfix/SNF_inst)

####
#The following exports the above as .csv files
####

write.csv(MOCE21_022018_inst,"MOCE21_022018_inst.csv")
write.csv(MOCE21_072219_inst,"MOCE21_072219_inst.csv")
write.csv(MOCE21_100319_inst,"MOCE21_100319_inst.csv")
write.csv(MOCE26_121219_inst,"MOCE26_121219_inst.csv")
write.csv(MOCE26_121719_inst,"MOCE26_121719_inst.csv")
write.csv(MOCE26_020820_inst,"MOCE26_020820_inst.csv")
write.csv(MOCE31_022118_inst,"MOCE31_022118_inst.csv")
write.csv(MOCE31_030619_inst,"MOCE31_030619_inst.csv")
write.csv(MOCE31_031519_inst,"MOCE31_031519_inst.csv")
write.csv(ALRU21_070118_inst,"ALRU21_070118_inst.csv")
write.csv(ALRU21_072218_inst,"ALRU21_072218_inst.csv")
write.csv(ALRU21_073118_inst,"ALRU21_073118_inst.csv")
write.csv(ALRU26_102519_inst,"ALRU26_102519_inst.csv")
write.csv(ALRU26_111519_inst,"ALRU26_111519_inst.csv")
write.csv(ALRU26_122019_inst,"ALRU26_122019_inst.csv")
write.csv(ALRU31_070818_inst,"ALRU31_070818_inst.csv")
write.csv(ALRU31_073018_inst,"ALRU31_073018_inst.csv")
write.csv(ALRU31_021618_inst,"ALRU31_021618_inst.csv")
write.csv(GLSE21_021518_inst,"GLSE21_021518_inst.csv")
write.csv(GLSE21_070218_inst,"GLSE21_070218_inst.csv")
write.csv(GLSE21_090219_inst,"GLSE21_090219_inst.csv")
write.csv(GLSE26_121419_inst,"GLSE26_121419_inst.csv")
write.csv(GLSE26_121819_inst,"GLSE26_121819_inst.csv")
write.csv(GLSE26_012120_inst,"GLSE26_012120_inst.csv")
write.csv(GLSE31_071118_inst,"GLSE31_071118_inst.csv")
write.csv(GLSE31_071618_inst,"GLSE31_071618_inst.csv")
write.csv(GLSE31_072418_inst,"GLSE31_072418_inst.csv")
write.csv(ROPS21_021318_inst,"ROPS21_021318_inst.csv")
write.csv(ROPS21_072318_inst,"ROPS21_072318_inst.csv")
write.csv(ROPS21_080218_inst,"ROPS21_080218_inst.csv")
write.csv(ROPS26_120319_inst,"ROPS26_120319_inst.csv")
write.csv(ROPS26_121619_inst,"ROPS26_121619_inst.csv")
write.csv(ROPS26_020620_inst,"ROPS26_020620_inst.csv")
write.csv(ROPS31_020918_inst,"ROPS31_020918_inst.csv")
write.csv(ROPS31_071418_inst,"ROPS31_071418_inst.csv")
write.csv(ROPS31_072618_inst,"ROPS31_072618_inst.csv")

####
#The following refits the modified beta function (Equation 5) to the simulated data generated above
####

####
#Read in data that were generated above
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

####
#Define functions
####

#Acclimation of all parameters (used for N-fixation with Morella, Alnus, and Robinia)
beta.lin.all <- function(ymax,a,b,c,d,e,f,T,Tgrow){
  y <- pmax(0,ymax*((e+f*Tgrow)-T)/((e+f*Tgrow)-(c+d*Tgrow))*(((T-(a+b*Tgrow))/((c+d*Tgrow)-(a+b*Tgrow)))^(((c+d*Tgrow)-(a+b*Tgrow))/((e+f*Tgrow)-(c+d*Tgrow)))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Tmin and Topt (used for N-fixation with Gliricidia)
beta.lin.Tmin.Topt <- function(ymax,Tmax,a,b,c,d,T,Tgrow){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-(c+d*Tgrow))*(((T-(a+b*Tgrow))/((c+d*Tgrow)-(a+b*Tgrow)))^(((c+d*Tgrow)-(a+b*Tgrow))/(Tmax-(c+d*Tgrow)))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

####
#Negative log-likelihood (NLL) functions
####

#NLL function for acclimation of all parameters (used for N-fixation with Morella, Alnus, and Robinia)
Nase_beta_linall_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                     ymax26a,ymax26b,ymax26c,
                                     ymax31a,ymax31b,ymax31c,
                                     a,b,c,d,e,f,
                                     T21a,T21b,T21c,
                                     T26a,T26b,T26c,
                                     T31a,T31b,T31c,
                                     Nasedat21a,Nasedat21b,Nasedat21c,
                                     Nasedat26a,Nasedat26b,Nasedat26c,
                                     Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.all(ymax21a,a,b,c,d,e,f,T21a,18.5)
  Nasemean21b <- beta.lin.all(ymax21b,a,b,c,d,e,f,T21b,18.5)
  Nasemean21c <- beta.lin.all(ymax21c,a,b,c,d,e,f,T21c,18.5)
  Nasemean26a <- beta.lin.all(ymax26a,a,b,c,d,e,f,T26a,23.5)
  Nasemean26b <- beta.lin.all(ymax26b,a,b,c,d,e,f,T26b,23.5)
  Nasemean26c <- beta.lin.all(ymax26c,a,b,c,d,e,f,T26c,23.5)
  Nasemean31a <- beta.lin.all(ymax31a,a,b,c,d,e,f,T31a,28.5)
  Nasemean31b <- beta.lin.all(ymax31b,a,b,c,d,e,f,T31b,28.5)
  Nasemean31c <- beta.lin.all(ymax31c,a,b,c,d,e,f,T31c,28.5)
  -(sum(dnorm(Nasedat21a,mean=Nasemean21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21b,mean=Nasemean21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21c,mean=Nasemean21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26a,mean=Nasemean26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26b,mean=Nasemean26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26c,mean=Nasemean26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Tmin and Topt (used for N-fixation with Gliricidia)
Nase_beta_linTminTopt_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                          ymax26a,ymax26b,ymax26c,
                                          ymax31a,ymax31b,ymax31c,
                                          Tmax,
                                          a,b,c,d,
                                          T21a,T21b,T21c,
                                          T26a,T26b,T26c,
                                          T31a,T31b,T31c,
                                          Nasedat21a,Nasedat21b,Nasedat21c,
                                          Nasedat26a,Nasedat26b,Nasedat26c,
                                          Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.Tmin.Topt(ymax21a,Tmax,a,b,c,d,T21a,18.5)
  Nasemean21b <- beta.lin.Tmin.Topt(ymax21b,Tmax,a,b,c,d,T21b,18.5)
  Nasemean21c <- beta.lin.Tmin.Topt(ymax21c,Tmax,a,b,c,d,T21c,18.5)
  Nasemean26a <- beta.lin.Tmin.Topt(ymax26a,Tmax,a,b,c,d,T26a,23.5)
  Nasemean26b <- beta.lin.Tmin.Topt(ymax26b,Tmax,a,b,c,d,T26b,23.5)
  Nasemean26c <- beta.lin.Tmin.Topt(ymax26c,Tmax,a,b,c,d,T26c,23.5)
  Nasemean31a <- beta.lin.Tmin.Topt(ymax31a,Tmax,a,b,c,d,T31a,28.5)
  Nasemean31b <- beta.lin.Tmin.Topt(ymax31b,Tmax,a,b,c,d,T31b,28.5)
  Nasemean31c <- beta.lin.Tmin.Topt(ymax31c,Tmax,a,b,c,d,T31c,28.5)
  -(sum(dnorm(Nasedat21a,mean=Nasemean21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21b,mean=Nasemean21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21c,mean=Nasemean21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26a,mean=Nasemean26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26b,mean=Nasemean26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26c,mean=Nasemean26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

####
#Maximum likelihood fits
####

#Morella
fit_Nase_beta_linall_MOCE_inst <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                           ymax26a=1,ymax26b=1,ymax26c=1,
                                                                           ymax31a=1,ymax31b=1,ymax31c=1,
                                                                           a=-20,b=0,c=19,d=0.67,e=44,f=0.03),
                                       data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                                 T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                                 T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                                 Nasedat21a=MOCE21a$Vmax.inst/9960.998949,Nasedat21b=MOCE21b$Vmax.inst/3833.121977,Nasedat21c=MOCE21c$Vmax.inst/3095.456985,
                                                 Nasedat26a=MOCE26a$Vmax.inst/6492.677546,Nasedat26b=MOCE26b$Vmax.inst/5974.076428,Nasedat26c=MOCE26c$Vmax.inst/7890.763063,
                                                 Nasedat31a=MOCE31a$Vmax.inst/8826.601393,Nasedat31b=MOCE31b$Vmax.inst/2371.923752,Nasedat31c=MOCE31c$Vmax.inst/629.7544222),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linall_MOCE_inst) #Summary
ft.m.inst <- coef(fit_Nase_beta_linall_MOCE_inst) #Parameter estimates
write.csv(ft.m.inst,"MOCE_inst_beta_fit.csv") #Exports parameter estimates

#Alnus
fit_Nase_beta_linall_ALRU_inst <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                           ymax26a=1.0,ymax26b=1,ymax26c=0.9,
                                                                           ymax31a=1.2,ymax31b=1.1,ymax31c=1.1,
                                                                           a=-20,b=0.49,c=31.2,d=0.064,e=41.9,f=0.027),
                                       data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                                 T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                                 T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                                 Nasedat21a=ALRU21a$Vmax.inst/10035.24751,Nasedat21b=ALRU21b$Vmax.inst/7990.765838,Nasedat21c=ALRU21c$Vmax.inst/4046.082544,
                                                 Nasedat26a=ALRU26a$Vmax.inst/2236.979133,Nasedat26b=ALRU26b$Vmax.inst/2744.252462,Nasedat26c=ALRU26c$Vmax.inst/1831.763172,
                                                 Nasedat31a=ALRU31a$Vmax.inst/406.3524389,Nasedat31b=ALRU31b$Vmax.inst/1832.723413,Nasedat31c=ALRU31c$Vmax.inst/1350.729381),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linall_ALRU_inst) #Summary
ft.a.inst <- coef(fit_Nase_beta_linall_ALRU_inst) #Parameter estimates
write.csv(ft.a.inst,"ALRU_inst_beta_fit.csv") #Exports parameter estimates

#Gliricidia
fit_Nase_beta_linTminTopt_GLSE_inst <- mle2(Nase_beta_linTminTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                     ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                     ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                     Tmax=45,
                                                                                     a=5,b=0,c=19,d=0.6),
                                            data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                                      T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                                      T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                                      Nasedat21a=GLSE21a$Vmax.inst/308.3421323,Nasedat21b=GLSE21b$Vmax.inst/2028.57371,Nasedat21c=GLSE21c$Vmax.inst/879.0395611,
                                                      Nasedat26a=GLSE26a$Vmax.inst/909.4624356,Nasedat26b=GLSE26b$Vmax.inst/1373.906634,Nasedat26c=GLSE26c$Vmax.inst/705.9735103,
                                                      Nasedat31a=GLSE31a$Vmax.inst/1196.373627,Nasedat31b=GLSE31b$Vmax.inst/312.4607561,Nasedat31c=GLSE31c$Vmax.inst/209.4577018),
                                            control=list(maxit=20000))
summary(fit_Nase_beta_linTminTopt_GLSE_inst) #Summary
ft.g.inst <- coef(fit_Nase_beta_linTminTopt_GLSE_inst) #Parameter estimates
write.csv(ft.g.inst,"GLSE_inst_beta_fit.csv") #Exports parameter estimates

#Robinia
fit_Nase_beta_linall_ROPS_inst <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=1,ymax21a=1,ymax21b=0.8,ymax21c=1.1,
                                                                           ymax26a=0.9,ymax26b=1,ymax26c=0.9,
                                                                           ymax31a=0.9,ymax31b=0.7,ymax31c=0.9,
                                                                           a=-13,b=0.8,c=32,d=0.03,e=47,f=0.03),
                                       data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                                 T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                                 T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                                 Nasedat21a=ROPS21a$Vmax.inst/1922.245351,Nasedat21b=ROPS21b$Vmax.inst/1690.867422,Nasedat21c=ROPS21c$Vmax.inst/450.8292229,
                                                 Nasedat26a=ROPS26a$Vmax.inst/2785.452785,Nasedat26b=ROPS26b$Vmax.inst/1572.439251,Nasedat26c=ROPS26c$Vmax.inst/1686.352,
                                                 Nasedat31a=ROPS31a$Vmax.inst/1646.769302,Nasedat31b=ROPS31b$Vmax.inst/920.7416901,Nasedat31c=ROPS31c$Vmax.inst/1978.187071),
                                       control=list(maxit=50000))
summary(fit_Nase_beta_linall_ROPS_inst) #Summary
ft.r.inst <- coef(fit_Nase_beta_linall_ROPS_inst) #Parameter estimates
write.csv(ft.r.inst,"ROPS_inst_beta_fit.csv") #Exports parameter estimates
