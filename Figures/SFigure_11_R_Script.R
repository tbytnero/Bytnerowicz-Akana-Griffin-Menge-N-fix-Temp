###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figure 11
###############################################################################################################
###############################################################################################################

####
#All figures were exported as PDFs from R Studio
#Dimensions for each figure are given with script specific to each figure
####

####
#Read N-fixation data in from "SNF_Temp" folder
####

M.og<-read.csv("MOCE_lin_beta_fit.csv")
A.og<-read.csv("ALRU_lin_beta_fit.csv")
G.og<-read.csv("GLSE_lin_beta_fit.csv")
R.og<-read.csv("ROPS_lin_beta_fit.csv")

####
#Read N-fixation data in from "SNF_temporal_decline" folder
####

M.inst<-read.csv("MOCE_inst_beta_fit.csv")
A.inst<-read.csv("ALRU_inst_beta_fit.csv")
G.inst<-read.csv("GLSE_inst_beta_fit.csv")
R.inst<-read.csv("ROPS_inst_beta_fit.csv")

####
#Extract growing temperature specific parameter estimates
####

###
#Measured parameter estimates
###

#Morella
Tmin21.M.o<-M.og[11,2]+M.og[12,2]*18.5
Tmin26.M.o<-M.og[11,2]+M.og[12,2]*23.5
Tmin31.M.o<-M.og[11,2]+M.og[12,2]*28.5
Topt21.M.o<-M.og[13,2]+M.og[14,2]*18.5
Topt26.M.o<-M.og[13,2]+M.og[14,2]*23.5
Topt31.M.o<-M.og[13,2]+M.og[14,2]*28.5
Tmax21.M.o<-M.og[15,2]+M.og[16,2]*18.5
Tmax26.M.o<-M.og[15,2]+M.og[16,2]*23.5
Tmax31.M.o<-M.og[15,2]+M.og[16,2]*28.5

#Alnus
Tmin21.A.o<-A.og[11,2]+A.og[12,2]*18.5
Tmin26.A.o<-A.og[11,2]+A.og[12,2]*23.5
Tmin31.A.o<-A.og[11,2]+A.og[12,2]*28.5
Topt21.A.o<-A.og[13,2]+A.og[14,2]*18.5
Topt26.A.o<-A.og[13,2]+A.og[14,2]*23.5
Topt31.A.o<-A.og[13,2]+A.og[14,2]*28.5
Tmax21.A.o<-A.og[15,2]+A.og[16,2]*18.5
Tmax26.A.o<-A.og[15,2]+A.og[16,2]*23.5
Tmax31.A.o<-A.og[15,2]+A.og[16,2]*28.5

#Gliricidia
Tmin21.G.o<-G.og[12,2]+G.og[13,2]*18.5
Tmin26.G.o<-G.og[12,2]+G.og[13,2]*23.5
Tmin31.G.o<-G.og[12,2]+G.og[13,2]*28.5
Topt21.G.o<-G.og[14,2]+G.og[15,2]*18.5
Topt26.G.o<-G.og[14,2]+G.og[15,2]*23.5
Topt31.G.o<-G.og[14,2]+G.og[15,2]*28.5
Tmax21.G.o<-G.og[11,2]
Tmax26.G.o<-G.og[11,2]
Tmax31.G.o<-G.og[11,2]

#Robinia
Tmin21.R.o<-R.og[11,2]+R.og[12,2]*18.5
Tmin26.R.o<-R.og[11,2]+R.og[12,2]*23.5
Tmin31.R.o<-R.og[11,2]+R.og[12,2]*28.5
Topt21.R.o<-R.og[13,2]+R.og[14,2]*18.5
Topt26.R.o<-R.og[13,2]+R.og[14,2]*23.5
Topt31.R.o<-R.og[13,2]+R.og[14,2]*28.5
Tmax21.R.o<-R.og[15,2]+R.og[16,2]*18.5
Tmax26.R.o<-R.og[15,2]+R.og[16,2]*23.5
Tmax31.R.o<-R.og[15,2]+R.og[16,2]*28.5

###
#Simulated parameter estimates minus temporal decline effect
###

#Morella
Tmax21.M.i<-M.inst[15,2]+M.inst[16,2]*18.5
Tmax26.M.i<-M.inst[15,2]+M.inst[16,2]*23.5
Tmax31.M.i<-M.inst[15,2]+M.inst[16,2]*28.5
Tmin21.M.i<-M.inst[11,2]+M.inst[12,2]*18.5
Tmin26.M.i<-M.inst[11,2]+M.inst[12,2]*23.5
Tmin31.M.i<-M.inst[11,2]+M.inst[12,2]*28.5
Topt21.M.i<-M.inst[13,2]+M.inst[14,2]*18.5
Topt26.M.i<-M.inst[13,2]+M.inst[14,2]*23.5
Topt31.M.i<-M.inst[13,2]+M.inst[14,2]*28.5

#Alnus
Tmax21.A.i<-A.inst[15,2]+A.inst[16,2]*18.5
Tmax26.A.i<-A.inst[15,2]+A.inst[16,2]*23.5
Tmax31.A.i<-A.inst[15,2]+A.inst[16,2]*28.5
Tmin21.A.i<-A.inst[11,2]+A.inst[12,2]*18.5
Tmin26.A.i<-A.inst[11,2]+A.inst[12,2]*23.5
Tmin31.A.i<-A.inst[11,2]+A.inst[12,2]*28.5
Topt21.A.i<-A.inst[13,2]+A.inst[14,2]*18.5
Topt26.A.i<-A.inst[13,2]+A.inst[14,2]*23.5
Topt31.A.i<-A.inst[13,2]+A.inst[14,2]*28.5

#Gliricidia
Tmin21.G.i<-G.inst[12,2]+G.inst[13,2]*18.5
Tmin26.G.i<-G.inst[12,2]+G.inst[13,2]*23.5
Tmin31.G.i<-G.inst[12,2]+G.inst[13,2]*28.5
Topt21.G.i<-G.inst[14,2]+G.inst[15,2]*18.5
Topt26.G.i<-G.inst[14,2]+G.inst[15,2]*23.5
Topt31.G.i<-G.inst[14,2]+G.inst[15,2]*28.5
Tmax21.G.i<-G.inst[11,2]
Tmax26.G.i<-G.inst[11,2]
Tmax31.G.i<-G.inst[11,2]

#Robinia
Tmax21.R.i<-R.inst[15,2]+R.inst[16,2]*18.5
Tmax26.R.i<-R.inst[15,2]+R.inst[16,2]*23.5
Tmax31.R.i<-R.inst[15,2]+R.inst[16,2]*28.5
Tmin21.R.i<-R.inst[11,2]+R.inst[12,2]*18.5
Tmin26.R.i<-R.inst[11,2]+R.inst[12,2]*23.5
Tmin31.R.i<-R.inst[11,2]+R.inst[12,2]*28.5
Topt21.R.i<-R.inst[13,2]+R.inst[14,2]*18.5
Topt26.R.i<-R.inst[13,2]+R.inst[14,2]*23.5
Topt31.R.i<-R.inst[13,2]+R.inst[14,2]*28.5

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
#Plot Supplementary Figure 11
####

#PDF dimension is 9x7 inches  

#Plotting Settings
par(pty="s")
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

##
#S11a
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.M.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin21.M.o,Topt21.M.o,Tmax21.M.o,25)
ins.25<-beta(1,Tmin21.M.i,Topt21.M.i,Tmax21.M.i,25)
ins.fit<-beta(1/ins.25,Tmin21.M.i,Topt21.M.i,Tmax21.M.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin21.M.o,Topt21.M.o,Tmax21.M.o,x),from=Tmin21.M.o,to=Tmax21.M.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[16:462]/x1.fit[251]~temp.sim[16:462],type="l",lwd=1.2,col="black",xlim=c(0,Tmax21.M.i),lty=1)
points(x2.slow.fit[16:462]/x2.slow.fit[251]~temp.sim[16:462],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax21.M.i))
points(x3.slow.fit[16:462]/x3.slow.fit[251]~temp.sim[16:462],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax21.M.i))
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)
legend("top",legend=c(expression('Measured (6 '*degree*'C hr'^-1*')'),
                          expression('Simulated (6 '*degree*'C hr'^-1*')'),
                          expression('Simulated (3 '*degree*'C hr'^-1*')'),
                          expression('Simulated (2 '*degree*'C hr'^-1*')')),
       col=c("black","black","royalblue1","yellowgreen"),lwd=c(1.2,1.2,1.2,1.2),lty=c(2,1,1,1),bty="n",
       y.intersp = 0.7,cex=1.1,seg.len=1.8,x.intersp = 0.5)

##
#S11b
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.A.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin21.A.o,Topt21.A.o,Tmax21.A.o,25)
ins.25<-beta(1,Tmin21.A.i,Topt21.A.i,Tmax21.A.i,25)
ins.fit<-beta(1/ins.25,Tmin21.A.i,Topt21.A.i,Tmax21.A.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin21.A.o,Topt21.A.o,Tmax21.A.o,x),from=Tmin21.A.o,to=Tmax21.A.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[1:427]/x1.fit[251]~temp.sim[1:427],type="l",lwd=1.2,col="black",xlim=c(0,Tmax21.A.i),lty=1)
points(x2.slow.fit[1:427]/x2.slow.fit[251]~temp.sim[1:427],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax21.A.i))
points(x3.slow.fit[1:427]/x3.slow.fit[251]~temp.sim[1:427],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax21.A.i))
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)

##
#S11c
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.G.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin21.G.o,Topt21.G.o,Tmax21.G.o,25)
ins.25<-beta(1,Tmin21.G.i,Topt21.G.i,Tmax21.G.i,25)
ins.fit<-beta(1/ins.25,Tmin21.G.i,Topt21.G.i,Tmax21.G.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin21.G.o,Topt21.G.o,Tmax21.G.o,x),from=Tmin21.G.o,to=Tmax21.G.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[22:468]/x1.fit[251]~temp.sim[22:468],type="l",lwd=1.2,col="black",xlim=c(0,Tmax21.G.i),lty=1)
points(x2.slow.fit[22:468]/x2.slow.fit[251]~temp.sim[22:468],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax21.G.i))
points(x3.slow.fit[22:468]/x3.slow.fit[251]~temp.sim[22:468],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax21.G.i))
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)

##
#S11d
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt21.R.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin21.R.o,Topt21.R.o,Tmax21.R.o,25)
ins.25<-beta(1,Tmin21.R.i,Topt21.R.i,Tmax21.R.i,25)ins.fit<-beta(1/ins.25,Tmin21.i,Topt21.i,Tmax21.i,temp.sim)
ins.fit<-beta(1/ins.25,Tmin21.R.i,Topt21.R.i,Tmax21.R.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin21.R.o,Topt21.R.o,Tmax21.R.o,x),from=Tmin21.R.o,to=Tmax21.R.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[25:473]/x1.fit[251]~temp.sim[25:473],type="l",lwd=1.2,col="black",xlim=c(0,Tmax21.R.i),lty=1)
points(x2.slow.fit[25:473]/x2.slow.fit[251]~temp.sim[25:473],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax21.R.i))
points(x3.slow.fit[25:473]/x3.slow.fit[251]~temp.sim[25:473],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax21.R.i))
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)

##
#S11e
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.M.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin26.M.o,Topt26.M.o,Tmax26.M.o,25)
ins.25<-beta(1,Tmin26.M.i,Topt26.M.i,Tmax26.M.i,25)
ins.fit<-beta(1/ins.25,Tmin26.M.i,Topt26.M.i,Tmax26.M.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin26.M.o,Topt26.M.o,Tmax26.M.o,x),from=Tmin26.M.o,to=Tmax26.M.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[77:459]/x1.fit[251]~temp.sim[77:459],type="l",lwd=1.2,col="black",xlim=c(0,Tmax26.M.i),lty=1)
points(x2.slow.fit[77:459]/x2.slow.fit[251]~temp.sim[77:459],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax26.M.i))
points(x3.slow.fit[77:459]/x3.slow.fit[251]~temp.sim[77:459],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax26.M.i))
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

##
#S11f
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.A.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin26.A.o,Topt26.A.o,Tmax26.A.o,25)
ins.25<-beta(1,Tmin26.A.i,Topt26.A.i,Tmax26.A.i,25)
ins.fit<-beta(1/ins.25,Tmin26.A.i,Topt26.A.i,Tmax26.A.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin26.A.o,Topt26.A.o,Tmax26.A.o,x),from=Tmin26.A.o,to=Tmax26.A.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[1:430]/x1.fit[251]~temp.sim[1:430],type="l",lwd=1.2,col="black",xlim=c(0,Tmax26.A.i),lty=1)
points(x2.slow.fit[1:430]/x2.slow.fit[251]~temp.sim[1:430],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax26.A.i))
points(x3.slow.fit[1:430]/x3.slow.fit[251]~temp.sim[1:430],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax26.A.i))
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

##
#S11g
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.G.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin26.G.o,Topt26.G.o,Tmax26.G.o,25)
ins.25<-beta(1,Tmin26.G.i,Topt26.G.i,Tmax26.G.i,25)
ins.fit<-beta(1/ins.25,Tmin26.G.i,Topt26.G.i,Tmax26.G.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin26.G.o,Topt26.G.o,Tmax26.G.o,x),from=Tmin26.G.o,to=Tmax26.G.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[77:468]/x1.fit[251]~temp.sim[77:468],type="l",lwd=1.2,col="black",xlim=c(0,Tmax26.G.i),lty=1)
points(x2.slow.fit[77:468]/x2.slow.fit[251]~temp.sim[77:468],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax26.G.i))
points(x3.slow.fit[77:468]/x3.slow.fit[251]~temp.sim[77:468],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax26.G.i))
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

##
#S11h
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt26.R.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin26.R.o,Topt26.R.o,Tmax26.R.o,25)
ins.25<-beta(1,Tmin26.R.i,Topt26.R.i,Tmax26.R.i,25)
ins.fit<-beta(1/ins.25,Tmin26.R.i,Topt26.R.i,Tmax26.R.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin26.R.o,Topt26.R.o,Tmax26.R.o,x),from=Tmin26.R.o,to=Tmax26.R.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[66:475]/x1.fit[251]~temp.sim[66:475],type="l",lwd=1.2,col="black",xlim=c(0,Tmax26.R.i),lty=1)
points(x2.slow.fit[66:475]/x2.slow.fit[251]~temp.sim[66:475],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax26.R.i))
points(x3.slow.fit[66:475]/x3.slow.fit[251]~temp.sim[66:475],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax26.R.i))
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

##
#S12i
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.M.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin31.M.o,Topt31.M.o,Tmax31.M.o,25)
ins.25<-beta(1,Tmin31.M.i,Topt31.M.i,Tmax31.M.i,25)
ins.fit<-beta(1/ins.25,Tmin31.M.i,Topt31.M.i,Tmax31.M.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=T,cex.axis=1.5,las=1)
curve(beta(1/og.25,Tmin31.M.o,Topt31.M.o,Tmax31.M.o,x),from=Tmin31.M.o,to=Tmax31.M.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[135:462]/x1.fit[251]~temp.sim[135:462],type="l",lwd=1.2,col="black",xlim=c(0,Tmax31.M.i))
points(x2.slow.fit[135:462]/x2.slow.fit[251]~temp.sim[135:462],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax31.M.i))
points(x3.slow.fit[135:462]/x3.slow.fit[251]~temp.sim[135:462],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax31.M.i))
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

##
#S11j
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.A.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin31.A.o,Topt31.A.o,Tmax31.A.o,25)
ins.25<-beta(1,Tmin31.A.i,Topt31.A.i,Tmax31.A.i,25)
ins.fit<-beta(1/ins.25,Tmin31.A.i,Topt31.A.i,Tmax31.A.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
curve(beta(1/og.25,Tmin31.A.o,Topt31.A.o,Tmax31.A.o,x),from=Tmin31.A.o,to=Tmax31.A.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[1:434]/x1.fit[251]~temp.sim[1:434],type="l",lwd=1.2,col="black",xlim=c(0,Tmax31.A.i))
points(x2.slow.fit[1:434]/x2.slow.fit[251]~temp.sim[1:434],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax31.A.i))
points(x3.slow.fit[1:434]/x3.slow.fit[251]~temp.sim[1:434],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax31.A.i))
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

##
#S11k
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.G.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin31.G.o,Topt31.G.o,Tmax31.G.o,25)
ins.25<-beta(1,Tmin31.G.i,Topt31.G.i,Tmax31.G.i,25)
ins.fit<-beta(1/ins.25,Tmin31.G.i,Topt31.G.i,Tmax31.G.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
curve(beta(1/og.25,Tmin31.G.o,Topt31.G.o,Tmax31.G.o,x),from=Tmin31.G.o,to=Tmax31.G.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[132:468]/x1.fit[251]~temp.sim[132:468],type="l",lwd=1.2,col="black",xlim=c(0,Tmax31.G.i))
points(x2.slow.fit[132:468]/x2.slow.fit[251]~temp.sim[132:468],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax31.G.i))
points(x3.slow.fit[132:468]/x3.slow.fit[251]~temp.sim[132:468],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax31.G.i))
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

##
#11l
##

#Set measured optimum temperature for
#reparameterizing temperature to tau in temporal-decline model (equation S7)
Topt <- Topt31.R.o
tau.sim <- seq(0,50,0.1) - Topt 

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

og.25<-beta(1,Tmin31.R.o,Topt31.R.o,Tmax31.R.o,25)
ins.25<-beta(1,Tmin31.R.i,Topt31.R.i,Tmax31.R.i,25)
ins.fit<-beta(1/ins.25,Tmin31.R.i,Topt31.R.i,Tmax31.R.i,temp.sim)

x1.fit<-ins.fit*SNF_realized.a
x2.slow.fit<-ins.fit*SNF_realized.b
x3.slow.fit<-ins.fit*SNF_realized.c

plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3),cex.axis=1.2,xlim=c(5,50),xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
curve(beta(1/og.25,Tmin31.R.o,Topt31.R.o,Tmax31.R.o,x),from=Tmin31.R.o,to=Tmax31.R.o,lwd=1.2,col="black",lty=2,add=TRUE)
points(x1.fit[109:476]/x1.fit[251]~temp.sim[109:476],type="l",lwd=1.2,col="black",xlim=c(0,Tmax31.R.i))
points(x2.slow.fit[109:476]/x2.slow.fit[251]~temp.sim[109:476],type="l",lwd=1.2,col="royalblue1",xlim=c(0,Tmax31.R.i))
points(x3.slow.fit[109:476]/x3.slow.fit[251]~temp.sim[109:476],type="l",lwd=1.2,col="yellowgreen",xlim=c(0,Tmax31.R.i))
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression('N-fixation (normalized to 1 at 25 '*degree*'C)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)