###############################################################################################################
###############################################################################################################
#This script fits piecewise functions (equation S4) to N-fixation data at constant temperatures
###############################################################################################################
###############################################################################################################

####
#Set Working Directory to "SNF_temporal_decline" folder
####

###
#Morella
###

##
#21:15 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
M21_25C<-read.csv("MOCE21_long25.csv")

#Plot data
plot(M21_25C$Time,M21_25C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE21; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point (only works here if I let max fall and max rise differ)
Break<-seq(0.59,2,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],p*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M21_25C[1:22800,],
             start=list(a=2600,b=0.053,c=0.79,d=0.87,p=2000,q=38))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.78
plot(sig~Break,type="l")

#
#30 deg. C
#

#Read in data
M21_30C<-read.csv("MOCE21_long30.csv")

#Plot data
plot(M21_30C$Time,M21_30C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE21; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.04,3.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M21_30C,
             start=list(a=2300,b=0.002,c=1.3,d=0.75,q=12))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.14
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
M21_35C<-read.csv("MOCE21_long35.csv")

#Plot data
plot(M21_35C$Time,M21_35C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE21; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.45,3,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M21_35C,
             start=list(a=3300,b=0.003,c=1.7,d=0.6,q=3))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #2.44
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
M21_40C<-read.csv("MOCE21_long40.csv")

#Plot data
plot(M21_40C$Time,M21_40C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE21; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,0.3,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M21_40C,
             start=list(a=1200,b=0.5,c=1.2,d=0.04,q=6))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.15
plot(sig~Break,type="l")

##
#26:20 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
M26_25C<-read.csv("MOCE26_long25.csv")

#Plot data
plot(M26_25C$Time,M26_25C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE26; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.49,1.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M26_25C[1200:length(M26_25C$Time),],
             start=list(a=6730,b=0.4,c=.53,d=0.81,q=6.4))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #1.31
plot(sig~Break,type="l")

#
#30 deg. C
#

#Read in data
M26_30C<-read.csv("MOCE26_long30.csv")

#Plot data
plot(M26_30C$Time,M26_30C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE26; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.41,0.8,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M26_30C[1350:length(M26_30C$Time),],
             start=list(a=7000,b=0.025,c=.84,d=0.52,q=2))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.78
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
M26_35C<-read.csv("MOCE26_long35.csv")

#Plot data
plot(M26_35C$Time,M26_35C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE26; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")
#No decline

#
#40 deg. C
#

#Read in data
M26_40C<-read.csv("MOCE26_long40.csv")

#Plot data
plot(M26_40C$Time,M26_40C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE26; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.2,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M26_40C[600:length(M26_40C$Time),],
             start=list(a=6050,b=0.25,c=.98,d=0.34,q=10))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.32
plot(sig~Break,type="l")

##
#31:25 deg. C growing temperature
##

#
#30 deg. C
#

#Read in data
M31_30C<-read.csv("MOCE31_long30.csv")

#Plot data
plot(M31_30C$Time,M31_30C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE31; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(M31_30C$Vmax)))
#No decline

#
#35 deg. C
#

#Read in data
M31_35C<-read.csv("MOCE31_long35.csv") 

#Plot data
plot(M31_35C$Time,M31_35C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE31; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(M31_35C$Vmax)))

#This finds best break point
Break<-c(seq(0.06,2.71,0.01),seq(2.8,4,0.01))
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M31_35C,
             start=list(a=2000,b=0.00001,c=2.2,d=0.83,q=24))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #3.9
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
M31_40C<-read.csv("MOCE31_long40.csv")

#Plot data
plot(M31_40C$Time,M31_40C$Vmax,pch=20,cex=0.1,col="blue",main="MOCE31; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(M31_40C$Vmax)))

#This finds best break point
Break<-c(seq(0.09,0.42,0.01),seq(0.44,0.97,0.01))
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=M31_40C,
             start=list(a=2420,b=0.0000001,c=4.3,d=0.97,q=70))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.34
plot(sig~Break,type="l")

###
#Alnus
###

##
#21:15 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
A21_25C<-read.csv("ALRU21_long25.csv")

#Plot data
plot(A21_25C$Time,A21_25C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU21; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")
#No decline

#
#30 deg. C
#

#Read in data
A21_30C<-read.csv("ALRU21_long30.csv")

#Plot data
plot(A21_30C$Time,A21_30C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU21; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-c(seq(0.3,1.54,0.01),seq(1.7,2.38,0.01),seq(2.43,3.5,0.01))
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A21_30C,
             start=list(a=6800,b=0.008,c=1,d=0.82,q=15))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #2.98
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
A21_35C<-read.csv("ALRU21_long35.csv")

#Plot data
plot(A21_35C$Time,A21_35C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU21; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,2,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A21_35C,
             start=list(a=1160,b=0.03,c=1.1,d=230/1160,q=9.1))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.51
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
A21_40C<-read.csv("ALRU21_long40.csv")

#Plot data
plot(A21_40C$Time,A21_40C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU21; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A21_40C,
             start=list(a=800,b=0.8,c=0.9,d=0.19,q=13))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.21
plot(sig~Break,type="l")

##
#26:20 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
A26_25C<-read.csv("ALRU26_long25.csv")

#Plot data
plot(A26_25C$Time,A26_25C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU26; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")
#No decline

#
#30 deg. C
#

#Read in data
A26_30C<-read.csv("ALRU26_long30.csv")

#Plot data
plot(A26_30C$Time,A26_30C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU26; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-c(seq(0.02,1.7,0.01),seq(1.8,3.2,0.01))
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A26_30C,
             start=list(a=2300,b=0.001,c=1.7,d=0.85,q=10))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #2.88
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
A26_35C<-read.csv("ALRU26_long35.csv")

#Plot data
plot(A26_35C$Time,A26_35C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU26; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,3,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A26_35C,
             start=list(a=170,b=0.005,c=1.3,d=0.84,q=26))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #1.83
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
A26_40C<-read.csv("ALRU26_long40.csv")

#Plot data
plot(A26_40C$Time,A26_40C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU26; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,0.91,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A26_40C,
             start=list(a=4300,b=0.16,c=3.5,d=0.09,q=25))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.07
plot(sig~Break,type="l")

##
#31:25 deg. C growing temperature
##

#
#30 deg. C
#

#Read in data
A31_30C<-read.csv("ALRU31_long30.csv")

#Plot data
plot(A31_30C$Time,A31_30C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU31; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")
#No decline

#
#35 deg. C
#

#Read in data
A31_35C<-read.csv("ALRU31_long35.csv")

#Plot data
plot(A31_35C$Time,A31_35C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU31; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,1.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A31_35C,
             start=list(a=3650,b=0.02,c=0.9,d=0.23,q=33))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.12
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
A31_40C<-read.csv("ALRU31_long40.csv")

#Plot data
plot(A31_40C$Time,A31_40C$Vmax,pch=20,cex=0.1,col="blue",main="ALRU31; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=A31_40C,
             start=list(a=2600,b=0.04,c=4,d=0.01,q=19))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.35
plot(sig~Break,type="l")

###
#Gliricidia
###

##
#21:15 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
G21_25C<-read.csv("GLSE21_long25.csv")

#Plot data
plot(G21_25C$Time,G21_25C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE21; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G21_25C$Vmax)))
#No decline

#
#30 deg. C
#

#Read in data
G21_30C<-read.csv("GLSE21_long30.csv")

#Plot data
plot(G21_30C$Time,G21_30C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE21; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G21_30C$Vmax)))

#This finds best break point
Break<-seq(0.02,2.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G21_30C,
             start=list(a=820,b=0.02,c=1,d=0.64,q=31))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #1.5
plot(sig~Break,type="l")

#
#30 deg. C
#

#Read in data
G21_35C<-read.csv("GLSE21_long35.csv")

#Plot data
plot(G21_35C$Time,G21_35C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE21; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G21_35C$Vmax)))
#No decline

#
#40 deg. C
#

#Read in data
G21_40C<-read.csv("GLSE21_long40.csv")

#Plot data
plot(G21_40C$Time,G21_40C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE21; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G21_40C$Vmax)))

#This finds best break point
Break<-seq(0.02,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G21_40C,
             start=list(a=1400,b=0.1,c=1,d=0.53,q=100))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.03
plot(sig~Break,type="l")

##
#26:20 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
G26_25C<-read.csv("GLSE26_long25.csv")

#Plot data
plot(G26_25C$Time,G26_25C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE26; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")
#No decline

#
#30 deg. C
#

#Read in data
G26_30C<-read.csv("GLSE26_long30.csv")

#Plot data
plot(G26_30C$Time,G26_30C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE26; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G26_30C,
             start=list(a=2400,b=0.067,c=.69,d=0.45,q=11))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.24
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
G26_35C<-read.csv("GLSE26_long35.csv")

#Plot data
plot(G26_35C$Time,G26_35C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE26; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.11,0.7,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G26_35C,
             start=list(a=2600,b=0.53,c=.86,d=0.25,q=11))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.34
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
G26_40C<-read.csv("GLSE26_long40.csv")

#Plot data
plot(G26_40C$Time,G26_40C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE26; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,0.48,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G26_40C,
             start=list(a=5100,b=0.61,c=1.2,d=0.015,q=11))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.07
plot(sig~Break,type="l")

##
#26:20 deg. C growing temperature
##

#
#30 deg. C
#

#Read in data
G31_30C<-read.csv("GLSE31_long30.csv")

#Plot data
plot(G31_30C$Time,G31_30C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE31; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G31_30C$Vmax)))
#No decrease

#
#35 deg. C
#

#Read in data
G31_35C<-read.csv("GLSE31_long35.csv")

#Plot data
plot(G31_35C$Time,G31_35C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE31; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G31_35C$Vmax)))

#This finds best break point
Break<-c(seq(1,2.56,0.01),seq(2.75,3.21,0.01),seq(3.23,4,0.01))
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G31_35C,
             start=list(a=1660,b=0.00001,c=2.2,d=0.83,q=23))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #3.9
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
G31_40C<-read.csv("GLSE31_long40.csv")

#Plot data
plot(G31_40C$Time,G31_40C$Vmax,pch=20,cex=0.1,col="blue",main="GLSE31; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(G31_40C$Vmax)))

#This finds best break point
Break<-seq(0.02,0.38,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=G31_40C,
             start=list(a=1240,b=0.26,c=0.68,d=0.52,q=127))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.18
plot(sig~Break,type="l")

###
#Robinia
###

##
#21:15 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
R21_25C<-read.csv("ROPS21_long25.csv")

#Plot data
plot(R21_25C$Time,R21_25C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS21; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R21_25C$Vmax)))

#This finds best break point
Break<-seq(0.02,4.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R21_25C,
             start=list(a=2900,b=0.000038,c=1.5,d=0.53,q=17))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #3.93
plot(sig~Break,type="l")

#
#30 deg. C
#

#Read in data
R21_30C<-read.csv("ROPS21_long30.csv")

#Plot data
plot(R21_30C$Time,R21_30C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS21; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R21_30C$Vmax)))

#This finds best break point
Break<-seq(0.02,1.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R21_30C,
             start=list(a=4500,b=0.02,c=1,d=0.29,q=7))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #1.13
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
R21_35C<-read.csv("ROPS21_long35.csv")

#Plot data
plot(R21_35C$Time,R21_35C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS21; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R21_35C$Vmax)))

#This finds best break point
Break<-seq(0.02,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R21_35C,
             start=list(a=5000,b=0.02,c=1.3,d=0.1,q=24))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.15
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
R21_40C<-read.csv("ROPS21_long40.csv")

#Plot data
plot(R21_40C$Time,R21_40C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS21; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R21_40C$Vmax)))

#This finds best break point
Break<-seq(0.02,0.15,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R21_40C,
             start=list(a=1920,b=0.009,c=30,d=0.08,q=28))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.08
plot(sig~Break,type="l")

##
#26:20 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
R26_25C<-read.csv("ROPS26_long25.csv")

#Plot data
plot(R26_25C$Time,R26_25C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS26; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")
#No decline

#
#30 deg. C
#

#Read in data
R26_30C<-read.csv("ROPS26_long30.csv")

#Plot data
plot(R26_30C$Time,R26_30C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS26; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.02,1.2,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R26_30C,
             start=list(a=2300,b=0.034,c=1.3,d=0.69,q=11))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.98
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
R26_35C<-read.csv("ROPS26_long35.csv")

#Plot data
plot(R26_35C$Time,R26_35C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS26; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.04,0.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R26_35C,
             start=list(a=3500,b=0.2,c=.9,d=0.35,q=11))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.04
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
R26_40C<-read.csv("ROPS26_long40.csv")

#Plot data
plot(R26_40C$Time,R26_40C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS26; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)")

#This finds best break point
Break<-seq(0.12,0.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R26_40C,
             start=list(a=2500,b=0.4,c=1,d=0.013,q=20))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.34
plot(sig~Break,type="l")

##
#31:25 deg. C growing temperature
##

#
#25 deg. C
#

#Read in data
R31_25C<-read.csv("ROPS31_long25.csv")

#Plot data
plot(R31_25C$Time,R31_25C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS31; 25C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R31_25C$Vmax)))

#This finds best break point
Break<-seq(0.1,1.18,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R31_25C,
             start=list(a=2600,b=6.3*10^-10,c=3.6,d=0.97,q=23))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.89
plot(sig~Break,type="l")

#
#30 deg. C
#

#Read in data
R31_30C<-read.csv("ROPS31_long30.csv")

#Plot data
plot(R31_30C$Time,R31_30C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS31; 30C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R31_30C$Vmax)))

#This finds best break point
Break<-c(seq(0.06,3.63,0.01),seq(3.7,4.51,0.01))
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R31_30C,
             start=list(a=1930,b=1.5*10^-9,c=3.5,d=0.96,q=35))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #4.49
plot(sig~Break,type="l")

#
#35 deg. C
#

#Read in data
R31_35C<-read.csv("ROPS31_long35.csv")

#Plot data
plot(R31_35C$Time,R31_35C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS31; 35C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R31_35C$Vmax)))

#This finds best break point
Break<-seq(0.02,1.5,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R31_35C,
             start=list(a=3700,b=0.06,c=1,d=0.26,q=19))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.15
plot(sig~Break,type="l")

#
#40 deg. C
#

#Read in data
R31_40C<-read.csv("ROPS31_long40.csv")

#Plot data
plot(R31_40C$Time,R31_40C$Vmax,pch=20,cex=0.1,col="blue",main="ROPS31; 40C",ylab="SNF (nmol C2H4 hr-1)",xlab="Time (hr)",ylim=c(0,max(R31_40C$Vmax)))

#This finds best break point
Break<-seq(0.02,1,0.01)
sig<-numeric(length(Break))
for(i in 1:length(Break)){
  model<-nls(Vmax~ifelse(Time<Break[i],a*(1-exp(-q*Time)),a*((1-d)/(1+b*exp(c*Time))+d)),
             data=R31_40C,
             start=list(a=1300,b=0.06,c=2,d=0.1,q=36))
  sig[i]<-summary(model)$sigma
}
Break[which.min(sig)] #0.08
plot(sig~Break,type="l")