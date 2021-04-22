###############################################################################################################
###############################################################################################################
#This script tests for the effect of growing temperature on leaf respiration in the light (RL),
#compares acclimation models for the effect of temperature on leaf respiration (Supplementary Table 10),
#and calculates the relative respiration rate at 25 deg. C for A-Ci curve calculations
###############################################################################################################
###############################################################################################################

#Load Necessary Package
library(bbmle)

####
#Read the data in from "Resp_Temp" folder
####

#RL data calculated using the Kok effect
RKok<-read.csv("RKok_Output.csv")

#Remove rows with NA values
RKok<-na.omit(RKok)

#Dark respiration (RD) data
Rd.data<-read.csv("Rd_Temp_data.csv")

#Subset RD data for each treatment (species x growing temperature)
Rd.A21<-Rd.data[Rd.data$Treatment=="ALRU21",]
Rd.A26<-Rd.data[Rd.data$Treatment=="ALRU26",]
Rd.A31<-Rd.data[Rd.data$Treatment=="ALRU31",]
Rd.M21<-Rd.data[Rd.data$Treatment=="MOCE21",]
Rd.M26<-Rd.data[Rd.data$Treatment=="MOCE26",]
Rd.M31<-Rd.data[Rd.data$Treatment=="MOCE31",]
Rd.G21<-Rd.data[Rd.data$Treatment=="GLSE21",]
Rd.G26<-Rd.data[Rd.data$Treatment=="GLSE26",]
Rd.G31<-Rd.data[Rd.data$Treatment=="GLSE31",]
Rd.R21<-Rd.data[Rd.data$Treatment=="ROPS21",]
Rd.R26<-Rd.data[Rd.data$Treatment=="ROPS26",]
Rd.R31<-Rd.data[Rd.data$Treatment=="ROPS31",]

#Subset RD data for each individual plant
Rd.A21a<-Rd.A21[Rd.A21$Replicate==1,]
Rd.A21b<-Rd.A21[Rd.A21$Replicate==2,]
Rd.A21c<-Rd.A21[Rd.A21$Replicate==3,]
Rd.A26a<-Rd.A26[Rd.A26$Replicate==1,]
Rd.A26b<-Rd.A26[Rd.A26$Replicate==2,]
Rd.A26c<-Rd.A26[Rd.A26$Replicate==3,]
Rd.A31a<-Rd.A31[Rd.A31$Replicate==1,]
Rd.A31b<-Rd.A31[Rd.A31$Replicate==2,]
Rd.A31c<-Rd.A31[Rd.A31$Replicate==3,]
Rd.G21a<-Rd.G21[Rd.G21$Replicate==1,]
Rd.G21b<-Rd.G21[Rd.G21$Replicate==2,]
Rd.G26a<-Rd.G26[Rd.G26$Replicate==1,]
Rd.G26b<-Rd.G26[Rd.G26$Replicate==2,]
Rd.G26c<-Rd.G26[Rd.G26$Replicate==3,]
Rd.G31a<-Rd.G31[Rd.G31$Replicate==1,]
Rd.G31b<-Rd.G31[Rd.G31$Replicate==2,]
Rd.G31c<-Rd.G31[Rd.G31$Replicate==3,]
Rd.M21a<-Rd.M21[Rd.M21$Replicate==1,]
Rd.M21b<-Rd.M21[Rd.M21$Replicate==2,]
Rd.M21c<-Rd.M21[Rd.M21$Replicate==3,]
Rd.M26a<-Rd.M26[Rd.M26$Replicate==1,]
Rd.M26b<-Rd.M26[Rd.M26$Replicate==2,]
Rd.M26c<-Rd.M26[Rd.M26$Replicate==3,]
Rd.M31a<-Rd.M31[Rd.M31$Replicate==1,]
Rd.M31b<-Rd.M31[Rd.M31$Replicate==2,]
Rd.M31c<-Rd.M31[Rd.M31$Replicate==3,]
Rd.R21a<-Rd.R21[Rd.R21$Replicate==1,]
Rd.R21b<-Rd.R21[Rd.R21$Replicate==2,]
Rd.R21c<-Rd.R21[Rd.R21$Replicate==3,]
Rd.R26a<-Rd.R26[Rd.R26$Replicate==1,]
Rd.R31a<-Rd.R31[Rd.R31$Replicate==1,]
Rd.R31b<-Rd.R31[Rd.R31$Replicate==2,]
Rd.R31c<-Rd.R31[Rd.R31$Replicate==3,]
Rd.R31d<-Rd.R31[Rd.R31$Replicate==4,]

###############################################################################################################
#Test for effect of growing temperature on RL
###############################################################################################################

##
#Morella
##

summary(lm(RKok$Rlight[18:26]~RKok$Growing.Temperature[18:26]))
#Significant; p=0.01445

#Calculate predicted RL at each growing temperature
-0.05244*18.5+1.72318 #0.75304
-0.05244*23.5+1.72318 #0.49084
-0.05244*28.5+1.72318 #0.22864

##
#Alnus
##

summary(lm(RKok$Rlight[1:9]~RKok$Growing.Temperature[1:9]))
#Not significant

#Calculate predicted RL across all growing temperatures
mean(RKok$Rlight[1:9]) #0.9382466

##
#Gliricidia
##

summary(lm(RKok$Rlight[10:17]~RKok$Growing.Temperature[10:17]))
#Not significant

#Calculate predicted RL across all growing temperatures
mean(RKok$Rlight[10:17]) #0.8394181

##
#Robinia
##

summary(lm(RKok$Rlight[27:34]~RKok$Growing.Temperature[27:34]))
#Not significant

#Calculate predicted RL across all growing temperatures
mean(RKok$Rlight[27:34]) #0.9832326

###############################################################################################################
#Compare acclimation models for the effect of temperature on leaf respiration
###############################################################################################################

####
#Define functions
####

#Acclimation of all parameters
norm.Topt.s.lin<-function(ymax,a,b,c,d,Tgrow,T){
  y<-ymax*exp(-(T-(a+b*Tgrow))^2/(2*(c+d*Tgrow)^2))
  y
}

#Acclimation of Topt
norm.Topt.lin<-function(ymax,a,b,s,Tgrow,T){
  y<-ymax*exp(-(T-(a+Tgrow))^2/(2*s^2))
  y
}

#Acclimation of s
norm.s.lin<-function(ymax,Topt,a,b,Tgrow,T){
  y<-ymax*exp(-(T-Topt)^2/(2*(a+b*Tgrow)^2))
  y
}

#No acclimation
norm<-function(ymax,Topt,s,T){
  y<-ymax*exp(-(T-Topt)^2/(2*s^2))
  y
}

####
#First model Morella and Alnus
#Gliricidia and Robinia are below
#This is because Gliricidia and Robinia are missing some data
#and therefore have negative log likelihood (NLL) functions
#that differ from Morella and Alnus
####

###
#NLL functions (for Morella and Alnus)
###

#Acclimation of all parameters
Resp_norm_Topt.s.lin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                  ymax26a,ymax26b,ymax26c,
                                  ymax31a,ymax31b,ymax31c,
                                  a,b,c,d,
                                  T21a,T21b,T21c,
                                  T26a,T26b,T26c,
                                  T31a,T31b,T31c,
                                  Respdat21a,Respdat21b,Respdat21c,
                                  Respdat26a,Respdat26b,Respdat26c,
                                  Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm.Topt.s.lin(ymax21a,a,b,c,d,18.5,T21a)
  Respmean21b <- norm.Topt.s.lin(ymax21b,a,b,c,d,18.5,T21b)
  Respmean21c <- norm.Topt.s.lin(ymax21c,a,b,c,d,18.5,T21c)
  Respmean26a <- norm.Topt.s.lin(ymax26a,a,b,c,d,23.5,T26a)
  Respmean26b <- norm.Topt.s.lin(ymax26b,a,b,c,d,23.5,T26b)
  Respmean26c <- norm.Topt.s.lin(ymax26c,a,b,c,d,23.5,T26c)
  Respmean31a <- norm.Topt.s.lin(ymax31a,a,b,c,d,28.5,T31a)
  Respmean31b <- norm.Topt.s.lin(ymax31b,a,b,c,d,28.5,T31b)
  Respmean31c <- norm.Topt.s.lin(ymax31c,a,b,c,d,28.5,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#Acclimation of Topt
Resp_norm_Topt.lin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                         ymax26a,ymax26b,ymax26c,
                                         ymax31a,ymax31b,ymax31c,
                                         a,b,s,
                                         T21a,T21b,T21c,
                                         T26a,T26b,T26c,
                                         T31a,T31b,T31c,
                                         Respdat21a,Respdat21b,Respdat21c,
                                         Respdat26a,Respdat26b,Respdat26c,
                                         Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm.Topt.lin(ymax21a,a,b,s,18.5,T21a)
  Respmean21b <- norm.Topt.lin(ymax21b,a,b,s,18.5,T21b)
  Respmean21c <- norm.Topt.lin(ymax21c,a,b,s,18.5,T21c)
  Respmean26a <- norm.Topt.lin(ymax26a,a,b,s,23.5,T26a)
  Respmean26b <- norm.Topt.lin(ymax26b,a,b,s,23.5,T26b)
  Respmean26c <- norm.Topt.lin(ymax26c,a,b,s,23.5,T26c)
  Respmean31a <- norm.Topt.lin(ymax31a,a,b,s,28.5,T31a)
  Respmean31b <- norm.Topt.lin(ymax31b,a,b,s,28.5,T31b)
  Respmean31c <- norm.Topt.lin(ymax31c,a,b,s,28.5,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#Acclimation of s
Resp_norm_s.lin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                       ymax26a,ymax26b,ymax26c,
                                       ymax31a,ymax31b,ymax31c,
                                       a,b,Topt,
                                       T21a,T21b,T21c,
                                       T26a,T26b,T26c,
                                       T31a,T31b,T31c,
                                       Respdat21a,Respdat21b,Respdat21c,
                                       Respdat26a,Respdat26b,Respdat26c,
                                       Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm.s.lin(ymax21a,Topt,a,b,18.5,T21a)
  Respmean21b <- norm.s.lin(ymax21b,Topt,a,b,18.5,T21b)
  Respmean21c <- norm.s.lin(ymax21c,Topt,a,b,18.5,T21c)
  Respmean26a <- norm.s.lin(ymax26a,Topt,a,b,23.5,T26a)
  Respmean26b <- norm.s.lin(ymax26b,Topt,a,b,23.5,T26b)
  Respmean26c <- norm.s.lin(ymax26c,Topt,a,b,23.5,T26c)
  Respmean31a <- norm.s.lin(ymax31a,Topt,a,b,28.5,T31a)
  Respmean31b <- norm.s.lin(ymax31b,Topt,a,b,28.5,T31b)
  Respmean31c <- norm.s.lin(ymax31c,Topt,a,b,28.5,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#NLL function without acclimation
Resp_norm_bin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                  ymax26a,ymax26b,ymax26c,
                                  ymax31a,ymax31b,ymax31c,
                                  Topt,s,
                                  T21a,T21b,T21c,
                                  T26a,T26b,T26c,
                                  T31a,T31b,T31c,
                                  Respdat21a,Respdat21b,Respdat21c,
                                  Respdat26a,Respdat26b,Respdat26c,
                                  Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm(ymax21a,Topt,s,T21a)
  Respmean21b <- norm(ymax21b,Topt,s,T21b)
  Respmean21c <- norm(ymax21c,Topt,s,T21c)
  Respmean26a <- norm(ymax26a,Topt,s,T26a)
  Respmean26b <- norm(ymax26b,Topt,s,T26b)
  Respmean26c <- norm(ymax26c,Topt,s,T26c)
  Respmean31a <- norm(ymax31a,Topt,s,T31a)
  Respmean31b <- norm(ymax31b,Topt,s,T31b)
  Respmean31c <- norm(ymax31c,Topt,s,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

###
#Maximum likelihood fits
###

##
#Morella
##

#Acclimation of all parameters
fit_Resp_norm_Topt.s.lin_MOCE <- mle2(Resp_norm_Topt.s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                ymax26a=8,ymax26b=8,ymax26c=8,
                                                                ymax31a=8,ymax31b=8,ymax31c=8,
                                                                a=65,b=0,
                                                                c=25,d=0),
                               data=list(T21a=Rd.M21a$Temperature,T21b=Rd.M21b$Temperature,T21c=Rd.M21c$Temperature,
                                         T26a=Rd.M26a$Temperature,T26b=Rd.M26b$Temperature,T26c=Rd.M26c$Temperature,
                                         T31a=Rd.M31a$Temperature,T31b=Rd.M31b$Temperature,T31c=Rd.M31c$Temperature,
                                         Respdat21a=Rd.M21a$Respiration,Respdat21b=Rd.M21b$Respiration,Respdat21c=Rd.M21c$Respiration,
                                         Respdat26a=Rd.M26a$Respiration,Respdat26b=Rd.M26b$Respiration,Respdat26c=Rd.M26c$Respiration,
                                         Respdat31a=Rd.M31a$Respiration,Respdat31b=Rd.M31b$Respiration,Respdat31c=Rd.M31c$Respiration),
                               control=list(maxit=20000))
summary(fit_Resp_norm_Topt.s.lin_MOCE)

#Acclimation of Topt
fit_Resp_norm_Topt.lin_MOCE <- mle2(Resp_norm_Topt.lin_normNLL,start=list(sdResp=-1,ymax21a=15,ymax21b=15,ymax21c=15,
                                                                              ymax26a=15,ymax26b=15,ymax26c=15,
                                                                              ymax31a=15,ymax31b=15,ymax31c=15,
                                                                              a=75,b=0,
                                                                              s=30),
                                      data=list(T21a=Rd.M21a$Temperature,T21b=Rd.M21b$Temperature,T21c=Rd.M21c$Temperature,
                                                T26a=Rd.M26a$Temperature,T26b=Rd.M26b$Temperature,T26c=Rd.M26c$Temperature,
                                                T31a=Rd.M31a$Temperature,T31b=Rd.M31b$Temperature,T31c=Rd.M31c$Temperature,
                                                Respdat21a=Rd.M21a$Respiration,Respdat21b=Rd.M21b$Respiration,Respdat21c=Rd.M21c$Respiration,
                                                Respdat26a=Rd.M26a$Respiration,Respdat26b=Rd.M26b$Respiration,Respdat26c=Rd.M26c$Respiration,
                                                Respdat31a=Rd.M31a$Respiration,Respdat31b=Rd.M31b$Respiration,Respdat31c=Rd.M31c$Respiration),
                                      control=list(maxit=20000))
summary(fit_Resp_norm_Topt.lin_MOCE)

#Acclimation of s
fit_Resp_norm_s.lin_MOCE <- mle2(Resp_norm_s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                          ymax26a=8,ymax26b=8,ymax26c=8,
                                                                          ymax31a=8,ymax31b=8,ymax31c=8,
                                                                          a=25,b=0,
                                                                          Topt=65),
                                    data=list(T21a=Rd.M21a$Temperature,T21b=Rd.M21b$Temperature,T21c=Rd.M21c$Temperature,
                                              T26a=Rd.M26a$Temperature,T26b=Rd.M26b$Temperature,T26c=Rd.M26c$Temperature,
                                              T31a=Rd.M31a$Temperature,T31b=Rd.M31b$Temperature,T31c=Rd.M31c$Temperature,
                                              Respdat21a=Rd.M21a$Respiration,Respdat21b=Rd.M21b$Respiration,Respdat21c=Rd.M21c$Respiration,
                                              Respdat26a=Rd.M26a$Respiration,Respdat26b=Rd.M26b$Respiration,Respdat26c=Rd.M26c$Respiration,
                                              Respdat31a=Rd.M31a$Respiration,Respdat31b=Rd.M31b$Respiration,Respdat31c=Rd.M31c$Respiration),
                                    control=list(maxit=20000))
summary(fit_Resp_norm_s.lin_MOCE)

#No acclimation
fit_Resp_norm_bin_MOCE <- mle2(Resp_norm_bin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                ymax26a=8,ymax26b=8,ymax26c=8,
                                                                ymax31a=8,ymax31b=8,ymax31c=8,
                                                                Topt=70,
                                                                s=25),
                               data=list(T21a=Rd.M21a$Temperature,T21b=Rd.M21b$Temperature,T21c=Rd.M21c$Temperature,
                                         T26a=Rd.M26a$Temperature,T26b=Rd.M26b$Temperature,T26c=Rd.M26c$Temperature,
                                         T31a=Rd.M31a$Temperature,T31b=Rd.M31b$Temperature,T31c=Rd.M31c$Temperature,
                                         Respdat21a=Rd.M21a$Respiration,Respdat21b=Rd.M21b$Respiration,Respdat21c=Rd.M21c$Respiration,
                                         Respdat26a=Rd.M26a$Respiration,Respdat26b=Rd.M26b$Respiration,Respdat26c=Rd.M26c$Respiration,
                                         Respdat31a=Rd.M31a$Respiration,Respdat31b=Rd.M31b$Respiration,Respdat31c=Rd.M31c$Respiration),
                               control=list(maxit=20000))
summary(fit_Resp_norm_bin_MOCE)

#Sample size
sum(c(length(Rd.M21$Respiration),length(Rd.M26$Respiration),length(Rd.M31$Respiration)))

#Delta AICc 
AICctab(fit_Resp_norm_Topt.s.lin_MOCE,fit_Resp_norm_Topt.lin_MOCE,fit_Resp_norm_s.lin_MOCE,fit_Resp_norm_bin_MOCE,nobs=277)

##
#Alnus
##

#Acclimation of all parameters
fit_Resp_norm_Topt.s.lin_ALRU <- mle2(Resp_norm_Topt.s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                              ymax26a=8,ymax26b=8,ymax26c=8,
                                                                              ymax31a=8,ymax31b=8,ymax31c=8,
                                                                              a=65,b=0,
                                                                              c=25,d=0),
                                      data=list(T21a=Rd.A21a$Temperature,T21b=Rd.A21b$Temperature,T21c=Rd.A21c$Temperature,
                                                T26a=Rd.A26a$Temperature,T26b=Rd.A26b$Temperature,T26c=Rd.A26c$Temperature,
                                                T31a=Rd.A31a$Temperature,T31b=Rd.A31b$Temperature,T31c=Rd.A31c$Temperature,
                                                Respdat21a=Rd.A21a$Respiration,Respdat21b=Rd.A21b$Respiration,Respdat21c=Rd.A21c$Respiration,
                                                Respdat26a=Rd.A26a$Respiration,Respdat26b=Rd.A26b$Respiration,Respdat26c=Rd.A26c$Respiration,
                                                Respdat31a=Rd.A31a$Respiration,Respdat31b=Rd.A31b$Respiration,Respdat31c=Rd.A31c$Respiration),
                                      control=list(maxit=20000))
summary(fit_Resp_norm_Topt.s.lin_ALRU)

#Acclimation of Topt
fit_Resp_norm_Topt.lin_ALRU <- mle2(Resp_norm_Topt.lin_normNLL,start=list(sdResp=-1,ymax21a=15,ymax21b=15,ymax21c=15,
                                                                          ymax26a=15,ymax26b=15,ymax26c=15,
                                                                          ymax31a=15,ymax31b=15,ymax31c=15,
                                                                          a=75,b=0,
                                                                          s=30),
                                    data=list(T21a=Rd.A21a$Temperature,T21b=Rd.A21b$Temperature,T21c=Rd.A21c$Temperature,
                                              T26a=Rd.A26a$Temperature,T26b=Rd.A26b$Temperature,T26c=Rd.A26c$Temperature,
                                              T31a=Rd.A31a$Temperature,T31b=Rd.A31b$Temperature,T31c=Rd.A31c$Temperature,
                                              Respdat21a=Rd.A21a$Respiration,Respdat21b=Rd.A21b$Respiration,Respdat21c=Rd.A21c$Respiration,
                                              Respdat26a=Rd.A26a$Respiration,Respdat26b=Rd.A26b$Respiration,Respdat26c=Rd.A26c$Respiration,
                                              Respdat31a=Rd.A31a$Respiration,Respdat31b=Rd.A31b$Respiration,Respdat31c=Rd.A31c$Respiration),
                                    control=list(maxit=20000))
summary(fit_Resp_norm_Topt.lin_ALRU)

#Acclimation of s
fit_Resp_norm_s.lin_ALRU <- mle2(Resp_norm_s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                    ymax26a=8,ymax26b=8,ymax26c=8,
                                                                    ymax31a=8,ymax31b=8,ymax31c=8,
                                                                    a=25,b=0,
                                                                    Topt=65),
                                 data=list(T21a=Rd.A21a$Temperature,T21b=Rd.A21b$Temperature,T21c=Rd.A21c$Temperature,
                                           T26a=Rd.A26a$Temperature,T26b=Rd.A26b$Temperature,T26c=Rd.A26c$Temperature,
                                           T31a=Rd.A31a$Temperature,T31b=Rd.A31b$Temperature,T31c=Rd.A31c$Temperature,
                                           Respdat21a=Rd.A21a$Respiration,Respdat21b=Rd.A21b$Respiration,Respdat21c=Rd.A21c$Respiration,
                                           Respdat26a=Rd.A26a$Respiration,Respdat26b=Rd.A26b$Respiration,Respdat26c=Rd.A26c$Respiration,
                                           Respdat31a=Rd.A31a$Respiration,Respdat31b=Rd.A31b$Respiration,Respdat31c=Rd.A31c$Respiration),
                                 control=list(maxit=20000))
summary(fit_Resp_norm_s.lin_ALRU)

#No acclimation
fit_Resp_norm_bin_ALRU <- mle2(Resp_norm_bin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                ymax26a=8,ymax26b=8,ymax26c=8,
                                                                ymax31a=8,ymax31b=8,ymax31c=8,
                                                                Topt=70,
                                                                s=25),
                               data=list(T21a=Rd.A21a$Temperature,T21b=Rd.A21b$Temperature,T21c=Rd.A21c$Temperature,
                                         T26a=Rd.A26a$Temperature,T26b=Rd.A26b$Temperature,T26c=Rd.A26c$Temperature,
                                         T31a=Rd.A31a$Temperature,T31b=Rd.A31b$Temperature,T31c=Rd.A31c$Temperature,
                                         Respdat21a=Rd.A21a$Respiration,Respdat21b=Rd.A21b$Respiration,Respdat21c=Rd.A21c$Respiration,
                                         Respdat26a=Rd.A26a$Respiration,Respdat26b=Rd.A26b$Respiration,Respdat26c=Rd.A26c$Respiration,
                                         Respdat31a=Rd.A31a$Respiration,Respdat31b=Rd.A31b$Respiration,Respdat31c=Rd.A31c$Respiration),
                               control=list(maxit=20000))
summary(fit_Resp_norm_bin_ALRU)

#Sample size
sum(c(length(Rd.A21$Respiration),length(Rd.A26$Respiration),length(Rd.A31$Respiration)))

#Delta AICc
AICctab(fit_Resp_norm_Topt.s.lin_ALRU,fit_Resp_norm_Topt.lin_ALRU,fit_Resp_norm_s.lin_ALRU,fit_Resp_norm_bin_ALRU,nobs=279)

###
#Gliricidia
###

##
#NLL functions
##

#Acclimation of all parameters
Resp_norm_Topt.s.lin_normNLL <- function(sdResp,ymax21a,ymax21b,
                                         ymax26a,ymax26b,ymax26c,
                                         ymax31a,ymax31b,ymax31c,
                                         a,b,c,d,
                                         T21a,T21b,
                                         T26a,T26b,T26c,
                                         T31a,T31b,T31c,
                                         Respdat21a,Respdat21b,
                                         Respdat26a,Respdat26b,Respdat26c,
                                         Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm.Topt.s.lin(ymax21a,a,b,c,d,18.5,T21a)
  Respmean21b <- norm.Topt.s.lin(ymax21b,a,b,c,d,18.5,T21b)
  Respmean26a <- norm.Topt.s.lin(ymax26a,a,b,c,d,23.5,T26a)
  Respmean26b <- norm.Topt.s.lin(ymax26b,a,b,c,d,23.5,T26b)
  Respmean26c <- norm.Topt.s.lin(ymax26c,a,b,c,d,23.5,T26c)
  Respmean31a <- norm.Topt.s.lin(ymax31a,a,b,c,d,28.5,T31a)
  Respmean31b <- norm.Topt.s.lin(ymax31b,a,b,c,d,28.5,T31b)
  Respmean31c <- norm.Topt.s.lin(ymax31c,a,b,c,d,28.5,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#Acclimation of Topt
Resp_norm_Topt.lin_normNLL <- function(sdResp,ymax21a,ymax21b,
                                       ymax26a,ymax26b,ymax26c,
                                       ymax31a,ymax31b,ymax31c,
                                       a,b,s,
                                       T21a,T21b,
                                       T26a,T26b,T26c,
                                       T31a,T31b,T31c,
                                       Respdat21a,Respdat21b,
                                       Respdat26a,Respdat26b,Respdat26c,
                                       Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm.Topt.lin(ymax21a,a,b,s,18.5,T21a)
  Respmean21b <- norm.Topt.lin(ymax21b,a,b,s,18.5,T21b)
  Respmean26a <- norm.Topt.lin(ymax26a,a,b,s,23.5,T26a)
  Respmean26b <- norm.Topt.lin(ymax26b,a,b,s,23.5,T26b)
  Respmean26c <- norm.Topt.lin(ymax26c,a,b,s,23.5,T26c)
  Respmean31a <- norm.Topt.lin(ymax31a,a,b,s,28.5,T31a)
  Respmean31b <- norm.Topt.lin(ymax31b,a,b,s,28.5,T31b)
  Respmean31c <- norm.Topt.lin(ymax31c,a,b,s,28.5,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#Acclimation of s
Resp_norm_s.lin_normNLL <- function(sdResp,ymax21a,ymax21b,
                                    ymax26a,ymax26b,ymax26c,
                                    ymax31a,ymax31b,ymax31c,
                                    a,b,Topt,
                                    T21a,T21b,
                                    T26a,T26b,T26c,
                                    T31a,T31b,T31c,
                                    Respdat21a,Respdat21b,
                                    Respdat26a,Respdat26b,Respdat26c,
                                    Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm.s.lin(ymax21a,Topt,a,b,18.5,T21a)
  Respmean21b <- norm.s.lin(ymax21b,Topt,a,b,18.5,T21b)
  Respmean26a <- norm.s.lin(ymax26a,Topt,a,b,23.5,T26a)
  Respmean26b <- norm.s.lin(ymax26b,Topt,a,b,23.5,T26b)
  Respmean26c <- norm.s.lin(ymax26c,Topt,a,b,23.5,T26c)
  Respmean31a <- norm.s.lin(ymax31a,Topt,a,b,28.5,T31a)
  Respmean31b <- norm.s.lin(ymax31b,Topt,a,b,28.5,T31b)
  Respmean31c <- norm.s.lin(ymax31c,Topt,a,b,28.5,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#No acclimation
Resp_norm_bin_normNLL <- function(sdResp,ymax21a,ymax21b,
                                  ymax26a,ymax26b,ymax26c,
                                  ymax31a,ymax31b,ymax31c,
                                  Topt,s,
                                  T21a,T21b,
                                  T26a,T26b,T26c,
                                  T31a,T31b,T31c,
                                  Respdat21a,Respdat21b,
                                  Respdat26a,Respdat26b,Respdat26c,
                                  Respdat31a,Respdat31b,Respdat31c){
  Respmean21a <- norm(ymax21a,Topt,s,T21a)
  Respmean21b <- norm(ymax21b,Topt,s,T21b)
  Respmean26a <- norm(ymax26a,Topt,s,T26a)
  Respmean26b <- norm(ymax26b,Topt,s,T26b)
  Respmean26c <- norm(ymax26c,Topt,s,T26c)
  Respmean31a <- norm(ymax31a,Topt,s,T31a)
  Respmean31b <- norm(ymax31b,Topt,s,T31b)
  Respmean31c <- norm(ymax31c,Topt,s,T31c)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26b,mean=Respmean26b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat26c,mean=Respmean26c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

##
#Maximum likelihood fits
##

#Acclimation of all parameters
fit_Resp_norm_Topt.s.lin_GLSE <- mle2(Resp_norm_Topt.s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,
                                                                              ymax26a=8,ymax26b=8,ymax26c=8,
                                                                              ymax31a=8,ymax31b=8,ymax31c=8,
                                                                              a=65,b=0,
                                                                              c=25,d=0),
                                      data=list(T21a=Rd.G21a$Temperature,T21b=Rd.G21b$Temperature,
                                                T26a=Rd.G26a$Temperature,T26b=Rd.G26b$Temperature,T26c=Rd.G26c$Temperature,
                                                T31a=Rd.G31a$Temperature,T31b=Rd.G31b$Temperature,T31c=Rd.G31c$Temperature,
                                                Respdat21a=Rd.G21a$Respiration,Respdat21b=Rd.G21b$Respiration,
                                                Respdat26a=Rd.G26a$Respiration,Respdat26b=Rd.G26b$Respiration,Respdat26c=Rd.G26c$Respiration,
                                                Respdat31a=Rd.G31a$Respiration,Respdat31b=Rd.G31b$Respiration,Respdat31c=Rd.G31c$Respiration),
                                      control=list(maxit=20000))
summary(fit_Resp_norm_Topt.s.lin_GLSE)

#Acclimation of Topt
fit_Resp_norm_Topt.lin_GLSE <- mle2(Resp_norm_Topt.lin_normNLL,start=list(sdResp=-1,ymax21a=6,ymax21b=6,
                                                                          ymax26a=6,ymax26b=6,ymax26c=6,
                                                                          ymax31a=6,ymax31b=6,ymax31c=6,
                                                                          a=60,b=-0.1,
                                                                          s=20),
                                    data=list(T21a=Rd.G21a$Temperature,T21b=Rd.G21b$Temperature,
                                              T26a=Rd.G26a$Temperature,T26b=Rd.G26b$Temperature,T26c=Rd.G26c$Temperature,
                                              T31a=Rd.G31a$Temperature,T31b=Rd.G31b$Temperature,T31c=Rd.G31c$Temperature,
                                              Respdat21a=Rd.G21a$Respiration,Respdat21b=Rd.G21b$Respiration,
                                              Respdat26a=Rd.G26a$Respiration,Respdat26b=Rd.G26b$Respiration,Respdat26c=Rd.G26c$Respiration,
                                              Respdat31a=Rd.G31a$Respiration,Respdat31b=Rd.G31b$Respiration,Respdat31c=Rd.G31c$Respiration),
                                    control=list(maxit=20000))
summary(fit_Resp_norm_Topt.lin_GLSE)

#Acclimation of s
fit_Resp_norm_s.lin_GLSE <- mle2(Resp_norm_s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,
                                                                    ymax26a=8,ymax26b=8,ymax26c=8,
                                                                    ymax31a=8,ymax31b=8,ymax31c=8,
                                                                    a=20,b=0.1,
                                                                    Topt=65),
                                 data=list(T21a=Rd.G21a$Temperature,T21b=Rd.G21b$Temperature,
                                           T26a=Rd.G26a$Temperature,T26b=Rd.G26b$Temperature,T26c=Rd.G26c$Temperature,
                                           T31a=Rd.G31a$Temperature,T31b=Rd.G31b$Temperature,T31c=Rd.G31c$Temperature,
                                           Respdat21a=Rd.G21a$Respiration,Respdat21b=Rd.G21b$Respiration,
                                           Respdat26a=Rd.G26a$Respiration,Respdat26b=Rd.G26b$Respiration,Respdat26c=Rd.G26c$Respiration,
                                           Respdat31a=Rd.G31a$Respiration,Respdat31b=Rd.G31b$Respiration,Respdat31c=Rd.G31c$Respiration),
                                 control=list(maxit=20000))
summary(fit_Resp_norm_s.lin_GLSE)

#No acclimation
fit_Resp_norm_bin_GLSE <- mle2(Resp_norm_bin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,
                                                                ymax26a=8,ymax26b=8,ymax26c=8,
                                                                ymax31a=8,ymax31b=8,ymax31c=8,
                                                                Topt=70,
                                                                s=25),
                               data=list(T21a=Rd.G21a$Temperature,T21b=Rd.G21b$Temperature,
                                         T26a=Rd.G26a$Temperature,T26b=Rd.G26b$Temperature,T26c=Rd.G26c$Temperature,
                                         T31a=Rd.G31a$Temperature,T31b=Rd.G31b$Temperature,T31c=Rd.G31c$Temperature,
                                         Respdat21a=Rd.G21a$Respiration,Respdat21b=Rd.G21b$Respiration,
                                         Respdat26a=Rd.G26a$Respiration,Respdat26b=Rd.G26b$Respiration,Respdat26c=Rd.G26c$Respiration,
                                         Respdat31a=Rd.G31a$Respiration,Respdat31b=Rd.G31b$Respiration,Respdat31c=Rd.G31c$Respiration),
                               control=list(maxit=20000))
summary(fit_Resp_norm_bin_GLSE)

#Sample size
sum(c(length(Rd.G21$Respiration),length(Rd.G26$Respiration),length(Rd.G31$Respiration)))

#Delta AICc
AICctab(fit_Resp_norm_Topt.s.lin_GLSE,fit_Resp_norm_Topt.lin_GLSE,fit_Resp_norm_s.lin_GLSE,fit_Resp_norm_bin_GLSE,nobs=248)

###
#Robinia
###

##
#NLL functions
##

#Acclimation of all parameters
Resp_norm_Topt.s.lin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                         ymax26a,
                                         ymax31a,ymax31b,ymax31c,ymax31d,
                                         a,b,c,d,
                                         T21a,T21b,T21c,
                                         T26a,
                                         T31a,T31b,T31c,T31d,
                                         Respdat21a,Respdat21b,Respdat21c,
                                         Respdat26a,
                                         Respdat31a,Respdat31b,Respdat31c,Respdat31d){
  Respmean21a <- norm.Topt.s.lin(ymax21a,a,b,c,d,18.5,T21a)
  Respmean21b <- norm.Topt.s.lin(ymax21b,a,b,c,d,18.5,T21b)
  Respmean21c <- norm.Topt.s.lin(ymax21c,a,b,c,d,18.5,T21c)
  Respmean26a <- norm.Topt.s.lin(ymax26a,a,b,c,d,23.5,T26a)
  Respmean31a <- norm.Topt.s.lin(ymax31a,a,b,c,d,28.5,T31a)
  Respmean31b <- norm.Topt.s.lin(ymax31b,a,b,c,d,28.5,T31b)
  Respmean31c <- norm.Topt.s.lin(ymax31c,a,b,c,d,28.5,T31c)
  Respmean31d <- norm.Topt.s.lin(ymax31d,a,b,c,d,28.5,T31d)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31d,mean=Respmean31d,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#Acclimation of Topt
Resp_norm_Topt.lin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                       ymax26a,
                                       ymax31a,ymax31b,ymax31c,ymax31d,
                                       a,b,s,
                                       T21a,T21b,T21c,
                                       T26a,
                                       T31a,T31b,T31c,T31d,
                                       Respdat21a,Respdat21b,Respdat21c,
                                       Respdat26a,
                                       Respdat31a,Respdat31b,Respdat31c,Respdat31d){
  Respmean21a <- norm.Topt.lin(ymax21a,a,b,s,18.5,T21a)
  Respmean21b <- norm.Topt.lin(ymax21b,a,b,s,18.5,T21b)
  Respmean21c <- norm.Topt.lin(ymax21c,a,b,s,18.5,T21c)
  Respmean26a <- norm.Topt.lin(ymax26a,a,b,s,23.5,T26a)
  Respmean31a <- norm.Topt.lin(ymax31a,a,b,s,28.5,T31a)
  Respmean31b <- norm.Topt.lin(ymax31b,a,b,s,28.5,T31b)
  Respmean31c <- norm.Topt.lin(ymax31c,a,b,s,28.5,T31c)
  Respmean31d <- norm.Topt.lin(ymax31d,a,b,s,28.5,T31d)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31d,mean=Respmean31d,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#Acclimation of s
Resp_norm_s.lin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                    ymax26a,
                                    ymax31a,ymax31b,ymax31c,ymax31d,
                                    a,b,Topt,
                                    T21a,T21b,T21c,
                                    T26a,
                                    T31a,T31b,T31c,
                                    Respdat21a,Respdat21b,Respdat21c,
                                    Respdat26a,
                                    Respdat31a,Respdat31b,Respdat31c,Respdat31d){
  Respmean21a <- norm.s.lin(ymax21a,Topt,a,b,18.5,T21a)
  Respmean21b <- norm.s.lin(ymax21b,Topt,a,b,18.5,T21b)
  Respmean21c <- norm.s.lin(ymax21c,Topt,a,b,18.5,T21c)
  Respmean26a <- norm.s.lin(ymax26a,Topt,a,b,23.5,T26a)
  Respmean31a <- norm.s.lin(ymax31a,Topt,a,b,28.5,T31a)
  Respmean31b <- norm.s.lin(ymax31b,Topt,a,b,28.5,T31b)
  Respmean31c <- norm.s.lin(ymax31c,Topt,a,b,28.5,T31c)
  Respmean31d <- norm.s.lin(ymax31d,Topt,a,b,28.5,T31d)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31d,mean=Respmean31d,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

#No acclimation
Resp_norm_bin_normNLL <- function(sdResp,ymax21a,ymax21b,ymax21c,
                                  ymax26a,
                                  ymax31a,ymax31b,ymax31c,ymax31d,
                                  Topt,s,
                                  T21a,T21b,T21c,
                                  T26a,
                                  T31a,T31b,T31c,T31d,
                                  Respdat21a,Respdat21b,Respdat21c,
                                  Respdat26a,
                                  Respdat31a,Respdat31b,Respdat31c,Respdat31d){
  Respmean21a <- norm(ymax21a,Topt,s,T21a)
  Respmean21b <- norm(ymax21b,Topt,s,T21b)
  Respmean21c <- norm(ymax21c,Topt,s,T21c)
  Respmean26a <- norm(ymax26a,Topt,s,T26a)
  Respmean31a <- norm(ymax31a,Topt,s,T31a)
  Respmean31b <- norm(ymax31b,Topt,s,T31b)
  Respmean31c <- norm(ymax31c,Topt,s,T31c)
  Respmean31d <- norm(ymax31d,Topt,s,T31d)
  -(sum(dnorm(Respdat21a,mean=Respmean21a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21b,mean=Respmean21b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat21c,mean=Respmean21c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat26a,mean=Respmean26a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31a,mean=Respmean31a,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31b,mean=Respmean31b,sd=exp(sdResp),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Respdat31c,mean=Respmean31c,sd=exp(sdResp),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Respdat31d,mean=Respmean31d,sd=exp(sdResp),log=TRUE),na.rm=TRUE))
}

##
#Maximum likelihood fits
##

#Acclimation of all parameters
fit_Resp_norm_Topt.s.lin_ROPS <- mle2(Resp_norm_Topt.s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                              ymax26a=8,
                                                                              ymax31a=8,ymax31b=8,ymax31c=8,ymax31d=8,
                                                                              a=60,b=0,
                                                                              c=20,d=0),
                                      data=list(T21a=Rd.R21a$Temperature,T21b=Rd.R21b$Temperature,T21c=Rd.R21c$Temperature,
                                                T26a=Rd.R26a$Temperature,
                                                T31a=Rd.R31a$Temperature,T31b=Rd.R31b$Temperature,T31c=Rd.R31c$Temperature,T31d=Rd.R31d$Temperature,
                                                Respdat21a=Rd.R21a$Respiration,Respdat21b=Rd.R21b$Respiration,Respdat21c=Rd.R21c$Respiration,
                                                Respdat26a=Rd.R26a$Respiration,
                                                Respdat31a=Rd.R31a$Respiration,Respdat31b=Rd.R31b$Respiration,Respdat31c=Rd.R31c$Respiration,Respdat31d=Rd.R31d$Respiration),
                                      control=list(maxit=20000))
summary(fit_Resp_norm_Topt.s.lin_ROPS)

#Acclimation of Topt
fit_Resp_norm_Topt.lin_ROPS <- mle2(Resp_norm_Topt.lin_normNLL,start=list(sdResp=-1,ymax21a=5,ymax21b=5,ymax21c=5,
                                                                          ymax26a=5,
                                                                          ymax31a=5,ymax31b=5,ymax31c=5,ymax31d=5,
                                                                          a=70,b=-0.5,
                                                                          s=20),
                                    data=list(T21a=Rd.R21a$Temperature,T21b=Rd.R21b$Temperature,T21c=Rd.R21c$Temperature,
                                              T26a=Rd.R26a$Temperature,
                                              T31a=Rd.R31a$Temperature,T31b=Rd.R31b$Temperature,T31c=Rd.R31c$Temperature,T31d=Rd.R31d$Temperature,
                                              Respdat21a=Rd.R21a$Respiration,Respdat21b=Rd.R21b$Respiration,Respdat21c=Rd.R21c$Respiration,
                                              Respdat26a=Rd.R26a$Respiration,
                                              Respdat31a=Rd.R31a$Respiration,Respdat31b=Rd.R31b$Respiration,Respdat31c=Rd.R31c$Respiration,Respdat31d=Rd.R31d$Respiration),
                                    control=list(maxit=20000))
summary(fit_Resp_norm_Topt.lin_ROPS)

#Acclimation of s
fit_Resp_norm_s.lin_ROPS <- mle2(Resp_norm_s.lin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                    ymax26a=8,
                                                                    ymax31a=8,ymax31b=8,ymax31c=8,ymax31d=8,
                                                                    a=25,b=-0.2,
                                                                    Topt=60),
                                 data=list(T21a=Rd.R21a$Temperature,T21b=Rd.R21b$Temperature,T21c=Rd.R21c$Temperature,
                                           T26a=Rd.R26a$Temperature,
                                           T31a=Rd.R31a$Temperature,T31b=Rd.R31b$Temperature,T31c=Rd.R31c$Temperature,T31d=Rd.R31d$Temperature,
                                           Respdat21a=Rd.R21a$Respiration,Respdat21b=Rd.R21b$Respiration,Respdat21c=Rd.R21c$Respiration,
                                           Respdat26a=Rd.R26a$Respiration,
                                           Respdat31a=Rd.R31a$Respiration,Respdat31b=Rd.R31b$Respiration,Respdat31c=Rd.R31c$Respiration,Respdat31d=Rd.R31d$Respiration),
                                 control=list(maxit=20000))
summary(fit_Resp_norm_s.lin_ROPS)

#No acclimation
fit_Resp_norm_bin_ROPS <- mle2(Resp_norm_bin_normNLL,start=list(sdResp=-1,ymax21a=8,ymax21b=8,ymax21c=8,
                                                                ymax26a=8,
                                                                ymax31a=8,ymax31b=8,ymax31c=8,ymax31d=8,
                                                                Topt=70,
                                                                s=25),
                               data=list(T21a=Rd.R21a$Temperature,T21b=Rd.R21b$Temperature,T21c=Rd.R21c$Temperature,
                                         T26a=Rd.R26a$Temperature,
                                         T31a=Rd.R31a$Temperature,T31b=Rd.R31b$Temperature,T31c=Rd.R31c$Temperature,T31d=Rd.R31d$Temperature,
                                         Respdat21a=Rd.R21a$Respiration,Respdat21b=Rd.R21b$Respiration,Respdat21c=Rd.R21c$Respiration,
                                         Respdat26a=Rd.R26a$Respiration,
                                         Respdat31a=Rd.R31a$Respiration,Respdat31b=Rd.R31b$Respiration,Respdat31c=Rd.R31c$Respiration,Respdat31d=Rd.R31d$Respiration),
                               control=list(maxit=20000))
summary(fit_Resp_norm_bin_ROPS)

#Sample size
sum(c(length(Rd.R21$Respiration),length(Rd.R26$Respiration),length(Rd.R31$Respiration)))

#Delta AICc
AICctab(fit_Resp_norm_Topt.s.lin_ROPS,fit_Resp_norm_Topt.lin_ROPS,fit_Resp_norm_s.lin_ROPS,fit_Resp_norm_bin_ROPS,nobs=248)

###############################################################################################################
#Calculate the relative respiration rate at 25 deg. C for A-Ci curve calculations
###############################################################################################################

#Coefficients
ft.Topt.s.lin.M<-coef(fit_Resp_norm_Topt.s.lin_MOCE)
ft.Topt.s.lin.A<-coef(fit_Resp_norm_Topt.s.lin_ALRU)
ft.Topt.s.lin.G<-coef(fit_Resp_norm_Topt.s.lin_GLSE)
ft.Topt.s.lin.R<-coef(fit_Resp_norm_Topt.s.lin_ROPS)

#Calculate relative respiration rate at 25 deg. C
Rd.25.M21<-norm.Topt.s.lin(1,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,25)
Rd.25.M26<-norm.Topt.s.lin(1,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,25)
Rd.25.M31<-norm.Topt.s.lin(1,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,25)
Rd.25.A21<-norm.Topt.s.lin(1,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,25)
Rd.25.A26<-norm.Topt.s.lin(1,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,25)
Rd.25.A31<-norm.Topt.s.lin(1,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,25)
Rd.25.G21<-norm.Topt.s.lin(1,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,25)
Rd.25.G26<-norm.Topt.s.lin(1,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,25)
Rd.25.G31<-norm.Topt.s.lin(1,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,25)
Rd.25.R21<-norm.Topt.s.lin(1,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,25)
Rd.25.R26<-norm.Topt.s.lin(1,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,25)
Rd.25.R31<-norm.Topt.s.lin(1,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,25)
