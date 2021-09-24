###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figures 5 and 6
###############################################################################################################
###############################################################################################################

####
#All figures were exported as PDFs from R Studio
#Dimensions for each figure are given with script specific to each figure
####

####
#Load Necessary Packages
####

library(bbmle)
library(MASS)

####
#Read in data
####

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
#Fit Topt and s acclimation model for the effect of temperature on leaf respiration
###############################################################################################################

####
#Define function
####

#Acclimation of all parameters
norm.Topt.s.lin<-function(ymax,a,b,c,d,Tgrow,T){
  y<-ymax*exp(-(T-(a+b*Tgrow))^2/(2*(c+d*Tgrow)^2))
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

####
#Calculate 95% CI
####

#Vector of measurement temperatures
Tsim<-seq(10,40,0.1)

##
#Morella
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m.r<-mvrnorm(1000,mu=ft.Topt.s.lin.M,Sigma=vcov(fit_Resp_norm_Topt.s.lin_MOCE))

#Respiration ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 6)
dist.m21.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m21.r<-rep(NA,length(Tsim))
high.m21.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m21.r[i,j]=norm.Topt.s.lin(1/Rd.25.M21,vmat.m.r[i,11],vmat.m.r[i,12],vmat.m.r[i,13],vmat.m.r[i,14],18.5,Tsim[j])
  }
  low.m21.r[j]<-quantile(na.omit(dist.m21.r[,j]),0.025)
  high.m21.r[j]<-quantile(na.omit(dist.m21.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 6)
dist.m26.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m26.r<-rep(NA,length(Tsim))
high.m26.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m26.r[i,j]=norm.Topt.s.lin(1/Rd.25.M26,vmat.m.r[i,11],vmat.m.r[i,12],vmat.m.r[i,13],vmat.m.r[i,14],23.5,Tsim[j])
  }
  low.m26.r[j]<-quantile(na.omit(dist.m26.r[,j]),0.025)
  high.m26.r[j]<-quantile(na.omit(dist.m26.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 6)
dist.m31.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m31.r<-rep(NA,length(Tsim))
high.m31.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m31.r[i,j]=norm.Topt.s.lin(1/Rd.25.M31,vmat.m.r[i,11],vmat.m.r[i,12],vmat.m.r[i,13],vmat.m.r[i,14],28.5,Tsim[j])
  }
  low.m31.r[j]<-quantile(na.omit(dist.m31.r[,j]),0.025)
  high.m31.r[j]<-quantile(na.omit(dist.m31.r[,j]),0.975)
}

##
#Alnus
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a.r<-mvrnorm(1000,mu=ft.Topt.s.lin.A,Sigma=vcov(fit_Resp_norm_Topt.s.lin_ALRU))

#Respiration ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 6)
dist.a21.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a21.r<-rep(NA,length(Tsim))
high.a21.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a21.r[i,j]=norm.Topt.s.lin(1/Rd.25.A21,vmat.a.r[i,11],vmat.a.r[i,12],vmat.a.r[i,13],vmat.a.r[i,14],18.5,Tsim[j])
  }
  low.a21.r[j]<-quantile(na.omit(dist.a21.r[,j]),0.025)
  high.a21.r[j]<-quantile(na.omit(dist.a21.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 6)
dist.a26.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a26.r<-rep(NA,length(Tsim))
high.a26.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a26.r[i,j]=norm.Topt.s.lin(1/Rd.25.A26,vmat.a.r[i,11],vmat.a.r[i,12],vmat.a.r[i,13],vmat.a.r[i,14],23.5,Tsim[j])
  }
  low.a26.r[j]<-quantile(na.omit(dist.a26.r[,j]),0.025)
  high.a26.r[j]<-quantile(na.omit(dist.a26.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 6)
dist.a31.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a31.r<-rep(NA,length(Tsim))
high.a31.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a31.r[i,j]=norm.Topt.s.lin(1/Rd.25.A31,vmat.a.r[i,11],vmat.a.r[i,12],vmat.a.r[i,13],vmat.a.r[i,14],28.5,Tsim[j])
  }
  low.a31.r[j]<-quantile(na.omit(dist.a31.r[,j]),0.025)
  high.a31.r[j]<-quantile(na.omit(dist.a31.r[,j]),0.975)
}

##
#Gliricidia
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g.r<-mvrnorm(1000,mu=ft.Topt.s.lin.G,Sigma=vcov(fit_Resp_norm_Topt.s.lin_GLSE))

#Respiration ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 6)
dist.g21.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g21.r<-rep(NA,length(Tsim))
high.g21.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g21.r[i,j]=norm.Topt.s.lin(1/Rd.25.G21,vmat.g.r[i,10],vmat.g.r[i,11],vmat.g.r[i,12],vmat.g.r[i,13],18.5,Tsim[j])
  }
  low.g21.r[j]<-quantile(na.omit(dist.g21.r[,j]),0.025)
  high.g21.r[j]<-quantile(na.omit(dist.g21.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 6)
dist.g26.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g26.r<-rep(NA,length(Tsim))
high.g26.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g26.r[i,j]=norm.Topt.s.lin(1/Rd.25.G26,vmat.g.r[i,10],vmat.g.r[i,11],vmat.g.r[i,12],vmat.g.r[i,13],23.5,Tsim[j])
  }
  low.g26.r[j]<-quantile(na.omit(dist.g26.r[,j]),0.025)
  high.g26.r[j]<-quantile(na.omit(dist.g26.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 6)
dist.g31.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g31.r<-rep(NA,length(Tsim))
high.g31.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g31.r[i,j]=norm.Topt.s.lin(1/Rd.25.G31,vmat.g.r[i,10],vmat.g.r[i,11],vmat.g.r[i,12],vmat.g.r[i,13],28.5,Tsim[j])
  }
  low.g31.r[j]<-quantile(na.omit(dist.g31.r[,j]),0.025)
  high.g31.r[j]<-quantile(na.omit(dist.g31.r[,j]),0.975)
}

##
#Robinia
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r.r<-mvrnorm(1000,mu=ft.Topt.s.lin.R,Sigma=vcov(fit_Resp_norm_Topt.s.lin_ROPS))

#Respiration ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 6)
dist.r21.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r21.r<-rep(NA,length(Tsim))
high.r21.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r21.r[i,j]=norm.Topt.s.lin(1/Rd.25.R21,vmat.r.r[i,10],vmat.r.r[i,11],vmat.r.r[i,12],vmat.r.r[i,13],18.5,Tsim[j])
  }
  low.r21.r[j]<-quantile(na.omit(dist.r21.r[,j]),0.025)
  high.r21.r[j]<-quantile(na.omit(dist.r21.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 6)
dist.r26.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r26.r<-rep(NA,length(Tsim))
high.r26.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r26.r[i,j]=norm.Topt.s.lin(1/Rd.25.R26,vmat.r.r[i,10],vmat.r.r[i,11],vmat.r.r[i,12],vmat.r.r[i,13],23.5,Tsim[j])
  }
  low.r26.r[j]<-quantile(na.omit(dist.r26.r[,j]),0.025)
  high.r26.r[j]<-quantile(na.omit(dist.r26.r[,j]),0.975)
}

#Respiration ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 6)
dist.r31.r=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r31.r<-rep(NA,length(Tsim))
high.r31.r<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r31.r[i,j]=norm.Topt.s.lin(1/Rd.25.R31,vmat.r.r[i,10],vmat.r.r[i,11],vmat.r.r[i,12],vmat.r.r[i,13],28.5,Tsim[j])
  }
  low.r31.r[j]<-quantile(na.omit(dist.r31.r[,j]),0.025)
  high.r31.r[j]<-quantile(na.omit(dist.r31.r[,j]),0.975)
}

###############################################################################################################
#Supplementary Figure 5
###############################################################################################################

#PDF dimension is 7x13 inches  

#Plotting Settings
par(mfrow=c(2,1))
par(mar=c(0.3,0.8,0.8,0.5))
par(oma=c(5,5.5,1,0.5))
par(pty="s")

#S5a
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(10,40),ylim=c(0,3.3),cex.lab=1.5,cex.axis=1.2,xaxt="n")
axis(1,at=c(10,20,30,40),labels=F,cex.axis=1.2)
curve(norm.Topt.s.lin(1/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,x),from=10,to=40,lty=1,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,x),from=10,to=40,lty=1,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,x),from=10,to=40,lty=1,lwd=3,col="orangered3",add=T)
curve(norm.Topt.s.lin(1/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,x),from=10,to=40,lty=2,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,x),from=10,to=40,lty=2,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,x),from=10,to=40,lty=2,lwd=3,col="orangered3",add=T)
curve(norm.Topt.s.lin(1/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,x),from=10,to=40,lty=3,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,x),from=10,to=40,lty=3,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,x),from=10,to=40,lty=3,lwd=3,col="orangered3",add=T)
curve(norm.Topt.s.lin(1/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,x),from=10,to=40,lty=4,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,x),from=10,to=40,lty=4,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(1/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,x),from=10,to=40,lty=4,lwd=3,col="orangered3",add=T)
legend(8,3.5,c(expression(underline(bolditalic("Morella"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
               expression(underline(bolditalic("Alnus"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
               expression(underline(bolditalic("Gliricidia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
               expression(underline(bolditalic("Robinia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C')),
       col=c(NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3"),
       lty=c(NA,1,1,1,NA,2,2,2,NA,3,3,3,NA,4,4,4),bty="n",lwd=3,y.intersp = 0.6,cex=1,seg.len=2,x.intersp = 0.5)
mtext(expression('Leaf respiration'),side=2,line=4.5,cex=1.5)
mtext(expression('(normalized to 1 at 25 '*degree*'C)'),side=2,line=3,cex=1.5)
mtext(text="a",side=3,cex=1.4,adj=0)

#S5b
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(10,40),ylim=c(0,3),cex.lab=1.5,cex.axis=1.2,xaxt="n")
axis(1,at=c(10,20,30,40),labels=c("10","20","30","40"),cex.axis=1.2)
curve(norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,x),from=10,to=40,lty=1,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,x),from=10,to=40,lty=1,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,x),from=10,to=40,lty=1,lwd=3,col="orangered3",add=T)
curve(norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,x),from=10,to=40,lty=2,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,x),from=10,to=40,lty=2,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,x),from=10,to=40,lty=2,lwd=3,col="orangered3",add=T)
curve(norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,x),from=10,to=40,lty=3,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,x),from=10,to=40,lty=3,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,x),from=10,to=40,lty=3,lwd=3,col="orangered3",add=T)
curve(norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,x),from=10,to=40,lty=4,lwd=3,col="dodgerblue1",add=T)
curve(norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,x),from=10,to=40,lty=4,lwd=3,col="gold1",add=T)
curve(norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,x),from=10,to=40,lty=4,lwd=3,col="orangered3",add=T)
mtext(expression(italic('R')[L]),side=2,line=4.5,cex=1.5)
mtext(expression('('*mu*'mol CO'[2]*' m'^-2*' s'^-1*')'),side=2,line=3,cex=1.5)
mtext(text="b",side=3,cex=1.4,adj=0)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=F)

###############################################################################################################
#Supplementary Figure 6
###############################################################################################################

#PDF dimension is 9x7 inches  

#Plotting Settings
par(pty="s")
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

#S6a
plot(Rd.M21a$Temperature,Rd.M21a$Respiration/(ft.Topt.s.lin.M[2]*Rd.25.M21),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m21.r,rev(high.m21.r)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(Rd.M21a$Temperature,Rd.M21a$Respiration/(ft.Topt.s.lin.M[2]*Rd.25.M21),pch=1,cex=1.5,col="black")
points(Rd.M21b$Temperature,Rd.M21b$Respiration/(ft.Topt.s.lin.M[3]*Rd.25.M21),pch=2,cex=1.5,col="black")
points(Rd.M21c$Temperature,Rd.M21c$Respiration/(ft.Topt.s.lin.M[4]*Rd.25.M21),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,x),
      from=10,to=40,lty=1,lwd=2,col="dodgerblue1",add=T)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)

#S6b
plot(Rd.A21a$Temperature,Rd.A21a$Respiration/(ft.Topt.s.lin.A[2]*Rd.25.A21),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a21.r,rev(high.a21.r)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(Rd.A21a$Temperature,Rd.A21a$Respiration/(ft.Topt.s.lin.A[2]*Rd.25.A21),pch=1,cex=1.5,col="black")
points(Rd.A21b$Temperature,Rd.A21b$Respiration/(ft.Topt.s.lin.A[3]*Rd.25.A21),pch=2,cex=1.5,col="black")
points(Rd.A21c$Temperature,Rd.A21c$Respiration/(ft.Topt.s.lin.A[4]*Rd.25.A21),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,x),
      from=10,to=40,lty=1,lwd=2,col="dodgerblue1",add=T)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)

#S6c
plot(Rd.G21a$Temperature,Rd.G21a$Respiration/(ft.Topt.s.lin.G[2]*Rd.25.G21),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g21.r,rev(high.g21.r)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(Rd.G21a$Temperature,Rd.G21a$Respiration/(ft.Topt.s.lin.G[2]*Rd.25.G21),pch=1,cex=1.5,col="black")
points(Rd.G21b$Temperature,Rd.G21b$Respiration/(ft.Topt.s.lin.G[3]*Rd.25.G21),pch=2,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,x),
      from=10,to=40,lty=1,lwd=2,col="dodgerblue1",add=T)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)

#S6d
plot(Rd.R21a$Temperature,Rd.R21a$Respiration/(ft.Topt.s.lin.R[2]*Rd.25.R21),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r21.r,rev(high.r21.r)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(Rd.R21a$Temperature,Rd.R21a$Respiration/(ft.Topt.s.lin.R[2]*Rd.25.R21),pch=1,cex=1.5,col="black")
points(Rd.R21b$Temperature,Rd.R21b$Respiration/(ft.Topt.s.lin.R[3]*Rd.25.R21),pch=2,cex=1.5,col="black")
points(Rd.R21c$Temperature,Rd.R21c$Respiration/(ft.Topt.s.lin.R[4]*Rd.25.R21),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,x),
      from=10,to=40,lty=1,lwd=2,col="dodgerblue1",add=T)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)

#S6e
plot(Rd.M26a$Temperature,Rd.M26a$Respiration/(ft.Topt.s.lin.M[5]*Rd.25.M26),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m26.r,rev(high.m26.r)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(Rd.M26a$Temperature,Rd.M26a$Respiration/(ft.Topt.s.lin.M[5]*Rd.25.M26),pch=1,cex=1.5,col="black")
points(Rd.M26b$Temperature,Rd.M26b$Respiration/(ft.Topt.s.lin.M[6]*Rd.25.M26),pch=2,cex=1.5,col="black")
points(Rd.M26c$Temperature,Rd.M26c$Respiration/(ft.Topt.s.lin.M[7]*Rd.25.M26),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,x),
      from=10,to=40,lty=1,lwd=2,col="gold1",add=T)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S6f
plot(Rd.A26a$Temperature,Rd.A26a$Respiration/(ft.Topt.s.lin.A[5]*Rd.25.A26),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a26.r,rev(high.a26.r)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(Rd.A26a$Temperature,Rd.A26a$Respiration/(ft.Topt.s.lin.A[5]*Rd.25.A26),pch=1,cex=1.5,col="black")
points(Rd.A26b$Temperature,Rd.A26b$Respiration/(ft.Topt.s.lin.A[6]*Rd.25.A26),pch=2,cex=1.5,col="black")
points(Rd.A26c$Temperature,Rd.A26c$Respiration/(ft.Topt.s.lin.A[7]*Rd.25.A26),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,x),
      from=10,to=40,lty=1,lwd=2,col="gold1",add=T)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S6g
plot(Rd.G26a$Temperature,Rd.G26a$Respiration/(ft.Topt.s.lin.G[4]*Rd.25.G26),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g26.r,rev(high.g26.r)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(Rd.G26a$Temperature,Rd.G26a$Respiration/(ft.Topt.s.lin.G[4]*Rd.25.G26),pch=1,cex=1.5,col="black")
points(Rd.G26b$Temperature,Rd.G26b$Respiration/(ft.Topt.s.lin.G[5]*Rd.25.G26),pch=2,cex=1.5,col="black")
points(Rd.G26c$Temperature,Rd.G26c$Respiration/(ft.Topt.s.lin.G[6]*Rd.25.G26),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,x),
      from=10,to=40,lty=1,lwd=2,col="gold1",add=T)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S6h
plot(Rd.R26a$Temperature,Rd.R26a$Respiration/(ft.Topt.s.lin.R[5]*Rd.25.R26),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r26.r,rev(high.r26.r)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(Rd.R26a$Temperature,Rd.R26a$Respiration/(ft.Topt.s.lin.R[5]*Rd.25.R26),pch=1,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,x),
      from=10,to=40,lty=1,lwd=2,col="gold1",add=T)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

#S6i
plot(Rd.M31a$Temperature,Rd.M31a$Respiration/(ft.Topt.s.lin.M[8]*Rd.25.M31),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.lab=1.5,cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m31.r,rev(high.m31.r)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(Rd.M31a$Temperature,Rd.M31a$Respiration/(ft.Topt.s.lin.M[8]*Rd.25.M31),pch=1,cex=1.5,col="black")
points(Rd.M31b$Temperature,Rd.M31b$Respiration/(ft.Topt.s.lin.M[9]*Rd.25.M31),pch=2,cex=1.5,col="black")
points(Rd.M31c$Temperature,Rd.M31c$Respiration/(ft.Topt.s.lin.M[10]*Rd.25.M31),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,x),
      from=10,to=40,lty=1,lwd=2,col="orangered3",add=T)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S6j
plot(Rd.A31a$Temperature,Rd.A31a$Respiration/(ft.Topt.s.lin.A[8]*Rd.25.A31),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a31.r,rev(high.a31.r)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(Rd.A31a$Temperature,Rd.A31a$Respiration/(ft.Topt.s.lin.A[8]*Rd.25.A31),pch=1,cex=1.5,col="black")
points(Rd.A31b$Temperature,Rd.A31b$Respiration/(ft.Topt.s.lin.A[9]*Rd.25.A31),pch=2,cex=1.5,col="black")
points(Rd.A31c$Temperature,Rd.A31c$Respiration/(ft.Topt.s.lin.A[10]*Rd.25.A31),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,x),
      from=10,to=40,lty=1,lwd=2,col="orangered3",add=T)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S6k
plot(Rd.G31a$Temperature,Rd.G31a$Respiration/(ft.Topt.s.lin.G[7]*Rd.25.G31),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g31.r,rev(high.g31.r)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(Rd.G31a$Temperature,Rd.G31a$Respiration/(ft.Topt.s.lin.G[7]*Rd.25.G31),pch=1,cex=1.5,col="black")
points(Rd.G31b$Temperature,Rd.G31b$Respiration/(ft.Topt.s.lin.G[8]*Rd.25.G31),pch=2,cex=1.5,col="black")
points(Rd.G31c$Temperature,Rd.G31c$Respiration/(ft.Topt.s.lin.G[9]*Rd.25.G31),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,x),
      from=10,to=40,lty=1,lwd=2,col="orangered3",add=T)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S6l
plot(Rd.R31a$Temperature,Rd.R31a$Respiration/(ft.Topt.s.lin.R[6]*Rd.25.R31),pch=1,cex=1.5,col="black",ylim=c(0,4),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,1,2,3,4),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r31.r,rev(high.r31.r)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(Rd.R31a$Temperature,Rd.R31a$Respiration/(ft.Topt.s.lin.R[6]*Rd.25.R31),pch=1,cex=1.5,col="black")
points(Rd.R31b$Temperature,Rd.R31b$Respiration/(ft.Topt.s.lin.R[7]*Rd.25.R31),pch=2,cex=1.5,col="black")
points(Rd.R31c$Temperature,Rd.R31c$Respiration/(ft.Topt.s.lin.R[8]*Rd.25.R31),pch=3,cex=1.5,col="black")
points(Rd.R31d$Temperature,Rd.R31d$Respiration/(ft.Topt.s.lin.R[9]*Rd.25.R31),pch=3,cex=1.5,col="black")
curve(norm.Topt.s.lin(1/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,x),
      from=10,to=40,lty=1,lwd=2,col="orangered3",add=T)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression('Leaf respiration (normalized to 1 at 25 '*degree*'C)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)