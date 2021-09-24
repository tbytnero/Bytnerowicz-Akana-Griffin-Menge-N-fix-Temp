###############################################################################################################
###############################################################################################################
#This script tests for a temperature effect on the background (no acetylene) production of ethylene by
#the plant-bacterial symbioses
###############################################################################################################
###############################################################################################################

#Load Necessary Package
library(bbmle)

#Read in data
BG.dat<-read.csv("C2H4_Background_C2H2_Temp.csv")

###
#Test for effect of temperature
###

#Define function
lin.func <- function(a,b,temp){
  y <- a + b*temp
  y
}

###
#First with species-specific intercepts
###

#Separate by species
ALRU.dat<-BG.dat[which(BG.dat$Species=="ALRU"),]
GLSE.dat<-BG.dat[which(BG.dat$Species=="GLSE"),]
MOCE.dat<-BG.dat[which(BG.dat$Species=="MOCE"),]
ROPS.dat<-BG.dat[which(BG.dat$Species=="ROPS"),]

#Define Log-likelihood function
BG_normNLL <- function(sdBG,aMOCE,aALRU,aGLSE,aROPS,b,
                       BGdatMOCE,BGdatALRU,BGdatGLSE,BGdatROPS,
                       tempMOCE,tempALRU,tempGLSE,tempROPS){
  BGmeanMOCE <- lin.func(aMOCE,b,tempMOCE)
  BGmeanALRU <- lin.func(aALRU,b,tempALRU)
  BGmeanGLSE <- lin.func(aGLSE,b,tempGLSE)
  BGmeanROPS <- lin.func(aROPS,b,tempROPS)
  -(sum(dnorm(BGdatMOCE,mean=BGmeanMOCE,sd=exp(sdBG),log=TRUE),na.dat=TRUE) +
      sum(dnorm(BGdatALRU,mean=BGmeanALRU,sd=exp(sdBG),log=TRUE),na.dat=TRUE) +
      sum(dnorm(BGdatGLSE,mean=BGmeanGLSE,sd=exp(sdBG),log=TRUE),na.dat=TRUE) +
      sum(dnorm(BGdatROPS,mean=BGmeanROPS,sd=exp(sdBG),log=TRUE),na.dat=TRUE))
}

##
#Maximum likelihood fit
##

#BG~Temperature (linear-linear)
fit_BG_lin <- mle2(BG_normNLL,start=list(sdBG=-1,aMOCE=0,aALRU=0,aGLSE=0,aROPS=0,b=0.1),
                   data=list(tempMOCE=MOCE.dat[,5],tempALRU=ALRU.dat[,5],tempGLSE=GLSE.dat[,5],tempROPS=ROPS.dat[,5],
                             BGdatMOCE=MOCE.dat[,8],BGdatALRU=ALRU.dat[,8],BGdatGLSE=GLSE.dat[,8],BGdatROPS=ROPS.dat[,8]),
                   control=list(maxit=20000))  

summary(fit_BG_lin) #Summary
confint(fit_BG_lin) #95% CI of slope overlaps zero (not significant)

###
#Now with growing temperature-specific intercepts
###

#Separate by growing temperature
dat21<-BG.dat[which(BG.dat$Gro_Temp=="21"),]
dat26<-BG.dat[which(BG.dat$Gro_Temp=="26"),]
dat31<-BG.dat[which(BG.dat$Gro_Temp=="31"),]

#Define Log-likelihood function
BG_normNLL2 <- function(sdBG,a21,a26,a31,b,
                       BGdat21,BGdat26,BGdat31,
                       temp21,temp26,temp31){
  BGmean21 <- lin.func(a21,b,temp21)
  BGmean26 <- lin.func(a26,b,temp26)
  BGmean31 <- lin.func(a31,b,temp31)
  -(sum(dnorm(BGdat21,mean=BGmean21,sd=exp(sdBG),log=TRUE),na.dat=TRUE) +
      sum(dnorm(BGdat26,mean=BGmean26,sd=exp(sdBG),log=TRUE),na.dat=TRUE) +
      sum(dnorm(BGdat31,mean=BGmean31,sd=exp(sdBG),log=TRUE),na.dat=TRUE))
}

#BG~Temperature (linear-linear)
fit_BG_lin2 <- mle2(BG_normNLL2,start=list(sdBG=-1,a21=0,a26=0,a31=0,b=0.1),
                   data=list(temp21=dat21[,5],temp26=dat26[,5],temp31=dat31[,5],
                             BGdat21=dat21[,8],BGdat26=dat26[,8],BGdat31=dat31[,8]),
                   control=list(maxit=20000))  

summary(fit_BG_lin2) #Summary
confint(fit_BG_lin2) #95% CI of slope overlaps zero (not significant)