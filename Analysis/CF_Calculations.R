###############################################################################################################
###############################################################################################################
#This script calculates the Conversion Factor (CF) values in Supplementary Table 13
#and examines the possible effect of temperature on CF values
###############################################################################################################
###############################################################################################################

#Load Necessary Package
library(bbmle)

#Load CF data from "Conversion_Factor" folder
CF.dat<-read.csv("CF_Temp_all.csv")

#Separate by species
ALRU.dat<-CF.dat[which(CF.dat$Species=="ALRU"),]
GLSE.dat<-CF.dat[which(CF.dat$Species=="GLSE"),]
MOCE.dat<-CF.dat[which(CF.dat$Species=="MOCE"),]
ROPS.dat<-CF.dat[which(CF.dat$Species=="ROPS"),]

#Remove CF values above 7 and below 1 (4 +/- 3)
ALRU.rm<-ALRU.dat[which(ALRU.dat[,23]<7&ALRU.dat[,23]>1),]
GLSE.rm<-GLSE.dat[which(GLSE.dat[,23]<7&GLSE.dat[,23]>1),]
MOCE.rm<-MOCE.dat[which(MOCE.dat[,23]<7&MOCE.dat[,23]>1),]
ROPS.rm<-ROPS.dat[which(ROPS.dat[,23]<7&ROPS.dat[,23]>1),]

###
#Test for effect of temperature
###

#Define function
lin.func <- function(a,b,temp){
  y <- a + b*temp
  y
}

#Define Log-likelihood function
CF_normNLL <- function(sdCF,aMOCE,aALRU,aGLSE,aROPS,b,
                    CFdatMOCE,CFdatALRU,CFdatGLSE,CFdatROPS,
                    tempMOCE,tempALRU,tempGLSE,tempROPS){
  CFmeanMOCE <- lin.func(aMOCE,b,tempMOCE)
  CFmeanALRU <- lin.func(aALRU,b,tempALRU)
  CFmeanGLSE <- lin.func(aGLSE,b,tempGLSE)
  CFmeanROPS <- lin.func(aROPS,b,tempROPS)
  -(sum(dnorm(CFdatMOCE,mean=CFmeanMOCE,sd=exp(sdCF),log=TRUE),na.rm=TRUE) +
      sum(dnorm(CFdatALRU,mean=CFmeanALRU,sd=exp(sdCF),log=TRUE),na.rm=TRUE) +
      sum(dnorm(CFdatGLSE,mean=CFmeanGLSE,sd=exp(sdCF),log=TRUE),na.rm=TRUE) +
      sum(dnorm(CFdatROPS,mean=CFmeanROPS,sd=exp(sdCF),log=TRUE),na.rm=TRUE))
}

##
#Maximum likelihood fits
##

#CF~Temperature (linear-linear)
fit_CF_lin <- mle2(CF_normNLL,start=list(sdCF=-1,aMOCE=4.0,aALRU=2.6,aGLSE=3.1,aROPS=3.4,b=0.03),
                   data=list(tempMOCE=MOCE.rm[,5],tempALRU=ALRU.rm[,5],tempGLSE=GLSE.rm[,5],tempROPS=ROPS.rm[,5],
                             CFdatMOCE=MOCE.rm[,23],CFdatALRU=ALRU.rm[,23],CFdatGLSE=GLSE.rm[,23],CFdatROPS=ROPS.rm[,23]),
                   control=list(maxit=20000))  

summary(fit_CF_lin) #Summary
confint(fit_CF_lin) #95% CI of slope overlaps zero (not significant)

#Log(CF)~Temperature (log-linear)
fit_CF_loglin <- mle2(CF_normNLL,start=list(sdCF=-1,aMOCE=1.3,aALRU=0.9,aGLSE=1.0,aROPS=1.2,b=0.008),
                   data=list(tempMOCE=MOCE.rm[,5],tempALRU=ALRU.rm[,5],tempGLSE=GLSE.rm[,5],tempROPS=ROPS.rm[,5],
                             CFdatMOCE=log(MOCE.rm[,23]),CFdatALRU=log(ALRU.rm[,23]),CFdatGLSE=log(GLSE.rm[,23]),CFdatROPS=log(ROPS.rm[,23])),
                   control=list(maxit=20000))  

summary(fit_CF_loglin) #Summary
confint(fit_CF_loglin) #95% CI of slope overlaps zero (not significant)

###
#Calculate values for Supplementary Table 13
###

#Calculate species-level mean CF
mean(ALRU.rm$CF_size_cor) #3.15
mean(MOCE.rm$CF_size_cor) #4.98
mean(GLSE.rm$CF_size_cor) #3.84
mean(ROPS.rm$CF_size_cor) #4.27

#Calculate standard error of species-level CF
sd(ALRU.rm$CF_size_cor)/sqrt(length(ALRU.rm$CF_size_cor)) #0.40
sd(MOCE.rm$CF_size_cor)/sqrt(length(MOCE.rm$CF_size_cor)) #0.39
sd(GLSE.rm$CF_size_cor)/sqrt(length(GLSE.rm$CF_size_cor)) #0.47
sd(ROPS.rm$CF_size_cor)/sqrt(length(ROPS.rm$CF_size_cor)) #0.33

#Sample size
length(ALRU.rm$CF_size_cor) #14
length(MOCE.rm$CF_size_cor) #17
length(GLSE.rm$CF_size_cor) #13
length(ROPS.rm$CF_size_cor) #15