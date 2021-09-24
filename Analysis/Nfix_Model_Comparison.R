###############################################################################################################
###############################################################################################################
#This script compares models for N-fixation (Supplementary Tables 6 and 8) and
#then generates parameter estimates and 95% CIs for N-fixation (Supplementary Tables 1 and 3) 
###############################################################################################################
###############################################################################################################

####
#Load Necessary Packages
####

library(bbmle)
library(MASS)

####
#Read the data
####

#Alnus rubra
ALRU21a <- read.csv("ALRU21_070118_Vmax_temp.csv")
ALRU21b <- read.csv("ALRU21_072218_Vmax_temp.csv")
ALRU21c <- read.csv("ALRU21_073118_Vmax_temp.csv")
ALRU26a <- read.csv("ALRU26_102519_Vmax_temp.csv")
ALRU26b <- read.csv("ALRU26_111519_Vmax_temp.csv")
ALRU26c <- read.csv("ALRU26_122019_Vmax_temp.csv")
ALRU31a <- read.csv("ALRU31_021618_Vmax_temp.csv")
ALRU31b <- read.csv("ALRU31_070818_Vmax_temp.csv")
ALRU31c <- read.csv("ALRU31_073018_Vmax_temp.csv")

#Gliricidia sepium
GLSE21a <- read.csv("GLSE21_021518_Vmax_temp.csv")
GLSE21b <- read.csv("GLSE21_070218_Vmax_temp.csv")
GLSE21c <- read.csv("GLSE21_090219_Vmax_temp.csv")
GLSE26a <- read.csv("GLSE26_121419_Vmax_temp.csv")
GLSE26b <- read.csv("GLSE26_121819_Vmax_temp.csv")
GLSE26c <- read.csv("GLSE26_012120_Vmax_temp.csv")
GLSE31a <- read.csv("GLSE31_071118_Vmax_temp.csv")
GLSE31b <- read.csv("GLSE31_071618_Vmax_temp.csv")
GLSE31c <- read.csv("GLSE31_072418_Vmax_temp.csv")

#Morella cerifera
MOCE21a <- read.csv("MOCE21_022018_Vmax_temp.csv")
MOCE21b <- read.csv("MOCE21_072219_Vmax_temp.csv")
MOCE21c <- read.csv("MOCE21_100319_Vmax_temp.csv")
MOCE26a <- read.csv("MOCE26_121219_Vmax_temp.csv")
MOCE26b <- read.csv("MOCE26_121719_Vmax_temp.csv")
MOCE26c <- read.csv("MOCE26_020820_Vmax_temp.csv")
MOCE31a <- read.csv("MOCE31_022118_Vmax_temp.csv")
MOCE31b <- read.csv("MOCE31_030619_Vmax_temp.csv")
MOCE31c <- read.csv("MOCE31_031519_Vmax_temp.csv")

#Robinia pseudoacacia
ROPS21a <- read.csv("ROPS21_021318_Vmax_temp.csv")
ROPS21b <- read.csv("ROPS21_072318_Vmax_temp.csv")
ROPS21c <- read.csv("ROPS21_080218_Vmax_temp.csv")
ROPS26a <- read.csv("ROPS26_120319_Vmax_temp.csv")
ROPS26b <- read.csv("ROPS26_121619_Vmax_temp.csv")
ROPS26c <- read.csv("ROPS26_020620_Vmax_temp.csv")
ROPS31a <- read.csv("ROPS31_020918_Vmax_temp.csv")
ROPS31b <- read.csv("ROPS31_071418_Vmax_temp.csv")
ROPS31c <- read.csv("ROPS31_072618_Vmax_temp.csv")

###############################################################################################################
#First we compare equations 2-5 (Supplementary Table 6)
###############################################################################################################

####
#Define functions
####

#Modified beta (equation 5)
beta <- function(ymax,Tmin,Topt,Tmax,T){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Peaked Arrhenius function (equation 4)
##Here, delta S is reparameterized as Hd/Topt+0.008314*log(Ea/(Hd-Ea)),
##following Medlyn et al. 2002 (reference 57)
##this helped with estimating starting parameter values and made it easier to fit the model
peak.ar <- function(k25,Ea,Topt,Hd,Tk){
  y <- k25*(Tk/298.15)*exp((Ea*(Tk - 298.15))/(298.15*0.008314*Tk)) * 
    (1+exp((298.15*(Hd/Topt+0.008314*log(Ea/(Hd-Ea))) - Hd)/(298.15*0.008314))) / 
    (1+exp((Tk*(Hd/Topt+0.008314*log(Ea/(Hd-Ea)))-Hd)/(Tk*0.008314)))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Quadratic (equation 3)
quad<-function(ymax,b,Topt,T){
  y <- pmax(0,ymax-b*(T-Topt)^2)
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Normal (equation 2)
norm<-function(ymax,Topt,s,T){
  y<-ymax*exp(-(T-Topt)^2/(2*s^2))
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

#NLL function for the modified beta (equation 5)
Nase_beta_normNLL_all <- function(sdNase,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM26a,ymaxM26b,ymaxM26c,ymaxM31a,ymaxM31b,ymaxM31c,
                                  ymaxA21a,ymaxA21b,ymaxA21c,ymaxA26a,ymaxA26b,ymaxA26c,ymaxA31a,ymaxA31b,ymaxA31c,
                                  ymaxG21a,ymaxG21b,ymaxG21c,ymaxG26a,ymaxG26b,ymaxG26c,ymaxG31a,ymaxG31b,ymaxG31c,
                                  ymaxR21a,ymaxR21b,ymaxR21c,ymaxR26a,ymaxR26b,ymaxR26c,ymaxR31a,ymaxR31b,ymaxR31c,
                                  TmaxM21,TmaxM26,TmaxM31,
                                  TmaxA21,TmaxA26,TmaxA31,
                                  TmaxG21,TmaxG26,TmaxG31,
                                  TmaxR21,TmaxR26,TmaxR31,
                                  TminM21,TminM26,TminM31,
                                  TminA21,TminA26,TminA31,
                                  TminG21,TminG26,TminG31,
                                  TminR21,TminR26,TminR31,
                                  ToptM21,ToptM26,ToptM31,
                                  ToptA21,ToptA26,ToptA31,
                                  ToptG21,ToptG26,ToptG31,
                                  ToptR21,ToptR26,ToptR31,
                                  TempM21a,TempM21b,TempM21c,TempM26a,TempM26b,TempM26c,TempM31a,TempM31b,TempM31c,
                                  TempA21a,TempA21b,TempA21c,TempA26a,TempA26b,TempA26c,TempA31a,TempA31b,TempA31c,
                                  TempG21a,TempG21b,TempG21c,TempG26a,TempG26b,TempG26c,TempG31a,TempG31b,TempG31c,
                                  TempR21a,TempR21b,TempR21c,TempR26a,TempR26b,TempR26c,TempR31a,TempR31b,TempR31c,
                                  NasedatM21a,NasedatM21b,NasedatM21c,NasedatM26a,NasedatM26b,NasedatM26c,NasedatM31a,NasedatM31b,NasedatM31c,
                                  NasedatA21a,NasedatA21b,NasedatA21c,NasedatA26a,NasedatA26b,NasedatA26c,NasedatA31a,NasedatA31b,NasedatA31c,
                                  NasedatG21a,NasedatG21b,NasedatG21c,NasedatG26a,NasedatG26b,NasedatG26c,NasedatG31a,NasedatG31b,NasedatG31c,
                                  NasedatR21a,NasedatR21b,NasedatR21c,NasedatR26a,NasedatR26b,NasedatR26c,NasedatR31a,NasedatR31b,NasedatR31c){
  NasemeanM21a <- beta(ymaxM21a,TminM21,ToptM21,TmaxM21,TempM21a)
  NasemeanM21b <- beta(ymaxM21b,TminM21,ToptM21,TmaxM21,TempM21b)
  NasemeanM21c <- beta(ymaxM21c,TminM21,ToptM21,TmaxM21,TempM21c)
  NasemeanM26a <- beta(ymaxM26a,TminM26,ToptM26,TmaxM26,TempM26a)
  NasemeanM26b <- beta(ymaxM26b,TminM26,ToptM26,TmaxM26,TempM26b)
  NasemeanM26c <- beta(ymaxM26c,TminM26,ToptM26,TmaxM26,TempM26c)
  NasemeanM31a <- beta(ymaxM31a,TminM31,ToptM31,TmaxM31,TempM31a)
  NasemeanM31b <- beta(ymaxM31b,TminM31,ToptM31,TmaxM31,TempM31b)
  NasemeanM31c <- beta(ymaxM31c,TminM31,ToptM31,TmaxM31,TempM31c)
  NasemeanA21a <- beta(ymaxA21a,TminA21,ToptA21,TmaxA21,TempA21a)
  NasemeanA21b <- beta(ymaxA21b,TminA21,ToptA21,TmaxA21,TempA21b)
  NasemeanA21c <- beta(ymaxA21c,TminA21,ToptA21,TmaxA21,TempA21c)
  NasemeanA26a <- beta(ymaxA26a,TminA26,ToptA26,TmaxA26,TempA26a)
  NasemeanA26b <- beta(ymaxA26b,TminA26,ToptA26,TmaxA26,TempA26b)
  NasemeanA26c <- beta(ymaxA26c,TminA26,ToptA26,TmaxA26,TempA26c)
  NasemeanA31a <- beta(ymaxA31a,TminA31,ToptA31,TmaxA31,TempA31a)
  NasemeanA31b <- beta(ymaxA31b,TminA31,ToptA31,TmaxA31,TempA31b)
  NasemeanA31c <- beta(ymaxA31c,TminA31,ToptA31,TmaxA31,TempA31c)
  NasemeanG21a <- beta(ymaxG21a,TminG21,ToptG21,TmaxG21,TempG21a)
  NasemeanG21b <- beta(ymaxG21b,TminG21,ToptG21,TmaxG21,TempG21b)
  NasemeanG21c <- beta(ymaxG21c,TminG21,ToptG21,TmaxG21,TempG21c)
  NasemeanG26a <- beta(ymaxG26a,TminG26,ToptG26,TmaxG26,TempG26a)
  NasemeanG26b <- beta(ymaxG26b,TminG26,ToptG26,TmaxG26,TempG26b)
  NasemeanG26c <- beta(ymaxG26c,TminG26,ToptG26,TmaxG26,TempG26c)
  NasemeanG31a <- beta(ymaxG31a,TminG31,ToptG31,TmaxG31,TempG31a)
  NasemeanG31b <- beta(ymaxG31b,TminG31,ToptG31,TmaxG31,TempG31b)
  NasemeanG31c <- beta(ymaxG31c,TminG31,ToptG31,TmaxG31,TempG31c)
  NasemeanR21a <- beta(ymaxR21a,TminR21,ToptR21,TmaxR21,TempR21a)
  NasemeanR21b <- beta(ymaxR21b,TminR21,ToptR21,TmaxR21,TempR21b)
  NasemeanR21c <- beta(ymaxR21c,TminR21,ToptR21,TmaxR21,TempR21c)
  NasemeanR26a <- beta(ymaxR26a,TminR26,ToptR26,TmaxR26,TempR26a)
  NasemeanR26b <- beta(ymaxR26b,TminR26,ToptR26,TmaxR26,TempR26b)
  NasemeanR26c <- beta(ymaxR26c,TminR26,ToptR26,TmaxR26,TempR26c)
  NasemeanR31a <- beta(ymaxR31a,TminR31,ToptR31,TmaxR31,TempR31a)
  NasemeanR31b <- beta(ymaxR31b,TminR31,ToptR31,TmaxR31,TempR31b)
  NasemeanR31c <- beta(ymaxR31c,TminR31,ToptR31,TmaxR31,TempR31c)
  -(sum(dnorm(NasedatM21a,mean=NasemeanM21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21b,mean=NasemeanM21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21c,mean=NasemeanM21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM26a,mean=NasemeanM26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26b,mean=NasemeanM26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26c,mean=NasemeanM26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM31a,mean=NasemeanM31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31b,mean=NasemeanM31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31c,mean=NasemeanM31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA21a,mean=NasemeanA21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21b,mean=NasemeanA21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21c,mean=NasemeanA21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA26a,mean=NasemeanA26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26b,mean=NasemeanA26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26c,mean=NasemeanA26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA31a,mean=NasemeanA31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31b,mean=NasemeanA31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31c,mean=NasemeanA31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG21a,mean=NasemeanG21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21b,mean=NasemeanG21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21c,mean=NasemeanG21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG26a,mean=NasemeanG26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26b,mean=NasemeanG26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26c,mean=NasemeanG26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG31a,mean=NasemeanG31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31b,mean=NasemeanG31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31c,mean=NasemeanG31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR21a,mean=NasemeanR21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21b,mean=NasemeanR21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21c,mean=NasemeanR21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR26a,mean=NasemeanR26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26b,mean=NasemeanR26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26c,mean=NasemeanR26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR31a,mean=NasemeanR31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31b,mean=NasemeanR31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31c,mean=NasemeanR31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for the peaked Arrhenius (equation 4)
Nase_peak.ar_normNLL_all <- function(sdNase,k25M21a,k25M21b,k25M21c,k25M26a,k25M26b,k25M26c,k25M31a,k25M31b,k25M31c,
                                      k25A21a,k25A21b,k25A21c,k25A26a,k25A26b,k25A26c,k25A31a,k25A31b,k25A31c,
                                      k25G21a,k25G21b,k25G21c,k25G26a,k25G26b,k25G26c,k25G31a,k25G31b,k25G31c,
                                      k25R21a,k25R21b,k25R21c,k25R26a,k25R26b,k25R26c,k25R31a,k25R31b,k25R31c,
                                      HdM21,HdM26,HdM31,
                                      HdA21,HdA26,HdA31,
                                      HdG21,HdG26,HdG31,
                                      HdR21,HdR26,HdR31,
                                      EaM21,EaM26,EaM31,
                                      EaA21,EaA26,EaA31,
                                      EaG21,EaG26,EaG31,
                                      EaR21,EaR26,EaR31,
                                      ToptM21,ToptM26,ToptM31,
                                      ToptA21,ToptA26,ToptA31,
                                      ToptG21,ToptG26,ToptG31,
                                      ToptR21,ToptR26,ToptR31,
                                      TkM21a,TkM21b,TkM21c,TkM26a,TkM26b,TkM26c,TkM31a,TkM31b,TkM31c,
                                      TkA21a,TkA21b,TkA21c,TkA26a,TkA26b,TkA26c,TkA31a,TkA31b,TkA31c,
                                      TkG21a,TkG21b,TkG21c,TkG26a,TkG26b,TkG26c,TkG31a,TkG31b,TkG31c,
                                      TkR21a,TkR21b,TkR21c,TkR26a,TkR26b,TkR26c,TkR31a,TkR31b,TkR31c,
                                      NasedatM21a,NasedatM21b,NasedatM21c,NasedatM26a,NasedatM26b,NasedatM26c,NasedatM31a,NasedatM31b,NasedatM31c,
                                      NasedatA21a,NasedatA21b,NasedatA21c,NasedatA26a,NasedatA26b,NasedatA26c,NasedatA31a,NasedatA31b,NasedatA31c,
                                      NasedatG21a,NasedatG21b,NasedatG21c,NasedatG26a,NasedatG26b,NasedatG26c,NasedatG31a,NasedatG31b,NasedatG31c,
                                      NasedatR21a,NasedatR21b,NasedatR21c,NasedatR26a,NasedatR26b,NasedatR26c,NasedatR31a,NasedatR31b,NasedatR31c){
  NasemeanM21a <- peak.ar(k25M21a,EaM21,ToptM21,HdM21,TkM21a)
  NasemeanM21b <- peak.ar(k25M21b,EaM21,ToptM21,HdM21,TkM21b)
  NasemeanM21c <- peak.ar(k25M21c,EaM21,ToptM21,HdM21,TkM21c)
  NasemeanM26a <- peak.ar(k25M26a,EaM26,ToptM26,HdM26,TkM26a)
  NasemeanM26b <- peak.ar(k25M26b,EaM26,ToptM26,HdM26,TkM26b)
  NasemeanM26c <- peak.ar(k25M26c,EaM26,ToptM26,HdM26,TkM26c)
  NasemeanM31a <- peak.ar(k25M31a,EaM31,ToptM31,HdM31,TkM31a)
  NasemeanM31b <- peak.ar(k25M31b,EaM31,ToptM31,HdM31,TkM31b)
  NasemeanM31c <- peak.ar(k25M31c,EaM31,ToptM31,HdM31,TkM31c)
  NasemeanA21a <- peak.ar(k25A21a,EaA21,ToptA21,HdA21,TkA21a)
  NasemeanA21b <- peak.ar(k25A21b,EaA21,ToptA21,HdA21,TkA21b)
  NasemeanA21c <- peak.ar(k25A21c,EaA21,ToptA21,HdA21,TkA21c)
  NasemeanA26a <- peak.ar(k25A26a,EaA26,ToptA26,HdA26,TkA26a)
  NasemeanA26b <- peak.ar(k25A26b,EaA26,ToptA26,HdA26,TkA26b)
  NasemeanA26c <- peak.ar(k25A26c,EaA26,ToptA26,HdA26,TkA26c)
  NasemeanA31a <- peak.ar(k25A31a,EaA31,ToptA31,HdA31,TkA31a)
  NasemeanA31b <- peak.ar(k25A31b,EaA31,ToptA31,HdA31,TkA31b)
  NasemeanA31c <- peak.ar(k25A31c,EaA31,ToptA31,HdA31,TkA31c)
  NasemeanG21a <- peak.ar(k25G21a,EaG21,ToptG21,HdG21,TkG21a)
  NasemeanG21b <- peak.ar(k25G21b,EaG21,ToptG21,HdG21,TkG21b)
  NasemeanG21c <- peak.ar(k25G21c,EaG21,ToptG21,HdG21,TkG21c)
  NasemeanG26a <- peak.ar(k25G26a,EaG26,ToptG26,HdG26,TkG26a)
  NasemeanG26b <- peak.ar(k25G26b,EaG26,ToptG26,HdG26,TkG26b)
  NasemeanG26c <- peak.ar(k25G26c,EaG26,ToptG26,HdG26,TkG26c)
  NasemeanG31a <- peak.ar(k25G31a,EaG31,ToptG31,HdG31,TkG31a)
  NasemeanG31b <- peak.ar(k25G31b,EaG31,ToptG31,HdG31,TkG31b)
  NasemeanG31c <- peak.ar(k25G31c,EaG31,ToptG31,HdG31,TkG31c)
  NasemeanR21a <- peak.ar(k25R21a,EaR21,ToptR21,HdR21,TkR21a)
  NasemeanR21b <- peak.ar(k25R21b,EaR21,ToptR21,HdR21,TkR21b)
  NasemeanR21c <- peak.ar(k25R21c,EaR21,ToptR21,HdR21,TkR21c)
  NasemeanR26a <- peak.ar(k25R26a,EaR26,ToptR26,HdR26,TkR26a)
  NasemeanR26b <- peak.ar(k25R26b,EaR26,ToptR26,HdR26,TkR26b)
  NasemeanR26c <- peak.ar(k25R26c,EaR26,ToptR26,HdR26,TkR26c)
  NasemeanR31a <- peak.ar(k25R31a,EaR31,ToptR31,HdR31,TkR31a)
  NasemeanR31b <- peak.ar(k25R31b,EaR31,ToptR31,HdR31,TkR31b)
  NasemeanR31c <- peak.ar(k25R31c,EaR31,ToptR31,HdR31,TkR31c)
  -(sum(dnorm(NasedatM21a,mean=NasemeanM21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21b,mean=NasemeanM21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21c,mean=NasemeanM21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM26a,mean=NasemeanM26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26b,mean=NasemeanM26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26c,mean=NasemeanM26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM31a,mean=NasemeanM31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31b,mean=NasemeanM31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31c,mean=NasemeanM31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA21a,mean=NasemeanA21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21b,mean=NasemeanA21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21c,mean=NasemeanA21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA26a,mean=NasemeanA26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26b,mean=NasemeanA26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26c,mean=NasemeanA26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA31a,mean=NasemeanA31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31b,mean=NasemeanA31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31c,mean=NasemeanA31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG21a,mean=NasemeanG21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21b,mean=NasemeanG21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21c,mean=NasemeanG21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG26a,mean=NasemeanG26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26b,mean=NasemeanG26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26c,mean=NasemeanG26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG31a,mean=NasemeanG31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31b,mean=NasemeanG31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31c,mean=NasemeanG31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR21a,mean=NasemeanR21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21b,mean=NasemeanR21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21c,mean=NasemeanR21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR26a,mean=NasemeanR26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26b,mean=NasemeanR26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26c,mean=NasemeanR26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR31a,mean=NasemeanR31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31b,mean=NasemeanR31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31c,mean=NasemeanR31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for the quadratic (equation 3)
Nase_quad_normNLL_all <- function(sdNase,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM26a,ymaxM26b,ymaxM26c,ymaxM31a,ymaxM31b,ymaxM31c,
                                  ymaxA21a,ymaxA21b,ymaxA21c,ymaxA26a,ymaxA26b,ymaxA26c,ymaxA31a,ymaxA31b,ymaxA31c,
                                  ymaxG21a,ymaxG21b,ymaxG21c,ymaxG26a,ymaxG26b,ymaxG26c,ymaxG31a,ymaxG31b,ymaxG31c,
                                  ymaxR21a,ymaxR21b,ymaxR21c,ymaxR26a,ymaxR26b,ymaxR26c,ymaxR31a,ymaxR31b,ymaxR31c,
                                  bM21,bM26,bM31,
                                  bA21,bA26,bA31,
                                  bG21,bG26,bG31,
                                  bR21,bR26,bR31,
                                  ToptM21,ToptM26,ToptM31,
                                  ToptA21,ToptA26,ToptA31,
                                  ToptG21,ToptG26,ToptG31,
                                  ToptR21,ToptR26,ToptR31,
                                  TempM21a,TempM21b,TempM21c,TempM26a,TempM26b,TempM26c,TempM31a,TempM31b,TempM31c,
                                  TempA21a,TempA21b,TempA21c,TempA26a,TempA26b,TempA26c,TempA31a,TempA31b,TempA31c,
                                  TempG21a,TempG21b,TempG21c,TempG26a,TempG26b,TempG26c,TempG31a,TempG31b,TempG31c,
                                  TempR21a,TempR21b,TempR21c,TempR26a,TempR26b,TempR26c,TempR31a,TempR31b,TempR31c,
                                  NasedatM21a,NasedatM21b,NasedatM21c,NasedatM26a,NasedatM26b,NasedatM26c,NasedatM31a,NasedatM31b,NasedatM31c,
                                  NasedatA21a,NasedatA21b,NasedatA21c,NasedatA26a,NasedatA26b,NasedatA26c,NasedatA31a,NasedatA31b,NasedatA31c,
                                  NasedatG21a,NasedatG21b,NasedatG21c,NasedatG26a,NasedatG26b,NasedatG26c,NasedatG31a,NasedatG31b,NasedatG31c,
                                  NasedatR21a,NasedatR21b,NasedatR21c,NasedatR26a,NasedatR26b,NasedatR26c,NasedatR31a,NasedatR31b,NasedatR31c){
  NasemeanM21a <- quad(ymaxM21a,bM21,ToptM21,TempM21a)
  NasemeanM21b <- quad(ymaxM21b,bM21,ToptM21,TempM21b)
  NasemeanM21c <- quad(ymaxM21c,bM21,ToptM21,TempM21c)
  NasemeanM26a <- quad(ymaxM26a,bM26,ToptM26,TempM26a)
  NasemeanM26b <- quad(ymaxM26b,bM26,ToptM26,TempM26b)
  NasemeanM26c <- quad(ymaxM26c,bM26,ToptM26,TempM26c)
  NasemeanM31a <- quad(ymaxM31a,bM31,ToptM31,TempM31a)
  NasemeanM31b <- quad(ymaxM31b,bM31,ToptM31,TempM31b)
  NasemeanM31c <- quad(ymaxM31c,bM31,ToptM31,TempM31c)
  NasemeanA21a <- quad(ymaxA21a,bA21,ToptA21,TempA21a)
  NasemeanA21b <- quad(ymaxA21b,bA21,ToptA21,TempA21b)
  NasemeanA21c <- quad(ymaxA21c,bA21,ToptA21,TempA21c)
  NasemeanA26a <- quad(ymaxA26a,bA26,ToptA26,TempA26a)
  NasemeanA26b <- quad(ymaxA26b,bA26,ToptA26,TempA26b)
  NasemeanA26c <- quad(ymaxA26c,bA26,ToptA26,TempA26c)
  NasemeanA31a <- quad(ymaxA31a,bA31,ToptA31,TempA31a)
  NasemeanA31b <- quad(ymaxA31b,bA31,ToptA31,TempA31b)
  NasemeanA31c <- quad(ymaxA31c,bA31,ToptA31,TempA31c)
  NasemeanG21a <- quad(ymaxG21a,bG21,ToptG21,TempG21a)
  NasemeanG21b <- quad(ymaxG21b,bG21,ToptG21,TempG21b)
  NasemeanG21c <- quad(ymaxG21c,bG21,ToptG21,TempG21c)
  NasemeanG26a <- quad(ymaxG26a,bG26,ToptG26,TempG26a)
  NasemeanG26b <- quad(ymaxG26b,bG26,ToptG26,TempG26b)
  NasemeanG26c <- quad(ymaxG26c,bG26,ToptG26,TempG26c)
  NasemeanG31a <- quad(ymaxG31a,bG31,ToptG31,TempG31a)
  NasemeanG31b <- quad(ymaxG31b,bG31,ToptG31,TempG31b)
  NasemeanG31c <- quad(ymaxG31c,bG31,ToptG31,TempG31c)
  NasemeanR21a <- quad(ymaxR21a,bR21,ToptR21,TempR21a)
  NasemeanR21b <- quad(ymaxR21b,bR21,ToptR21,TempR21b)
  NasemeanR21c <- quad(ymaxR21c,bR21,ToptR21,TempR21c)
  NasemeanR26a <- quad(ymaxR26a,bR26,ToptR26,TempR26a)
  NasemeanR26b <- quad(ymaxR26b,bR26,ToptR26,TempR26b)
  NasemeanR26c <- quad(ymaxR26c,bR26,ToptR26,TempR26c)
  NasemeanR31a <- quad(ymaxR31a,bR31,ToptR31,TempR31a)
  NasemeanR31b <- quad(ymaxR31b,bR31,ToptR31,TempR31b)
  NasemeanR31c <- quad(ymaxR31c,bR31,ToptR31,TempR31c)
  -(sum(dnorm(NasedatM21a,mean=NasemeanM21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21b,mean=NasemeanM21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21c,mean=NasemeanM21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM26a,mean=NasemeanM26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26b,mean=NasemeanM26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26c,mean=NasemeanM26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM31a,mean=NasemeanM31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31b,mean=NasemeanM31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31c,mean=NasemeanM31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA21a,mean=NasemeanA21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21b,mean=NasemeanA21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21c,mean=NasemeanA21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA26a,mean=NasemeanA26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26b,mean=NasemeanA26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26c,mean=NasemeanA26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA31a,mean=NasemeanA31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31b,mean=NasemeanA31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31c,mean=NasemeanA31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG21a,mean=NasemeanG21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21b,mean=NasemeanG21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21c,mean=NasemeanG21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG26a,mean=NasemeanG26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26b,mean=NasemeanG26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26c,mean=NasemeanG26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG31a,mean=NasemeanG31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31b,mean=NasemeanG31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31c,mean=NasemeanG31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR21a,mean=NasemeanR21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21b,mean=NasemeanR21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21c,mean=NasemeanR21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR26a,mean=NasemeanR26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26b,mean=NasemeanR26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26c,mean=NasemeanR26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR31a,mean=NasemeanR31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31b,mean=NasemeanR31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31c,mean=NasemeanR31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for the normal (equation 2)
Nase_norm_normNLL_all <- function(sdNase,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM26a,ymaxM26b,ymaxM26c,ymaxM31a,ymaxM31b,ymaxM31c,
                                  ymaxA21a,ymaxA21b,ymaxA21c,ymaxA26a,ymaxA26b,ymaxA26c,ymaxA31a,ymaxA31b,ymaxA31c,
                                  ymaxG21a,ymaxG21b,ymaxG21c,ymaxG26a,ymaxG26b,ymaxG26c,ymaxG31a,ymaxG31b,ymaxG31c,
                                  ymaxR21a,ymaxR21b,ymaxR21c,ymaxR26a,ymaxR26b,ymaxR26c,ymaxR31a,ymaxR31b,ymaxR31c,
                                  ToptM21,ToptM26,ToptM31,
                                  ToptA21,ToptA26,ToptA31,
                                  ToptG21,ToptG26,ToptG31,
                                  ToptR21,ToptR26,ToptR31,
                                  sM21,sM26,sM31,
                                  sA21,sA26,sA31,
                                  sG21,sG26,sG31,
                                  sR21,sR26,sR31,
                                  TempM21a,TempM21b,TempM21c,TempM26a,TempM26b,TempM26c,TempM31a,TempM31b,TempM31c,
                                  TempA21a,TempA21b,TempA21c,TempA26a,TempA26b,TempA26c,TempA31a,TempA31b,TempA31c,
                                  TempG21a,TempG21b,TempG21c,TempG26a,TempG26b,TempG26c,TempG31a,TempG31b,TempG31c,
                                  TempR21a,TempR21b,TempR21c,TempR26a,TempR26b,TempR26c,TempR31a,TempR31b,TempR31c,
                                  NasedatM21a,NasedatM21b,NasedatM21c,NasedatM26a,NasedatM26b,NasedatM26c,NasedatM31a,NasedatM31b,NasedatM31c,
                                  NasedatA21a,NasedatA21b,NasedatA21c,NasedatA26a,NasedatA26b,NasedatA26c,NasedatA31a,NasedatA31b,NasedatA31c,
                                  NasedatG21a,NasedatG21b,NasedatG21c,NasedatG26a,NasedatG26b,NasedatG26c,NasedatG31a,NasedatG31b,NasedatG31c,
                                  NasedatR21a,NasedatR21b,NasedatR21c,NasedatR26a,NasedatR26b,NasedatR26c,NasedatR31a,NasedatR31b,NasedatR31c){
  NasemeanM21a <- norm(ymaxM21a,ToptM21,sM21,TempM21a)
  NasemeanM21b <- norm(ymaxM21b,ToptM21,sM21,TempM21b)
  NasemeanM21c <- norm(ymaxM21c,ToptM21,sM21,TempM21c)
  NasemeanM26a <- norm(ymaxM26a,ToptM26,sM26,TempM26a)
  NasemeanM26b <- norm(ymaxM26b,ToptM26,sM26,TempM26b)
  NasemeanM26c <- norm(ymaxM26c,ToptM26,sM26,TempM26c)
  NasemeanM31a <- norm(ymaxM31a,ToptM31,sM31,TempM31a)
  NasemeanM31b <- norm(ymaxM31b,ToptM31,sM31,TempM31b)
  NasemeanM31c <- norm(ymaxM31c,ToptM31,sM31,TempM31c)
  NasemeanA21a <- norm(ymaxA21a,ToptA21,sA21,TempA21a)
  NasemeanA21b <- norm(ymaxA21b,ToptA21,sA21,TempA21b)
  NasemeanA21c <- norm(ymaxA21c,ToptA21,sA21,TempA21c)
  NasemeanA26a <- norm(ymaxA26a,ToptA26,sA26,TempA26a)
  NasemeanA26b <- norm(ymaxA26b,ToptA26,sA26,TempA26b)
  NasemeanA26c <- norm(ymaxA26c,ToptA26,sA26,TempA26c)
  NasemeanA31a <- norm(ymaxA31a,ToptA31,sA31,TempA31a)
  NasemeanA31b <- norm(ymaxA31b,ToptA31,sA31,TempA31b)
  NasemeanA31c <- norm(ymaxA31c,ToptA31,sA31,TempA31c)
  NasemeanG21a <- norm(ymaxG21a,ToptG21,sG21,TempG21a)
  NasemeanG21b <- norm(ymaxG21b,ToptG21,sG21,TempG21b)
  NasemeanG21c <- norm(ymaxG21c,ToptG21,sG21,TempG21c)
  NasemeanG26a <- norm(ymaxG26a,ToptG26,sG26,TempG26a)
  NasemeanG26b <- norm(ymaxG26b,ToptG26,sG26,TempG26b)
  NasemeanG26c <- norm(ymaxG26c,ToptG26,sG26,TempG26c)
  NasemeanG31a <- norm(ymaxG31a,ToptG31,sG31,TempG31a)
  NasemeanG31b <- norm(ymaxG31b,ToptG31,sG31,TempG31b)
  NasemeanG31c <- norm(ymaxG31c,ToptG31,sG31,TempG31c)
  NasemeanR21a <- norm(ymaxR21a,ToptR21,sR21,TempR21a)
  NasemeanR21b <- norm(ymaxR21b,ToptR21,sR21,TempR21b)
  NasemeanR21c <- norm(ymaxR21c,ToptR21,sR21,TempR21c)
  NasemeanR26a <- norm(ymaxR26a,ToptR26,sR26,TempR26a)
  NasemeanR26b <- norm(ymaxR26b,ToptR26,sR26,TempR26b)
  NasemeanR26c <- norm(ymaxR26c,ToptR26,sR26,TempR26c)
  NasemeanR31a <- norm(ymaxR31a,ToptR31,sR31,TempR31a)
  NasemeanR31b <- norm(ymaxR31b,ToptR31,sR31,TempR31b)
  NasemeanR31c <- norm(ymaxR31c,ToptR31,sR31,TempR31c)
  -(sum(dnorm(NasedatM21a,mean=NasemeanM21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21b,mean=NasemeanM21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM21c,mean=NasemeanM21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM26a,mean=NasemeanM26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26b,mean=NasemeanM26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM26c,mean=NasemeanM26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatM31a,mean=NasemeanM31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31b,mean=NasemeanM31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatM31c,mean=NasemeanM31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA21a,mean=NasemeanA21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21b,mean=NasemeanA21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA21c,mean=NasemeanA21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA26a,mean=NasemeanA26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26b,mean=NasemeanA26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA26c,mean=NasemeanA26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatA31a,mean=NasemeanA31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31b,mean=NasemeanA31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatA31c,mean=NasemeanA31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG21a,mean=NasemeanG21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21b,mean=NasemeanG21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG21c,mean=NasemeanG21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG26a,mean=NasemeanG26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26b,mean=NasemeanG26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG26c,mean=NasemeanG26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatG31a,mean=NasemeanG31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31b,mean=NasemeanG31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatG31c,mean=NasemeanG31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR21a,mean=NasemeanR21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21b,mean=NasemeanR21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR21c,mean=NasemeanR21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR26a,mean=NasemeanR26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26b,mean=NasemeanR26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR26c,mean=NasemeanR26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(NasedatR31a,mean=NasemeanR31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31b,mean=NasemeanR31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(NasedatR31c,mean=NasemeanR31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

####
#Now use maximum likelihood to fit functions
####

#For the modified beta (equation 5)
fit_Nase_beta_all <- mle2(Nase_beta_normNLL_all,start=list(sdNase=-1,ymaxM21a=1,ymaxM21b=1,ymaxM21c=1,ymaxM26a=1,ymaxM26b=1,ymaxM26c=1,ymaxM31a=1,ymaxM31b=1,ymaxM31c=1,
                                                           ymaxA21a=1,ymaxA21b=1,ymaxA21c=1,ymaxA26a=1,ymaxA26b=1,ymaxA26c=1,ymaxA31a=1,ymaxA31b=1,ymaxA31c=1,
                                                           ymaxG21a=1,ymaxG21b=1,ymaxG21c=1,ymaxG26a=1,ymaxG26b=1,ymaxG26c=1,ymaxG31a=1,ymaxG31b=1,ymaxG31c=1,
                                                           ymaxR21a=1,ymaxR21b=1,ymaxR21c=1,ymaxR26a=1,ymaxR26b=1,ymaxR26c=1,ymaxR31a=1,ymaxR31b=1,ymaxR31c=1,
                                                          TmaxM21=44,TmaxM26=44,TmaxM31=44,
                                                          TmaxA21=44,TmaxA26=44,TmaxA31=44,
                                                          TmaxG21=44,TmaxG26=44,TmaxG31=44,
                                                          TmaxR21=44,TmaxR26=44,TmaxR31=44,
                                                          TminM21=-2,TminM26=-2,TminM31=-2,
                                                          TminA21=-2,TminA26=-2,TminA31=-2,
                                                          TminG21=-2,TminG26=-2,TminG31=-2,
                                                          TminR21=-2,TminR26=-2,TminR31=-2,
                                                          ToptM21=33,ToptM26=33,ToptM31=33,
                                                          ToptA21=33,ToptA26=33,ToptA31=33,
                                                          ToptG21=33,ToptG26=33,ToptG31=33,
                                                          ToptR21=33,ToptR26=33,ToptR31=33),
                                  data=list(TempM21a=MOCE21a$Temperature,TempM21b=MOCE21b$Temperature,TempM21c=MOCE21c$Temperature,
                                            TempM26a=MOCE26a$Temperature,TempM26b=MOCE26b$Temperature,TempM26c=MOCE26c$Temperature,
                                            TempM31a=MOCE31a$Temperature,TempM31b=MOCE31b$Temperature,TempM31c=MOCE31c$Temperature,
                                            TempA21a=ALRU21a$Temperature,TempA21b=ALRU21b$Temperature,TempA21c=ALRU21c$Temperature,
                                            TempA26a=ALRU26a$Temperature,TempA26b=ALRU26b$Temperature,TempA26c=ALRU26c$Temperature,
                                            TempA31a=ALRU31a$Temperature,TempA31b=ALRU31b$Temperature,TempA31c=ALRU31c$Temperature,
                                            TempG21a=GLSE21a$Temperature,TempG21b=GLSE21b$Temperature,TempG21c=GLSE21c$Temperature,
                                            TempG26a=GLSE26a$Temperature,TempG26b=GLSE26b$Temperature,TempG26c=GLSE26c$Temperature,
                                            TempG31a=GLSE31a$Temperature,TempG31b=GLSE31b$Temperature,TempG31c=GLSE31c$Temperature,
                                            TempR21a=ROPS21a$Temperature,TempR21b=ROPS21b$Temperature,TempR21c=ROPS21c$Temperature,
                                            TempR26a=ROPS26a$Temperature,TempR26b=ROPS26b$Temperature,TempR26c=ROPS26c$Temperature,
                                            TempR31a=ROPS31a$Temperature,TempR31b=ROPS31b$Temperature,TempR31c=ROPS31c$Temperature,
                                            NasedatM21a=MOCE21a$Vmax/9960.998949,NasedatM21b=MOCE21b$Vmax/3833.121977,NasedatM21c=MOCE21c$Vmax/3095.456985,
                                            NasedatM26a=MOCE26a$Vmax/6492.677546,NasedatM26b=MOCE26b$Vmax/5974.076428,NasedatM26c=MOCE26c$Vmax/7890.763063,
                                            NasedatM31a=MOCE31a$Vmax/8826.601393,NasedatM31b=MOCE31b$Vmax/2371.923752,NasedatM31c=MOCE31c$Vmax/629.7544222,
                                            NasedatA21a=ALRU21a$Vmax/10035.24751,NasedatA21b=ALRU21b$Vmax/7990.765838,NasedatA21c=ALRU21c$Vmax/4046.082544,
                                            NasedatA26a=ALRU26a$Vmax/2236.979133,NasedatA26b=ALRU26b$Vmax/2744.252462,NasedatA26c=ALRU26c$Vmax/1831.763172,
                                            NasedatA31a=ALRU31a$Vmax/406.3524389,NasedatA31b=ALRU31b$Vmax/1832.723413,NasedatA31c=ALRU31c$Vmax/1350.729381,
                                            NasedatG21a=GLSE21a$Vmax/308.3421323,NasedatG21b=GLSE21b$Vmax/2028.57371,NasedatG21c=GLSE21c$Vmax/879.0395611,
                                            NasedatG26a=GLSE26a$Vmax/909.4624356,NasedatG26b=GLSE26b$Vmax/1373.906634,NasedatG26c=GLSE26c$Vmax/705.9735103,
                                            NasedatG31a=GLSE31a$Vmax/1196.373627,NasedatG31b=GLSE31b$Vmax/312.4607561,NasedatG31c=GLSE31c$Vmax/209.4577018,
                                            NasedatR21a=ROPS21a$Vmax/1922.245351,NasedatR21b=ROPS21b$Vmax/1690.867422,NasedatR21c=ROPS21c$Vmax/450.8292229,
                                            NasedatR26a=ROPS26a$Vmax/2785.452785,NasedatR26b=ROPS26b$Vmax/1572.439251,NasedatR26c=ROPS26c$Vmax/1686.352,
                                            NasedatR31a=ROPS31a$Vmax/1646.769302,NasedatR31b=ROPS31b$Vmax/920.7416901,NasedatR31c=ROPS31c$Vmax/1978.187071),
                                  control=list(maxit=20000))
summary(fit_Nase_beta_all)

#For the peaked Arrhenius (equation 4)
fit_Nase_peak.ar_all <- mle2(Nase_peak.ar_normNLL_all,start=list(sdNase=-1,k25M21a=1,k25M21b=1,k25M21c=1,k25M26a=1,k25M26b=1,k25M26c=1,k25M31a=.4,k25M31b=.4,k25M31c=.2,
                                                                  k25A21a=1,k25A21b=1,k25A21c=1,k25A26a=1,k25A26b=1,k25A26c=1,k25A31a=1,k25A31b=1,k25A31c=1,
                                                                  k25G21a=1,k25G21b=1,k25G21c=1,k25G26a=1,k25G26b=1,k25G26c=1,k25G31a=1,k25G31b=1,k25G31c=1,
                                                                  k25R21a=1,k25R21b=1,k25R21c=1,k25R26a=1,k25R26b=1,k25R26c=1,k25R31a=1,k25R31b=1,k25R31c=1,
                                                                  HdM21=200,HdM26=200,HdM31=400,
                                                                  HdA21=200,HdA26=200,HdA31=200,
                                                                  HdG21=200,HdG26=200,HdG31=200,
                                                                  HdR21=200,HdR26=200,HdR31=200,
                                                                  EaM21=45,EaM26=45,EaM31=65,
                                                                  EaA21=45,EaA26=45,EaA31=45,
                                                                  EaG21=45,EaG26=45,EaG31=45,
                                                                  EaR21=45,EaR26=45,EaR31=45,
                                                                  ToptM21=303,ToptM26=306,ToptM31=310,
                                                                  ToptA21=305,ToptA26=307,ToptA31=305,
                                                                  ToptG21=303,ToptG26=306,ToptG31=308,
                                                                  ToptR21=305,ToptR26=305,ToptR31=305),
                              data=list(TkM21a=MOCE21a$Temperature+273.15,TkM21b=MOCE21b$Temperature+273.15,TkM21c=MOCE21c$Temperature+273.15,
                                        TkM26a=MOCE26a$Temperature+273.15,TkM26b=MOCE26b$Temperature+273.15,TkM26c=MOCE26c$Temperature+273.15,
                                        TkM31a=MOCE31a$Temperature+273.15,TkM31b=MOCE31b$Temperature+273.15,TkM31c=MOCE31c$Temperature+273.15,
                                        TkA21a=ALRU21a$Temperature+273.15,TkA21b=ALRU21b$Temperature+273.15,TkA21c=ALRU21c$Temperature+273.15,
                                        TkA26a=ALRU26a$Temperature+273.15,TkA26b=ALRU26b$Temperature+273.15,TkA26c=ALRU26c$Temperature+273.15,
                                        TkA31a=ALRU31a$Temperature+273.15,TkA31b=ALRU31b$Temperature+273.15,TkA31c=ALRU31c$Temperature+273.15,
                                        TkG21a=GLSE21a$Temperature+273.15,TkG21b=GLSE21b$Temperature+273.15,TkG21c=GLSE21c$Temperature+273.15,
                                        TkG26a=GLSE26a$Temperature+273.15,TkG26b=GLSE26b$Temperature+273.15,TkG26c=GLSE26c$Temperature+273.15,
                                        TkG31a=GLSE31a$Temperature+273.15,TkG31b=GLSE31b$Temperature+273.15,TkG31c=GLSE31c$Temperature+273.15,
                                        TkR21a=ROPS21a$Temperature+273.15,TkR21b=ROPS21b$Temperature+273.15,TkR21c=ROPS21c$Temperature+273.15,
                                        TkR26a=ROPS26a$Temperature+273.15,TkR26b=ROPS26b$Temperature+273.15,TkR26c=ROPS26c$Temperature+273.15,
                                        TkR31a=ROPS31a$Temperature+273.15,TkR31b=ROPS31b$Temperature+273.15,TkR31c=ROPS31c$Temperature+273.15,
                                        NasedatM21a=MOCE21a$Vmax/9960.998949,NasedatM21b=MOCE21b$Vmax/3833.121977,NasedatM21c=MOCE21c$Vmax/3095.456985,
                                        NasedatM26a=MOCE26a$Vmax/6492.677546,NasedatM26b=MOCE26b$Vmax/5974.076428,NasedatM26c=MOCE26c$Vmax/7890.763063,
                                        NasedatM31a=MOCE31a$Vmax/8826.601393,NasedatM31b=MOCE31b$Vmax/2371.7,NasedatM31c=MOCE31c$Vmax/627,
                                        NasedatA21a=ALRU21a$Vmax/10035.24751,NasedatA21b=ALRU21b$Vmax/7990.765838,NasedatA21c=ALRU21c$Vmax/4046.082544,
                                        NasedatA26a=ALRU26a$Vmax/2236.979133,NasedatA26b=ALRU26b$Vmax/2744.252462,NasedatA26c=ALRU26c$Vmax/1831.763172,
                                        NasedatA31a=ALRU31a$Vmax/406.3524389,NasedatA31b=ALRU31b$Vmax/1832.723413,NasedatA31c=ALRU31c$Vmax/1350.729381,
                                        NasedatG21a=GLSE21a$Vmax/308.3421323,NasedatG21b=GLSE21b$Vmax/2028.57371,NasedatG21c=GLSE21c$Vmax/879.0395611,
                                        NasedatG26a=GLSE26a$Vmax/909.4624356,NasedatG26b=GLSE26b$Vmax/1373.906634,NasedatG26c=GLSE26c$Vmax/705.9735103,
                                        NasedatG31a=GLSE31a$Vmax/1196.373627,NasedatG31b=GLSE31b$Vmax/312.4607561,NasedatG31c=GLSE31c$Vmax/209.4577018,
                                        NasedatR21a=ROPS21a$Vmax/1922.245351,NasedatR21b=ROPS21b$Vmax/1690.867422,NasedatR21c=ROPS21c$Vmax/450.8292229,
                                        NasedatR26a=ROPS26a$Vmax/2785.452785,NasedatR26b=ROPS26b$Vmax/1572.439251,NasedatR26c=ROPS26c$Vmax/1686.352,
                                        NasedatR31a=ROPS31a$Vmax/1646.769302,NasedatR31b=ROPS31b$Vmax/920.7416901,NasedatR31c=ROPS31c$Vmax/1978.187071),
                              control=list(maxit=20000))
summary(fit_Nase_peak.ar_all)

#For the quadratic (equation 3)
fit_Nase_quad_all <- mle2(Nase_quad_normNLL_all,start=list(sdNase=-1,ymaxM21a=1,ymaxM21b=1,ymaxM21c=1,ymaxM26a=1,ymaxM26b=1,ymaxM26c=1,ymaxM31a=1,ymaxM31b=1,ymaxM31c=1,
                                                           ymaxA21a=1,ymaxA21b=1,ymaxA21c=1,ymaxA26a=1,ymaxA26b=1,ymaxA26c=1,ymaxA31a=1,ymaxA31b=1,ymaxA31c=1,
                                                           ymaxG21a=1,ymaxG21b=1,ymaxG21c=1,ymaxG26a=1,ymaxG26b=1,ymaxG26c=1,ymaxG31a=1,ymaxG31b=1,ymaxG31c=1,
                                                           ymaxR21a=1,ymaxR21b=1,ymaxR21c=1,ymaxR26a=1,ymaxR26b=1,ymaxR26c=1,ymaxR31a=1,ymaxR31b=1,ymaxR31c=1,
                                                           bM21=0.01,bM26=0.01,bM31=0.01,
                                                           bA21=0.01,bA26=0.01,bA31=0.01,
                                                           bG21=0.01,bG26=0.01,bG31=0.01,
                                                           bR21=0.01,bR26=0.01,bR31=0.01,
                                                           ToptM21=29,ToptM26=34,ToptM31=37,
                                                           ToptA21=32,ToptA26=34,ToptA31=33,
                                                           ToptG21=30,ToptG26=34,ToptG31=36,
                                                           ToptR21=32,ToptR26=32,ToptR31=32),
                          data=list(TempM21a=MOCE21a$Temperature,TempM21b=MOCE21b$Temperature,TempM21c=MOCE21c$Temperature,
                                    TempM26a=MOCE26a$Temperature,TempM26b=MOCE26b$Temperature,TempM26c=MOCE26c$Temperature,
                                    TempM31a=MOCE31a$Temperature,TempM31b=MOCE31b$Temperature,TempM31c=MOCE31c$Temperature,
                                    TempA21a=ALRU21a$Temperature,TempA21b=ALRU21b$Temperature,TempA21c=ALRU21c$Temperature,
                                    TempA26a=ALRU26a$Temperature,TempA26b=ALRU26b$Temperature,TempA26c=ALRU26c$Temperature,
                                    TempA31a=ALRU31a$Temperature,TempA31b=ALRU31b$Temperature,TempA31c=ALRU31c$Temperature,
                                    TempG21a=GLSE21a$Temperature,TempG21b=GLSE21b$Temperature,TempG21c=GLSE21c$Temperature,
                                    TempG26a=GLSE26a$Temperature,TempG26b=GLSE26b$Temperature,TempG26c=GLSE26c$Temperature,
                                    TempG31a=GLSE31a$Temperature,TempG31b=GLSE31b$Temperature,TempG31c=GLSE31c$Temperature,
                                    TempR21a=ROPS21a$Temperature,TempR21b=ROPS21b$Temperature,TempR21c=ROPS21c$Temperature,
                                    TempR26a=ROPS26a$Temperature,TempR26b=ROPS26b$Temperature,TempR26c=ROPS26c$Temperature,
                                    TempR31a=ROPS31a$Temperature,TempR31b=ROPS31b$Temperature,TempR31c=ROPS31c$Temperature,
                                    NasedatM21a=MOCE21a$Vmax/9960.998949,NasedatM21b=MOCE21b$Vmax/3833.121977,NasedatM21c=MOCE21c$Vmax/3095.456985,
                                    NasedatM26a=MOCE26a$Vmax/6492.677546,NasedatM26b=MOCE26b$Vmax/5974.076428,NasedatM26c=MOCE26c$Vmax/7890.763063,
                                    NasedatM31a=MOCE31a$Vmax/8826.601393,NasedatM31b=MOCE31b$Vmax/2371.923752,NasedatM31c=MOCE31c$Vmax/629.7544222,
                                    NasedatA21a=ALRU21a$Vmax/10035.24751,NasedatA21b=ALRU21b$Vmax/7990.765838,NasedatA21c=ALRU21c$Vmax/4046.082544,
                                    NasedatA26a=ALRU26a$Vmax/2236.979133,NasedatA26b=ALRU26b$Vmax/2744.252462,NasedatA26c=ALRU26c$Vmax/1831.763172,
                                    NasedatA31a=ALRU31a$Vmax/406.3524389,NasedatA31b=ALRU31b$Vmax/1832.723413,NasedatA31c=ALRU31c$Vmax/1350.729381,
                                    NasedatG21a=GLSE21a$Vmax/308.3421323,NasedatG21b=GLSE21b$Vmax/2028.57371,NasedatG21c=GLSE21c$Vmax/879.0395611,
                                    NasedatG26a=GLSE26a$Vmax/909.4624356,NasedatG26b=GLSE26b$Vmax/1373.906634,NasedatG26c=GLSE26c$Vmax/705.9735103,
                                    NasedatG31a=GLSE31a$Vmax/1196.373627,NasedatG31b=GLSE31b$Vmax/312.4607561,NasedatG31c=GLSE31c$Vmax/209.4577018,
                                    NasedatR21a=ROPS21a$Vmax/1922.245351,NasedatR21b=ROPS21b$Vmax/1690.867422,NasedatR21c=ROPS21c$Vmax/450.8292229,
                                    NasedatR26a=ROPS26a$Vmax/2785.452785,NasedatR26b=ROPS26b$Vmax/1572.439251,NasedatR26c=ROPS26c$Vmax/1686.352,
                                    NasedatR31a=ROPS31a$Vmax/1646.769302,NasedatR31b=ROPS31b$Vmax/920.7416901,NasedatR31c=ROPS31c$Vmax/1978.187071),
                          control=list(maxit=20000))
summary(fit_Nase_quad_all)

#For the normal (equation 2)
fit_Nase_norm_all <- mle2(Nase_norm_normNLL_all,start=list(sdNase=-1,ymaxM21a=1,ymaxM21b=1,ymaxM21c=1,ymaxM26a=1,ymaxM26b=1,ymaxM26c=1,ymaxM31a=1,ymaxM31b=1,ymaxM31c=1,
                                                           ymaxA21a=1,ymaxA21b=1,ymaxA21c=1,ymaxA26a=1,ymaxA26b=1,ymaxA26c=1,ymaxA31a=1,ymaxA31b=1,ymaxA31c=1,
                                                           ymaxG21a=1,ymaxG21b=1,ymaxG21c=1,ymaxG26a=1,ymaxG26b=1,ymaxG26c=1,ymaxG31a=1,ymaxG31b=1,ymaxG31c=1,
                                                           ymaxR21a=1,ymaxR21b=1,ymaxR21c=1,ymaxR26a=1,ymaxR26b=1,ymaxR26c=1,ymaxR31a=1,ymaxR31b=1,ymaxR31c=1,
                                                           sM21=8,sM26=8,sM31=8,
                                                           sA21=8,sA26=8,sA31=8,
                                                           sG21=8,sG26=8,sG31=8,
                                                           sR21=8,sR26=8,sR31=8,
                                                           ToptM21=33,ToptM26=33,ToptM31=33,
                                                           ToptA21=33,ToptA26=33,ToptA31=33,
                                                           ToptG21=33,ToptG26=33,ToptG31=33,
                                                           ToptR21=33,ToptR26=33,ToptR31=33),
                          data=list(TempM21a=MOCE21a$Temperature,TempM21b=MOCE21b$Temperature,TempM21c=MOCE21c$Temperature,
                                    TempM26a=MOCE26a$Temperature,TempM26b=MOCE26b$Temperature,TempM26c=MOCE26c$Temperature,
                                    TempM31a=MOCE31a$Temperature,TempM31b=MOCE31b$Temperature,TempM31c=MOCE31c$Temperature,
                                    TempA21a=ALRU21a$Temperature,TempA21b=ALRU21b$Temperature,TempA21c=ALRU21c$Temperature,
                                    TempA26a=ALRU26a$Temperature,TempA26b=ALRU26b$Temperature,TempA26c=ALRU26c$Temperature,
                                    TempA31a=ALRU31a$Temperature,TempA31b=ALRU31b$Temperature,TempA31c=ALRU31c$Temperature,
                                    TempG21a=GLSE21a$Temperature,TempG21b=GLSE21b$Temperature,TempG21c=GLSE21c$Temperature,
                                    TempG26a=GLSE26a$Temperature,TempG26b=GLSE26b$Temperature,TempG26c=GLSE26c$Temperature,
                                    TempG31a=GLSE31a$Temperature,TempG31b=GLSE31b$Temperature,TempG31c=GLSE31c$Temperature,
                                    TempR21a=ROPS21a$Temperature,TempR21b=ROPS21b$Temperature,TempR21c=ROPS21c$Temperature,
                                    TempR26a=ROPS26a$Temperature,TempR26b=ROPS26b$Temperature,TempR26c=ROPS26c$Temperature,
                                    TempR31a=ROPS31a$Temperature,TempR31b=ROPS31b$Temperature,TempR31c=ROPS31c$Temperature,
                                    NasedatM21a=MOCE21a$Vmax/9960.998949,NasedatM21b=MOCE21b$Vmax/3833.121977,NasedatM21c=MOCE21c$Vmax/3095.456985,
                                    NasedatM26a=MOCE26a$Vmax/6492.677546,NasedatM26b=MOCE26b$Vmax/5974.076428,NasedatM26c=MOCE26c$Vmax/7890.763063,
                                    NasedatM31a=MOCE31a$Vmax/8826.601393,NasedatM31b=MOCE31b$Vmax/2371.923752,NasedatM31c=MOCE31c$Vmax/629.7544222,
                                    NasedatA21a=ALRU21a$Vmax/10035.24751,NasedatA21b=ALRU21b$Vmax/7990.765838,NasedatA21c=ALRU21c$Vmax/4046.082544,
                                    NasedatA26a=ALRU26a$Vmax/2236.979133,NasedatA26b=ALRU26b$Vmax/2744.252462,NasedatA26c=ALRU26c$Vmax/1831.763172,
                                    NasedatA31a=ALRU31a$Vmax/406.3524389,NasedatA31b=ALRU31b$Vmax/1832.723413,NasedatA31c=ALRU31c$Vmax/1350.729381,
                                    NasedatG21a=GLSE21a$Vmax/308.3421323,NasedatG21b=GLSE21b$Vmax/2028.57371,NasedatG21c=GLSE21c$Vmax/879.0395611,
                                    NasedatG26a=GLSE26a$Vmax/909.4624356,NasedatG26b=GLSE26b$Vmax/1373.906634,NasedatG26c=GLSE26c$Vmax/705.9735103,
                                    NasedatG31a=GLSE31a$Vmax/1196.373627,NasedatG31b=GLSE31b$Vmax/312.4607561,NasedatG31c=GLSE31c$Vmax/209.4577018,
                                    NasedatR21a=ROPS21a$Vmax/1922.245351,NasedatR21b=ROPS21b$Vmax/1690.867422,NasedatR21c=ROPS21c$Vmax/450.8292229,
                                    NasedatR26a=ROPS26a$Vmax/2785.452785,NasedatR26b=ROPS26b$Vmax/1572.439251,NasedatR26c=ROPS26c$Vmax/1686.352,
                                    NasedatR31a=ROPS31a$Vmax/1646.769302,NasedatR31b=ROPS31b$Vmax/920.7416901,NasedatR31c=ROPS31c$Vmax/1978.187071),
                          control=list(maxit=20000))
summary(fit_Nase_norm_all)

#Calculate delta AIC values
AICtab(fit_Nase_beta_all,fit_Nase_peak.ar_all,fit_Nase_quad_all,fit_Nase_norm_all)

#Calculate sample size for Supplementary Table 6
length(c(MOCE21a$Vmax,MOCE21b$Vmax,MOCE21c$Vmax,MOCE26a$Vmax,MOCE26b$Vmax,MOCE26c$Vmax,MOCE31a$Vmax,MOCE31b$Vmax,MOCE31c$Vmax))
length(c(ALRU21a$Vmax,ALRU21b$Vmax,ALRU21c$Vmax,ALRU26a$Vmax,ALRU26b$Vmax,ALRU26c$Vmax,ALRU31a$Vmax,ALRU31b$Vmax,ALRU31c$Vmax))
length(c(GLSE21a$Vmax,GLSE21b$Vmax,GLSE21c$Vmax,GLSE26a$Vmax,GLSE26b$Vmax,GLSE26c$Vmax,GLSE31a$Vmax,GLSE31b$Vmax,GLSE31c$Vmax))
length(c(ROPS21a$Vmax,ROPS21b$Vmax,ROPS21c$Vmax,ROPS26a$Vmax,ROPS26b$Vmax,ROPS26c$Vmax,ROPS31a$Vmax,ROPS31b$Vmax,ROPS31c$Vmax))

###############################################################################################################
#Model comparison for linear acclimation of modified beta function parameters (Supplementary Table 8)
###############################################################################################################

####
#Define functions
####

#Acclimation of all parameters
beta.lin.all <- function(ymax,a,b,c,d,e,f,T,Tgrow){
  y <- pmax(0,ymax*((e+f*Tgrow)-T)/((e+f*Tgrow)-(c+d*Tgrow))*(((T-(a+b*Tgrow))/((c+d*Tgrow)-(a+b*Tgrow)))^(((c+d*Tgrow)-(a+b*Tgrow))/((e+f*Tgrow)-(c+d*Tgrow)))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Tmin and Topt
beta.lin.Tmin.Topt <- function(ymax,Tmax,a,b,c,d,T,Tgrow){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-(c+d*Tgrow))*(((T-(a+b*Tgrow))/((c+d*Tgrow)-(a+b*Tgrow)))^(((c+d*Tgrow)-(a+b*Tgrow))/(Tmax-(c+d*Tgrow)))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Tmin and Tmax
beta.lin.Tmin.Tmax <- function(ymax,Topt,a,b,c,d,T,Tgrow){
  y <- pmax(0,ymax*((c+d*Tgrow)-T)/((c+d*Tgrow)-Topt)*(((T-(a+b*Tgrow))/(Topt-(a+b*Tgrow)))^((Topt-(a+b*Tgrow))/((c+d*Tgrow)-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Topt and Tmax
beta.lin.Topt.Tmax <- function(ymax,Tmin,a,b,c,d,T,Tgrow){
  y <- pmax(0,ymax*((c+d*Tgrow)-T)/((c+d*Tgrow)-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/((c+d*Tgrow)-(a+b*Tgrow)))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Tmin
beta.lin.Tmin <- function(ymax,Tmax,Topt,a,b,T,Tgrow){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-Topt)*(((T-(a+b*Tgrow))/(Topt-(a+b*Tgrow)))^((Topt-(a+b*Tgrow))/(Tmax-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Topt
beta.lin.Topt <- function(ymax,Tmax,Tmin,a,b,T,Tgrow){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/(Tmax-(a+b*Tgrow)))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Tmax
beta.lin.Tmax <- function(ymax,Tmin,Topt,a,b,T,Tgrow){
  y <- pmax(0,ymax*((a+b*Tgrow)-T)/((a+b*Tgrow)-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/((a+b*Tgrow)-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#No acclimation
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
#NLL functions
####

#NLL function for acclimation of all parameters
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

#NLL function for acclimation of Tmin and Topt
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

#NLL function for acclimation of Tmin and Tmax
Nase_beta_linTminTmax_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                          ymax26a,ymax26b,ymax26c,
                                          ymax31a,ymax31b,ymax31c,
                                          Topt,
                                          a,b,c,d,
                                          T21a,T21b,T21c,
                                          T26a,T26b,T26c,
                                          T31a,T31b,T31c,
                                          Nasedat21a,Nasedat21b,Nasedat21c,
                                          Nasedat26a,Nasedat26b,Nasedat26c,
                                          Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.Tmin.Tmax(ymax21a,Topt,a,b,c,d,T21a,18.5)
  Nasemean21b <- beta.lin.Tmin.Tmax(ymax21b,Topt,a,b,c,d,T21b,18.5)
  Nasemean21c <- beta.lin.Tmin.Tmax(ymax21c,Topt,a,b,c,d,T21c,18.5)
  Nasemean26a <- beta.lin.Tmin.Tmax(ymax26a,Topt,a,b,c,d,T26a,23.5)
  Nasemean26b <- beta.lin.Tmin.Tmax(ymax26b,Topt,a,b,c,d,T26b,23.5)
  Nasemean26c <- beta.lin.Tmin.Tmax(ymax26c,Topt,a,b,c,d,T26c,23.5)
  Nasemean31a <- beta.lin.Tmin.Tmax(ymax31a,Topt,a,b,c,d,T31a,28.5)
  Nasemean31b <- beta.lin.Tmin.Tmax(ymax31b,Topt,a,b,c,d,T31b,28.5)
  Nasemean31c <- beta.lin.Tmin.Tmax(ymax31c,Topt,a,b,c,d,T31c,28.5)
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

#NLL function for acclimation of Topt and Tmax
Nase_beta_linToptTmax_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                          ymax26a,ymax26b,ymax26c,
                                          ymax31a,ymax31b,ymax31c,
                                          Tmin,
                                          a,b,c,d,
                                          T21a,T21b,T21c,
                                          T26a,T26b,T26c,
                                          T31a,T31b,T31c,
                                          Nasedat21a,Nasedat21b,Nasedat21c,
                                          Nasedat26a,Nasedat26b,Nasedat26c,
                                          Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.Topt.Tmax(ymax21a,Tmin,a,b,c,d,T21a,18.5)
  Nasemean21b <- beta.lin.Topt.Tmax(ymax21b,Tmin,a,b,c,d,T21b,18.5)
  Nasemean21c <- beta.lin.Topt.Tmax(ymax21c,Tmin,a,b,c,d,T21c,18.5)
  Nasemean26a <- beta.lin.Topt.Tmax(ymax26a,Tmin,a,b,c,d,T26a,23.5)
  Nasemean26b <- beta.lin.Topt.Tmax(ymax26b,Tmin,a,b,c,d,T26b,23.5)
  Nasemean26c <- beta.lin.Topt.Tmax(ymax26c,Tmin,a,b,c,d,T26c,23.5)
  Nasemean31a <- beta.lin.Topt.Tmax(ymax31a,Tmin,a,b,c,d,T31a,28.5)
  Nasemean31b <- beta.lin.Topt.Tmax(ymax31b,Tmin,a,b,c,d,T31b,28.5)
  Nasemean31c <- beta.lin.Topt.Tmax(ymax31c,Tmin,a,b,c,d,T31c,28.5)
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

#NLL function for acclimation of Tmin
Nase_beta_linTmin_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                      ymax26a,ymax26b,ymax26c,
                                      ymax31a,ymax31b,ymax31c,
                                      Tmax,
                                      Topt,
                                      a,b,
                                      T21a,T21b,T21c,
                                      T26a,T26b,T26c,
                                      T31a,T31b,T31c,
                                      Nasedat21a,Nasedat21b,Nasedat21c,
                                      Nasedat26a,Nasedat26b,Nasedat26c,
                                      Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.Tmin(ymax21a,Tmax,Topt,a,b,T21a,18.5)
  Nasemean21b <- beta.lin.Tmin(ymax21b,Tmax,Topt,a,b,T21b,18.5)
  Nasemean21c <- beta.lin.Tmin(ymax21c,Tmax,Topt,a,b,T21c,18.5)
  Nasemean26a <- beta.lin.Tmin(ymax26a,Tmax,Topt,a,b,T26a,23.5)
  Nasemean26b <- beta.lin.Tmin(ymax26b,Tmax,Topt,a,b,T26b,23.5)
  Nasemean26c <- beta.lin.Tmin(ymax26c,Tmax,Topt,a,b,T26c,23.5)
  Nasemean31a <- beta.lin.Tmin(ymax31a,Tmax,Topt,a,b,T31a,28.5)
  Nasemean31b <- beta.lin.Tmin(ymax31b,Tmax,Topt,a,b,T31b,28.5)
  Nasemean31c <- beta.lin.Tmin(ymax31c,Tmax,Topt,a,b,T31c,28.5)
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

#NLL function for acclimation of Topt
Nase_beta_linTopt_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                      ymax26a,ymax26b,ymax26c,
                                      ymax31a,ymax31b,ymax31c,
                                      Tmax,
                                      Tmin,
                                      a,b,
                                      T21a,T21b,T21c,
                                      T26a,T26b,T26c,
                                      T31a,T31b,T31c,
                                      Nasedat21a,Nasedat21b,Nasedat21c,
                                      Nasedat26a,Nasedat26b,Nasedat26c,
                                      Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.Topt(ymax21a,Tmax,Tmin,a,b,T21a,18.5)
  Nasemean21b <- beta.lin.Topt(ymax21b,Tmax,Tmin,a,b,T21b,18.5)
  Nasemean21c <- beta.lin.Topt(ymax21c,Tmax,Tmin,a,b,T21c,18.5)
  Nasemean26a <- beta.lin.Topt(ymax26a,Tmax,Tmin,a,b,T26a,23.5)
  Nasemean26b <- beta.lin.Topt(ymax26b,Tmax,Tmin,a,b,T26b,23.5)
  Nasemean26c <- beta.lin.Topt(ymax26c,Tmax,Tmin,a,b,T26c,23.5)
  Nasemean31a <- beta.lin.Topt(ymax31a,Tmax,Tmin,a,b,T31a,28.5)
  Nasemean31b <- beta.lin.Topt(ymax31b,Tmax,Tmin,a,b,T31b,28.5)
  Nasemean31c <- beta.lin.Topt(ymax31c,Tmax,Tmin,a,b,T31c,28.5)
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

#NLL function for acclimation of Tmax
Nase_beta_linTmax_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                      ymax26a,ymax26b,ymax26c,
                                      ymax31a,ymax31b,ymax31c,
                                      Tmin,
                                      Topt,
                                      a,b,
                                      T21a,T21b,T21c,
                                      T26a,T26b,T26c,
                                      T31a,T31b,T31c,
                                      Nasedat21a,Nasedat21b,Nasedat21c,
                                      Nasedat26a,Nasedat26b,Nasedat26c,
                                      Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta.lin.Tmax(ymax21a,Tmin,Topt,a,b,T21a,18.5)
  Nasemean21b <- beta.lin.Tmax(ymax21b,Tmin,Topt,a,b,T21b,18.5)
  Nasemean21c <- beta.lin.Tmax(ymax21c,Tmin,Topt,a,b,T21c,18.5)
  Nasemean26a <- beta.lin.Tmax(ymax26a,Tmin,Topt,a,b,T26a,23.5)
  Nasemean26b <- beta.lin.Tmax(ymax26b,Tmin,Topt,a,b,T26b,23.5)
  Nasemean26c <- beta.lin.Tmax(ymax26c,Tmin,Topt,a,b,T26c,23.5)
  Nasemean31a <- beta.lin.Tmax(ymax31a,Tmin,Topt,a,b,T31a,28.5)
  Nasemean31b <- beta.lin.Tmax(ymax31b,Tmin,Topt,a,b,T31b,28.5)
  Nasemean31c <- beta.lin.Tmax(ymax31c,Tmin,Topt,a,b,T31c,28.5)
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

#NLL function for modified beta without acclimation
Nase_beta_bin_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,
                                  ymax26a,ymax26b,ymax26c,
                                  ymax31a,ymax31b,ymax31c,
                                  Tmin,
                                  Topt,
                                  Tmax,
                                  T21a,T21b,T21c,
                                  T26a,T26b,T26c,
                                  T31a,T31b,T31c,
                                  Nasedat21a,Nasedat21b,Nasedat21c,
                                  Nasedat26a,Nasedat26b,Nasedat26c,
                                  Nasedat31a,Nasedat31b,Nasedat31c){
  Nasemean21a <- beta(ymax21a,Tmin,Topt,Tmax,T21a)
  Nasemean21b <- beta(ymax21b,Tmin,Topt,Tmax,T21b)
  Nasemean21c <- beta(ymax21c,Tmin,Topt,Tmax,T21c)
  Nasemean26a <- beta(ymax26a,Tmin,Topt,Tmax,T26a)
  Nasemean26b <- beta(ymax26b,Tmin,Topt,Tmax,T26b)
  Nasemean26c <- beta(ymax26c,Tmin,Topt,Tmax,T26c)
  Nasemean31a <- beta(ymax31a,Tmin,Topt,Tmax,T31a)
  Nasemean31b <- beta(ymax31b,Tmin,Topt,Tmax,T31b)
  Nasemean31c <- beta(ymax31c,Tmin,Topt,Tmax,T31c)
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

###
#Morella
###

#Acclimation of all parameters
fit_Nase_beta_linall_MOCE <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                      ymax26a=1,ymax26b=1,ymax26c=1,
                                                                      ymax31a=1,ymax31b=1,ymax31c=1,
                                                                      a=-20,b=0,c=19,d=0.67,e=44,f=0.03),
                                  data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                            T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                            T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                            Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                            Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                            Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                  control=list(maxit=20000))
summary(fit_Nase_beta_linall_MOCE)

#Acclimation of Tmin and Topt
fit_Nase_beta_linTminTopt_MOCE <- mle2(Nase_beta_linTminTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmax=45,
                                                                                a=-20,b=0,c=19,d=0.67),
                                       data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                                 T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                                 T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                                 Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                                 Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                                 Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTopt_MOCE)

#Acclimation of Tmin and Tmax
fit_Nase_beta_linTminTmax_MOCE <- mle2(Nase_beta_linTminTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Topt=34,
                                                                                a=-20,b=0,c=44,d=0.03),
                                       data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                                 T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                                 T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                                 Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                                 Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                                 Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTmax_MOCE)

#Acclimation of Topt and Tmax
fit_Nase_beta_linToptTmax_MOCE <- mle2(Nase_beta_linToptTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmin=-2,
                                                                                a=19,b=0.67,c=44,d=0.03),
                                       data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                                 T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                                 T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                                 Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                                 Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                                 Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linToptTmax_MOCE)

#Acclimation of Tmin
fit_Nase_beta_linTmin_MOCE <- mle2(Nase_beta_linTmin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                         ymax26a=1,ymax26b=1,ymax26c=1,
                                                                         ymax31a=1,ymax31b=1,ymax31c=1,
                                                                         Topt=34,Tmax=45,
                                                                         a=-20,b=1),
                                    data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                              T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                              T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                              Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                              Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                              Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                    control=list(maxit=20000))
summary(fit_Nase_beta_linTmin_MOCE)

#Acclimation of Topt
fit_Nase_beta_linTopt_MOCE <- mle2(Nase_beta_linTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmax=45,
                                                                        Tmin=-2,
                                                                        a=19,b=0.67),
                                   data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                             T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                             T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                             Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                             Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                             Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTopt_MOCE)

#Acclimation of Tmax
fit_Nase_beta_linTmax_MOCE <- mle2(Nase_beta_linTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmin=0,Topt=34,
                                                                        a=44,b=0.03),
                                   data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                             T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                             T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                             Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                             Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                             Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmax_MOCE)

#No acclimation
fit_Nase_beta_bin_MOCE <- mle2(Nase_beta_bin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                Tmax=47,
                                                                Tmin=-2,
                                                                Topt=30),
                               data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                         T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                         T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                         Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                         Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                         Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                               control=list(maxit=20000))
summary(fit_Nase_beta_bin_MOCE)

#Delta AIC 
AICtab(fit_Nase_beta_linall_MOCE,fit_Nase_beta_linTminTopt_MOCE,fit_Nase_beta_linTminTmax_MOCE,
       fit_Nase_beta_linToptTmax_MOCE,fit_Nase_beta_linTmin_MOCE,fit_Nase_beta_linTopt_MOCE,
       fit_Nase_beta_linTmax_MOCE,fit_Nase_beta_bin_MOCE)

###
#Alnus
###

#Acclimation of all parameters
fit_Nase_beta_linall_ALRU <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-2.40096442,ymax21a=0.99710396,ymax21b=0.97485499,ymax21c=0.97033957,
                                                                       ymax26a=1.00553305,ymax26b=1.00558194,ymax26c=0.91760223,
                                                                       ymax31a=1.17671495,ymax31b=1.08051574,ymax31c=1.08032221,
                                                                       a=-20.60264188,b=0.51241949,c=31.20446094,d=0.06552585,e=41.84794232,f=0.02993945),
                                   data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                             T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                             T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                             Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                             Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                             Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linall_ALRU)

#Acclimation of Tmin and Topt
fit_Nase_beta_linTminTopt_ALRU <- mle2(Nase_beta_linTminTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmax=45,
                                                                                a=-25,b=0.8,c=31,d=0.06),
                                       data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                                 T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                                 T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                                 Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                                 Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                                 Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTopt_ALRU)

#Acclimation of Tmin and Tmax
fit_Nase_beta_linTminTmax_ALRU <- mle2(Nase_beta_linTminTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                 ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                 ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                 Topt=32.7,
                                                                                 a=-21,b=0.5,c=41.5,d=0.01),
                                        data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                                  T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                                  T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                                  Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                                  Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                                  Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                        control=list(maxit=20000))
summary(fit_Nase_beta_linTminTmax_ALRU)

#Acclimation of Topt and Tmax
fit_Nase_beta_linToptTmax_ALRU <- mle2(Nase_beta_linToptTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmin=-2,
                                                                                a=31,b=0.06,c=41,d=0.07),
                                       data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                                 T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                                 T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                                 Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                                 Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                                 Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linToptTmax_ALRU)

#Acclimation of Tmin
fit_Nase_beta_linTmin_ALRU <- mle2(Nase_beta_linTmin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Topt=32,Tmax=45,
                                                                        a=-25,b=0.8),
                                   data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                             T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                             T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                             Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                             Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                             Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmin_ALRU)

#Acclimation of Topt
fit_Nase_beta_linTopt_ALRU <- mle2(Nase_beta_linTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmin=-2,Tmax=45,
                                                                        a=31,b=0.06),
                                   data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                             T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                             T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                             Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                             Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                             Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTopt_ALRU)

#Acclimation of Tmax
fit_Nase_beta_linTmax_ALRU <- mle2(Nase_beta_linTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmin=-9,Topt=32.8,
                                                                        a=41,b=0),
                                   data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                             T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                             T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                             Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                             Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                             Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmax_ALRU)

#No acclimation
fit_Nase_beta_all_bin_ALRU <- mle2(Nase_beta_bin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                    ymax26a=1,ymax26b=1,ymax26c=1,
                                                                    ymax31a=1,ymax31b=1,ymax31c=1,
                                                                    Tmax=47,
                                                                    Tmin=-2,
                                                                    Topt=30),
                                   data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                             T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                             T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                             Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                             Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                             Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_all_bin_ALRU) 

#Delta AIC 
AICtab(fit_Nase_beta_linall_ALRU,fit_Nase_beta_linTminTopt_ALRU,fit_Nase_beta_linTminTmax_ALRU,
       fit_Nase_beta_linToptTmax_ALRU,fit_Nase_beta_linTmin_ALRU,fit_Nase_beta_linTopt_ALRU,
       fit_Nase_beta_linTmax_ALRU,fit_Nase_beta_all_bin_ALRU)

###
#Gliricidia
###

#Acclimation of all parameters
fit_Nase_beta_linall_GLSE <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                      ymax26a=1,ymax26b=1,ymax26c=1,
                                                                      ymax31a=1,ymax31b=1,ymax31c=1,
                                                                      a=5,b=0,c=19,d=0.6,e=45,f=0),
                                  data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                            T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                            T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                            Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                            Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                            Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                  control=list(maxit=20000))
summary(fit_Nase_beta_linall_GLSE)

#Acclimation of Tmin and Topt
fit_Nase_beta_linTminTopt_GLSE <- mle2(Nase_beta_linTminTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmax=45,
                                                                                a=5,b=0,c=19,d=0.6),
                                       data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                                 T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                                 T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                                 Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                                 Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                                 Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTopt_GLSE)

#Acclimation of Tmin and Tmax
fit_Nase_beta_linTminTmax_GLSE <- mle2(Nase_beta_linTminTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Topt=33,
                                                                                a=5,b=0,c=45,d=0),
                                       data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                                 T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                                 T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                                 Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                                 Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                                 Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTmax_GLSE)

#Acclimation of Topt and Tmax
fit_Nase_beta_linToptTmax_GLSE <- mle2(Nase_beta_linToptTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmin=5,
                                                                                a=19,b=0.6,c=45,d=0),
                                       data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                                 T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                                 T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                                 Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                                 Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                                 Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linToptTmax_GLSE)

#Acclimation of Tmin
fit_Nase_beta_linTmin_GLSE <- mle2(Nase_beta_linTmin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Topt=33,Tmax=45,
                                                                        a=5,b=0),
                                   data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                             T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                             T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                             Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                             Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                             Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmin_GLSE)

#Acclimation of Topt
fit_Nase_beta_linTopt_GLSE <- mle2(Nase_beta_linTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmax=45,
                                                                        Tmin=4,
                                                                        a=19,b=0.6),
                                   data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                             T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                             T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                             Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                             Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                             Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTopt_GLSE)

#Acclimation of Tmax
fit_Nase_beta_linTmax_GLSE <- mle2(Nase_beta_linTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmin=5,Topt=33,
                                                                        a=45,b=0),
                                   data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                             T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                             T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                             Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                             Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                             Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmax_GLSE)

#No acclimation
fit_Nase_beta_bin_GLSE <- mle2(Nase_beta_bin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                Tmax=47,
                                                                Tmin=-2,
                                                                Topt=30),
                               data=list(T21a=GLSE21a$Temperature,T21b=GLSE21b$Temperature,T21c=GLSE21c$Temperature,
                                         T26a=GLSE26a$Temperature,T26b=GLSE26b$Temperature,T26c=GLSE26c$Temperature,
                                         T31a=GLSE31a$Temperature,T31b=GLSE31b$Temperature,T31c=GLSE31c$Temperature,
                                         Nasedat21a=GLSE21a$Vmax/308.3421323,Nasedat21b=GLSE21b$Vmax/2028.57371,Nasedat21c=GLSE21c$Vmax/879.0395611,
                                         Nasedat26a=GLSE26a$Vmax/909.4624356,Nasedat26b=GLSE26b$Vmax/1373.906634,Nasedat26c=GLSE26c$Vmax/705.9735103,
                                         Nasedat31a=GLSE31a$Vmax/1196.373627,Nasedat31b=GLSE31b$Vmax/312.4607561,Nasedat31c=GLSE31c$Vmax/209.4577018),
                               control=list(maxit=20000))
summary(fit_Nase_beta_bin_GLSE)

#Delta AIC
AICtab(fit_Nase_beta_linall_GLSE,fit_Nase_beta_linTminTopt_GLSE,fit_Nase_beta_linTminTmax_GLSE,
       fit_Nase_beta_linToptTmax_GLSE,fit_Nase_beta_linTmin_GLSE,fit_Nase_beta_linTopt_GLSE,
       fit_Nase_beta_linTmax_GLSE,fit_Nase_beta_bin_GLSE)

###
#Robinia
###

#Acclimation of all parameters
fit_Nase_beta_linall_ROPS <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-2.10736306,ymax21a=0.95624295,ymax21b=0.81793514,ymax21c=1.08048702,
                                                                      ymax26a=0.89691111,ymax26b=0.96189080,ymax26c=0.88615555,
                                                                      ymax31a=0.88263265,ymax31b=0.70436753,ymax31c=0.88864834,
                                                                      a=-15.83976396,b=0.93637792,c=31.44431832,d=0.02429493,e=44.70312824,f=0.05161931),
                                  data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                            T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                            T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                            Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                            Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                            Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                  control=list(maxit=50000))
summary(fit_Nase_beta_linall_ROPS)

#Acclimation of Tmin and Topt
fit_Nase_beta_linTminTopt_ROPS <- mle2(Nase_beta_linTminTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmax=44,
                                                                                a=-16,b=1,c=31,d=0.01),
                                       data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                                 T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                                 T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                                 Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                                 Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                                 Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTopt_ROPS)

#Acclimation of Tmin and Tmax
fit_Nase_beta_linTminTmax_ROPS <- mle2(Nase_beta_linTminTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Topt=32,
                                                                                a=-16,b=1,c=43,d=0.1),
                                       data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                                 T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                                 T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                                 Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                                 Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                                 Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTmax_ROPS)

#Acclimation of Topt and Tmax
fit_Nase_beta_linToptTmax_ROPS <- mle2(Nase_beta_linToptTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                                Tmin=2,
                                                                                a=33,b=0,c=42,d=0.15),
                                       data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                                 T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                                 T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                                 Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                                 Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                                 Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linToptTmax_ROPS)

#Acclimation of Tmin
fit_Nase_beta_linTmin_ROPS <- mle2(Nase_beta_linTmin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Topt=32,Tmax=45,
                                                                        a=-61,b=2.6),
                                   data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                             T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                             T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                             Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                             Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                             Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmin_ROPS)

#Acclimation of Topt
fit_Nase_beta_linTopt_ROPS <- mle2(Nase_beta_linTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmin=2,Tmax=45,
                                                                        a=33,b=0),
                                   data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                             T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                             T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                             Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                             Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                             Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTopt_ROPS)

#Acclimation of Tmax
fit_Nase_beta_linTmax_ROPS <- mle2(Nase_beta_linTmax_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                        ymax26a=1,ymax26b=1,ymax26c=1,
                                                                        ymax31a=1,ymax31b=1,ymax31c=1,
                                                                        Tmin=2,Topt=32,
                                                                        a=42,b=0.15),
                                   data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                             T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                             T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                             Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                             Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                             Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                   control=list(maxit=20000))
summary(fit_Nase_beta_linTmax_ROPS)

#No acclimation
fit_Nase_beta_bin_ROPS <- mle2(Nase_beta_bin_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                ymax26a=1,ymax26b=1,ymax26c=1,
                                                                ymax31a=1,ymax31b=1,ymax31c=1,
                                                                Tmax=47,
                                                                Tmin=-2,
                                                                Topt=30),
                               data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                         T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                         T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                         Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                         Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                         Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                               control=list(maxit=20000))
summary(fit_Nase_beta_bin_ROPS)

#Delta AIC
AICtab(fit_Nase_beta_linall_ROPS,fit_Nase_beta_linTminTopt_ROPS,fit_Nase_beta_linTminTmax_ROPS,
       fit_Nase_beta_linToptTmax_ROPS,fit_Nase_beta_linTmin_ROPS,fit_Nase_beta_linTopt_ROPS,
       fit_Nase_beta_linTmax_ROPS,fit_Nase_beta_bin_ROPS)

###############################################################################################################
#Generate parameter estimates and calculate 95% CI for Supplementary Tables 1 and 3
###############################################################################################################

####
#Start with Supplementary Table 1
####

#Morella
coef(fit_Nase_beta_linall_MOCE)
fit_Nase_beta_linall_MOCE.pr<-profile(fit_Nase_beta_linall_MOCE)
confint(fit_Nase_beta_linall_MOCE.pr)
write.csv(coef(fit_Nase_beta_linall_MOCE),"MOCE_lin_beta_fit.csv") #Export parameter estimates

#Alnus
coef(fit_Nase_beta_linall_ALRU)
fit_Nase_beta_linall_ALRU.pr<-profile(fit_Nase_beta_linall_ALRU)
confint(fit_Nase_beta_linall_ALRU.pr)
write.csv(coef(fit_Nase_beta_linall_ALRU),"ALRU_lin_beta_fit.csv") #Export parameter estimates

#Gliricidia
coef(fit_Nase_beta_linTminTopt_GLSE)
fit_Nase_beta_linTminTopt_GLSE.pr<-profile(fit_Nase_beta_linTminTopt_GLSE)
confint(fit_Nase_beta_linTminTopt_GLSE.pr)
write.csv(coef(fit_Nase_beta_linTminTopt_GLSE),"GLSE_lin_beta_fit.csv") #Export parameter estimates

#Robinia
coef(fit_Nase_beta_linall_ROPS)
fit_Nase_beta_linall_ROPS.pr<-profile(fit_Nase_beta_linall_ROPS)
confint(fit_Nase_beta_linall_ROPS.pr)
write.csv(coef(fit_Nase_beta_linall_ROPS),"ROPS_lin_beta_fit.csv") #Export parameter estimates

###
#Next fit tropical and temperate temperature response functions
###

###
#NLL functions
###

#NLL function for acclimation of all parameters
Nase_beta_linall_2sp_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,ymax21d,ymax21e,ymax21f,
                                         ymax26a,ymax26b,ymax26c,ymax26d,ymax26e,ymax26f,
                                         ymax31a,ymax31b,ymax31c,ymax31d,ymax31e,ymax31f,
                                         a,b,c,d,e,f,
                                         T21a,T21b,T21c,T21d,T21e,T21f,
                                         T26a,T26b,T26c,T26d,T26e,T26f,
                                         T31a,T31b,T31c,T31d,T31e,T31f,
                                         Nasedat21a,Nasedat21b,Nasedat21c,Nasedat21d,Nasedat21e,Nasedat21f,
                                         Nasedat26a,Nasedat26b,Nasedat26c,Nasedat26d,Nasedat26e,Nasedat26f,
                                         Nasedat31a,Nasedat31b,Nasedat31c,Nasedat31d,Nasedat31e,Nasedat31f){
  Nasemean21a <- beta.lin.all(ymax21a,a,b,c,d,e,f,T21a,18.5)
  Nasemean21b <- beta.lin.all(ymax21b,a,b,c,d,e,f,T21b,18.5)
  Nasemean21c <- beta.lin.all(ymax21c,a,b,c,d,e,f,T21c,18.5)
  Nasemean21d <- beta.lin.all(ymax21d,a,b,c,d,e,f,T21d,18.5)
  Nasemean21e <- beta.lin.all(ymax21e,a,b,c,d,e,f,T21e,18.5)
  Nasemean21f <- beta.lin.all(ymax21f,a,b,c,d,e,f,T21f,18.5)
  Nasemean26a <- beta.lin.all(ymax26a,a,b,c,d,e,f,T26a,23.5)
  Nasemean26b <- beta.lin.all(ymax26b,a,b,c,d,e,f,T26b,23.5)
  Nasemean26c <- beta.lin.all(ymax26c,a,b,c,d,e,f,T26c,23.5)
  Nasemean26d <- beta.lin.all(ymax26d,a,b,c,d,e,f,T26d,23.5)
  Nasemean26e <- beta.lin.all(ymax26e,a,b,c,d,e,f,T26e,23.5)
  Nasemean26f <- beta.lin.all(ymax26f,a,b,c,d,e,f,T26f,23.5)
  Nasemean31a <- beta.lin.all(ymax31a,a,b,c,d,e,f,T31a,28.5)
  Nasemean31b <- beta.lin.all(ymax31b,a,b,c,d,e,f,T31b,28.5)
  Nasemean31c <- beta.lin.all(ymax31c,a,b,c,d,e,f,T31c,28.5)
  Nasemean31d <- beta.lin.all(ymax31d,a,b,c,d,e,f,T31d,28.5)
  Nasemean31e <- beta.lin.all(ymax31e,a,b,c,d,e,f,T31e,28.5)
  Nasemean31f <- beta.lin.all(ymax31f,a,b,c,d,e,f,T31f,28.5)
  -(sum(dnorm(Nasedat21a,mean=Nasemean21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21b,mean=Nasemean21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21c,mean=Nasemean21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21d,mean=Nasemean21d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21e,mean=Nasemean21e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21f,mean=Nasemean21f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26a,mean=Nasemean26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26b,mean=Nasemean26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26c,mean=Nasemean26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat26d,mean=Nasemean26d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26e,mean=Nasemean26e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26f,mean=Nasemean26f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31d,mean=Nasemean31d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31e,mean=Nasemean31e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31f,mean=Nasemean31f,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Tmin and Topt
Nase_beta_lin.Tmin.Topt_2sp_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,ymax21d,ymax21e,ymax21f,
                                                ymax26a,ymax26b,ymax26c,ymax26d,ymax26e,ymax26f,
                                                ymax31a,ymax31b,ymax31c,ymax31d,ymax31e,ymax31f,
                                                Tmax,a,b,c,d,
                                                T21a,T21b,T21c,T21d,T21e,T21f,
                                                T26a,T26b,T26c,T26d,T26e,T26f,
                                                T31a,T31b,T31c,T31d,T31e,T31f,
                                                Nasedat21a,Nasedat21b,Nasedat21c,Nasedat21d,Nasedat21e,Nasedat21f,
                                                Nasedat26a,Nasedat26b,Nasedat26c,Nasedat26d,Nasedat26e,Nasedat26f,
                                                Nasedat31a,Nasedat31b,Nasedat31c,Nasedat31d,Nasedat31e,Nasedat31f){
  Nasemean21a <- beta.lin.Tmin.Topt(ymax21a,Tmax,a,b,c,d,T21a,18.5)
  Nasemean21b <- beta.lin.Tmin.Topt(ymax21b,Tmax,a,b,c,d,T21b,18.5)
  Nasemean21c <- beta.lin.Tmin.Topt(ymax21c,Tmax,a,b,c,d,T21c,18.5)
  Nasemean21d <- beta.lin.Tmin.Topt(ymax21d,Tmax,a,b,c,d,T21d,18.5)
  Nasemean21e <- beta.lin.Tmin.Topt(ymax21e,Tmax,a,b,c,d,T21e,18.5)
  Nasemean21f <- beta.lin.Tmin.Topt(ymax21f,Tmax,a,b,c,d,T21f,18.5)
  Nasemean26a <- beta.lin.Tmin.Topt(ymax26a,Tmax,a,b,c,d,T26a,23.5)
  Nasemean26b <- beta.lin.Tmin.Topt(ymax26b,Tmax,a,b,c,d,T26b,23.5)
  Nasemean26c <- beta.lin.Tmin.Topt(ymax26c,Tmax,a,b,c,d,T26c,23.5)
  Nasemean26d <- beta.lin.Tmin.Topt(ymax26d,Tmax,a,b,c,d,T26d,23.5)
  Nasemean26e <- beta.lin.Tmin.Topt(ymax26e,Tmax,a,b,c,d,T26e,23.5)
  Nasemean26f <- beta.lin.Tmin.Topt(ymax26f,Tmax,a,b,c,d,T26f,23.5)
  Nasemean31a <- beta.lin.Tmin.Topt(ymax31a,Tmax,a,b,c,d,T31a,28.5)
  Nasemean31b <- beta.lin.Tmin.Topt(ymax31b,Tmax,a,b,c,d,T31b,28.5)
  Nasemean31c <- beta.lin.Tmin.Topt(ymax31c,Tmax,a,b,c,d,T31c,28.5)
  Nasemean31d <- beta.lin.Tmin.Topt(ymax31d,Tmax,a,b,c,d,T31d,28.5)
  Nasemean31e <- beta.lin.Tmin.Topt(ymax31e,Tmax,a,b,c,d,T31e,28.5)
  Nasemean31f <- beta.lin.Tmin.Topt(ymax31f,Tmax,a,b,c,d,T31f,28.5)
  -(sum(dnorm(Nasedat21a,mean=Nasemean21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21b,mean=Nasemean21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21c,mean=Nasemean21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21d,mean=Nasemean21d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21e,mean=Nasemean21e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21f,mean=Nasemean21f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26a,mean=Nasemean26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26b,mean=Nasemean26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26c,mean=Nasemean26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat26d,mean=Nasemean26d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26e,mean=Nasemean26e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26f,mean=Nasemean26f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31d,mean=Nasemean31d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31e,mean=Nasemean31e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31f,mean=Nasemean31f,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

###
#Maximum likelihood fits
###

#Tropical species (Morella and Gliricidia)
fit_Nase_beta_linTminTopt_Trop <- mle2(Nase_beta_lin.Tmin.Topt_2sp_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,ymax21d=1,ymax21e=1,ymax21f=1,
                                                                                      ymax26a=1,ymax26b=1,ymax26c=1,ymax26d=1,ymax26e=1,ymax26f=1,
                                                                                      ymax31a=1,ymax31b=1,ymax31c=1,ymax31d=1,ymax31e=1,ymax31f=1,
                                                                                      Tmax=45.35,a=-14.58,b=0.94,c=18.39,d=0.63),
                                       data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                                 T21d=GLSE21a$Temperature,T21e=GLSE21b$Temperature,T21f=GLSE21c$Temperature,
                                                 T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                                 T26d=GLSE26a$Temperature,T26e=GLSE26b$Temperature,T26f=GLSE26c$Temperature,
                                                 T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                                 T31d=GLSE31a$Temperature,T31e=GLSE31b$Temperature,T31f=GLSE31c$Temperature,
                                                 Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                                 Nasedat21d=GLSE21a$Vmax/308.3421323,Nasedat21e=GLSE21b$Vmax/2028.57371,Nasedat21f=GLSE21c$Vmax/879.0395611,
                                                 Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                                 Nasedat26d=GLSE26a$Vmax/909.4624356,Nasedat26e=GLSE26b$Vmax/1373.906634,Nasedat26f=GLSE26c$Vmax/705.9735103,
                                                 Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222,
                                                 Nasedat31d=GLSE31a$Vmax/1196.373627,Nasedat31e=GLSE31b$Vmax/312.4607561,Nasedat31f=GLSE31c$Vmax/209.4577018),
                                       control=list(maxit=20000))
summary(fit_Nase_beta_linTminTopt_Trop)

#Temperate species (Alnus and Robinia)
fit_Nase_beta_linall_Temp <- mle2(Nase_beta_linall_2sp_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,ymax21d=1,ymax21e=1,ymax21f=1,
                                                                          ymax26a=1,ymax26b=1,ymax26c=1,ymax26d=1,ymax26e=1,ymax26f=1,
                                                                          ymax31a=1,ymax31b=1,ymax31c=1,ymax31d=1,ymax31e=1,ymax31f=1,
                                                                          a=-17.93,b=0.71,c=31.34,d=0.04,e=43.33,f=0.04),
                                  data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                            T21d=ROPS21a$Temperature,T21e=ROPS21b$Temperature,T21f=ROPS21c$Temperature,
                                            T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                            T26d=ROPS26a$Temperature,T26e=ROPS26b$Temperature,T26f=ROPS26c$Temperature,
                                            T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                            T31d=ROPS31a$Temperature,T31e=ROPS31b$Temperature,T31f=ROPS31c$Temperature,
                                            Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                            Nasedat21d=ROPS21a$Vmax/1922.245351,Nasedat21e=ROPS21b$Vmax/1690.867422,Nasedat21f=ROPS21c$Vmax/450.8292229,
                                            Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                            Nasedat26d=ROPS26a$Vmax/2785.452785,Nasedat26e=ROPS26b$Vmax/1572.439251,Nasedat26f=ROPS26c$Vmax/1686.352,
                                            Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381,
                                            Nasedat31d=ROPS31a$Vmax/1646.769302,Nasedat31e=ROPS31b$Vmax/920.7416901,Nasedat31f=ROPS31c$Vmax/1978.187071),
                                  control=list(maxit=20000))
summary(fit_Nase_beta_linall_Temp)

###
#Parameter estimates and 95% CI for tropical and temperate species
###

#Tropical
coef(fit_Nase_beta_linTminTopt_Trop)
fit_Nase_beta_linTminTopt_Trop.pr<-profile(fit_Nase_beta_linTminTopt_Trop)
confint(fit_Nase_beta_linTminTopt_Trop.pr)

#Temperate
fit_Nase_beta_linall_Temp.pr<-profile(fit_Nase_beta_linall_Temp)
confint(fit_Nase_beta_linall_Temp.pr)
coef(fit_Nase_beta_linall_Temp)

####
#Now calculate parameters and 95% CI for Supplementary Table 3
####

###
#Morella
###

##
#Parameter estimates
##

#Tmin
coef(fit_Nase_beta_linall_MOCE)[11]+coef(fit_Nase_beta_linall_MOCE)[12]*23.5

#Topt
coef(fit_Nase_beta_linall_MOCE)[13]+coef(fit_Nase_beta_linall_MOCE)[14]*23.5

#Tmax
coef(fit_Nase_beta_linall_MOCE)[15]+coef(fit_Nase_beta_linall_MOCE)[16]*23.5

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_MOCE),Sigma=vcov(fit_Nase_beta_linall_MOCE))

#Tmin
dist.mTmin26<-rep(NA,1000)
for(i in 1:1000){
  dist.mTmin26[i]=vmat.m[i,11]+vmat.m[i,12]*23.5
}
quantile(dist.mTmin26,0.025)
quantile(dist.mTmin26,0.975)

#Topt
dist.mTopt26<-rep(NA,1000)
for(i in 1:1000){
  dist.mTopt26[i]=vmat.m[i,13]+vmat.m[i,14]*23.5
}
quantile(dist.mTopt26,0.025)
quantile(dist.mTopt26,0.975)

#Tmax
dist.mTmax26<-rep(NA,1000)
for(i in 1:1000){
  dist.mTmax26[i]=vmat.m[i,15]+vmat.m[i,16]*23.5
}
quantile(dist.mTmax26,0.025)
quantile(dist.mTmax26,0.975)

###
#Alnus
###

##
#Parameter estimates
##

#Tmin
coef(fit_Nase_beta_linall_ALRU)[11]+coef(fit_Nase_beta_linall_ALRU)[12]*18.5

#Topt
coef(fit_Nase_beta_linall_ALRU)[13]+coef(fit_Nase_beta_linall_ALRU)[14]*18.5

#Tmax
coef(fit_Nase_beta_linall_ALRU)[15]+coef(fit_Nase_beta_linall_ALRU)[16]*18.5

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_ALRU),Sigma=vcov(fit_Nase_beta_linall_ALRU))

#Tmin
dist.aTmin21<-rep(NA,1000)
for(i in 1:1000){
  dist.aTmin21[i]=vmat.a[i,11]+vmat.a[i,12]*18.5
}
quantile(dist.aTmin21,0.025)
quantile(dist.aTmin21,0.975)

#Topt
dist.aTopt21<-rep(NA,1000)
for(i in 1:1000){
  dist.aTopt21[i]=vmat.a[i,13]+vmat.a[i,14]*18.5
}
quantile(dist.aTopt21,0.025)
quantile(dist.aTopt21,0.975)

#Tmax
dist.aTmax21<-rep(NA,1000)
for(i in 1:1000){
  dist.aTmax21[i]=vmat.a[i,15]+vmat.a[i,16]*18.5
}
quantile(dist.aTmax21,0.025)
quantile(dist.aTmax21,0.975)

###
#Gliricidia
###

##
#Parameter estimates
##

#Tmin
coef(fit_Nase_beta_linTminTopt_GLSE)[12]+coef(fit_Nase_beta_linTminTopt_GLSE)[13]*23.5

#Topt
coef(fit_Nase_beta_linTminTopt_GLSE)[14]+coef(fit_Nase_beta_linTminTopt_GLSE)[15]*23.5

#Tmax estimate is the same as in Supplementary Table 1

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g=mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_GLSE),Sigma=vcov(fit_Nase_beta_linTminTopt_GLSE))

#Tmin
dist.gTmin26<-rep(NA,1000)
for(i in 1:1000){
  dist.gTmin26[i]=vmat.g[i,12]+vmat.g[i,13]*23.5
}
quantile(dist.gTmin26,0.025)
quantile(dist.gTmin26,0.975)

#Topt
dist.gTopt26<-rep(NA,1000)
for(i in 1:1000){
  dist.gTopt26[i]=vmat.g[i,14]+vmat.g[i,15]*23.5
}
quantile(dist.gTopt26,0.025)
quantile(dist.gTopt26,0.975)

#Tmax 95% CI is the same as in Supplementary Table 1

###
#Robinia
###

##
#Parameter estimates
##

#Tmin
coef(fit_Nase_beta_linall_ROPS)[11]+coef(fit_Nase_beta_linall_ROPS)[12]*18.5

#Topt
coef(fit_Nase_beta_linall_ROPS)[13]+coef(fit_Nase_beta_linall_ROPS)[14]*18.5

#Tmax
coef(fit_Nase_beta_linall_ROPS)[15]+coef(fit_Nase_beta_linall_ROPS)[16]*18.5

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_ROPS),Sigma=vcov(fit_Nase_beta_linall_ROPS))

#Tmin
dist.rTmin21<-rep(NA,1000)
for(i in 1:1000){
  dist.rTmin21[i]=vmat.r[i,11]+vmat.r[i,12]*18.5
}
quantile(dist.rTmin21,0.025)
quantile(dist.rTmin21,0.975)

#Topt
dist.rTopt21<-rep(NA,1000)
for(i in 1:1000){
  dist.rTopt21[i]=vmat.r[i,13]+vmat.r[i,14]*18.5
}
quantile(dist.rTopt21,0.025)
quantile(dist.rTopt21,0.975)

#Tmax
dist.rTmax21<-rep(NA,1000)
for(i in 1:1000){
  dist.rTmax21[i]=vmat.r[i,15]+vmat.r[i,16]*18.5
}
quantile(dist.rTmax21,0.025)
quantile(dist.rTmax21,0.975)

###
#Tropical
###

##
#Parameter estimates
##

#Tmin
coef(fit_Nase_beta_linTminTopt_Trop)[21]+coef(fit_Nase_beta_linTminTopt_Trop)[22]*23.5

#Topt
coef(fit_Nase_beta_linTminTopt_Trop)[23]+coef(fit_Nase_beta_linTminTopt_Trop)[24]*23.5

#Tmax estimate is the same as in Supplementary Table 1

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.Trop<-mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_Trop),Sigma=vcov(fit_Nase_beta_linTminTopt_Trop))

#Tmin
dist.TropTmin26<-rep(NA,1000)
for(i in 1:1000){
  dist.TropTmin26[i]=vmat.Trop[i,21]+vmat.Trop[i,22]*23.5
}
quantile(dist.TropTmin26,0.025)
quantile(dist.TropTmin26,0.975)

#Topt
dist.TropTopt26<-rep(NA,1000)
for(i in 1:1000){
  dist.TropTopt26[i]=vmat.Trop[i,23]+vmat.Trop[i,24]*23.5
}
quantile(dist.TropTopt26,0.025)
quantile(dist.TropTopt26,0.975)

#Tmax 95% CI is the same as in Supplementary Table 1

###
#Temperate
###

##
#Parameter estimates
##

#Tmin
coef(fit_Nase_beta_linall_Temp)[20]+coef(fit_Nase_beta_linall_Temp)[21]*18.5

#Topt
coef(fit_Nase_beta_linall_Temp)[22]+coef(fit_Nase_beta_linall_Temp)[23]*18.5

#Tmax
coef(fit_Nase_beta_linall_Temp)[24]+coef(fit_Nase_beta_linall_Temp)[25]*18.5

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.Temp<-mvrnorm(1000,mu=coef(fit_Nase_beta_linall_Temp),Sigma=vcov(fit_Nase_beta_linall_Temp))

#Tmin
dist.TempTmin21<-rep(NA,1000)
for(i in 1:1000){
  dist.TempTmin21[i]=vmat.Temp[i,20]+vmat.Temp[i,21]*18.5
}
quantile(dist.TempTmin21,0.025)
quantile(dist.TempTmin21,0.975)

#Topt
dist.TempTopt21<-rep(NA,1000)
for(i in 1:1000){
  dist.TempTopt21[i]=vmat.Temp[i,22]+vmat.Temp[i,23]*18.5
}
quantile(dist.TempTopt21,0.025)
quantile(dist.TempTopt21,0.975)

#Tmax
dist.TempTmax21<-rep(NA,1000)
for(i in 1:1000){
  dist.TempTmax21[i]=vmat.Temp[i,24]+vmat.Temp[i,25]*18.5
}
quantile(dist.TempTmax21,0.025)
quantile(dist.TempTmax21,0.975)

###
#For overall function, fit maximum likelihood model (below)
###

#NLL function
Nase_beta_bin_4sp_normNLL <- function(sdNase,ymax21a,ymax21b,ymax21c,ymax21d,ymax21e,ymax21f,
                                      ymax21g,ymax21h,ymax21i,ymax21j,ymax21k,ymax21l,
                                      ymax26a,ymax26b,ymax26c,ymax26d,ymax26e,ymax26f,
                                      ymax26g,ymax26h,ymax26i,ymax26j,ymax26k,ymax26l,
                                      ymax31a,ymax31b,ymax31c,ymax31d,ymax31e,ymax31f,
                                      ymax31g,ymax31h,ymax31i,ymax31j,ymax31k,ymax31l,
                                      Tmax,
                                      Tmin,
                                      Topt,
                                      T21a,T21b,T21c,T21d,T21e,T21f,
                                      T21g,T21h,T21i,T21j,T21k,T21l,
                                      T26a,T26b,T26c,T26d,T26e,T26f,
                                      T26g,T26h,T26i,T26j,T26k,T26l,
                                      T31a,T31b,T31c,T31d,T31e,T31f,
                                      T31g,T31h,T31i,T31j,T31k,T31l,
                                      Nasedat21a,Nasedat21b,Nasedat21c,Nasedat21d,Nasedat21e,Nasedat21f,
                                      Nasedat21g,Nasedat21h,Nasedat21i,Nasedat21j,Nasedat21k,Nasedat21l,
                                      Nasedat26a,Nasedat26b,Nasedat26c,Nasedat26d,Nasedat26e,Nasedat26f,
                                      Nasedat26g,Nasedat26h,Nasedat26i,Nasedat26j,Nasedat26k,Nasedat26l,
                                      Nasedat31a,Nasedat31b,Nasedat31c,Nasedat31d,Nasedat31e,Nasedat31f,
                                      Nasedat31g,Nasedat31h,Nasedat31i,Nasedat31j,Nasedat31k,Nasedat31l){
  Nasemean21a <- beta(ymax21a,Tmin,Topt,Tmax,T21a)
  Nasemean21b <- beta(ymax21b,Tmin,Topt,Tmax,T21b)
  Nasemean21c <- beta(ymax21c,Tmin,Topt,Tmax,T21c)
  Nasemean21d <- beta(ymax21d,Tmin,Topt,Tmax,T21d)
  Nasemean21e <- beta(ymax21e,Tmin,Topt,Tmax,T21e)
  Nasemean21f <- beta(ymax21f,Tmin,Topt,Tmax,T21f)
  Nasemean21g <- beta(ymax21g,Tmin,Topt,Tmax,T21g)
  Nasemean21h <- beta(ymax21h,Tmin,Topt,Tmax,T21h)
  Nasemean21i <- beta(ymax21i,Tmin,Topt,Tmax,T21i)
  Nasemean21j <- beta(ymax21j,Tmin,Topt,Tmax,T21j)
  Nasemean21k <- beta(ymax21k,Tmin,Topt,Tmax,T21k)
  Nasemean21l <- beta(ymax21l,Tmin,Topt,Tmax,T21l)
  Nasemean26a <- beta(ymax26a,Tmin,Topt,Tmax,T26a)
  Nasemean26b <- beta(ymax26b,Tmin,Topt,Tmax,T26b)
  Nasemean26c <- beta(ymax26c,Tmin,Topt,Tmax,T26c)
  Nasemean26d <- beta(ymax26d,Tmin,Topt,Tmax,T26d)
  Nasemean26e <- beta(ymax26e,Tmin,Topt,Tmax,T26e)
  Nasemean26f <- beta(ymax26f,Tmin,Topt,Tmax,T26f)
  Nasemean26g <- beta(ymax26g,Tmin,Topt,Tmax,T26g)
  Nasemean26h <- beta(ymax26h,Tmin,Topt,Tmax,T26h)
  Nasemean26i <- beta(ymax26i,Tmin,Topt,Tmax,T26i)
  Nasemean26j <- beta(ymax26j,Tmin,Topt,Tmax,T26j)
  Nasemean26k <- beta(ymax26k,Tmin,Topt,Tmax,T26k)
  Nasemean26l <- beta(ymax26l,Tmin,Topt,Tmax,T26l)
  Nasemean31a <- beta(ymax31a,Tmin,Topt,Tmax,T31a)
  Nasemean31b <- beta(ymax31b,Tmin,Topt,Tmax,T31b)
  Nasemean31c <- beta(ymax31c,Tmin,Topt,Tmax,T31c)
  Nasemean31d <- beta(ymax31d,Tmin,Topt,Tmax,T31d)
  Nasemean31e <- beta(ymax31e,Tmin,Topt,Tmax,T31e)
  Nasemean31f <- beta(ymax31f,Tmin,Topt,Tmax,T31f)
  Nasemean31g <- beta(ymax31g,Tmin,Topt,Tmax,T31g)
  Nasemean31h <- beta(ymax31h,Tmin,Topt,Tmax,T31h)
  Nasemean31i <- beta(ymax31i,Tmin,Topt,Tmax,T31i)
  Nasemean31j <- beta(ymax31j,Tmin,Topt,Tmax,T31j)
  Nasemean31k <- beta(ymax31k,Tmin,Topt,Tmax,T31k)
  Nasemean31l <- beta(ymax31l,Tmin,Topt,Tmax,T31l)
  -(sum(dnorm(Nasedat21a,mean=Nasemean21a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21b,mean=Nasemean21b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21c,mean=Nasemean21c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21d,mean=Nasemean21d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21e,mean=Nasemean21e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21f,mean=Nasemean21f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21g,mean=Nasemean21g,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21h,mean=Nasemean21h,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21i,mean=Nasemean21i,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21j,mean=Nasemean21j,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21k,mean=Nasemean21k,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat21l,mean=Nasemean21l,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26a,mean=Nasemean26a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26b,mean=Nasemean26b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26c,mean=Nasemean26c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat26d,mean=Nasemean26d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26e,mean=Nasemean26e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26f,mean=Nasemean26f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat26g,mean=Nasemean26g,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26h,mean=Nasemean26h,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26i,mean=Nasemean26i,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26j,mean=Nasemean26j,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26k,mean=Nasemean26k,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat26l,mean=Nasemean26l,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31d,mean=Nasemean31d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31e,mean=Nasemean31e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31f,mean=Nasemean31f,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31g,mean=Nasemean31g,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31h,mean=Nasemean31h,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31i,mean=Nasemean31i,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31j,mean=Nasemean31j,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31k,mean=Nasemean31k,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31l,mean=Nasemean31l,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#Maximum likelihood fit
fit_Nase_beta_4sp_bin <- mle2(Nase_beta_bin_4sp_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,ymax21d=1,ymax21e=1,ymax21f=1,
                                                                   ymax21g=1,ymax21h=1,ymax21i=1,ymax21j=1,ymax21k=1,ymax21l=1,
                                                                   ymax26a=1,ymax26b=1,ymax26c=1,ymax26d=1,ymax26e=1,ymax26f=1,
                                                                   ymax26g=1,ymax26h=1,ymax26i=1,ymax26j=1,ymax26k=1,ymax26l=1,
                                                                   ymax31a=1,ymax31b=1,ymax31c=1,ymax31d=1,ymax31e=1,ymax31f=1,
                                                                   ymax31g=1,ymax31h=1,ymax31i=1,ymax31j=1,ymax31k=1,ymax31l=1,
                                                                   Tmax=44,Tmin=-2,Topt=32),
                              data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                        T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                        T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                        T21d=ALRU21a$Temperature,T21e=ALRU21b$Temperature,T21f=ALRU21c$Temperature,
                                        T26d=ALRU26a$Temperature,T26e=ALRU26b$Temperature,T26f=ALRU26c$Temperature,
                                        T31d=ALRU31a$Temperature,T31e=ALRU31b$Temperature,T31f=ALRU31c$Temperature,
                                        T21g=GLSE21a$Temperature,T21h=GLSE21b$Temperature,T21i=GLSE21c$Temperature,
                                        T26g=GLSE26a$Temperature,T26h=GLSE26b$Temperature,T26i=GLSE26c$Temperature,
                                        T31g=GLSE31a$Temperature,T31h=GLSE31b$Temperature,T31i=GLSE31c$Temperature,
                                        T21j=ROPS21a$Temperature,T21k=ROPS21b$Temperature,T21l=ROPS21c$Temperature,
                                        T26j=ROPS26a$Temperature,T26k=ROPS26b$Temperature,T26l=ROPS26c$Temperature,
                                        T31j=ROPS31a$Temperature,T31k=ROPS31b$Temperature,T31l=ROPS31c$Temperature,
                                        Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                        Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                        Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222,
                                        Nasedat21d=ALRU21a$Vmax/10035.24751,Nasedat21e=ALRU21b$Vmax/7990.765838,Nasedat21f=ALRU21c$Vmax/4046.082544,
                                        Nasedat26d=ALRU26a$Vmax/2236.979133,Nasedat26e=ALRU26b$Vmax/2744.252462,Nasedat26f=ALRU26c$Vmax/1831.763172,
                                        Nasedat31d=ALRU31a$Vmax/406.3524389,Nasedat31e=ALRU31b$Vmax/1832.723413,Nasedat31f=ALRU31c$Vmax/1350.729381,
                                        Nasedat21g=GLSE21a$Vmax/308.3421323,Nasedat21h=GLSE21b$Vmax/2028.57371,Nasedat21i=GLSE21c$Vmax/879.0395611,
                                        Nasedat26g=GLSE26a$Vmax/909.4624356,Nasedat26h=GLSE26b$Vmax/1373.906634,Nasedat26i=GLSE26c$Vmax/705.9735103,
                                        Nasedat31g=GLSE31a$Vmax/1196.373627,Nasedat31h=GLSE31b$Vmax/312.4607561,Nasedat31i=GLSE31c$Vmax/209.4577018,
                                        Nasedat21j=ROPS21a$Vmax/1922.245351,Nasedat21k=ROPS21b$Vmax/1690.867422,Nasedat21l=ROPS21c$Vmax/450.8292229,
                                        Nasedat26j=ROPS26a$Vmax/2785.452785,Nasedat26k=ROPS26b$Vmax/1572.439251,Nasedat26l=ROPS26c$Vmax/1686.352,
                                        Nasedat31j=ROPS31a$Vmax/1646.769302,Nasedat31k=ROPS31b$Vmax/920.7416901,Nasedat31l=ROPS31c$Vmax/1978.187071),
                              control=list(maxit=20000))
summary(fit_Nase_beta_4sp_bin)

#Parameter estimates and 95% CI
coef(fit_Nase_beta_4sp_bin)
fit_Nase_beta_4sp_bin.pr<-profile(fit_Nase_beta_4sp_bin)
confint(fit_Nase_beta_4sp_bin.pr)
