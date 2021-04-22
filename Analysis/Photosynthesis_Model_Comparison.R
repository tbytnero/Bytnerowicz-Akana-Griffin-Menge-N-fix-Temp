###############################################################################################################
###############################################################################################################
#This script compares models for photosynthesis (both A275 and Asat; Supplementary Tables 7 and 9) and
#then generates parameter estimates and 95% CIs for photosynthesis (Supplementary Table 2) 
###############################################################################################################
###############################################################################################################

####
#Load Necessary Packages
####

library(bbmle)
library(MASS)

####
#Read the data in from "Photo_Temp" folder
####

###
#A275
###

#Alnus rubra
A21.A275.dat<-read.csv("ACi.dat.A21.csv")
A26.A275.dat<-read.csv("ACi.dat.A26.csv")
A31.A275.dat<-read.csv("ACi.dat.A31.csv")

#Morella cerifera
M21.A275.dat<-read.csv("ACi.dat.M21.csv")
M26.A275.dat<-read.csv("ACi.dat.M26.csv")
M31.A275.dat<-read.csv("ACi.dat.M31.csv")

#Gliricidia sepium
G21.A275.dat<-read.csv("ACi.dat.G21.csv")
G26.A275.dat<-read.csv("ACi.dat.G26.csv")
G31.A275.dat<-read.csv("ACi.dat.G31.csv")

#Robinia pseudoacacia
R21.A275.dat<-read.csv("ACi.dat.R21.csv")
R26.A275.dat<-read.csv("ACi.dat.R26.csv")
R31.A275.dat<-read.csv("ACi.dat.R31.csv")

#Calculate temperature in Kelvin
A21.A275.dat$TsK<-A21.A275.dat$Temp+273.15
A26.A275.dat$TsK<-A26.A275.dat$Temp+273.15
A31.A275.dat$TsK<-A31.A275.dat$Temp+273.15
M21.A275.dat$TsK<-M21.A275.dat$Temp+273.15
M26.A275.dat$TsK<-M26.A275.dat$Temp+273.15
M31.A275.dat$TsK<-M31.A275.dat$Temp+273.15
G21.A275.dat$TsK<-G21.A275.dat$Temp+273.15
G26.A275.dat$TsK<-G26.A275.dat$Temp+273.15
G31.A275.dat$TsK<-G31.A275.dat$Temp+273.15
R21.A275.dat$TsK<-R21.A275.dat$Temp+273.15
R26.A275.dat$TsK<-R26.A275.dat$Temp+273.15
R31.A275.dat$TsK<-R31.A275.dat$Temp+273.15

###
#Asat
###

A400.data<-read.csv("A400_data.csv")

#Alnus rubra
A21.Asat.dat<-A400.data[A400.data$Tr=="A21",]
A26.Asat.dat<-A400.data[A400.data$Tr=="A26",]
A31.Asat.dat<-A400.data[A400.data$Tr=="A31",]

#Gliricidia sepium
G21.Asat.dat<-A400.data[A400.data$Tr=="G21",]
G26.Asat.dat<-A400.data[A400.data$Tr=="G26",]
G31.Asat.dat<-A400.data[A400.data$Tr=="G31",]

#Morella cerifera
M21.Asat.dat<-A400.data[A400.data$Tr=="M21",]
M26.Asat.dat<-A400.data[A400.data$Tr=="M26",]
M31.Asat.dat<-A400.data[A400.data$Tr=="M31",]

#Robinia pseudoacacia
R21.Asat.dat<-A400.data[A400.data$Tr=="R21",]
R26.Asat.dat<-A400.data[A400.data$Tr=="R26",]
R31.Asat.dat<-A400.data[A400.data$Tr=="R31",]

#Calculate temperature in Kelvin
A21.Asat.dat$TsK<-A21.Asat.dat$Temp+273.15
A26.Asat.dat$TsK<-A26.Asat.dat$Temp+273.15
A31.Asat.dat$TsK<-A31.Asat.dat$Temp+273.15
M21.Asat.dat$TsK<-M21.Asat.dat$Temp+273.15
M26.Asat.dat$TsK<-M26.Asat.dat$Temp+273.15
M31.Asat.dat$TsK<-M31.Asat.dat$Temp+273.15
G21.Asat.dat$TsK<-G21.Asat.dat$Temp+273.15
G26.Asat.dat$TsK<-G26.Asat.dat$Temp+273.15
G31.Asat.dat$TsK<-G31.Asat.dat$Temp+273.15
R21.Asat.dat$TsK<-R21.Asat.dat$Temp+273.15
R26.Asat.dat$TsK<-R26.Asat.dat$Temp+273.15
R31.Asat.dat$TsK<-R31.Asat.dat$Temp+273.15

###############################################################################################################
#First we compare equations 2-5 (Supplementary Table 7)
###############################################################################################################

####
#Define functions
####

#Modified beta (equation 5)
beta <- function(ymax,Tmin,Topt,Tmax,T){
  y <- ymax*(Tmax-T)/(Tmax-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)))
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
  y <- ymax-b*(T-Topt)^2
  y
}

#Normal (equation 2)
norm<-function(ymax,Topt,s,T){
  y<-ymax*exp(-(T-Topt)^2/(2*s^2))
  y
}

####
#Negative log-likelihood (NLL) functions
####

#NLL function for the modified beta (equation 5)
Photo_beta_normNLL_all <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                  ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                  ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                  ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                  ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                  ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                  ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                  ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                  ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                  ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                  ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                  ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
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
                                  TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                  TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                  TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                  TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                  TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                  TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                  TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                  TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                  TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                  TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                  TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                  TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                  PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                  PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                  PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                  PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                  PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                  PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                  PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                  PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                  PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                  PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                  PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                  PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta(ymaxM21a,TminM21,ToptM21,TmaxM21,TempM21a)
  PhotomeanM21b <- beta(ymaxM21b,TminM21,ToptM21,TmaxM21,TempM21b)
  PhotomeanM21c <- beta(ymaxM21c,TminM21,ToptM21,TmaxM21,TempM21c)
  PhotomeanM21d <- beta(ymaxM21d,TminM21,ToptM21,TmaxM21,TempM21d)
  PhotomeanM21e <- beta(ymaxM21e,TminM21,ToptM21,TmaxM21,TempM21e)
  PhotomeanM21f <- beta(ymaxM21f,TminM21,ToptM21,TmaxM21,TempM21f)
  PhotomeanM26a <- beta(ymaxM26a,TminM26,ToptM26,TmaxM26,TempM26a)
  PhotomeanM26b <- beta(ymaxM26b,TminM26,ToptM26,TmaxM26,TempM26b)
  PhotomeanM26c <- beta(ymaxM26c,TminM26,ToptM26,TmaxM26,TempM26c)
  PhotomeanM26d <- beta(ymaxM26d,TminM26,ToptM26,TmaxM26,TempM26d)
  PhotomeanM26e <- beta(ymaxM26e,TminM26,ToptM26,TmaxM26,TempM26e)
  PhotomeanM26f <- beta(ymaxM26f,TminM26,ToptM26,TmaxM26,TempM26f)
  PhotomeanM31a <- beta(ymaxM31a,TminM31,ToptM31,TmaxM31,TempM31a)
  PhotomeanM31b <- beta(ymaxM31b,TminM31,ToptM31,TmaxM31,TempM31b)
  PhotomeanM31c <- beta(ymaxM31c,TminM31,ToptM31,TmaxM31,TempM31c)
  PhotomeanM31d <- beta(ymaxM31d,TminM31,ToptM31,TmaxM31,TempM31d)
  PhotomeanM31e <- beta(ymaxM31e,TminM31,ToptM31,TmaxM31,TempM31e)
  PhotomeanM31f <- beta(ymaxM31f,TminM31,ToptM31,TmaxM31,TempM31f)
  PhotomeanA21a <- beta(ymaxA21a,TminA21,ToptA21,TmaxA21,TempA21a)
  PhotomeanA21b <- beta(ymaxA21b,TminA21,ToptA21,TmaxA21,TempA21b)
  PhotomeanA21c <- beta(ymaxA21c,TminA21,ToptA21,TmaxA21,TempA21c)
  PhotomeanA21d <- beta(ymaxA21d,TminA21,ToptA21,TmaxA21,TempA21d)
  PhotomeanA21e <- beta(ymaxA21e,TminA21,ToptA21,TmaxA21,TempA21e)
  PhotomeanA21f <- beta(ymaxA21f,TminA21,ToptA21,TmaxA21,TempA21f)
  PhotomeanA26a <- beta(ymaxA26a,TminA26,ToptA26,TmaxA26,TempA26a)
  PhotomeanA26b <- beta(ymaxA26b,TminA26,ToptA26,TmaxA26,TempA26b)
  PhotomeanA26c <- beta(ymaxA26c,TminA26,ToptA26,TmaxA26,TempA26c)
  PhotomeanA26d <- beta(ymaxA26d,TminA26,ToptA26,TmaxA26,TempA26d)
  PhotomeanA26e <- beta(ymaxA26e,TminA26,ToptA26,TmaxA26,TempA26e)
  PhotomeanA26f <- beta(ymaxA26f,TminA26,ToptA26,TmaxA26,TempA26f)
  PhotomeanA31a <- beta(ymaxA31a,TminA31,ToptA31,TmaxA31,TempA31a)
  PhotomeanA31b <- beta(ymaxA31b,TminA31,ToptA31,TmaxA31,TempA31b)
  PhotomeanA31c <- beta(ymaxA31c,TminA31,ToptA31,TmaxA31,TempA31c)
  PhotomeanA31d <- beta(ymaxA31d,TminA31,ToptA31,TmaxA31,TempA31d)
  PhotomeanA31e <- beta(ymaxA31e,TminA31,ToptA31,TmaxA31,TempA31e)
  PhotomeanA31f <- beta(ymaxA31f,TminA31,ToptA31,TmaxA31,TempA31f)
  PhotomeanG21a <- beta(ymaxG21a,TminG21,ToptG21,TmaxG21,TempG21a)
  PhotomeanG21b <- beta(ymaxG21b,TminG21,ToptG21,TmaxG21,TempG21b)
  PhotomeanG21c <- beta(ymaxG21c,TminG21,ToptG21,TmaxG21,TempG21c)
  PhotomeanG21d <- beta(ymaxG21d,TminG21,ToptG21,TmaxG21,TempG21d)
  PhotomeanG21e <- beta(ymaxG21e,TminG21,ToptG21,TmaxG21,TempG21e)
  PhotomeanG21f <- beta(ymaxG21f,TminG21,ToptG21,TmaxG21,TempG21f)
  PhotomeanG26a <- beta(ymaxG26a,TminG26,ToptG26,TmaxG26,TempG26a)
  PhotomeanG26b <- beta(ymaxG26b,TminG26,ToptG26,TmaxG26,TempG26b)
  PhotomeanG26c <- beta(ymaxG26c,TminG26,ToptG26,TmaxG26,TempG26c)
  PhotomeanG26d <- beta(ymaxG26d,TminG26,ToptG26,TmaxG26,TempG26d)
  PhotomeanG26e <- beta(ymaxG26e,TminG26,ToptG26,TmaxG26,TempG26e)
  PhotomeanG26f <- beta(ymaxG26f,TminG26,ToptG26,TmaxG26,TempG26f)
  PhotomeanG26g <- beta(ymaxG26g,TminG26,ToptG26,TmaxG26,TempG26g)
  PhotomeanG31a <- beta(ymaxG31a,TminG31,ToptG31,TmaxG31,TempG31a)
  PhotomeanG31b <- beta(ymaxG31b,TminG31,ToptG31,TmaxG31,TempG31b)
  PhotomeanG31c <- beta(ymaxG31c,TminG31,ToptG31,TmaxG31,TempG31c)
  PhotomeanG31d <- beta(ymaxG31d,TminG31,ToptG31,TmaxG31,TempG31d)
  PhotomeanG31e <- beta(ymaxG31e,TminG31,ToptG31,TmaxG31,TempG31e)
  PhotomeanG31f <- beta(ymaxG31f,TminG31,ToptG31,TmaxG31,TempG31f)
  PhotomeanR21a <- beta(ymaxR21a,TminR21,ToptR21,TmaxR21,TempR21a)
  PhotomeanR21b <- beta(ymaxR21b,TminR21,ToptR21,TmaxR21,TempR21b)
  PhotomeanR21c <- beta(ymaxR21c,TminR21,ToptR21,TmaxR21,TempR21c)
  PhotomeanR21d <- beta(ymaxR21d,TminR21,ToptR21,TmaxR21,TempR21d)
  PhotomeanR21e <- beta(ymaxR21e,TminR21,ToptR21,TmaxR21,TempR21e)
  PhotomeanR21f <- beta(ymaxR21f,TminR21,ToptR21,TmaxR21,TempR21f)
  PhotomeanR21g <- beta(ymaxR21g,TminR21,ToptR21,TmaxR21,TempR21g)
  PhotomeanR21h <- beta(ymaxR21h,TminR21,ToptR21,TmaxR21,TempR21h)
  PhotomeanR26a <- beta(ymaxR26a,TminR26,ToptR26,TmaxR26,TempR26a)
  PhotomeanR26b <- beta(ymaxR26b,TminR26,ToptR26,TmaxR26,TempR26b)
  PhotomeanR26c <- beta(ymaxR26c,TminR26,ToptR26,TmaxR26,TempR26c)
  PhotomeanR26d <- beta(ymaxR26d,TminR26,ToptR26,TmaxR26,TempR26d)
  PhotomeanR26e <- beta(ymaxR26e,TminR26,ToptR26,TmaxR26,TempR26e)
  PhotomeanR26f <- beta(ymaxR26f,TminR26,ToptR26,TmaxR26,TempR26f)
  PhotomeanR31a <- beta(ymaxR31a,TminR31,ToptR31,TmaxR31,TempR31a)
  PhotomeanR31b <- beta(ymaxR31b,TminR31,ToptR31,TmaxR31,TempR31b)
  PhotomeanR31c <- beta(ymaxR31c,TminR31,ToptR31,TmaxR31,TempR31c)
  PhotomeanR31d <- beta(ymaxR31d,TminR31,ToptR31,TmaxR31,TempR31d)
  PhotomeanR31e <- beta(ymaxR31e,TminR31,ToptR31,TmaxR31,TempR31e)
  PhotomeanR31f <- beta(ymaxR31f,TminR31,ToptR31,TmaxR31,TempR31f)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for the peaked Arrhenius (equation 4)
Photo_peak.ar_normNLL_all <- function(sdPhoto,k25M21a,k25M21b,k25M21c,k25M21d,k25M21e,k25M21f,
                                      k25M26a,k25M26b,k25M26c,k25M26d,k25M26e,k25M26f,
                                      k25M31a,k25M31b,k25M31c,k25M31d,k25M31e,k25M31f,
                                      k25A21a,k25A21b,k25A21c,k25A21d,k25A21e,k25A21f,
                                      k25A26a,k25A26b,k25A26c,k25A26d,k25A26e,k25A26f,
                                      k25A31a,k25A31b,k25A31c,k25A31d,k25A31e,k25A31f,
                                      k25G21a,k25G21b,k25G21c,k25G21d,k25G21e,k25G21f,
                                      k25G26a,k25G26b,k25G26c,k25G26d,k25G26e,k25G26f,k25G26g,
                                      k25G31a,k25G31b,k25G31c,k25G31d,k25G31e,k25G31f,
                                      k25R21a,k25R21b,k25R21c,k25R21d,k25R21e,k25R21f,k25R21g,k25R21h,
                                      k25R26a,k25R26b,k25R26c,k25R26d,k25R26e,k25R26f,
                                      k25R31a,k25R31b,k25R31c,k25R31d,k25R31e,k25R31f,
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
                                      TkM21a,TkM21b,TkM21c,TkM21d,TkM21e,TkM21f,
                                      TkM26a,TkM26b,TkM26c,TkM26d,TkM26e,TkM26f,
                                      TkM31a,TkM31b,TkM31c,TkM31d,TkM31e,TkM31f,
                                      TkA21a,TkA21b,TkA21c,TkA21d,TkA21e,TkA21f,
                                      TkA26a,TkA26b,TkA26c,TkA26d,TkA26e,TkA26f,
                                      TkA31a,TkA31b,TkA31c,TkA31d,TkA31e,TkA31f,
                                      TkG21a,TkG21b,TkG21c,TkG21d,TkG21e,TkG21f,
                                      TkG26a,TkG26b,TkG26c,TkG26d,TkG26e,TkG26f,TkG26g,
                                      TkG31a,TkG31b,TkG31c,TkG31d,TkG31e,TkG31f,
                                      TkR21a,TkR21b,TkR21c,TkR21d,TkR21e,TkR21f,TkR21g,TkR21h,
                                      TkR26a,TkR26b,TkR26c,TkR26d,TkR26e,TkR26f,
                                      TkR31a,TkR31b,TkR31c,TkR31d,TkR31e,TkR31f,
                                      PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                      PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                      PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                      PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                      PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                      PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                      PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                      PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                      PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                      PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                      PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                      PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- peak.ar(k25M21a,EaM21,ToptM21,HdM21,TkM21a)
  PhotomeanM21b <- peak.ar(k25M21b,EaM21,ToptM21,HdM21,TkM21b)
  PhotomeanM21c <- peak.ar(k25M21c,EaM21,ToptM21,HdM21,TkM21c)
  PhotomeanM21d <- peak.ar(k25M21d,EaM21,ToptM21,HdM21,TkM21d)
  PhotomeanM21e <- peak.ar(k25M21e,EaM21,ToptM21,HdM21,TkM21e)
  PhotomeanM21f <- peak.ar(k25M21f,EaM21,ToptM21,HdM21,TkM21f)
  PhotomeanM26a <- peak.ar(k25M26a,EaM26,ToptM26,HdM26,TkM26a)
  PhotomeanM26b <- peak.ar(k25M26b,EaM26,ToptM26,HdM26,TkM26b)
  PhotomeanM26c <- peak.ar(k25M26c,EaM26,ToptM26,HdM26,TkM26c)
  PhotomeanM26d <- peak.ar(k25M26d,EaM26,ToptM26,HdM26,TkM26d)
  PhotomeanM26e <- peak.ar(k25M26e,EaM26,ToptM26,HdM26,TkM26e)
  PhotomeanM26f <- peak.ar(k25M26f,EaM26,ToptM26,HdM26,TkM26f)
  PhotomeanM31a <- peak.ar(k25M31a,EaM31,ToptM31,HdM31,TkM31a)
  PhotomeanM31b <- peak.ar(k25M31b,EaM31,ToptM31,HdM31,TkM31b)
  PhotomeanM31c <- peak.ar(k25M31c,EaM31,ToptM31,HdM31,TkM31c)
  PhotomeanM31d <- peak.ar(k25M31d,EaM31,ToptM31,HdM31,TkM31d)
  PhotomeanM31e <- peak.ar(k25M31e,EaM31,ToptM31,HdM31,TkM31e)
  PhotomeanM31f <- peak.ar(k25M31f,EaM31,ToptM31,HdM31,TkM31f)
  PhotomeanA21a <- peak.ar(k25A21a,EaA21,ToptA21,HdA21,TkA21a)
  PhotomeanA21b <- peak.ar(k25A21b,EaA21,ToptA21,HdA21,TkA21b)
  PhotomeanA21c <- peak.ar(k25A21c,EaA21,ToptA21,HdA21,TkA21c)
  PhotomeanA21d <- peak.ar(k25A21d,EaA21,ToptA21,HdA21,TkA21d)
  PhotomeanA21e <- peak.ar(k25A21e,EaA21,ToptA21,HdA21,TkA21e)
  PhotomeanA21f <- peak.ar(k25A21f,EaA21,ToptA21,HdA21,TkA21f)
  PhotomeanA26a <- peak.ar(k25A26a,EaA26,ToptA26,HdA26,TkA26a)
  PhotomeanA26b <- peak.ar(k25A26b,EaA26,ToptA26,HdA26,TkA26b)
  PhotomeanA26c <- peak.ar(k25A26c,EaA26,ToptA26,HdA26,TkA26c)
  PhotomeanA26d <- peak.ar(k25A26d,EaA26,ToptA26,HdA26,TkA26d)
  PhotomeanA26e <- peak.ar(k25A26e,EaA26,ToptA26,HdA26,TkA26e)
  PhotomeanA26f <- peak.ar(k25A26f,EaA26,ToptA26,HdA26,TkA26f)
  PhotomeanA31a <- peak.ar(k25A31a,EaA31,ToptA31,HdA31,TkA31a)
  PhotomeanA31b <- peak.ar(k25A31b,EaA31,ToptA31,HdA31,TkA31b)
  PhotomeanA31c <- peak.ar(k25A31c,EaA31,ToptA31,HdA31,TkA31c)
  PhotomeanA31d <- peak.ar(k25A31d,EaA31,ToptA31,HdA31,TkA31d)
  PhotomeanA31e <- peak.ar(k25A31e,EaA31,ToptA31,HdA31,TkA31e)
  PhotomeanA31f <- peak.ar(k25A31f,EaA31,ToptA31,HdA31,TkA31f)
  PhotomeanG21a <- peak.ar(k25G21a,EaG21,ToptG21,HdG21,TkG21a)
  PhotomeanG21b <- peak.ar(k25G21b,EaG21,ToptG21,HdG21,TkG21b)
  PhotomeanG21c <- peak.ar(k25G21c,EaG21,ToptG21,HdG21,TkG21c)
  PhotomeanG21d <- peak.ar(k25G21d,EaG21,ToptG21,HdG21,TkG21d)
  PhotomeanG21e <- peak.ar(k25G21e,EaG21,ToptG21,HdG21,TkG21e)
  PhotomeanG21f <- peak.ar(k25G21f,EaG21,ToptG21,HdG21,TkG21f)
  PhotomeanG26a <- peak.ar(k25G26a,EaG26,ToptG26,HdG26,TkG26a)
  PhotomeanG26b <- peak.ar(k25G26b,EaG26,ToptG26,HdG26,TkG26b)
  PhotomeanG26c <- peak.ar(k25G26c,EaG26,ToptG26,HdG26,TkG26c)
  PhotomeanG26d <- peak.ar(k25G26d,EaG26,ToptG26,HdG26,TkG26d)
  PhotomeanG26e <- peak.ar(k25G26e,EaG26,ToptG26,HdG26,TkG26e)
  PhotomeanG26f <- peak.ar(k25G26f,EaG26,ToptG26,HdG26,TkG26f)
  PhotomeanG26g <- peak.ar(k25G26g,EaG26,ToptG26,HdG26,TkG26g)
  PhotomeanG31a <- peak.ar(k25G31a,EaG31,ToptG31,HdG31,TkG31a)
  PhotomeanG31b <- peak.ar(k25G31b,EaG31,ToptG31,HdG31,TkG31b)
  PhotomeanG31c <- peak.ar(k25G31c,EaG31,ToptG31,HdG31,TkG31c)
  PhotomeanG31d <- peak.ar(k25G31d,EaG31,ToptG31,HdG31,TkG31d)
  PhotomeanG31e <- peak.ar(k25G31e,EaG31,ToptG31,HdG31,TkG31e)
  PhotomeanG31f <- peak.ar(k25G31f,EaG31,ToptG31,HdG31,TkG31f)
  PhotomeanR21a <- peak.ar(k25R21a,EaR21,ToptR21,HdR21,TkR21a)
  PhotomeanR21b <- peak.ar(k25R21b,EaR21,ToptR21,HdR21,TkR21b)
  PhotomeanR21c <- peak.ar(k25R21c,EaR21,ToptR21,HdR21,TkR21c)
  PhotomeanR21d <- peak.ar(k25R21d,EaR21,ToptR21,HdR21,TkR21d)
  PhotomeanR21e <- peak.ar(k25R21e,EaR21,ToptR21,HdR21,TkR21e)
  PhotomeanR21f <- peak.ar(k25R21f,EaR21,ToptR21,HdR21,TkR21f)
  PhotomeanR21g <- peak.ar(k25R21g,EaR21,ToptR21,HdR21,TkR21g)
  PhotomeanR21h <- peak.ar(k25R21h,EaR21,ToptR21,HdR21,TkR21h)
  PhotomeanR26a <- peak.ar(k25R26a,EaR26,ToptR26,HdR26,TkR26a)
  PhotomeanR26b <- peak.ar(k25R26b,EaR26,ToptR26,HdR26,TkR26b)
  PhotomeanR26c <- peak.ar(k25R26c,EaR26,ToptR26,HdR26,TkR26c)
  PhotomeanR26d <- peak.ar(k25R26d,EaR26,ToptR26,HdR26,TkR26d)
  PhotomeanR26e <- peak.ar(k25R26e,EaR26,ToptR26,HdR26,TkR26e)
  PhotomeanR26f <- peak.ar(k25R26f,EaR26,ToptR26,HdR26,TkR26f)
  PhotomeanR31a <- peak.ar(k25R31a,EaR31,ToptR31,HdR31,TkR31a)
  PhotomeanR31b <- peak.ar(k25R31b,EaR31,ToptR31,HdR31,TkR31b)
  PhotomeanR31c <- peak.ar(k25R31c,EaR31,ToptR31,HdR31,TkR31c)
  PhotomeanR31d <- peak.ar(k25R31d,EaR31,ToptR31,HdR31,TkR31d)
  PhotomeanR31e <- peak.ar(k25R31e,EaR31,ToptR31,HdR31,TkR31e)
  PhotomeanR31f <- peak.ar(k25R31f,EaR31,ToptR31,HdR31,TkR31f)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for the quadratic (equation 3)
Photo_quad_normNLL_all <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                  ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                  ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                  ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                  ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                  ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                  ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                  ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                  ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                  ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                  ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                  ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                  bM21,bM26,bM31,
                                  bA21,bA26,bA31,
                                  bG21,bG26,bG31,
                                  bR21,bR26,bR31,
                                  ToptM21,ToptM26,ToptM31,
                                  ToptA21,ToptA26,ToptA31,
                                  ToptG21,ToptG26,ToptG31,
                                  ToptR21,ToptR26,ToptR31,
                                  TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                  TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                  TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                  TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                  TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                  TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                  TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                  TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                  TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                  TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                  TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                  TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                  PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                  PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                  PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                  PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                  PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                  PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                  PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                  PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                  PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                  PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                  PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                  PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- quad(ymaxM21a,bM21,ToptM21,TempM21a)
  PhotomeanM21b <- quad(ymaxM21b,bM21,ToptM21,TempM21b)
  PhotomeanM21c <- quad(ymaxM21c,bM21,ToptM21,TempM21c)
  PhotomeanM21d <- quad(ymaxM21d,bM21,ToptM21,TempM21d)
  PhotomeanM21e <- quad(ymaxM21e,bM21,ToptM21,TempM21e)
  PhotomeanM21f <- quad(ymaxM21f,bM21,ToptM21,TempM21f)
  PhotomeanM26a <- quad(ymaxM26a,bM26,ToptM26,TempM26a)
  PhotomeanM26b <- quad(ymaxM26b,bM26,ToptM26,TempM26b)
  PhotomeanM26c <- quad(ymaxM26c,bM26,ToptM26,TempM26c)
  PhotomeanM26d <- quad(ymaxM26d,bM26,ToptM26,TempM26d)
  PhotomeanM26e <- quad(ymaxM26e,bM26,ToptM26,TempM26e)
  PhotomeanM26f <- quad(ymaxM26f,bM26,ToptM26,TempM26f)
  PhotomeanM31a <- quad(ymaxM31a,bM31,ToptM31,TempM31a)
  PhotomeanM31b <- quad(ymaxM31b,bM31,ToptM31,TempM31b)
  PhotomeanM31c <- quad(ymaxM31c,bM31,ToptM31,TempM31c)
  PhotomeanM31d <- quad(ymaxM31d,bM31,ToptM31,TempM31d)
  PhotomeanM31e <- quad(ymaxM31e,bM31,ToptM31,TempM31e)
  PhotomeanM31f <- quad(ymaxM31f,bM31,ToptM31,TempM31f)
  PhotomeanA21a <- quad(ymaxA21a,bA21,ToptA21,TempA21a)
  PhotomeanA21b <- quad(ymaxA21b,bA21,ToptA21,TempA21b)
  PhotomeanA21c <- quad(ymaxA21c,bA21,ToptA21,TempA21c)
  PhotomeanA21d <- quad(ymaxA21d,bA21,ToptA21,TempA21d)
  PhotomeanA21e <- quad(ymaxA21e,bA21,ToptA21,TempA21e)
  PhotomeanA21f <- quad(ymaxA21f,bA21,ToptA21,TempA21f)
  PhotomeanA26a <- quad(ymaxA26a,bA26,ToptA26,TempA26a)
  PhotomeanA26b <- quad(ymaxA26b,bA26,ToptA26,TempA26b)
  PhotomeanA26c <- quad(ymaxA26c,bA26,ToptA26,TempA26c)
  PhotomeanA26d <- quad(ymaxA26d,bA26,ToptA26,TempA26d)
  PhotomeanA26e <- quad(ymaxA26e,bA26,ToptA26,TempA26e)
  PhotomeanA26f <- quad(ymaxA26f,bA26,ToptA26,TempA26f)
  PhotomeanA31a <- quad(ymaxA31a,bA31,ToptA31,TempA31a)
  PhotomeanA31b <- quad(ymaxA31b,bA31,ToptA31,TempA31b)
  PhotomeanA31c <- quad(ymaxA31c,bA31,ToptA31,TempA31c)
  PhotomeanA31d <- quad(ymaxA31d,bA31,ToptA31,TempA31d)
  PhotomeanA31e <- quad(ymaxA31e,bA31,ToptA31,TempA31e)
  PhotomeanA31f <- quad(ymaxA31f,bA31,ToptA31,TempA31f)
  PhotomeanG21a <- quad(ymaxG21a,bG21,ToptG21,TempG21a)
  PhotomeanG21b <- quad(ymaxG21b,bG21,ToptG21,TempG21b)
  PhotomeanG21c <- quad(ymaxG21c,bG21,ToptG21,TempG21c)
  PhotomeanG21d <- quad(ymaxG21d,bG21,ToptG21,TempG21d)
  PhotomeanG21e <- quad(ymaxG21e,bG21,ToptG21,TempG21e)
  PhotomeanG21f <- quad(ymaxG21f,bG21,ToptG21,TempG21f)
  PhotomeanG26a <- quad(ymaxG26a,bG26,ToptG26,TempG26a)
  PhotomeanG26b <- quad(ymaxG26b,bG26,ToptG26,TempG26b)
  PhotomeanG26c <- quad(ymaxG26c,bG26,ToptG26,TempG26c)
  PhotomeanG26d <- quad(ymaxG26d,bG26,ToptG26,TempG26d)
  PhotomeanG26e <- quad(ymaxG26e,bG26,ToptG26,TempG26e)
  PhotomeanG26f <- quad(ymaxG26f,bG26,ToptG26,TempG26f)
  PhotomeanG26g <- quad(ymaxG26g,bG26,ToptG26,TempG26g)
  PhotomeanG31a <- quad(ymaxG31a,bG31,ToptG31,TempG31a)
  PhotomeanG31b <- quad(ymaxG31b,bG31,ToptG31,TempG31b)
  PhotomeanG31c <- quad(ymaxG31c,bG31,ToptG31,TempG31c)
  PhotomeanG31d <- quad(ymaxG31d,bG31,ToptG31,TempG31d)
  PhotomeanG31e <- quad(ymaxG31e,bG31,ToptG31,TempG31e)
  PhotomeanG31f <- quad(ymaxG31f,bG31,ToptG31,TempG31f)
  PhotomeanR21a <- quad(ymaxR21a,bR21,ToptR21,TempR21a)
  PhotomeanR21b <- quad(ymaxR21b,bR21,ToptR21,TempR21b)
  PhotomeanR21c <- quad(ymaxR21c,bR21,ToptR21,TempR21c)
  PhotomeanR21d <- quad(ymaxR21d,bR21,ToptR21,TempR21d)
  PhotomeanR21e <- quad(ymaxR21e,bR21,ToptR21,TempR21e)
  PhotomeanR21f <- quad(ymaxR21f,bR21,ToptR21,TempR21f)
  PhotomeanR21g <- quad(ymaxR21g,bR21,ToptR21,TempR21g)
  PhotomeanR21h <- quad(ymaxR21h,bR21,ToptR21,TempR21h)
  PhotomeanR26a <- quad(ymaxR26a,bR26,ToptR26,TempR26a)
  PhotomeanR26b <- quad(ymaxR26b,bR26,ToptR26,TempR26b)
  PhotomeanR26c <- quad(ymaxR26c,bR26,ToptR26,TempR26c)
  PhotomeanR26d <- quad(ymaxR26d,bR26,ToptR26,TempR26d)
  PhotomeanR26e <- quad(ymaxR26e,bR26,ToptR26,TempR26e)
  PhotomeanR26f <- quad(ymaxR26f,bR26,ToptR26,TempR26f)
  PhotomeanR31a <- quad(ymaxR31a,bR31,ToptR31,TempR31a)
  PhotomeanR31b <- quad(ymaxR31b,bR31,ToptR31,TempR31b)
  PhotomeanR31c <- quad(ymaxR31c,bR31,ToptR31,TempR31c)
  PhotomeanR31d <- quad(ymaxR31d,bR31,ToptR31,TempR31d)
  PhotomeanR31e <- quad(ymaxR31e,bR31,ToptR31,TempR31e)
  PhotomeanR31f <- quad(ymaxR31f,bR31,ToptR31,TempR31f)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for the normal (equation 2)
Photo_norm_normNLL_all <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                  ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                  ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                  ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                  ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                  ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                  ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                  ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                  ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                  ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                  ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                  ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                  ToptM21,ToptM26,ToptM31,
                                  ToptA21,ToptA26,ToptA31,
                                  ToptG21,ToptG26,ToptG31,
                                  ToptR21,ToptR26,ToptR31,
                                  sM21,sM26,sM31,
                                  sA21,sA26,sA31,
                                  sG21,sG26,sG31,
                                  sR21,sR26,sR31,
                                  TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                  TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                  TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                  TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                  TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                  TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                  TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                  TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                  TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                  TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                  TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                  TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                  PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                  PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                  PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                  PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                  PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                  PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                  PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                  PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                  PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                  PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                  PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                  PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- norm(ymaxM21a,ToptM21,sM21,TempM21a)
  PhotomeanM21b <- norm(ymaxM21b,ToptM21,sM21,TempM21b)
  PhotomeanM21c <- norm(ymaxM21c,ToptM21,sM21,TempM21c)
  PhotomeanM21d <- norm(ymaxM21d,ToptM21,sM21,TempM21d)
  PhotomeanM21e <- norm(ymaxM21e,ToptM21,sM21,TempM21e)
  PhotomeanM21f <- norm(ymaxM21f,ToptM21,sM21,TempM21f)
  PhotomeanM26a <- norm(ymaxM26a,ToptM26,sM26,TempM26a)
  PhotomeanM26b <- norm(ymaxM26b,ToptM26,sM26,TempM26b)
  PhotomeanM26c <- norm(ymaxM26c,ToptM26,sM26,TempM26c)
  PhotomeanM26d <- norm(ymaxM26d,ToptM26,sM26,TempM26d)
  PhotomeanM26e <- norm(ymaxM26e,ToptM26,sM26,TempM26e)
  PhotomeanM26f <- norm(ymaxM26f,ToptM26,sM26,TempM26f)
  PhotomeanM31a <- norm(ymaxM31a,ToptM31,sM31,TempM31a)
  PhotomeanM31b <- norm(ymaxM31b,ToptM31,sM31,TempM31b)
  PhotomeanM31c <- norm(ymaxM31c,ToptM31,sM31,TempM31c)
  PhotomeanM31d <- norm(ymaxM31d,ToptM31,sM31,TempM31d)
  PhotomeanM31e <- norm(ymaxM31e,ToptM31,sM31,TempM31e)
  PhotomeanM31f <- norm(ymaxM31f,ToptM31,sM31,TempM31f)
  PhotomeanA21a <- norm(ymaxA21a,ToptA21,sA21,TempA21a)
  PhotomeanA21b <- norm(ymaxA21b,ToptA21,sA21,TempA21b)
  PhotomeanA21c <- norm(ymaxA21c,ToptA21,sA21,TempA21c)
  PhotomeanA21d <- norm(ymaxA21d,ToptA21,sA21,TempA21d)
  PhotomeanA21e <- norm(ymaxA21e,ToptA21,sA21,TempA21e)
  PhotomeanA21f <- norm(ymaxA21f,ToptA21,sA21,TempA21f)
  PhotomeanA26a <- norm(ymaxA26a,ToptA26,sA26,TempA26a)
  PhotomeanA26b <- norm(ymaxA26b,ToptA26,sA26,TempA26b)
  PhotomeanA26c <- norm(ymaxA26c,ToptA26,sA26,TempA26c)
  PhotomeanA26d <- norm(ymaxA26d,ToptA26,sA26,TempA26d)
  PhotomeanA26e <- norm(ymaxA26e,ToptA26,sA26,TempA26e)
  PhotomeanA26f <- norm(ymaxA26f,ToptA26,sA26,TempA26f)
  PhotomeanA31a <- norm(ymaxA31a,ToptA31,sA31,TempA31a)
  PhotomeanA31b <- norm(ymaxA31b,ToptA31,sA31,TempA31b)
  PhotomeanA31c <- norm(ymaxA31c,ToptA31,sA31,TempA31c)
  PhotomeanA31d <- norm(ymaxA31d,ToptA31,sA31,TempA31d)
  PhotomeanA31e <- norm(ymaxA31e,ToptA31,sA31,TempA31e)
  PhotomeanA31f <- norm(ymaxA31f,ToptA31,sA31,TempA31f)
  PhotomeanG21a <- norm(ymaxG21a,ToptG21,sG21,TempG21a)
  PhotomeanG21b <- norm(ymaxG21b,ToptG21,sG21,TempG21b)
  PhotomeanG21c <- norm(ymaxG21c,ToptG21,sG21,TempG21c)
  PhotomeanG21d <- norm(ymaxG21d,ToptG21,sG21,TempG21d)
  PhotomeanG21e <- norm(ymaxG21e,ToptG21,sG21,TempG21e)
  PhotomeanG21f <- norm(ymaxG21f,ToptG21,sG21,TempG21f)
  PhotomeanG26a <- norm(ymaxG26a,ToptG26,sG26,TempG26a)
  PhotomeanG26b <- norm(ymaxG26b,ToptG26,sG26,TempG26b)
  PhotomeanG26c <- norm(ymaxG26c,ToptG26,sG26,TempG26c)
  PhotomeanG26d <- norm(ymaxG26d,ToptG26,sG26,TempG26d)
  PhotomeanG26e <- norm(ymaxG26e,ToptG26,sG26,TempG26e)
  PhotomeanG26f <- norm(ymaxG26f,ToptG26,sG26,TempG26f)
  PhotomeanG26g <- norm(ymaxG26g,ToptG26,sG26,TempG26g)
  PhotomeanG31a <- norm(ymaxG31a,ToptG31,sG31,TempG31a)
  PhotomeanG31b <- norm(ymaxG31b,ToptG31,sG31,TempG31b)
  PhotomeanG31c <- norm(ymaxG31c,ToptG31,sG31,TempG31c)
  PhotomeanG31d <- norm(ymaxG31d,ToptG31,sG31,TempG31d)
  PhotomeanG31e <- norm(ymaxG31e,ToptG31,sG31,TempG31e)
  PhotomeanG31f <- norm(ymaxG31f,ToptG31,sG31,TempG31f)
  PhotomeanR21a <- norm(ymaxR21a,ToptR21,sR21,TempR21a)
  PhotomeanR21b <- norm(ymaxR21b,ToptR21,sR21,TempR21b)
  PhotomeanR21c <- norm(ymaxR21c,ToptR21,sR21,TempR21c)
  PhotomeanR21d <- norm(ymaxR21d,ToptR21,sR21,TempR21d)
  PhotomeanR21e <- norm(ymaxR21e,ToptR21,sR21,TempR21e)
  PhotomeanR21f <- norm(ymaxR21f,ToptR21,sR21,TempR21f)
  PhotomeanR21g <- norm(ymaxR21g,ToptR21,sR21,TempR21g)
  PhotomeanR21h <- norm(ymaxR21h,ToptR21,sR21,TempR21h)
  PhotomeanR26a <- norm(ymaxR26a,ToptR26,sR26,TempR26a)
  PhotomeanR26b <- norm(ymaxR26b,ToptR26,sR26,TempR26b)
  PhotomeanR26c <- norm(ymaxR26c,ToptR26,sR26,TempR26c)
  PhotomeanR26d <- norm(ymaxR26d,ToptR26,sR26,TempR26d)
  PhotomeanR26e <- norm(ymaxR26e,ToptR26,sR26,TempR26e)
  PhotomeanR26f <- norm(ymaxR26f,ToptR26,sR26,TempR26f)
  PhotomeanR31a <- norm(ymaxR31a,ToptR31,sR31,TempR31a)
  PhotomeanR31b <- norm(ymaxR31b,ToptR31,sR31,TempR31b)
  PhotomeanR31c <- norm(ymaxR31c,ToptR31,sR31,TempR31c)
  PhotomeanR31d <- norm(ymaxR31d,ToptR31,sR31,TempR31d)
  PhotomeanR31e <- norm(ymaxR31e,ToptR31,sR31,TempR31e)
  PhotomeanR31f <- norm(ymaxR31f,ToptR31,sR31,TempR31f)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

####
#Now use maximum likelihood to fit functions
####

###
#A275 first
###

#For the modified beta (equation 5)
fit_Photo_beta_all_A275 <- mle2(Photo_beta_normNLL_all,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=5,ymaxM21e=7,ymaxM21f=9,
                                                           ymaxM26a=8,ymaxM26b=6,ymaxM26c=5,ymaxM26d=5,ymaxM26e=10,ymaxM26f=6,
                                                           ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=5,ymaxM31e=16,ymaxM31f=12,
                                                           ymaxA21a=14,ymaxA21b=14,ymaxA21c=12,ymaxA21d=13,ymaxA21e=13,ymaxA21f=18,
                                                           ymaxA26a=12,ymaxA26b=11,ymaxA26c=7,ymaxA26d=6,ymaxA26e=5,ymaxA26f=5,
                                                           ymaxA31a=6,ymaxA31b=7,ymaxA31c=11,ymaxA31d=12,ymaxA31e=13,ymaxA31f=13,
                                                           ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                           ymaxG26a=9,ymaxG26b=8,ymaxG26c=6,ymaxG26d=9,ymaxG26e=5,ymaxG26f=8,ymaxG26g=6,
                                                           ymaxG31a=10,ymaxG31b=15,ymaxG31c=5,ymaxG31d=4,ymaxG31e=6,ymaxG31f=5,
                                                           ymaxR21a=7,ymaxR21b=10,ymaxR21c=8,ymaxR21d=5,ymaxR21e=8,ymaxR21f=5,ymaxR21g=7,ymaxR21h=12,
                                                           ymaxR26a=9,ymaxR26b=8,ymaxR26c=9,ymaxR26d=12,ymaxR26e=6,ymaxR26f=7,
                                                           ymaxR31a=15,ymaxR31b=13,ymaxR31c=7,ymaxR31d=6,ymaxR31e=7,ymaxR31f=10,
                                                           TmaxM21=40,TmaxM26=40,TmaxM31=40,
                                                           TmaxA21=42,TmaxA26=40,TmaxA31=40,
                                                           TmaxG21=43,TmaxG26=40,TmaxG31=43,
                                                           TmaxR21=40,TmaxR26=40,TmaxR31=40,
                                                           TminM21=5,TminM26=5,TminM31=5,
                                                           TminA21=1,TminA26=5,TminA31=5,
                                                           TminG21=11,TminG26=5,TminG31=7,
                                                           TminR21=5,TminR26=5,TminR31=5,
                                                           ToptM21=27,ToptM26=27,ToptM31=27,
                                                           ToptA21=25,ToptA26=27,ToptA31=27,
                                                           ToptG21=25,ToptG26=27,ToptG31=32,
                                                           ToptR21=27,ToptR26=27,ToptR31=27),
                          data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                    TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                    TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                    TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                    TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                    TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                    TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                    TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                    TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                    TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                    TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                    TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                    PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                    PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                    PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                    PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                    PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                    PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                    PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                    PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                    PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                    PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                    PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                    PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                          control=list(maxit=20000))
summary(fit_Photo_beta_all_A275)

#For the peaked Arrhenius (equation 4)
fit_Photo_peak.ar_all_A275 <- mle2(Photo_peak.ar_normNLL_all,start=list(sdPhoto=-1,k25M21a=5.5,k25M21b=5.1,k25M21c=6,k25M21d=5.6,k25M21e=7.6,k25M21f=9.7,
                                                                        k25M26a=8.9,k25M26b=6.6,k25M26c=5.7,k25M26d=5.6,k25M26e=11.1,k25M26f=6.2,
                                                                        k25M31a=7.9*.9,k25M31b=5.8*.9,k25M31c=5.9*.9,k25M31d=3.4*.9,k25M31e=16.6*.9,k25M31f=12.5*.9,
                                                                        k25A21a=15.4,k25A21b=14.3,k25A21c=12.8,k25A21d=13.9,k25A21e=13.5,k25A21f=19.5,
                                                                        k25A26a=12.3,k25A26b=11.3,k25A26c=7.7,k25A26d=6.5,k25A26e=3.4,k25A26f=3.3,
                                                                        k25A31a=6.5*.92,k25A31b=7.8*.92,k25A31c=11.4*.92,k25A31d=12.2*.92,k25A31e=13.3*.92,k25A31f=13.1*.92,
                                                                        k25G21a=16.6,k25G21b=8,k25G21c=11.2,k25G21d=16.1,k25G21e=8.3,k25G21f=10.9,
                                                                        k25G26a=10.4,k25G26b=8.3,k25G26c=6.9,k25G26d=10,k25G26e=5.8,k25G26f=8.3,k25G26g=6.2,
                                                                        k25G31a=10.7*.74,k25G31b=15.8*.74,k25G31c=5.3*.74,k25G31d=3.7*.74,k25G31e=6.7*.74,k25G31f=5.5*.74,
                                                                        k25R21a=7.3,k25R21b=10.5,k25R21c=8.8,k25R21d=5,k25R21e=8.2,k25R21f=5.3,k25R21g=7.5,k25R21h=13.3,
                                                                        k25R26a=9.6*.87,k25R26b=8.6*.87,k25R26c=9.7*.87,k25R26d=13.7*.87,k25R26e=6.2*.87,k25R26f=7.2*.87,
                                                                        k25R31a=15.9*.89,k25R31b=13.8*.89,k25R31c=7*.89,k25R31d=6*.89,k25R31e=7.7*.89,k25R31f=10.7*.89,
                                                                        HdM21=200,HdM26=200,HdM31=200,
                                                                        HdA21=200,HdA26=200,HdA31=200,
                                                                        HdG21=200,HdG26=200,HdG31=200,
                                                                        HdR21=200,HdR26=200,HdR31=200,
                                                                        EaM21=45,EaM26=45,EaM31=45,
                                                                        EaA21=45,EaA26=45,EaA31=45,
                                                                        EaG21=45,EaG26=45,EaG31=60,
                                                                        EaR21=45,EaR26=45,EaR31=45,
                                                                        ToptM21=300,ToptM26=300,ToptM31=300,
                                                                        ToptA21=300,ToptA26=300,ToptA31=300,
                                                                        ToptG21=300,ToptG26=300,ToptG31=304,
                                                                        ToptR21=300,ToptR26=300,ToptR31=300),
                                   data=list(TkM21a=M21.A275.dat$TsK[1:3],TkM21b=M21.A275.dat$TsK[4:8],TkM21c=M21.A275.dat$TsK[9:11],TkM21d=M21.A275.dat$TsK[12:16],TkM21e=M21.A275.dat$TsK[17:19],TkM21f=M21.A275.dat$TsK[20:24],
                                             TkM26a=M26.A275.dat$TsK[1:4],TkM26b=M26.A275.dat$TsK[5:8],TkM26c=M26.A275.dat$TsK[9:12],TkM26d=M26.A275.dat$TsK[13:15],TkM26e=M26.A275.dat$TsK[16:19],TkM26f=M26.A275.dat$TsK[20:23],
                                             TkM31a=M31.A275.dat$TsK[1:5],TkM31b=M31.A275.dat$TsK[6:8],TkM31c=M31.A275.dat$TsK[9:12],TkM31d=M31.A275.dat$TsK[13:15],TkM31e=M31.A275.dat$TsK[16:20],TkM31f=M31.A275.dat$TsK[21:23],
                                             TkA21a=A21.A275.dat$TsK[1:3],TkA21b=A21.A275.dat$TsK[4:8],TkA21c=A21.A275.dat$TsK[9:11],TkA21d=A21.A275.dat$TsK[12:16],TkA21e=A21.A275.dat$TsK[17:19],TkA21f=A21.A275.dat$TsK[20:24],
                                             TkA26a=A26.A275.dat$TsK[1:4],TkA26b=A26.A275.dat$TsK[5:8],TkA26c=A26.A275.dat$TsK[9:12],TkA26d=A26.A275.dat$TsK[13:15],TkA26e=A26.A275.dat$TsK[16:19],TkA26f=A26.A275.dat$TsK[20:22],
                                             TkA31a=A31.A275.dat$TsK[1:5],TkA31b=A31.A275.dat$TsK[6:7],TkA31c=A31.A275.dat$TsK[8:12],TkA31d=A31.A275.dat$TsK[13:15],TkA31e=A31.A275.dat$TsK[16:20],TkA31f=A31.A275.dat$TsK[21:23],
                                             TkG21a=G21.A275.dat$TsK[1:3],TkG21b=G21.A275.dat$TsK[4:8],TkG21c=G21.A275.dat$TsK[9:11],TkG21d=G21.A275.dat$TsK[12:16],TkG21e=G21.A275.dat$TsK[17:19],TkG21f=G21.A275.dat$TsK[20:24],
                                             TkG26a=G26.A275.dat$TsK[1:4],TkG26b=G26.A275.dat$TsK[5:8],TkG26c=G26.A275.dat$TsK[9:12],TkG26d=G26.A275.dat$TsK[13:16],TkG26e=G26.A275.dat$TsK[17:20],TkG26f=G26.A275.dat$TsK[21:24],TkG26g=G26.A275.dat$TsK[25:27],
                                             TkG31a=G31.A275.dat$TsK[1:4],TkG31b=G31.A275.dat$TsK[5:7],TkG31c=G31.A275.dat$TsK[8:12],TkG31d=G31.A275.dat$TsK[13:14],TkG31e=G31.A275.dat$TsK[15:20],TkG31f=G31.A275.dat$TsK[21:23],
                                             TkR21a=R21.A275.dat$TsK[1:2],TkR21b=R21.A275.dat$TsK[3:6],TkR21c=R21.A275.dat$TsK[7:9],TkR21d=R21.A275.dat$TsK[10:14],TkR21e=R21.A275.dat$TsK[15:17],TkR21f=R21.A275.dat$TsK[18:22],TkR21g=R21.A275.dat$TsK[23:25],TkR21h=R21.A275.dat$TsK[26:30],
                                             TkR26a=R26.A275.dat$TsK[1:4],TkR26b=R26.A275.dat$TsK[5:8],TkR26c=R26.A275.dat$TsK[9:13],TkR26d=R26.A275.dat$TsK[14:17],TkR26e=R26.A275.dat$TsK[18:21],TkR26f=R26.A275.dat$TsK[22:25],
                                             TkR31a=R31.A275.dat$TsK[1:5],TkR31b=R31.A275.dat$TsK[6:8],TkR31c=R31.A275.dat$TsK[9:13],TkR31d=R31.A275.dat$TsK[14:16],TkR31e=R31.A275.dat$TsK[17:21],TkR31f=R31.A275.dat$TsK[22:24],
                                             PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                             PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                             PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                             PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                             PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                             PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                             PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                             PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                             PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                             PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                             PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                             PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                   control=list(maxit=20000))
summary(fit_Photo_peak.ar_all_A275)

#For the quadratic (equation 3)
fit_Photo_quad_all_A275 <- mle2(Photo_quad_normNLL_all,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=5,ymaxM21e=7,ymaxM21f=9,
                                                           ymaxM26a=8,ymaxM26b=6,ymaxM26c=5,ymaxM26d=5,ymaxM26e=10,ymaxM26f=6,
                                                           ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=5,ymaxM31e=16,ymaxM31f=12,
                                                           ymaxA21a=14,ymaxA21b=14,ymaxA21c=12,ymaxA21d=13,ymaxA21e=13,ymaxA21f=18,
                                                           ymaxA26a=12,ymaxA26b=11,ymaxA26c=7,ymaxA26d=6,ymaxA26e=5,ymaxA26f=5,
                                                           ymaxA31a=6,ymaxA31b=7,ymaxA31c=11,ymaxA31d=12,ymaxA31e=13,ymaxA31f=13,
                                                           ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                           ymaxG26a=9,ymaxG26b=8,ymaxG26c=6,ymaxG26d=9,ymaxG26e=5,ymaxG26f=8,ymaxG26g=6,
                                                           ymaxG31a=10,ymaxG31b=15,ymaxG31c=5,ymaxG31d=5,ymaxG31e=6,ymaxG31f=5,
                                                           ymaxR21a=7,ymaxR21b=10,ymaxR21c=8,ymaxR21d=5,ymaxR21e=8,ymaxR21f=5,ymaxR21g=7,ymaxR21h=12,
                                                           ymaxR26a=9,ymaxR26b=8,ymaxR26c=9,ymaxR26d=12,ymaxR26e=6,ymaxR26f=7,
                                                           ymaxR31a=15,ymaxR31b=13,ymaxR31c=7,ymaxR31d=6,ymaxR31e=7,ymaxR31f=10,
                                                           bM21=0.04,bM26=0.04,bM31=0.04,
                                                           bA21=0.04,bA26=0.04,bA31=0.04,
                                                           bG21=0.04,bG26=0.04,bG31=0.04,
                                                           bR21=0.04,bR26=0.04,bR31=0.04,
                                                           ToptM21=27,ToptM26=27,ToptM31=27,
                                                           ToptA21=27,ToptA26=27,ToptA31=27,
                                                           ToptG21=27,ToptG26=27,ToptG31=27,
                                                           ToptR21=27,ToptR26=27,ToptR31=27),
                          data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                    TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                    TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                    TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                    TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                    TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                    TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                    TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                    TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                    TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                    TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                    TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                    PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                    PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                    PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                    PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                    PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                    PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                    PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                    PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                    PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                    PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                    PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                    PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                          control=list(maxit=20000))
summary(fit_Photo_quad_all_A275)

#For the normal (equation 2)
fit_Photo_norm_all_A275 <- mle2(Photo_norm_normNLL_all,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=5,ymaxM21e=7,ymaxM21f=9,
                                                           ymaxM26a=8,ymaxM26b=6,ymaxM26c=5,ymaxM26d=5,ymaxM26e=10,ymaxM26f=6,
                                                           ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=5,ymaxM31e=16,ymaxM31f=12,
                                                           ymaxA21a=14,ymaxA21b=14,ymaxA21c=12,ymaxA21d=13,ymaxA21e=13,ymaxA21f=18,
                                                           ymaxA26a=12,ymaxA26b=11,ymaxA26c=7,ymaxA26d=6,ymaxA26e=5,ymaxA26f=5,
                                                           ymaxA31a=6,ymaxA31b=7,ymaxA31c=11,ymaxA31d=12,ymaxA31e=13,ymaxA31f=13,
                                                           ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                           ymaxG26a=9,ymaxG26b=8,ymaxG26c=6,ymaxG26d=9,ymaxG26e=5,ymaxG26f=8,ymaxG26g=6,
                                                           ymaxG31a=10,ymaxG31b=15,ymaxG31c=5,ymaxG31d=5,ymaxG31e=6,ymaxG31f=5,
                                                           ymaxR21a=7,ymaxR21b=10,ymaxR21c=8,ymaxR21d=5,ymaxR21e=8,ymaxR21f=5,ymaxR21g=7,ymaxR21h=12,
                                                           ymaxR26a=9,ymaxR26b=8,ymaxR26c=9,ymaxR26d=12,ymaxR26e=6,ymaxR26f=7,
                                                           ymaxR31a=15,ymaxR31b=13,ymaxR31c=7,ymaxR31d=6,ymaxR31e=7,ymaxR31f=10,
                                                           sM21=10,sM26=10,sM31=10,
                                                           sA21=10,sA26=10,sA31=10,
                                                           sG21=10,sG26=10,sG31=10,
                                                           sR21=10,sR26=10,sR31=10,
                                                           ToptM21=27,ToptM26=27,ToptM31=27,
                                                           ToptA21=27,ToptA26=27,ToptA31=27,
                                                           ToptG21=27,ToptG26=27,ToptG31=27,
                                                           ToptR21=27,ToptR26=27,ToptR31=27),
                          data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                    TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                    TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                    TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                    TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                    TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                    TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                    TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                    TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                    TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                    TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                    TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                    PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                    PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                    PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                    PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                    PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                    PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                    PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                    PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                    PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                    PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                    PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                    PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                          control=list(maxit=20000))
summary(fit_Photo_norm_all_A275)

#Calculate delta AIC values
AICctab(fit_Photo_beta_all_A275,fit_Photo_peak.ar_all_A275,fit_Photo_quad_all_A275,fit_Photo_norm_all_A275,nobs=292)

#Calculate sample size
length(c(M21.A275.dat$A275,M26.A275.dat$A275,M31.A275.dat$A275,A21.A275.dat$A275,A26.A275.dat$A275,A31.A275.dat$A275,
         G21.A275.dat$A275,G26.A275.dat$A275,G31.A275.dat$A275,R21.A275.dat$A275,R26.A275.dat$A275,R31.A275.dat$A275))

###
#Asat
###

#For the modified beta (equation 5)
fit_Photo_beta_all_Asat <- mle2(Photo_beta_normNLL_all,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                TmaxM21=43,TmaxM26=41,TmaxM31=42,
                                                                TmaxA21=43,TmaxA26=43,TmaxA31=43,
                                                                TmaxG21=45,TmaxG26=41,TmaxG31=42,
                                                                TmaxR21=45,TmaxR26=42,TmaxR31=41,
                                                                TminM21=4,TminM26=3,TminM31=-11,
                                                                TminA21=-1,TminA26=-3,TminA31=-13,
                                                                TminG21=8,TminG26=8,TminG31=7,
                                                                TminR21=2,TminR26=7,TminR31=5,
                                                                ToptM21=28,ToptM26=29,ToptM31=31,
                                                                ToptA21=26,ToptA26=30,ToptA31=31,
                                                                ToptG21=24,ToptG26=25,ToptG31=31,
                                                                ToptR21=28,ToptR26=26,ToptR31=29),
                               data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                         TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                         TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                         TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                         TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                         TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                         TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                         TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                         TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                         TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                         TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                         TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                         PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                         PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                         PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                         PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                         PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                         PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                         PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                         PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                         PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                         PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                         PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                         PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                               control=list(maxit=20000))
summary(fit_Photo_beta_all_Asat)

#For the peaked Arrhenius (equation 4)
fit_Photo_peak.ar_all_Asat <- mle2(Photo_peak.ar_normNLL_all,start=list(sdPhoto=-1,k25M21a=6.5*.95,k25M21b=6*.95,k25M21c=6.8*.95,k25M21d=6.4*.95,k25M21e=8.8*.95,k25M21f=10.7*.95,
                                                                        k25M26a=10.4*.88,k25M26b=7.8*.88,k25M26c=7.1*.88,k25M26d=6.3*.88,k25M26e=12.3*.88,k25M26f=7.2*.88,
                                                                        k25M31a=8.1*.81,k25M31b=6.1*.81,k25M31c=7.5*.81,k25M31d=4.6*.81,k25M31e=18.4*.81,k25M31f=14.6*.81,
                                                                        k25A21a=15.6*.98,k25A21b=14.6*.98,k25A21c=13.4*.98,k25A21d=13.7*.98,k25A21e=13.8*.98,k25A21f=21.1*.98,
                                                                        k25A26a=15.1*.86,k25A26b=13.9*.86,k25A26c=9.8*.86,k25A26d=7.5*.86,k25A26e=4.6*.86,k25A26f=4.2*.86,
                                                                        k25A31a=7.8*.84,k25A31b=8.7*.84,k25A31c=13.2*.84,k25A31d=13.6*.84,k25A31e=15.2*.84,k25A31f=14.3*.84,
                                                                        k25G21a=16.6,k25G21b=7.8,k25G21c=10.2,k25G21d=13.9,k25G21e=6.3,k25G21f=9.5,
                                                                        k25G26a=8.3,k25G26b=7.9,k25G26c=5.7,k25G26d=9.1,k25G26e=3.9,k25G26f=5.2,k25G26g=5.5,
                                                                        k25G31a=11*.72,k25G31b=14*.72,k25G31c=6*.72,k25G31d=4*.72,k25G31e=7*.72,k25G31f=7*.72,
                                                                        k25R21a=6.4*.96,k25R21b=7.2*.96,k25R21c=8.4*.96,k25R21d=4.7*.96,k25R21e=7.6*.96,k25R21f=6.1*.96,k25R21g=7.5*.96,k25R21h=13.9*.96,
                                                                        k25R26a=8.3*.97,k25R26b=7.7*.97,k25R26c=8.9*.97,k25R26d=12.1*.97,k25R26e=6.6*.97,k25R26f=8.1*.97,
                                                                        k25R31a=15*.88,k25R31b=11.6*.88,k25R31c=8.3*.88,k25R31d=7.1*.88,k25R31e=7.9*.88,k25R31f=9.6*.88,
                                                                        HdM21=183,HdM26=230,HdM31=213,
                                                                        HdA21=189,HdA26=213,HdA31=199,
                                                                        HdG21=172,HdG26=214,HdG31=255,
                                                                        HdR21=156,HdR26=209,HdR31=229,
                                                                        EaM21=70,EaM26=75,EaM31=63,
                                                                        EaA21=43,EaA26=60,EaA31=51,
                                                                        EaG21=111,EaG26=82,EaG31=86,
                                                                        EaR21=63,EaR26=76,EaR31=78,
                                                                        ToptM21=300,ToptM26=300,ToptM31=302,
                                                                        ToptA21=299,ToptA26=301,ToptA31=302,
                                                                        ToptG21=298,ToptG26=300,ToptG31=304,
                                                                        ToptR21=300,ToptR26=300,ToptR31=302),
                                   data=list(TkM21a=M21.Asat.dat$TsK[1:3],TkM21b=M21.Asat.dat$TsK[4:8],TkM21c=M21.Asat.dat$TsK[9:11],TkM21d=M21.Asat.dat$TsK[12:16],TkM21e=M21.Asat.dat$TsK[17:19],TkM21f=M21.Asat.dat$TsK[20:24],
                                             TkM26a=M26.Asat.dat$TsK[1:4],TkM26b=M26.Asat.dat$TsK[5:8],TkM26c=M26.Asat.dat$TsK[9:12],TkM26d=M26.Asat.dat$TsK[13:16],TkM26e=M26.Asat.dat$TsK[17:20],TkM26f=M26.Asat.dat$TsK[21:24],
                                             TkM31a=M31.Asat.dat$TsK[1:5],TkM31b=M31.Asat.dat$TsK[6:8],TkM31c=M31.Asat.dat$TsK[9:12],TkM31d=M31.Asat.dat$TsK[13:15],TkM31e=M31.Asat.dat$TsK[16:20],TkM31f=M31.Asat.dat$TsK[21:23],
                                             TkA21a=A21.Asat.dat$TsK[1:3],TkA21b=A21.Asat.dat$TsK[4:8],TkA21c=A21.Asat.dat$TsK[9:11],TkA21d=A21.Asat.dat$TsK[12:16],TkA21e=A21.Asat.dat$TsK[17:19],TkA21f=A21.Asat.dat$TsK[20:24],
                                             TkA26a=A26.Asat.dat$TsK[1:4],TkA26b=A26.Asat.dat$TsK[5:8],TkA26c=A26.Asat.dat$TsK[9:12],TkA26d=A26.Asat.dat$TsK[13:16],TkA26e=A26.Asat.dat$TsK[17:20],TkA26f=A26.Asat.dat$TsK[21:24],
                                             TkA31a=A31.Asat.dat$TsK[1:5],TkA31b=A31.Asat.dat$TsK[6:8],TkA31c=A31.Asat.dat$TsK[9:13],TkA31d=A31.Asat.dat$TsK[14:16],TkA31e=A31.Asat.dat$TsK[17:21],TkA31f=A31.Asat.dat$TsK[22:24],
                                             TkG21a=G21.Asat.dat$TsK[1:3],TkG21b=G21.Asat.dat$TsK[4:8],TkG21c=G21.Asat.dat$TsK[9:11],TkG21d=G21.Asat.dat$TsK[12:16],TkG21e=G21.Asat.dat$TsK[17:19],TkG21f=G21.Asat.dat$TsK[20:24],
                                             TkG26a=G26.Asat.dat$TsK[1:4],TkG26b=G26.Asat.dat$TsK[5:8],TkG26c=G26.Asat.dat$TsK[9:12],TkG26d=G26.Asat.dat$TsK[13:16],TkG26e=G26.Asat.dat$TsK[17:20],TkG26f=G26.Asat.dat$TsK[21:24],TkG26g=G26.Asat.dat$TsK[25:28],
                                             TkG31a=G31.Asat.dat$TsK[1:4],TkG31b=G31.Asat.dat$TsK[5:7],TkG31c=G31.Asat.dat$TsK[8:12],TkG31d=G31.Asat.dat$TsK[13:15],TkG31e=G31.Asat.dat$TsK[16:21],TkG31f=G31.Asat.dat$TsK[22:24],
                                             TkR21a=R21.Asat.dat$TsK[1:3],TkR21b=R21.Asat.dat$TsK[4:8],TkR21c=R21.Asat.dat$TsK[9:11],TkR21d=R21.Asat.dat$TsK[12:16],TkR21e=R21.Asat.dat$TsK[17:19],TkR21f=R21.Asat.dat$TsK[20:24],TkR21g=R21.Asat.dat$TsK[25:27],TkR21h=R21.Asat.dat$TsK[28:32],
                                             TkR26a=R26.Asat.dat$TsK[1:4],TkR26b=R26.Asat.dat$TsK[5:8],TkR26c=R26.Asat.dat$TsK[9:13],TkR26d=R26.Asat.dat$TsK[14:17],TkR26e=R26.Asat.dat$TsK[18:21],TkR26f=R26.Asat.dat$TsK[22:25],
                                             TkR31a=R31.Asat.dat$TsK[1:5],TkR31b=R31.Asat.dat$TsK[6:8],TkR31c=R31.Asat.dat$TsK[9:13],TkR31d=R31.Asat.dat$TsK[14:16],TkR31e=R31.Asat.dat$TsK[17:21],TkR31f=R31.Asat.dat$TsK[22:24],
                                             PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                             PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                             PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                             PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                             PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                             PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                             PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                             PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                             PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                             PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                             PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                             PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                   control=list(maxit=20000))
summary(fit_Photo_peak.ar_all_Asat)

#For the quadratic (equation 3)
fit_Photo_quad_all_Asat <- mle2(Photo_quad_normNLL_all,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                bM21=0.04,bM26=0.04,bM31=0.04,
                                                                bA21=0.04,bA26=0.04,bA31=0.04,
                                                                bG21=0.04,bG26=0.04,bG31=0.04,
                                                                bR21=0.04,bR26=0.04,bR31=0.04,
                                                                ToptM21=28,ToptM26=29,ToptM31=31,
                                                                ToptA21=26,ToptA26=30,ToptA31=31,
                                                                ToptG21=25,ToptG26=25,ToptG31=31,
                                                                ToptR21=28,ToptR26=26,ToptR31=29),
                               data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                         TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                         TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                         TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                         TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                         TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                         TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                         TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                         TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                         TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                         TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                         TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                         PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                         PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                         PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                         PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                         PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                         PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                         PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                         PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                         PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                         PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                         PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                         PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                               control=list(maxit=20000))
summary(fit_Photo_quad_all_Asat)

#For the normal (equation 2)
fit_Photo_norm_all_Asat <- mle2(Photo_norm_normNLL_all,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                sM21=10,sM26=10,sM31=10,
                                                                sA21=10,sA26=10,sA31=10,
                                                                sG21=10,sG26=10,sG31=10,
                                                                sR21=10,sR26=10,sR31=10,
                                                                ToptM21=28,ToptM26=29,ToptM31=31,
                                                                ToptA21=26,ToptA26=30,ToptA31=31,
                                                                ToptG21=25,ToptG26=25,ToptG31=31,
                                                                ToptR21=28,ToptR26=26,ToptR31=29),
                               data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                         TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                         TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                         TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                         TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                         TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                         TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                         TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                         TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                         TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                         TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                         TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                         PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                         PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                         PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                         PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                         PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                         PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                         PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                         PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                         PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                         PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                         PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                         PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                               control=list(maxit=20000))
summary(fit_Photo_norm_all_Asat)

#Calculate delta AIC values
AICctab(fit_Photo_beta_all_Asat,fit_Photo_peak.ar_all_Asat,fit_Photo_quad_all_Asat,fit_Photo_norm_all_Asat,nobs=300)

#Calculate sample size
length(c(M21.Asat.dat$A400,M26.Asat.dat$A400,M31.Asat.dat$A400,A21.Asat.dat$A400,A26.Asat.dat$A400,A31.Asat.dat$A400,
         G21.Asat.dat$A400,G26.Asat.dat$A400,G31.Asat.dat$A400,R21.Asat.dat$A400,R26.Asat.dat$A400,R31.Asat.dat$A400))

###############################################################################################################
#Model comparison for linear acclimation of modified beta function parameters (Supplementary Table 9)
###############################################################################################################

####
#Define functions
####

#Acclimation of all parameters
beta.lin.all <- function(ymax,a,b,c,d,e,f,T,Tgrow){
  y <- ymax*((e+f*Tgrow)-T)/((e+f*Tgrow)-(c+d*Tgrow))*(((T-(a+b*Tgrow))/((c+d*Tgrow)-(a+b*Tgrow)))^(((c+d*Tgrow)-(a+b*Tgrow))/((e+f*Tgrow)-(c+d*Tgrow))))
  y
}

#Acclimation of Tmin and Topt
beta.Tmin.Topt.lin <- function(ymax,a,b,c,d,Tmax,T,Tgrow){
  y <- ymax*(Tmax-T)/(Tmax-(c+d*Tgrow))*(((T-(a+b*Tgrow))/((c+d*Tgrow)-(a+b*Tgrow)))^(((c+d*Tgrow)-(a+b*Tgrow))/(Tmax-(c+d*Tgrow))))
  y
}

#Acclimation of Tmin and Tmax
beta.Tmin.Tmax.lin <- function(ymax,a,b,Topt,c,d,T,Tgrow){
  y <- ymax*((c+d*Tgrow)-T)/((c+d*Tgrow)-Topt)*(((T-(a+b*Tgrow))/(Topt-(a+b*Tgrow)))^((Topt-(a+b*Tgrow))/((c+d*Tgrow)-Topt)))
  y
}

#Acclimation of Topt and Tmax
beta.Topt.Tmax.lin <- function(ymax,Tmin,a,b,c,d,T,Tgrow){
  y <- ymax*((c+d*Tgrow)-T)/((c+d*Tgrow)-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/((c+d*Tgrow)-(a+b*Tgrow))))
  y
}

#Acclimation of Tmin
beta.Tmin.lin <- function(ymax,a,b,Topt,Tmax,T,Tgrow){
  y <- ymax*(Tmax-T)/(Tmax-Topt)*(((T-(a+b*Tgrow))/(Topt-(a+b*Tgrow)))^((Topt-(a+b*Tgrow))/(Tmax-Topt)))
  y
}

#Acclimation of Topt
beta.Topt.lin <- function(ymax,Tmin,a,b,Tmax,T,Tgrow){
  y <- ymax*(Tmax-T)/(Tmax-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/(Tmax-(a+b*Tgrow))))
  y
}

#Acclimation of Tmax
beta.Tmax.lin <- function(ymax,Tmin,Topt,a,b,T,Tgrow){
  y <- ymax*((a+b*Tgrow)-T)/((a+b*Tgrow)-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/((a+b*Tgrow)-Topt)))
  y
}

#No acclimation
beta <- function(ymax,Tmin,Topt,Tmax,T){
  y <- ymax*(Tmax-T)/(Tmax-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt)))
  y
}

####
#NLL functions
####

#NLL function for acclimation of all parameters
Photo_beta_normNLL_all_lin.all <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                           ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                           ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                           ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                           ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                           ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                           ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                           ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                           ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                           ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                           ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                           ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                           aM,
                                           aA,
                                           aG,
                                           aR,
                                           bM,
                                           bA,
                                           bG,
                                           bR,
                                           cM,
                                           cA,
                                           cG,
                                           cR,
                                           dM,
                                           dA,
                                           dG,
                                           dR,
                                           eM,
                                           eA,
                                           eG,
                                           eR,
                                           fM,
                                           fA,
                                           fG,
                                           fR,
                                           TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                           TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                           TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                           TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                           TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                           TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                           TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                           TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                           TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                           TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                           TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                           TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                           PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                           PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                           PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                           PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                           PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                           PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                           PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                           PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                           PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                           PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                           PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                           PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.lin.all(ymaxM21a,aM,bM,cM,dM,eM,fM,TempM21a,18.5)
  PhotomeanM21b <- beta.lin.all(ymaxM21b,aM,bM,cM,dM,eM,fM,TempM21b,18.5)
  PhotomeanM21c <- beta.lin.all(ymaxM21c,aM,bM,cM,dM,eM,fM,TempM21c,18.5)
  PhotomeanM21d <- beta.lin.all(ymaxM21d,aM,bM,cM,dM,eM,fM,TempM21d,18.5)
  PhotomeanM21e <- beta.lin.all(ymaxM21e,aM,bM,cM,dM,eM,fM,TempM21e,18.5)
  PhotomeanM21f <- beta.lin.all(ymaxM21f,aM,bM,cM,dM,eM,fM,TempM21f,18.5)
  PhotomeanM26a <- beta.lin.all(ymaxM26a,aM,bM,cM,dM,eM,fM,TempM26a,23.5)
  PhotomeanM26b <- beta.lin.all(ymaxM26b,aM,bM,cM,dM,eM,fM,TempM26b,23.5)
  PhotomeanM26c <- beta.lin.all(ymaxM26c,aM,bM,cM,dM,eM,fM,TempM26c,23.5)
  PhotomeanM26d <- beta.lin.all(ymaxM26d,aM,bM,cM,dM,eM,fM,TempM26d,23.5)
  PhotomeanM26e <- beta.lin.all(ymaxM26e,aM,bM,cM,dM,eM,fM,TempM26e,23.5)
  PhotomeanM26f <- beta.lin.all(ymaxM26f,aM,bM,cM,dM,eM,fM,TempM26f,23.5)
  PhotomeanM31a <- beta.lin.all(ymaxM31a,aM,bM,cM,dM,eM,fM,TempM31a,28.5)
  PhotomeanM31b <- beta.lin.all(ymaxM31b,aM,bM,cM,dM,eM,fM,TempM31b,28.5)
  PhotomeanM31c <- beta.lin.all(ymaxM31c,aM,bM,cM,dM,eM,fM,TempM31c,28.5)
  PhotomeanM31d <- beta.lin.all(ymaxM31d,aM,bM,cM,dM,eM,fM,TempM31d,28.5)
  PhotomeanM31e <- beta.lin.all(ymaxM31e,aM,bM,cM,dM,eM,fM,TempM31e,28.5)
  PhotomeanM31f <- beta.lin.all(ymaxM31f,aM,bM,cM,dM,eM,fM,TempM31f,28.5)
  PhotomeanA21a <- beta.lin.all(ymaxA21a,aA,bA,cA,dA,eA,fA,TempA21a,18.5)
  PhotomeanA21b <- beta.lin.all(ymaxA21b,aA,bA,cA,dA,eA,fA,TempA21b,18.5)
  PhotomeanA21c <- beta.lin.all(ymaxA21c,aA,bA,cA,dA,eA,fA,TempA21c,18.5)
  PhotomeanA21d <- beta.lin.all(ymaxA21d,aA,bA,cA,dA,eA,fA,TempA21d,18.5)
  PhotomeanA21e <- beta.lin.all(ymaxA21e,aA,bA,cA,dA,eA,fA,TempA21e,18.5)
  PhotomeanA21f <- beta.lin.all(ymaxA21f,aA,bA,cA,dA,eA,fA,TempA21f,18.5)
  PhotomeanA26a <- beta.lin.all(ymaxA26a,aA,bA,cA,dA,eA,fA,TempA26a,23.5)
  PhotomeanA26b <- beta.lin.all(ymaxA26b,aA,bA,cA,dA,eA,fA,TempA26b,23.5)
  PhotomeanA26c <- beta.lin.all(ymaxA26c,aA,bA,cA,dA,eA,fA,TempA26c,23.5)
  PhotomeanA26d <- beta.lin.all(ymaxA26d,aA,bA,cA,dA,eA,fA,TempA26d,23.5)
  PhotomeanA26e <- beta.lin.all(ymaxA26e,aA,bA,cA,dA,eA,fA,TempA26e,23.5)
  PhotomeanA26f <- beta.lin.all(ymaxA26f,aA,bA,cA,dA,eA,fA,TempA26f,23.5)
  PhotomeanA31a <- beta.lin.all(ymaxA31a,aA,bA,cA,dA,eA,fA,TempA31a,28.5)
  PhotomeanA31b <- beta.lin.all(ymaxA31b,aA,bA,cA,dA,eA,fA,TempA31b,28.5)
  PhotomeanA31c <- beta.lin.all(ymaxA31c,aA,bA,cA,dA,eA,fA,TempA31c,28.5)
  PhotomeanA31d <- beta.lin.all(ymaxA31d,aA,bA,cA,dA,eA,fA,TempA31d,28.5)
  PhotomeanA31e <- beta.lin.all(ymaxA31e,aA,bA,cA,dA,eA,fA,TempA31e,28.5)
  PhotomeanA31f <- beta.lin.all(ymaxA31f,aA,bA,cA,dA,eA,fA,TempA31f,28.5)
  PhotomeanG21a <- beta.lin.all(ymaxG21a,aG,bG,cG,dG,eG,fG,TempG21a,18.5)
  PhotomeanG21b <- beta.lin.all(ymaxG21b,aG,bG,cG,dG,eG,fG,TempG21b,18.5)
  PhotomeanG21c <- beta.lin.all(ymaxG21c,aG,bG,cG,dG,eG,fG,TempG21c,18.5)
  PhotomeanG21d <- beta.lin.all(ymaxG21d,aG,bG,cG,dG,eG,fG,TempG21d,18.5)
  PhotomeanG21e <- beta.lin.all(ymaxG21e,aG,bG,cG,dG,eG,fG,TempG21e,18.5)
  PhotomeanG21f <- beta.lin.all(ymaxG21f,aG,bG,cG,dG,eG,fG,TempG21f,18.5)
  PhotomeanG26a <- beta.lin.all(ymaxG26a,aG,bG,cG,dG,eG,fG,TempG26a,23.5)
  PhotomeanG26b <- beta.lin.all(ymaxG26b,aG,bG,cG,dG,eG,fG,TempG26b,23.5)
  PhotomeanG26c <- beta.lin.all(ymaxG26c,aG,bG,cG,dG,eG,fG,TempG26c,23.5)
  PhotomeanG26d <- beta.lin.all(ymaxG26d,aG,bG,cG,dG,eG,fG,TempG26d,23.5)
  PhotomeanG26e <- beta.lin.all(ymaxG26e,aG,bG,cG,dG,eG,fG,TempG26e,23.5)
  PhotomeanG26f <- beta.lin.all(ymaxG26f,aG,bG,cG,dG,eG,fG,TempG26f,23.5)
  PhotomeanG26g <- beta.lin.all(ymaxG26g,aG,bG,cG,dG,eG,fG,TempG26g,23.5)
  PhotomeanG31a <- beta.lin.all(ymaxG31a,aG,bG,cG,dG,eG,fG,TempG31a,28.5)
  PhotomeanG31b <- beta.lin.all(ymaxG31b,aG,bG,cG,dG,eG,fG,TempG31b,28.5)
  PhotomeanG31c <- beta.lin.all(ymaxG31c,aG,bG,cG,dG,eG,fG,TempG31c,28.5)
  PhotomeanG31d <- beta.lin.all(ymaxG31d,aG,bG,cG,dG,eG,fG,TempG31d,28.5)
  PhotomeanG31e <- beta.lin.all(ymaxG31e,aG,bG,cG,dG,eG,fG,TempG31e,28.5)
  PhotomeanG31f <- beta.lin.all(ymaxG31f,aG,bG,cG,dG,eG,fG,TempG31f,28.5)
  PhotomeanR21a <- beta.lin.all(ymaxR21a,aR,bR,cR,dR,eR,fR,TempR21a,18.5)
  PhotomeanR21b <- beta.lin.all(ymaxR21b,aR,bR,cR,dR,eR,fR,TempR21b,18.5)
  PhotomeanR21c <- beta.lin.all(ymaxR21c,aR,bR,cR,dR,eR,fR,TempR21c,18.5)
  PhotomeanR21d <- beta.lin.all(ymaxR21d,aR,bR,cR,dR,eR,fR,TempR21d,18.5)
  PhotomeanR21e <- beta.lin.all(ymaxR21e,aR,bR,cR,dR,eR,fR,TempR21e,18.5)
  PhotomeanR21f <- beta.lin.all(ymaxR21f,aR,bR,cR,dR,eR,fR,TempR21f,18.5)
  PhotomeanR21g <- beta.lin.all(ymaxR21g,aR,bR,cR,dR,eR,fR,TempR21g,18.5)
  PhotomeanR21h <- beta.lin.all(ymaxR21h,aR,bR,cR,dR,eR,fR,TempR21h,18.5)
  PhotomeanR26a <- beta.lin.all(ymaxR26a,aR,bR,cR,dR,eR,fR,TempR26a,23.5)
  PhotomeanR26b <- beta.lin.all(ymaxR26b,aR,bR,cR,dR,eR,fR,TempR26b,23.5)
  PhotomeanR26c <- beta.lin.all(ymaxR26c,aR,bR,cR,dR,eR,fR,TempR26c,23.5)
  PhotomeanR26d <- beta.lin.all(ymaxR26d,aR,bR,cR,dR,eR,fR,TempR26d,23.5)
  PhotomeanR26e <- beta.lin.all(ymaxR26e,aR,bR,cR,dR,eR,fR,TempR26e,23.5)
  PhotomeanR26f <- beta.lin.all(ymaxR26f,aR,bR,cR,dR,eR,fR,TempR26f,23.5)
  PhotomeanR31a <- beta.lin.all(ymaxR31a,aR,bR,cR,dR,eR,fR,TempR31a,28.5)
  PhotomeanR31b <- beta.lin.all(ymaxR31b,aR,bR,cR,dR,eR,fR,TempR31b,28.5)
  PhotomeanR31c <- beta.lin.all(ymaxR31c,aR,bR,cR,dR,eR,fR,TempR31c,28.5)
  PhotomeanR31d <- beta.lin.all(ymaxR31d,aR,bR,cR,dR,eR,fR,TempR31d,28.5)
  PhotomeanR31e <- beta.lin.all(ymaxR31e,aR,bR,cR,dR,eR,fR,TempR31e,28.5)
  PhotomeanR31f <- beta.lin.all(ymaxR31f,aR,bR,cR,dR,eR,fR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Tmin and Topt
Photo_beta_normNLL_all_Tmin.Topt.lin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                                 ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                                 ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                                 ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                                 ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                                 ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                                 ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                                 ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                                 ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                                 ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                                 ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                                 ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                                 TmaxM,
                                                 TmaxA,
                                                 TmaxG,
                                                 TmaxR,
                                                 aM,
                                                 aA,
                                                 aG,
                                                 aR,
                                                 bM,
                                                 bA,
                                                 bG,
                                                 bR,
                                                 cM,
                                                 cA,
                                                 cG,
                                                 cR,
                                                 dM,
                                                 dA,
                                                 dG,
                                                 dR,
                                                 TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                                 TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                                 TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                                 TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                                 TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                                 TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                                 TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                                 TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                                 TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                                 TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                                 TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                                 TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                                 PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                                 PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                                 PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                                 PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                                 PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                                 PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                                 PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                                 PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                                 PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                                 PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                                 PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                                 PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.Tmin.Topt.lin(ymaxM21a,aM,bM,cM,dM,TmaxM,TempM21a,18.5)
  PhotomeanM21b <- beta.Tmin.Topt.lin(ymaxM21b,aM,bM,cM,dM,TmaxM,TempM21b,18.5)
  PhotomeanM21c <- beta.Tmin.Topt.lin(ymaxM21c,aM,bM,cM,dM,TmaxM,TempM21c,18.5)
  PhotomeanM21d <- beta.Tmin.Topt.lin(ymaxM21d,aM,bM,cM,dM,TmaxM,TempM21d,18.5)
  PhotomeanM21e <- beta.Tmin.Topt.lin(ymaxM21e,aM,bM,cM,dM,TmaxM,TempM21e,18.5)
  PhotomeanM21f <- beta.Tmin.Topt.lin(ymaxM21f,aM,bM,cM,dM,TmaxM,TempM21f,18.5)
  PhotomeanM26a <- beta.Tmin.Topt.lin(ymaxM26a,aM,bM,cM,dM,TmaxM,TempM26a,23.5)
  PhotomeanM26b <- beta.Tmin.Topt.lin(ymaxM26b,aM,bM,cM,dM,TmaxM,TempM26b,23.5)
  PhotomeanM26c <- beta.Tmin.Topt.lin(ymaxM26c,aM,bM,cM,dM,TmaxM,TempM26c,23.5)
  PhotomeanM26d <- beta.Tmin.Topt.lin(ymaxM26d,aM,bM,cM,dM,TmaxM,TempM26d,23.5)
  PhotomeanM26e <- beta.Tmin.Topt.lin(ymaxM26e,aM,bM,cM,dM,TmaxM,TempM26e,23.5)
  PhotomeanM26f <- beta.Tmin.Topt.lin(ymaxM26f,aM,bM,cM,dM,TmaxM,TempM26f,23.5)
  PhotomeanM31a <- beta.Tmin.Topt.lin(ymaxM31a,aM,bM,cM,dM,TmaxM,TempM31a,28.5)
  PhotomeanM31b <- beta.Tmin.Topt.lin(ymaxM31b,aM,bM,cM,dM,TmaxM,TempM31b,28.5)
  PhotomeanM31c <- beta.Tmin.Topt.lin(ymaxM31c,aM,bM,cM,dM,TmaxM,TempM31c,28.5)
  PhotomeanM31d <- beta.Tmin.Topt.lin(ymaxM31d,aM,bM,cM,dM,TmaxM,TempM31d,28.5)
  PhotomeanM31e <- beta.Tmin.Topt.lin(ymaxM31e,aM,bM,cM,dM,TmaxM,TempM31e,28.5)
  PhotomeanM31f <- beta.Tmin.Topt.lin(ymaxM31f,aM,bM,cM,dM,TmaxM,TempM31f,28.5)
  PhotomeanA21a <- beta.Tmin.Topt.lin(ymaxA21a,aA,bA,cA,dA,TmaxA,TempA21a,18.5)
  PhotomeanA21b <- beta.Tmin.Topt.lin(ymaxA21b,aA,bA,cA,dA,TmaxA,TempA21b,18.5)
  PhotomeanA21c <- beta.Tmin.Topt.lin(ymaxA21c,aA,bA,cA,dA,TmaxA,TempA21c,18.5)
  PhotomeanA21d <- beta.Tmin.Topt.lin(ymaxA21d,aA,bA,cA,dA,TmaxA,TempA21d,18.5)
  PhotomeanA21e <- beta.Tmin.Topt.lin(ymaxA21e,aA,bA,cA,dA,TmaxA,TempA21e,18.5)
  PhotomeanA21f <- beta.Tmin.Topt.lin(ymaxA21f,aA,bA,cA,dA,TmaxA,TempA21f,18.5)
  PhotomeanA26a <- beta.Tmin.Topt.lin(ymaxA26a,aA,bA,cA,dA,TmaxA,TempA26a,23.5)
  PhotomeanA26b <- beta.Tmin.Topt.lin(ymaxA26b,aA,bA,cA,dA,TmaxA,TempA26b,23.5)
  PhotomeanA26c <- beta.Tmin.Topt.lin(ymaxA26c,aA,bA,cA,dA,TmaxA,TempA26c,23.5)
  PhotomeanA26d <- beta.Tmin.Topt.lin(ymaxA26d,aA,bA,cA,dA,TmaxA,TempA26d,23.5)
  PhotomeanA26e <- beta.Tmin.Topt.lin(ymaxA26e,aA,bA,cA,dA,TmaxA,TempA26e,23.5)
  PhotomeanA26f <- beta.Tmin.Topt.lin(ymaxA26f,aA,bA,cA,dA,TmaxA,TempA26f,23.5)
  PhotomeanA31a <- beta.Tmin.Topt.lin(ymaxA31a,aA,bA,cA,dA,TmaxA,TempA31a,28.5)
  PhotomeanA31b <- beta.Tmin.Topt.lin(ymaxA31b,aA,bA,cA,dA,TmaxA,TempA31b,28.5)
  PhotomeanA31c <- beta.Tmin.Topt.lin(ymaxA31c,aA,bA,cA,dA,TmaxA,TempA31c,28.5)
  PhotomeanA31d <- beta.Tmin.Topt.lin(ymaxA31d,aA,bA,cA,dA,TmaxA,TempA31d,28.5)
  PhotomeanA31e <- beta.Tmin.Topt.lin(ymaxA31e,aA,bA,cA,dA,TmaxA,TempA31e,28.5)
  PhotomeanA31f <- beta.Tmin.Topt.lin(ymaxA31f,aA,bA,cA,dA,TmaxA,TempA31f,28.5)
  PhotomeanG21a <- beta.Tmin.Topt.lin(ymaxG21a,aG,bG,cG,dG,TmaxG,TempG21a,18.5)
  PhotomeanG21b <- beta.Tmin.Topt.lin(ymaxG21b,aG,bG,cG,dG,TmaxG,TempG21b,18.5)
  PhotomeanG21c <- beta.Tmin.Topt.lin(ymaxG21c,aG,bG,cG,dG,TmaxG,TempG21c,18.5)
  PhotomeanG21d <- beta.Tmin.Topt.lin(ymaxG21d,aG,bG,cG,dG,TmaxG,TempG21d,18.5)
  PhotomeanG21e <- beta.Tmin.Topt.lin(ymaxG21e,aG,bG,cG,dG,TmaxG,TempG21e,18.5)
  PhotomeanG21f <- beta.Tmin.Topt.lin(ymaxG21f,aG,bG,cG,dG,TmaxG,TempG21f,18.5)
  PhotomeanG26a <- beta.Tmin.Topt.lin(ymaxG26a,aG,bG,cG,dG,TmaxG,TempG26a,23.5)
  PhotomeanG26b <- beta.Tmin.Topt.lin(ymaxG26b,aG,bG,cG,dG,TmaxG,TempG26b,23.5)
  PhotomeanG26c <- beta.Tmin.Topt.lin(ymaxG26c,aG,bG,cG,dG,TmaxG,TempG26c,23.5)
  PhotomeanG26d <- beta.Tmin.Topt.lin(ymaxG26d,aG,bG,cG,dG,TmaxG,TempG26d,23.5)
  PhotomeanG26e <- beta.Tmin.Topt.lin(ymaxG26e,aG,bG,cG,dG,TmaxG,TempG26e,23.5)
  PhotomeanG26f <- beta.Tmin.Topt.lin(ymaxG26f,aG,bG,cG,dG,TmaxG,TempG26f,23.5)
  PhotomeanG26g <- beta.Tmin.Topt.lin(ymaxG26g,aG,bG,cG,dG,TmaxG,TempG26g,23.5)
  PhotomeanG31a <- beta.Tmin.Topt.lin(ymaxG31a,aG,bG,cG,dG,TmaxG,TempG31a,28.5)
  PhotomeanG31b <- beta.Tmin.Topt.lin(ymaxG31b,aG,bG,cG,dG,TmaxG,TempG31b,28.5)
  PhotomeanG31c <- beta.Tmin.Topt.lin(ymaxG31c,aG,bG,cG,dG,TmaxG,TempG31c,28.5)
  PhotomeanG31d <- beta.Tmin.Topt.lin(ymaxG31d,aG,bG,cG,dG,TmaxG,TempG31d,28.5)
  PhotomeanG31e <- beta.Tmin.Topt.lin(ymaxG31e,aG,bG,cG,dG,TmaxG,TempG31e,28.5)
  PhotomeanG31f <- beta.Tmin.Topt.lin(ymaxG31f,aG,bG,cG,dG,TmaxG,TempG31f,28.5)
  PhotomeanR21a <- beta.Tmin.Topt.lin(ymaxR21a,aR,bR,cR,dR,TmaxR,TempR21a,18.5)
  PhotomeanR21b <- beta.Tmin.Topt.lin(ymaxR21b,aR,bR,cR,dR,TmaxR,TempR21b,18.5)
  PhotomeanR21c <- beta.Tmin.Topt.lin(ymaxR21c,aR,bR,cR,dR,TmaxR,TempR21c,18.5)
  PhotomeanR21d <- beta.Tmin.Topt.lin(ymaxR21d,aR,bR,cR,dR,TmaxR,TempR21d,18.5)
  PhotomeanR21e <- beta.Tmin.Topt.lin(ymaxR21e,aR,bR,cR,dR,TmaxR,TempR21e,18.5)
  PhotomeanR21f <- beta.Tmin.Topt.lin(ymaxR21f,aR,bR,cR,dR,TmaxR,TempR21f,18.5)
  PhotomeanR21g <- beta.Tmin.Topt.lin(ymaxR21g,aR,bR,cR,dR,TmaxR,TempR21g,18.5)
  PhotomeanR21h <- beta.Tmin.Topt.lin(ymaxR21h,aR,bR,cR,dR,TmaxR,TempR21h,18.5)
  PhotomeanR26a <- beta.Tmin.Topt.lin(ymaxR26a,aR,bR,cR,dR,TmaxR,TempR26a,23.5)
  PhotomeanR26b <- beta.Tmin.Topt.lin(ymaxR26b,aR,bR,cR,dR,TmaxR,TempR26b,23.5)
  PhotomeanR26c <- beta.Tmin.Topt.lin(ymaxR26c,aR,bR,cR,dR,TmaxR,TempR26c,23.5)
  PhotomeanR26d <- beta.Tmin.Topt.lin(ymaxR26d,aR,bR,cR,dR,TmaxR,TempR26d,23.5)
  PhotomeanR26e <- beta.Tmin.Topt.lin(ymaxR26e,aR,bR,cR,dR,TmaxR,TempR26e,23.5)
  PhotomeanR26f <- beta.Tmin.Topt.lin(ymaxR26f,aR,bR,cR,dR,TmaxR,TempR26f,23.5)
  PhotomeanR31a <- beta.Tmin.Topt.lin(ymaxR31a,aR,bR,cR,dR,TmaxR,TempR31a,28.5)
  PhotomeanR31b <- beta.Tmin.Topt.lin(ymaxR31b,aR,bR,cR,dR,TmaxR,TempR31b,28.5)
  PhotomeanR31c <- beta.Tmin.Topt.lin(ymaxR31c,aR,bR,cR,dR,TmaxR,TempR31c,28.5)
  PhotomeanR31d <- beta.Tmin.Topt.lin(ymaxR31d,aR,bR,cR,dR,TmaxR,TempR31d,28.5)
  PhotomeanR31e <- beta.Tmin.Topt.lin(ymaxR31e,aR,bR,cR,dR,TmaxR,TempR31e,28.5)
  PhotomeanR31f <- beta.Tmin.Topt.lin(ymaxR31f,aR,bR,cR,dR,TmaxR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Tmin and Tmax
Photo_beta_normNLL_all_Tmin.Tmax.lin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                                 ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                                 ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                                 ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                                 ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                                 ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                                 ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                                 ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                                 ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                                 ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                                 ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                                 ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                                 aM,
                                                 aA,
                                                 aG,
                                                 aR,
                                                 bM,
                                                 bA,
                                                 bG,
                                                 bR,
                                                 ToptM,
                                                 ToptA,
                                                 ToptG,
                                                 ToptR,
                                                 cM,
                                                 cA,
                                                 cG,
                                                 cR,
                                                 dM,
                                                 dA,
                                                 dG,
                                                 dR,
                                                 TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                                 TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                                 TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                                 TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                                 TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                                 TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                                 TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                                 TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                                 TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                                 TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                                 TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                                 TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                                 PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                                 PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                                 PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                                 PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                                 PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                                 PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                                 PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                                 PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                                 PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                                 PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                                 PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                                 PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.Tmin.Tmax.lin(ymaxM21a,aM,bM,ToptM,cM,dM,TempM21a,18.5)
  PhotomeanM21b <- beta.Tmin.Tmax.lin(ymaxM21b,aM,bM,ToptM,cM,dM,TempM21b,18.5)
  PhotomeanM21c <- beta.Tmin.Tmax.lin(ymaxM21c,aM,bM,ToptM,cM,dM,TempM21c,18.5)
  PhotomeanM21d <- beta.Tmin.Tmax.lin(ymaxM21d,aM,bM,ToptM,cM,dM,TempM21d,18.5)
  PhotomeanM21e <- beta.Tmin.Tmax.lin(ymaxM21e,aM,bM,ToptM,cM,dM,TempM21e,18.5)
  PhotomeanM21f <- beta.Tmin.Tmax.lin(ymaxM21f,aM,bM,ToptM,cM,dM,TempM21f,18.5)
  PhotomeanM26a <- beta.Tmin.Tmax.lin(ymaxM26a,aM,bM,ToptM,cM,dM,TempM26a,23.5)
  PhotomeanM26b <- beta.Tmin.Tmax.lin(ymaxM26b,aM,bM,ToptM,cM,dM,TempM26b,23.5)
  PhotomeanM26c <- beta.Tmin.Tmax.lin(ymaxM26c,aM,bM,ToptM,cM,dM,TempM26c,23.5)
  PhotomeanM26d <- beta.Tmin.Tmax.lin(ymaxM26d,aM,bM,ToptM,cM,dM,TempM26d,23.5)
  PhotomeanM26e <- beta.Tmin.Tmax.lin(ymaxM26e,aM,bM,ToptM,cM,dM,TempM26e,23.5)
  PhotomeanM26f <- beta.Tmin.Tmax.lin(ymaxM26f,aM,bM,ToptM,cM,dM,TempM26f,23.5)
  PhotomeanM31a <- beta.Tmin.Tmax.lin(ymaxM31a,aM,bM,ToptM,cM,dM,TempM31a,28.5)
  PhotomeanM31b <- beta.Tmin.Tmax.lin(ymaxM31b,aM,bM,ToptM,cM,dM,TempM31b,28.5)
  PhotomeanM31c <- beta.Tmin.Tmax.lin(ymaxM31c,aM,bM,ToptM,cM,dM,TempM31c,28.5)
  PhotomeanM31d <- beta.Tmin.Tmax.lin(ymaxM31d,aM,bM,ToptM,cM,dM,TempM31d,28.5)
  PhotomeanM31e <- beta.Tmin.Tmax.lin(ymaxM31e,aM,bM,ToptM,cM,dM,TempM31e,28.5)
  PhotomeanM31f <- beta.Tmin.Tmax.lin(ymaxM31f,aM,bM,ToptM,cM,dM,TempM31f,28.5)
  PhotomeanA21a <- beta.Tmin.Tmax.lin(ymaxA21a,aA,bA,ToptA,cA,dA,TempA21a,18.5)
  PhotomeanA21b <- beta.Tmin.Tmax.lin(ymaxA21b,aA,bA,ToptA,cA,dA,TempA21b,18.5)
  PhotomeanA21c <- beta.Tmin.Tmax.lin(ymaxA21c,aA,bA,ToptA,cA,dA,TempA21c,18.5)
  PhotomeanA21d <- beta.Tmin.Tmax.lin(ymaxA21d,aA,bA,ToptA,cA,dA,TempA21d,18.5)
  PhotomeanA21e <- beta.Tmin.Tmax.lin(ymaxA21e,aA,bA,ToptA,cA,dA,TempA21e,18.5)
  PhotomeanA21f <- beta.Tmin.Tmax.lin(ymaxA21f,aA,bA,ToptA,cA,dA,TempA21f,18.5)
  PhotomeanA26a <- beta.Tmin.Tmax.lin(ymaxA26a,aA,bA,ToptA,cA,dA,TempA26a,23.5)
  PhotomeanA26b <- beta.Tmin.Tmax.lin(ymaxA26b,aA,bA,ToptA,cA,dA,TempA26b,23.5)
  PhotomeanA26c <- beta.Tmin.Tmax.lin(ymaxA26c,aA,bA,ToptA,cA,dA,TempA26c,23.5)
  PhotomeanA26d <- beta.Tmin.Tmax.lin(ymaxA26d,aA,bA,ToptA,cA,dA,TempA26d,23.5)
  PhotomeanA26e <- beta.Tmin.Tmax.lin(ymaxA26e,aA,bA,ToptA,cA,dA,TempA26e,23.5)
  PhotomeanA26f <- beta.Tmin.Tmax.lin(ymaxA26f,aA,bA,ToptA,cA,dA,TempA26f,23.5)
  PhotomeanA31a <- beta.Tmin.Tmax.lin(ymaxA31a,aA,bA,ToptA,cA,dA,TempA31a,28.5)
  PhotomeanA31b <- beta.Tmin.Tmax.lin(ymaxA31b,aA,bA,ToptA,cA,dA,TempA31b,28.5)
  PhotomeanA31c <- beta.Tmin.Tmax.lin(ymaxA31c,aA,bA,ToptA,cA,dA,TempA31c,28.5)
  PhotomeanA31d <- beta.Tmin.Tmax.lin(ymaxA31d,aA,bA,ToptA,cA,dA,TempA31d,28.5)
  PhotomeanA31e <- beta.Tmin.Tmax.lin(ymaxA31e,aA,bA,ToptA,cA,dA,TempA31e,28.5)
  PhotomeanA31f <- beta.Tmin.Tmax.lin(ymaxA31f,aA,bA,ToptA,cA,dA,TempA31f,28.5)
  PhotomeanG21a <- beta.Tmin.Tmax.lin(ymaxG21a,aG,bG,ToptG,cG,dG,TempG21a,18.5)
  PhotomeanG21b <- beta.Tmin.Tmax.lin(ymaxG21b,aG,bG,ToptG,cG,dG,TempG21b,18.5)
  PhotomeanG21c <- beta.Tmin.Tmax.lin(ymaxG21c,aG,bG,ToptG,cG,dG,TempG21c,18.5)
  PhotomeanG21d <- beta.Tmin.Tmax.lin(ymaxG21d,aG,bG,ToptG,cG,dG,TempG21d,18.5)
  PhotomeanG21e <- beta.Tmin.Tmax.lin(ymaxG21e,aG,bG,ToptG,cG,dG,TempG21e,18.5)
  PhotomeanG21f <- beta.Tmin.Tmax.lin(ymaxG21f,aG,bG,ToptG,cG,dG,TempG21f,18.5)
  PhotomeanG26a <- beta.Tmin.Tmax.lin(ymaxG26a,aG,bG,ToptG,cG,dG,TempG26a,23.5)
  PhotomeanG26b <- beta.Tmin.Tmax.lin(ymaxG26b,aG,bG,ToptG,cG,dG,TempG26b,23.5)
  PhotomeanG26c <- beta.Tmin.Tmax.lin(ymaxG26c,aG,bG,ToptG,cG,dG,TempG26c,23.5)
  PhotomeanG26d <- beta.Tmin.Tmax.lin(ymaxG26d,aG,bG,ToptG,cG,dG,TempG26d,23.5)
  PhotomeanG26e <- beta.Tmin.Tmax.lin(ymaxG26e,aG,bG,ToptG,cG,dG,TempG26e,23.5)
  PhotomeanG26f <- beta.Tmin.Tmax.lin(ymaxG26f,aG,bG,ToptG,cG,dG,TempG26f,23.5)
  PhotomeanG26g <- beta.Tmin.Tmax.lin(ymaxG26g,aG,bG,ToptG,cG,dG,TempG26g,23.5)
  PhotomeanG31a <- beta.Tmin.Tmax.lin(ymaxG31a,aG,bG,ToptG,cG,dG,TempG31a,28.5)
  PhotomeanG31b <- beta.Tmin.Tmax.lin(ymaxG31b,aG,bG,ToptG,cG,dG,TempG31b,28.5)
  PhotomeanG31c <- beta.Tmin.Tmax.lin(ymaxG31c,aG,bG,ToptG,cG,dG,TempG31c,28.5)
  PhotomeanG31d <- beta.Tmin.Tmax.lin(ymaxG31d,aG,bG,ToptG,cG,dG,TempG31d,28.5)
  PhotomeanG31e <- beta.Tmin.Tmax.lin(ymaxG31e,aG,bG,ToptG,cG,dG,TempG31e,28.5)
  PhotomeanG31f <- beta.Tmin.Tmax.lin(ymaxG31f,aG,bG,ToptG,cG,dG,TempG31f,28.5)
  PhotomeanR21a <- beta.Tmin.Tmax.lin(ymaxR21a,aR,bR,ToptR,cR,dR,TempR21a,18.5)
  PhotomeanR21b <- beta.Tmin.Tmax.lin(ymaxR21b,aR,bR,ToptR,cR,dR,TempR21b,18.5)
  PhotomeanR21c <- beta.Tmin.Tmax.lin(ymaxR21c,aR,bR,ToptR,cR,dR,TempR21c,18.5)
  PhotomeanR21d <- beta.Tmin.Tmax.lin(ymaxR21d,aR,bR,ToptR,cR,dR,TempR21d,18.5)
  PhotomeanR21e <- beta.Tmin.Tmax.lin(ymaxR21e,aR,bR,ToptR,cR,dR,TempR21e,18.5)
  PhotomeanR21f <- beta.Tmin.Tmax.lin(ymaxR21f,aR,bR,ToptR,cR,dR,TempR21f,18.5)
  PhotomeanR21g <- beta.Tmin.Tmax.lin(ymaxR21g,aR,bR,ToptR,cR,dR,TempR21g,18.5)
  PhotomeanR21h <- beta.Tmin.Tmax.lin(ymaxR21h,aR,bR,ToptR,cR,dR,TempR21h,18.5)
  PhotomeanR26a <- beta.Tmin.Tmax.lin(ymaxR26a,aR,bR,ToptR,cR,dR,TempR26a,23.5)
  PhotomeanR26b <- beta.Tmin.Tmax.lin(ymaxR26b,aR,bR,ToptR,cR,dR,TempR26b,23.5)
  PhotomeanR26c <- beta.Tmin.Tmax.lin(ymaxR26c,aR,bR,ToptR,cR,dR,TempR26c,23.5)
  PhotomeanR26d <- beta.Tmin.Tmax.lin(ymaxR26d,aR,bR,ToptR,cR,dR,TempR26d,23.5)
  PhotomeanR26e <- beta.Tmin.Tmax.lin(ymaxR26e,aR,bR,ToptR,cR,dR,TempR26e,23.5)
  PhotomeanR26f <- beta.Tmin.Tmax.lin(ymaxR26f,aR,bR,ToptR,cR,dR,TempR26f,23.5)
  PhotomeanR31a <- beta.Tmin.Tmax.lin(ymaxR31a,aR,bR,ToptR,cR,dR,TempR31a,28.5)
  PhotomeanR31b <- beta.Tmin.Tmax.lin(ymaxR31b,aR,bR,ToptR,cR,dR,TempR31b,28.5)
  PhotomeanR31c <- beta.Tmin.Tmax.lin(ymaxR31c,aR,bR,ToptR,cR,dR,TempR31c,28.5)
  PhotomeanR31d <- beta.Tmin.Tmax.lin(ymaxR31d,aR,bR,ToptR,cR,dR,TempR31d,28.5)
  PhotomeanR31e <- beta.Tmin.Tmax.lin(ymaxR31e,aR,bR,ToptR,cR,dR,TempR31e,28.5)
  PhotomeanR31f <- beta.Tmin.Tmax.lin(ymaxR31f,aR,bR,ToptR,cR,dR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Topt and Tmax
Photo_beta_normNLL_all_Topt.Tmax.lin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                                 ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                                 ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                                 ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                                 ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                                 ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                                 ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                                 ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                                 ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                                 ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                                 ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                                 ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                                 TminM,
                                                 TminA,
                                                 TminG,
                                                 TminR,
                                                 aM,
                                                 aA,
                                                 aG,
                                                 aR,
                                                 bM,
                                                 bA,
                                                 bG,
                                                 bR,
                                                 cM,
                                                 cA,
                                                 cG,
                                                 cR,
                                                 dM,
                                                 dA,
                                                 dG,
                                                 dR,
                                                 TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                                 TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                                 TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                                 TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                                 TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                                 TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                                 TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                                 TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                                 TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                                 TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                                 TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                                 TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                                 PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                                 PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                                 PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                                 PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                                 PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                                 PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                                 PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                                 PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                                 PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                                 PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                                 PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                                 PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.Topt.Tmax.lin(ymaxM21a,TminM,aM,bM,cM,dM,TempM21a,18.5)
  PhotomeanM21b <- beta.Topt.Tmax.lin(ymaxM21b,TminM,aM,bM,cM,dM,TempM21b,18.5)
  PhotomeanM21c <- beta.Topt.Tmax.lin(ymaxM21c,TminM,aM,bM,cM,dM,TempM21c,18.5)
  PhotomeanM21d <- beta.Topt.Tmax.lin(ymaxM21d,TminM,aM,bM,cM,dM,TempM21d,18.5)
  PhotomeanM21e <- beta.Topt.Tmax.lin(ymaxM21e,TminM,aM,bM,cM,dM,TempM21e,18.5)
  PhotomeanM21f <- beta.Topt.Tmax.lin(ymaxM21f,TminM,aM,bM,cM,dM,TempM21f,18.5)
  PhotomeanM26a <- beta.Topt.Tmax.lin(ymaxM26a,TminM,aM,bM,cM,dM,TempM26a,23.5)
  PhotomeanM26b <- beta.Topt.Tmax.lin(ymaxM26b,TminM,aM,bM,cM,dM,TempM26b,23.5)
  PhotomeanM26c <- beta.Topt.Tmax.lin(ymaxM26c,TminM,aM,bM,cM,dM,TempM26c,23.5)
  PhotomeanM26d <- beta.Topt.Tmax.lin(ymaxM26d,TminM,aM,bM,cM,dM,TempM26d,23.5)
  PhotomeanM26e <- beta.Topt.Tmax.lin(ymaxM26e,TminM,aM,bM,cM,dM,TempM26e,23.5)
  PhotomeanM26f <- beta.Topt.Tmax.lin(ymaxM26f,TminM,aM,bM,cM,dM,TempM26f,23.5)
  PhotomeanM31a <- beta.Topt.Tmax.lin(ymaxM31a,TminM,aM,bM,cM,dM,TempM31a,28.5)
  PhotomeanM31b <- beta.Topt.Tmax.lin(ymaxM31b,TminM,aM,bM,cM,dM,TempM31b,28.5)
  PhotomeanM31c <- beta.Topt.Tmax.lin(ymaxM31c,TminM,aM,bM,cM,dM,TempM31c,28.5)
  PhotomeanM31d <- beta.Topt.Tmax.lin(ymaxM31d,TminM,aM,bM,cM,dM,TempM31d,28.5)
  PhotomeanM31e <- beta.Topt.Tmax.lin(ymaxM31e,TminM,aM,bM,cM,dM,TempM31e,28.5)
  PhotomeanM31f <- beta.Topt.Tmax.lin(ymaxM31f,TminM,aM,bM,cM,dM,TempM31f,28.5)
  PhotomeanA21a <- beta.Topt.Tmax.lin(ymaxA21a,TminA,aA,bA,cA,dA,TempA21a,18.5)
  PhotomeanA21b <- beta.Topt.Tmax.lin(ymaxA21b,TminA,aA,bA,cA,dA,TempA21b,18.5)
  PhotomeanA21c <- beta.Topt.Tmax.lin(ymaxA21c,TminA,aA,bA,cA,dA,TempA21c,18.5)
  PhotomeanA21d <- beta.Topt.Tmax.lin(ymaxA21d,TminA,aA,bA,cA,dA,TempA21d,18.5)
  PhotomeanA21e <- beta.Topt.Tmax.lin(ymaxA21e,TminA,aA,bA,cA,dA,TempA21e,18.5)
  PhotomeanA21f <- beta.Topt.Tmax.lin(ymaxA21f,TminA,aA,bA,cA,dA,TempA21f,18.5)
  PhotomeanA26a <- beta.Topt.Tmax.lin(ymaxA26a,TminA,aA,bA,cA,dA,TempA26a,23.5)
  PhotomeanA26b <- beta.Topt.Tmax.lin(ymaxA26b,TminA,aA,bA,cA,dA,TempA26b,23.5)
  PhotomeanA26c <- beta.Topt.Tmax.lin(ymaxA26c,TminA,aA,bA,cA,dA,TempA26c,23.5)
  PhotomeanA26d <- beta.Topt.Tmax.lin(ymaxA26d,TminA,aA,bA,cA,dA,TempA26d,23.5)
  PhotomeanA26e <- beta.Topt.Tmax.lin(ymaxA26e,TminA,aA,bA,cA,dA,TempA26e,23.5)
  PhotomeanA26f <- beta.Topt.Tmax.lin(ymaxA26f,TminA,aA,bA,cA,dA,TempA26f,23.5)
  PhotomeanA31a <- beta.Topt.Tmax.lin(ymaxA31a,TminA,aA,bA,cA,dA,TempA31a,28.5)
  PhotomeanA31b <- beta.Topt.Tmax.lin(ymaxA31b,TminA,aA,bA,cA,dA,TempA31b,28.5)
  PhotomeanA31c <- beta.Topt.Tmax.lin(ymaxA31c,TminA,aA,bA,cA,dA,TempA31c,28.5)
  PhotomeanA31d <- beta.Topt.Tmax.lin(ymaxA31d,TminA,aA,bA,cA,dA,TempA31d,28.5)
  PhotomeanA31e <- beta.Topt.Tmax.lin(ymaxA31e,TminA,aA,bA,cA,dA,TempA31e,28.5)
  PhotomeanA31f <- beta.Topt.Tmax.lin(ymaxA31f,TminA,aA,bA,cA,dA,TempA31f,28.5)
  PhotomeanG21a <- beta.Topt.Tmax.lin(ymaxG21a,TminG,aG,bG,cG,dG,TempG21a,18.5)
  PhotomeanG21b <- beta.Topt.Tmax.lin(ymaxG21b,TminG,aG,bG,cG,dG,TempG21b,18.5)
  PhotomeanG21c <- beta.Topt.Tmax.lin(ymaxG21c,TminG,aG,bG,cG,dG,TempG21c,18.5)
  PhotomeanG21d <- beta.Topt.Tmax.lin(ymaxG21d,TminG,aG,bG,cG,dG,TempG21d,18.5)
  PhotomeanG21e <- beta.Topt.Tmax.lin(ymaxG21e,TminG,aG,bG,cG,dG,TempG21e,18.5)
  PhotomeanG21f <- beta.Topt.Tmax.lin(ymaxG21f,TminG,aG,bG,cG,dG,TempG21f,18.5)
  PhotomeanG26a <- beta.Topt.Tmax.lin(ymaxG26a,TminG,aG,bG,cG,dG,TempG26a,23.5)
  PhotomeanG26b <- beta.Topt.Tmax.lin(ymaxG26b,TminG,aG,bG,cG,dG,TempG26b,23.5)
  PhotomeanG26c <- beta.Topt.Tmax.lin(ymaxG26c,TminG,aG,bG,cG,dG,TempG26c,23.5)
  PhotomeanG26d <- beta.Topt.Tmax.lin(ymaxG26d,TminG,aG,bG,cG,dG,TempG26d,23.5)
  PhotomeanG26e <- beta.Topt.Tmax.lin(ymaxG26e,TminG,aG,bG,cG,dG,TempG26e,23.5)
  PhotomeanG26f <- beta.Topt.Tmax.lin(ymaxG26f,TminG,aG,bG,cG,dG,TempG26f,23.5)
  PhotomeanG26g <- beta.Topt.Tmax.lin(ymaxG26g,TminG,aG,bG,cG,dG,TempG26g,23.5)
  PhotomeanG31a <- beta.Topt.Tmax.lin(ymaxG31a,TminG,aG,bG,cG,dG,TempG31a,28.5)
  PhotomeanG31b <- beta.Topt.Tmax.lin(ymaxG31b,TminG,aG,bG,cG,dG,TempG31b,28.5)
  PhotomeanG31c <- beta.Topt.Tmax.lin(ymaxG31c,TminG,aG,bG,cG,dG,TempG31c,28.5)
  PhotomeanG31d <- beta.Topt.Tmax.lin(ymaxG31d,TminG,aG,bG,cG,dG,TempG31d,28.5)
  PhotomeanG31e <- beta.Topt.Tmax.lin(ymaxG31e,TminG,aG,bG,cG,dG,TempG31e,28.5)
  PhotomeanG31f <- beta.Topt.Tmax.lin(ymaxG31f,TminG,aG,bG,cG,dG,TempG31f,28.5)
  PhotomeanR21a <- beta.Topt.Tmax.lin(ymaxR21a,TminR,aR,bR,cR,dR,TempR21a,18.5)
  PhotomeanR21b <- beta.Topt.Tmax.lin(ymaxR21b,TminR,aR,bR,cR,dR,TempR21b,18.5)
  PhotomeanR21c <- beta.Topt.Tmax.lin(ymaxR21c,TminR,aR,bR,cR,dR,TempR21c,18.5)
  PhotomeanR21d <- beta.Topt.Tmax.lin(ymaxR21d,TminR,aR,bR,cR,dR,TempR21d,18.5)
  PhotomeanR21e <- beta.Topt.Tmax.lin(ymaxR21e,TminR,aR,bR,cR,dR,TempR21e,18.5)
  PhotomeanR21f <- beta.Topt.Tmax.lin(ymaxR21f,TminR,aR,bR,cR,dR,TempR21f,18.5)
  PhotomeanR21g <- beta.Topt.Tmax.lin(ymaxR21g,TminR,aR,bR,cR,dR,TempR21g,18.5)
  PhotomeanR21h <- beta.Topt.Tmax.lin(ymaxR21h,TminR,aR,bR,cR,dR,TempR21h,18.5)
  PhotomeanR26a <- beta.Topt.Tmax.lin(ymaxR26a,TminR,aR,bR,cR,dR,TempR26a,23.5)
  PhotomeanR26b <- beta.Topt.Tmax.lin(ymaxR26b,TminR,aR,bR,cR,dR,TempR26b,23.5)
  PhotomeanR26c <- beta.Topt.Tmax.lin(ymaxR26c,TminR,aR,bR,cR,dR,TempR26c,23.5)
  PhotomeanR26d <- beta.Topt.Tmax.lin(ymaxR26d,TminR,aR,bR,cR,dR,TempR26d,23.5)
  PhotomeanR26e <- beta.Topt.Tmax.lin(ymaxR26e,TminR,aR,bR,cR,dR,TempR26e,23.5)
  PhotomeanR26f <- beta.Topt.Tmax.lin(ymaxR26f,TminR,aR,bR,cR,dR,TempR26f,23.5)
  PhotomeanR31a <- beta.Topt.Tmax.lin(ymaxR31a,TminR,aR,bR,cR,dR,TempR31a,28.5)
  PhotomeanR31b <- beta.Topt.Tmax.lin(ymaxR31b,TminR,aR,bR,cR,dR,TempR31b,28.5)
  PhotomeanR31c <- beta.Topt.Tmax.lin(ymaxR31c,TminR,aR,bR,cR,dR,TempR31c,28.5)
  PhotomeanR31d <- beta.Topt.Tmax.lin(ymaxR31d,TminR,aR,bR,cR,dR,TempR31d,28.5)
  PhotomeanR31e <- beta.Topt.Tmax.lin(ymaxR31e,TminR,aR,bR,cR,dR,TempR31e,28.5)
  PhotomeanR31f <- beta.Topt.Tmax.lin(ymaxR31f,TminR,aR,bR,cR,dR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Tmin
Photo_beta_normNLL_all_Tmin.lin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                            ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                            ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                            ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                            ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                            ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                            ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                            ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                            ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                            ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                            ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                            ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                            TmaxM,
                                            TmaxA,
                                            TmaxG,
                                            TmaxR,
                                            aM,
                                            aA,
                                            aG,
                                            aR,
                                            bM,
                                            bA,
                                            bG,
                                            bR,
                                            ToptM,
                                            ToptA,
                                            ToptG,
                                            ToptR,
                                            TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                            TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                            TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                            TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                            TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                            TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                            TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                            TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                            TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                            TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                            TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                            TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                            PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                            PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                            PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                            PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                            PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                            PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                            PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                            PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                            PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                            PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                            PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                            PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.Tmin.lin(ymaxM21a,aM,bM,ToptM,TmaxM,TempM21a,18.5)
  PhotomeanM21b <- beta.Tmin.lin(ymaxM21b,aM,bM,ToptM,TmaxM,TempM21b,18.5)
  PhotomeanM21c <- beta.Tmin.lin(ymaxM21c,aM,bM,ToptM,TmaxM,TempM21c,18.5)
  PhotomeanM21d <- beta.Tmin.lin(ymaxM21d,aM,bM,ToptM,TmaxM,TempM21d,18.5)
  PhotomeanM21e <- beta.Tmin.lin(ymaxM21e,aM,bM,ToptM,TmaxM,TempM21e,18.5)
  PhotomeanM21f <- beta.Tmin.lin(ymaxM21f,aM,bM,ToptM,TmaxM,TempM21f,18.5)
  PhotomeanM26a <- beta.Tmin.lin(ymaxM26a,aM,bM,ToptM,TmaxM,TempM26a,23.5)
  PhotomeanM26b <- beta.Tmin.lin(ymaxM26b,aM,bM,ToptM,TmaxM,TempM26b,23.5)
  PhotomeanM26c <- beta.Tmin.lin(ymaxM26c,aM,bM,ToptM,TmaxM,TempM26c,23.5)
  PhotomeanM26d <- beta.Tmin.lin(ymaxM26d,aM,bM,ToptM,TmaxM,TempM26d,23.5)
  PhotomeanM26e <- beta.Tmin.lin(ymaxM26e,aM,bM,ToptM,TmaxM,TempM26e,23.5)
  PhotomeanM26f <- beta.Tmin.lin(ymaxM26f,aM,bM,ToptM,TmaxM,TempM26f,23.5)
  PhotomeanM31a <- beta.Tmin.lin(ymaxM31a,aM,bM,ToptM,TmaxM,TempM31a,28.5)
  PhotomeanM31b <- beta.Tmin.lin(ymaxM31b,aM,bM,ToptM,TmaxM,TempM31b,28.5)
  PhotomeanM31c <- beta.Tmin.lin(ymaxM31c,aM,bM,ToptM,TmaxM,TempM31c,28.5)
  PhotomeanM31d <- beta.Tmin.lin(ymaxM31d,aM,bM,ToptM,TmaxM,TempM31d,28.5)
  PhotomeanM31e <- beta.Tmin.lin(ymaxM31e,aM,bM,ToptM,TmaxM,TempM31e,28.5)
  PhotomeanM31f <- beta.Tmin.lin(ymaxM31f,aM,bM,ToptM,TmaxM,TempM31f,28.5)
  PhotomeanA21a <- beta.Tmin.lin(ymaxA21a,aA,bA,ToptA,TmaxA,TempA21a,18.5)
  PhotomeanA21b <- beta.Tmin.lin(ymaxA21b,aA,bA,ToptA,TmaxA,TempA21b,18.5)
  PhotomeanA21c <- beta.Tmin.lin(ymaxA21c,aA,bA,ToptA,TmaxA,TempA21c,18.5)
  PhotomeanA21d <- beta.Tmin.lin(ymaxA21d,aA,bA,ToptA,TmaxA,TempA21d,18.5)
  PhotomeanA21e <- beta.Tmin.lin(ymaxA21e,aA,bA,ToptA,TmaxA,TempA21e,18.5)
  PhotomeanA21f <- beta.Tmin.lin(ymaxA21f,aA,bA,ToptA,TmaxA,TempA21f,18.5)
  PhotomeanA26a <- beta.Tmin.lin(ymaxA26a,aA,bA,ToptA,TmaxA,TempA26a,23.5)
  PhotomeanA26b <- beta.Tmin.lin(ymaxA26b,aA,bA,ToptA,TmaxA,TempA26b,23.5)
  PhotomeanA26c <- beta.Tmin.lin(ymaxA26c,aA,bA,ToptA,TmaxA,TempA26c,23.5)
  PhotomeanA26d <- beta.Tmin.lin(ymaxA26d,aA,bA,ToptA,TmaxA,TempA26d,23.5)
  PhotomeanA26e <- beta.Tmin.lin(ymaxA26e,aA,bA,ToptA,TmaxA,TempA26e,23.5)
  PhotomeanA26f <- beta.Tmin.lin(ymaxA26f,aA,bA,ToptA,TmaxA,TempA26f,23.5)
  PhotomeanA31a <- beta.Tmin.lin(ymaxA31a,aA,bA,ToptA,TmaxA,TempA31a,28.5)
  PhotomeanA31b <- beta.Tmin.lin(ymaxA31b,aA,bA,ToptA,TmaxA,TempA31b,28.5)
  PhotomeanA31c <- beta.Tmin.lin(ymaxA31c,aA,bA,ToptA,TmaxA,TempA31c,28.5)
  PhotomeanA31d <- beta.Tmin.lin(ymaxA31d,aA,bA,ToptA,TmaxA,TempA31d,28.5)
  PhotomeanA31e <- beta.Tmin.lin(ymaxA31e,aA,bA,ToptA,TmaxA,TempA31e,28.5)
  PhotomeanA31f <- beta.Tmin.lin(ymaxA31f,aA,bA,ToptA,TmaxA,TempA31f,28.5)
  PhotomeanG21a <- beta.Tmin.lin(ymaxG21a,aG,bG,ToptG,TmaxG,TempG21a,18.5)
  PhotomeanG21b <- beta.Tmin.lin(ymaxG21b,aG,bG,ToptG,TmaxG,TempG21b,18.5)
  PhotomeanG21c <- beta.Tmin.lin(ymaxG21c,aG,bG,ToptG,TmaxG,TempG21c,18.5)
  PhotomeanG21d <- beta.Tmin.lin(ymaxG21d,aG,bG,ToptG,TmaxG,TempG21d,18.5)
  PhotomeanG21e <- beta.Tmin.lin(ymaxG21e,aG,bG,ToptG,TmaxG,TempG21e,18.5)
  PhotomeanG21f <- beta.Tmin.lin(ymaxG21f,aG,bG,ToptG,TmaxG,TempG21f,18.5)
  PhotomeanG26a <- beta.Tmin.lin(ymaxG26a,aG,bG,ToptG,TmaxG,TempG26a,23.5)
  PhotomeanG26b <- beta.Tmin.lin(ymaxG26b,aG,bG,ToptG,TmaxG,TempG26b,23.5)
  PhotomeanG26c <- beta.Tmin.lin(ymaxG26c,aG,bG,ToptG,TmaxG,TempG26c,23.5)
  PhotomeanG26d <- beta.Tmin.lin(ymaxG26d,aG,bG,ToptG,TmaxG,TempG26d,23.5)
  PhotomeanG26e <- beta.Tmin.lin(ymaxG26e,aG,bG,ToptG,TmaxG,TempG26e,23.5)
  PhotomeanG26f <- beta.Tmin.lin(ymaxG26f,aG,bG,ToptG,TmaxG,TempG26f,23.5)
  PhotomeanG26g <- beta.Tmin.lin(ymaxG26g,aG,bG,ToptG,TmaxG,TempG26g,23.5)
  PhotomeanG31a <- beta.Tmin.lin(ymaxG31a,aG,bG,ToptG,TmaxG,TempG31a,28.5)
  PhotomeanG31b <- beta.Tmin.lin(ymaxG31b,aG,bG,ToptG,TmaxG,TempG31b,28.5)
  PhotomeanG31c <- beta.Tmin.lin(ymaxG31c,aG,bG,ToptG,TmaxG,TempG31c,28.5)
  PhotomeanG31d <- beta.Tmin.lin(ymaxG31d,aG,bG,ToptG,TmaxG,TempG31d,28.5)
  PhotomeanG31e <- beta.Tmin.lin(ymaxG31e,aG,bG,ToptG,TmaxG,TempG31e,28.5)
  PhotomeanG31f <- beta.Tmin.lin(ymaxG31f,aG,bG,ToptG,TmaxG,TempG31f,28.5)
  PhotomeanR21a <- beta.Tmin.lin(ymaxR21a,aR,bR,ToptR,TmaxR,TempR21a,18.5)
  PhotomeanR21b <- beta.Tmin.lin(ymaxR21b,aR,bR,ToptR,TmaxR,TempR21b,18.5)
  PhotomeanR21c <- beta.Tmin.lin(ymaxR21c,aR,bR,ToptR,TmaxR,TempR21c,18.5)
  PhotomeanR21d <- beta.Tmin.lin(ymaxR21d,aR,bR,ToptR,TmaxR,TempR21d,18.5)
  PhotomeanR21e <- beta.Tmin.lin(ymaxR21e,aR,bR,ToptR,TmaxR,TempR21e,18.5)
  PhotomeanR21f <- beta.Tmin.lin(ymaxR21f,aR,bR,ToptR,TmaxR,TempR21f,18.5)
  PhotomeanR21g <- beta.Tmin.lin(ymaxR21g,aR,bR,ToptR,TmaxR,TempR21g,18.5)
  PhotomeanR21h <- beta.Tmin.lin(ymaxR21h,aR,bR,ToptR,TmaxR,TempR21h,18.5)
  PhotomeanR26a <- beta.Tmin.lin(ymaxR26a,aR,bR,ToptR,TmaxR,TempR26a,23.5)
  PhotomeanR26b <- beta.Tmin.lin(ymaxR26b,aR,bR,ToptR,TmaxR,TempR26b,23.5)
  PhotomeanR26c <- beta.Tmin.lin(ymaxR26c,aR,bR,ToptR,TmaxR,TempR26c,23.5)
  PhotomeanR26d <- beta.Tmin.lin(ymaxR26d,aR,bR,ToptR,TmaxR,TempR26d,23.5)
  PhotomeanR26e <- beta.Tmin.lin(ymaxR26e,aR,bR,ToptR,TmaxR,TempR26e,23.5)
  PhotomeanR26f <- beta.Tmin.lin(ymaxR26f,aR,bR,ToptR,TmaxR,TempR26f,23.5)
  PhotomeanR31a <- beta.Tmin.lin(ymaxR31a,aR,bR,ToptR,TmaxR,TempR31a,28.5)
  PhotomeanR31b <- beta.Tmin.lin(ymaxR31b,aR,bR,ToptR,TmaxR,TempR31b,28.5)
  PhotomeanR31c <- beta.Tmin.lin(ymaxR31c,aR,bR,ToptR,TmaxR,TempR31c,28.5)
  PhotomeanR31d <- beta.Tmin.lin(ymaxR31d,aR,bR,ToptR,TmaxR,TempR31d,28.5)
  PhotomeanR31e <- beta.Tmin.lin(ymaxR31e,aR,bR,ToptR,TmaxR,TempR31e,28.5)
  PhotomeanR31f <- beta.Tmin.lin(ymaxR31f,aR,bR,ToptR,TmaxR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Topt
Photo_beta_normNLL_all_Topt.lin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                           ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                           ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                           ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                           ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                           ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                           ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                           ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                           ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                           ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                           ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                           ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                           TmaxM,
                                           TmaxA,
                                           TmaxG,
                                           TmaxR,
                                           TminM,
                                           TminA,
                                           TminG,
                                           TminR,
                                           aM,
                                           aA,
                                           aG,
                                           aR,
                                           bM,
                                           bA,
                                           bG,
                                           bR,
                                           TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                           TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                           TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                           TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                           TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                           TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                           TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                           TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                           TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                           TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                           TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                           TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                           PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                           PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                           PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                           PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                           PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                           PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                           PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                           PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                           PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                           PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                           PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                           PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.Topt.lin(ymaxM21a,TminM,aM,bM,TmaxM,TempM21a,18.5)
  PhotomeanM21b <- beta.Topt.lin(ymaxM21b,TminM,aM,bM,TmaxM,TempM21b,18.5)
  PhotomeanM21c <- beta.Topt.lin(ymaxM21c,TminM,aM,bM,TmaxM,TempM21c,18.5)
  PhotomeanM21d <- beta.Topt.lin(ymaxM21d,TminM,aM,bM,TmaxM,TempM21d,18.5)
  PhotomeanM21e <- beta.Topt.lin(ymaxM21e,TminM,aM,bM,TmaxM,TempM21e,18.5)
  PhotomeanM21f <- beta.Topt.lin(ymaxM21f,TminM,aM,bM,TmaxM,TempM21f,18.5)
  PhotomeanM26a <- beta.Topt.lin(ymaxM26a,TminM,aM,bM,TmaxM,TempM26a,23.5)
  PhotomeanM26b <- beta.Topt.lin(ymaxM26b,TminM,aM,bM,TmaxM,TempM26b,23.5)
  PhotomeanM26c <- beta.Topt.lin(ymaxM26c,TminM,aM,bM,TmaxM,TempM26c,23.5)
  PhotomeanM26d <- beta.Topt.lin(ymaxM26d,TminM,aM,bM,TmaxM,TempM26d,23.5)
  PhotomeanM26e <- beta.Topt.lin(ymaxM26e,TminM,aM,bM,TmaxM,TempM26e,23.5)
  PhotomeanM26f <- beta.Topt.lin(ymaxM26f,TminM,aM,bM,TmaxM,TempM26f,23.5)
  PhotomeanM31a <- beta.Topt.lin(ymaxM31a,TminM,aM,bM,TmaxM,TempM31a,28.5)
  PhotomeanM31b <- beta.Topt.lin(ymaxM31b,TminM,aM,bM,TmaxM,TempM31b,28.5)
  PhotomeanM31c <- beta.Topt.lin(ymaxM31c,TminM,aM,bM,TmaxM,TempM31c,28.5)
  PhotomeanM31d <- beta.Topt.lin(ymaxM31d,TminM,aM,bM,TmaxM,TempM31d,28.5)
  PhotomeanM31e <- beta.Topt.lin(ymaxM31e,TminM,aM,bM,TmaxM,TempM31e,28.5)
  PhotomeanM31f <- beta.Topt.lin(ymaxM31f,TminM,aM,bM,TmaxM,TempM31f,28.5)
  PhotomeanA21a <- beta.Topt.lin(ymaxA21a,TminA,aA,bA,TmaxA,TempA21a,18.5)
  PhotomeanA21b <- beta.Topt.lin(ymaxA21b,TminA,aA,bA,TmaxA,TempA21b,18.5)
  PhotomeanA21c <- beta.Topt.lin(ymaxA21c,TminA,aA,bA,TmaxA,TempA21c,18.5)
  PhotomeanA21d <- beta.Topt.lin(ymaxA21d,TminA,aA,bA,TmaxA,TempA21d,18.5)
  PhotomeanA21e <- beta.Topt.lin(ymaxA21e,TminA,aA,bA,TmaxA,TempA21e,18.5)
  PhotomeanA21f <- beta.Topt.lin(ymaxA21f,TminA,aA,bA,TmaxA,TempA21f,18.5)
  PhotomeanA26a <- beta.Topt.lin(ymaxA26a,TminA,aA,bA,TmaxA,TempA26a,23.5)
  PhotomeanA26b <- beta.Topt.lin(ymaxA26b,TminA,aA,bA,TmaxA,TempA26b,23.5)
  PhotomeanA26c <- beta.Topt.lin(ymaxA26c,TminA,aA,bA,TmaxA,TempA26c,23.5)
  PhotomeanA26d <- beta.Topt.lin(ymaxA26d,TminA,aA,bA,TmaxA,TempA26d,23.5)
  PhotomeanA26e <- beta.Topt.lin(ymaxA26e,TminA,aA,bA,TmaxA,TempA26e,23.5)
  PhotomeanA26f <- beta.Topt.lin(ymaxA26f,TminA,aA,bA,TmaxA,TempA26f,23.5)
  PhotomeanA31a <- beta.Topt.lin(ymaxA31a,TminA,aA,bA,TmaxA,TempA31a,28.5)
  PhotomeanA31b <- beta.Topt.lin(ymaxA31b,TminA,aA,bA,TmaxA,TempA31b,28.5)
  PhotomeanA31c <- beta.Topt.lin(ymaxA31c,TminA,aA,bA,TmaxA,TempA31c,28.5)
  PhotomeanA31d <- beta.Topt.lin(ymaxA31d,TminA,aA,bA,TmaxA,TempA31d,28.5)
  PhotomeanA31e <- beta.Topt.lin(ymaxA31e,TminA,aA,bA,TmaxA,TempA31e,28.5)
  PhotomeanA31f <- beta.Topt.lin(ymaxA31f,TminA,aA,bA,TmaxA,TempA31f,28.5)
  PhotomeanG21a <- beta.Topt.lin(ymaxG21a,TminG,aG,bG,TmaxG,TempG21a,18.5)
  PhotomeanG21b <- beta.Topt.lin(ymaxG21b,TminG,aG,bG,TmaxG,TempG21b,18.5)
  PhotomeanG21c <- beta.Topt.lin(ymaxG21c,TminG,aG,bG,TmaxG,TempG21c,18.5)
  PhotomeanG21d <- beta.Topt.lin(ymaxG21d,TminG,aG,bG,TmaxG,TempG21d,18.5)
  PhotomeanG21e <- beta.Topt.lin(ymaxG21e,TminG,aG,bG,TmaxG,TempG21e,18.5)
  PhotomeanG21f <- beta.Topt.lin(ymaxG21f,TminG,aG,bG,TmaxG,TempG21f,18.5)
  PhotomeanG26a <- beta.Topt.lin(ymaxG26a,TminG,aG,bG,TmaxG,TempG26a,23.5)
  PhotomeanG26b <- beta.Topt.lin(ymaxG26b,TminG,aG,bG,TmaxG,TempG26b,23.5)
  PhotomeanG26c <- beta.Topt.lin(ymaxG26c,TminG,aG,bG,TmaxG,TempG26c,23.5)
  PhotomeanG26d <- beta.Topt.lin(ymaxG26d,TminG,aG,bG,TmaxG,TempG26d,23.5)
  PhotomeanG26e <- beta.Topt.lin(ymaxG26e,TminG,aG,bG,TmaxG,TempG26e,23.5)
  PhotomeanG26f <- beta.Topt.lin(ymaxG26f,TminG,aG,bG,TmaxG,TempG26f,23.5)
  PhotomeanG26g <- beta.Topt.lin(ymaxG26g,TminG,aG,bG,TmaxG,TempG26g,23.5)
  PhotomeanG31a <- beta.Topt.lin(ymaxG31a,TminG,aG,bG,TmaxG,TempG31a,28.5)
  PhotomeanG31b <- beta.Topt.lin(ymaxG31b,TminG,aG,bG,TmaxG,TempG31b,28.5)
  PhotomeanG31c <- beta.Topt.lin(ymaxG31c,TminG,aG,bG,TmaxG,TempG31c,28.5)
  PhotomeanG31d <- beta.Topt.lin(ymaxG31d,TminG,aG,bG,TmaxG,TempG31d,28.5)
  PhotomeanG31e <- beta.Topt.lin(ymaxG31e,TminG,aG,bG,TmaxG,TempG31e,28.5)
  PhotomeanG31f <- beta.Topt.lin(ymaxG31f,TminG,aG,bG,TmaxG,TempG31f,28.5)
  PhotomeanR21a <- beta.Topt.lin(ymaxR21a,TminR,aR,bR,TmaxR,TempR21a,18.5)
  PhotomeanR21b <- beta.Topt.lin(ymaxR21b,TminR,aR,bR,TmaxR,TempR21b,18.5)
  PhotomeanR21c <- beta.Topt.lin(ymaxR21c,TminR,aR,bR,TmaxR,TempR21c,18.5)
  PhotomeanR21d <- beta.Topt.lin(ymaxR21d,TminR,aR,bR,TmaxR,TempR21d,18.5)
  PhotomeanR21e <- beta.Topt.lin(ymaxR21e,TminR,aR,bR,TmaxR,TempR21e,18.5)
  PhotomeanR21f <- beta.Topt.lin(ymaxR21f,TminR,aR,bR,TmaxR,TempR21f,18.5)
  PhotomeanR21g <- beta.Topt.lin(ymaxR21g,TminR,aR,bR,TmaxR,TempR21g,18.5)
  PhotomeanR21h <- beta.Topt.lin(ymaxR21h,TminR,aR,bR,TmaxR,TempR21h,18.5)
  PhotomeanR26a <- beta.Topt.lin(ymaxR26a,TminR,aR,bR,TmaxR,TempR26a,23.5)
  PhotomeanR26b <- beta.Topt.lin(ymaxR26b,TminR,aR,bR,TmaxR,TempR26b,23.5)
  PhotomeanR26c <- beta.Topt.lin(ymaxR26c,TminR,aR,bR,TmaxR,TempR26c,23.5)
  PhotomeanR26d <- beta.Topt.lin(ymaxR26d,TminR,aR,bR,TmaxR,TempR26d,23.5)
  PhotomeanR26e <- beta.Topt.lin(ymaxR26e,TminR,aR,bR,TmaxR,TempR26e,23.5)
  PhotomeanR26f <- beta.Topt.lin(ymaxR26f,TminR,aR,bR,TmaxR,TempR26f,23.5)
  PhotomeanR31a <- beta.Topt.lin(ymaxR31a,TminR,aR,bR,TmaxR,TempR31a,28.5)
  PhotomeanR31b <- beta.Topt.lin(ymaxR31b,TminR,aR,bR,TmaxR,TempR31b,28.5)
  PhotomeanR31c <- beta.Topt.lin(ymaxR31c,TminR,aR,bR,TmaxR,TempR31c,28.5)
  PhotomeanR31d <- beta.Topt.lin(ymaxR31d,TminR,aR,bR,TmaxR,TempR31d,28.5)
  PhotomeanR31e <- beta.Topt.lin(ymaxR31e,TminR,aR,bR,TmaxR,TempR31e,28.5)
  PhotomeanR31f <- beta.Topt.lin(ymaxR31f,TminR,aR,bR,TmaxR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for acclimation of Tmax
Photo_beta_normNLL_all_Tmax.lin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                            ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                            ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                            ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                            ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                            ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                            ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                            ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                            ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                            ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                            ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                            ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                            aM,
                                            aA,
                                            aG,
                                            aR,
                                            bM,
                                            bA,
                                            bG,
                                            bR,
                                            TminM,
                                            TminA,
                                            TminG,
                                            TminR,
                                            ToptM,
                                            ToptA,
                                            ToptG,
                                            ToptR,
                                            TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                            TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                            TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                            TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                            TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                            TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                            TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                            TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                            TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                            TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                            TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                            TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                            PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                            PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                            PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                            PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                            PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                            PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                            PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                            PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                            PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                            PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                            PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                            PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta.Tmax.lin(ymaxM21a,TminM,ToptM,aM,bM,TempM21a,18.5)
  PhotomeanM21b <- beta.Tmax.lin(ymaxM21b,TminM,ToptM,aM,bM,TempM21b,18.5)
  PhotomeanM21c <- beta.Tmax.lin(ymaxM21c,TminM,ToptM,aM,bM,TempM21c,18.5)
  PhotomeanM21d <- beta.Tmax.lin(ymaxM21d,TminM,ToptM,aM,bM,TempM21d,18.5)
  PhotomeanM21e <- beta.Tmax.lin(ymaxM21e,TminM,ToptM,aM,bM,TempM21e,18.5)
  PhotomeanM21f <- beta.Tmax.lin(ymaxM21f,TminM,ToptM,aM,bM,TempM21f,18.5)
  PhotomeanM26a <- beta.Tmax.lin(ymaxM26a,TminM,ToptM,aM,bM,TempM26a,23.5)
  PhotomeanM26b <- beta.Tmax.lin(ymaxM26b,TminM,ToptM,aM,bM,TempM26b,23.5)
  PhotomeanM26c <- beta.Tmax.lin(ymaxM26c,TminM,ToptM,aM,bM,TempM26c,23.5)
  PhotomeanM26d <- beta.Tmax.lin(ymaxM26d,TminM,ToptM,aM,bM,TempM26d,23.5)
  PhotomeanM26e <- beta.Tmax.lin(ymaxM26e,TminM,ToptM,aM,bM,TempM26e,23.5)
  PhotomeanM26f <- beta.Tmax.lin(ymaxM26f,TminM,ToptM,aM,bM,TempM26f,23.5)
  PhotomeanM31a <- beta.Tmax.lin(ymaxM31a,TminM,ToptM,aM,bM,TempM31a,28.5)
  PhotomeanM31b <- beta.Tmax.lin(ymaxM31b,TminM,ToptM,aM,bM,TempM31b,28.5)
  PhotomeanM31c <- beta.Tmax.lin(ymaxM31c,TminM,ToptM,aM,bM,TempM31c,28.5)
  PhotomeanM31d <- beta.Tmax.lin(ymaxM31d,TminM,ToptM,aM,bM,TempM31d,28.5)
  PhotomeanM31e <- beta.Tmax.lin(ymaxM31e,TminM,ToptM,aM,bM,TempM31e,28.5)
  PhotomeanM31f <- beta.Tmax.lin(ymaxM31f,TminM,ToptM,aM,bM,TempM31f,28.5)
  PhotomeanA21a <- beta.Tmax.lin(ymaxA21a,TminA,ToptA,aA,bA,TempA21a,18.5)
  PhotomeanA21b <- beta.Tmax.lin(ymaxA21b,TminA,ToptA,aA,bA,TempA21b,18.5)
  PhotomeanA21c <- beta.Tmax.lin(ymaxA21c,TminA,ToptA,aA,bA,TempA21c,18.5)
  PhotomeanA21d <- beta.Tmax.lin(ymaxA21d,TminA,ToptA,aA,bA,TempA21d,18.5)
  PhotomeanA21e <- beta.Tmax.lin(ymaxA21e,TminA,ToptA,aA,bA,TempA21e,18.5)
  PhotomeanA21f <- beta.Tmax.lin(ymaxA21f,TminA,ToptA,aA,bA,TempA21f,18.5)
  PhotomeanA26a <- beta.Tmax.lin(ymaxA26a,TminA,ToptA,aA,bA,TempA26a,23.5)
  PhotomeanA26b <- beta.Tmax.lin(ymaxA26b,TminA,ToptA,aA,bA,TempA26b,23.5)
  PhotomeanA26c <- beta.Tmax.lin(ymaxA26c,TminA,ToptA,aA,bA,TempA26c,23.5)
  PhotomeanA26d <- beta.Tmax.lin(ymaxA26d,TminA,ToptA,aA,bA,TempA26d,23.5)
  PhotomeanA26e <- beta.Tmax.lin(ymaxA26e,TminA,ToptA,aA,bA,TempA26e,23.5)
  PhotomeanA26f <- beta.Tmax.lin(ymaxA26f,TminA,ToptA,aA,bA,TempA26f,23.5)
  PhotomeanA31a <- beta.Tmax.lin(ymaxA31a,TminA,ToptA,aA,bA,TempA31a,28.5)
  PhotomeanA31b <- beta.Tmax.lin(ymaxA31b,TminA,ToptA,aA,bA,TempA31b,28.5)
  PhotomeanA31c <- beta.Tmax.lin(ymaxA31c,TminA,ToptA,aA,bA,TempA31c,28.5)
  PhotomeanA31d <- beta.Tmax.lin(ymaxA31d,TminA,ToptA,aA,bA,TempA31d,28.5)
  PhotomeanA31e <- beta.Tmax.lin(ymaxA31e,TminA,ToptA,aA,bA,TempA31e,28.5)
  PhotomeanA31f <- beta.Tmax.lin(ymaxA31f,TminA,ToptA,aA,bA,TempA31f,28.5)
  PhotomeanG21a <- beta.Tmax.lin(ymaxG21a,TminG,ToptG,aG,bG,TempG21a,18.5)
  PhotomeanG21b <- beta.Tmax.lin(ymaxG21b,TminG,ToptG,aG,bG,TempG21b,18.5)
  PhotomeanG21c <- beta.Tmax.lin(ymaxG21c,TminG,ToptG,aG,bG,TempG21c,18.5)
  PhotomeanG21d <- beta.Tmax.lin(ymaxG21d,TminG,ToptG,aG,bG,TempG21d,18.5)
  PhotomeanG21e <- beta.Tmax.lin(ymaxG21e,TminG,ToptG,aG,bG,TempG21e,18.5)
  PhotomeanG21f <- beta.Tmax.lin(ymaxG21f,TminG,ToptG,aG,bG,TempG21f,18.5)
  PhotomeanG26a <- beta.Tmax.lin(ymaxG26a,TminG,ToptG,aG,bG,TempG26a,23.5)
  PhotomeanG26b <- beta.Tmax.lin(ymaxG26b,TminG,ToptG,aG,bG,TempG26b,23.5)
  PhotomeanG26c <- beta.Tmax.lin(ymaxG26c,TminG,ToptG,aG,bG,TempG26c,23.5)
  PhotomeanG26d <- beta.Tmax.lin(ymaxG26d,TminG,ToptG,aG,bG,TempG26d,23.5)
  PhotomeanG26e <- beta.Tmax.lin(ymaxG26e,TminG,ToptG,aG,bG,TempG26e,23.5)
  PhotomeanG26f <- beta.Tmax.lin(ymaxG26f,TminG,ToptG,aG,bG,TempG26f,23.5)
  PhotomeanG26g <- beta.Tmax.lin(ymaxG26g,TminG,ToptG,aG,bG,TempG26g,23.5)
  PhotomeanG31a <- beta.Tmax.lin(ymaxG31a,TminG,ToptG,aG,bG,TempG31a,28.5)
  PhotomeanG31b <- beta.Tmax.lin(ymaxG31b,TminG,ToptG,aG,bG,TempG31b,28.5)
  PhotomeanG31c <- beta.Tmax.lin(ymaxG31c,TminG,ToptG,aG,bG,TempG31c,28.5)
  PhotomeanG31d <- beta.Tmax.lin(ymaxG31d,TminG,ToptG,aG,bG,TempG31d,28.5)
  PhotomeanG31e <- beta.Tmax.lin(ymaxG31e,TminG,ToptG,aG,bG,TempG31e,28.5)
  PhotomeanG31f <- beta.Tmax.lin(ymaxG31f,TminG,ToptG,aG,bG,TempG31f,28.5)
  PhotomeanR21a <- beta.Tmax.lin(ymaxR21a,TminR,ToptR,aR,bR,TempR21a,18.5)
  PhotomeanR21b <- beta.Tmax.lin(ymaxR21b,TminR,ToptR,aR,bR,TempR21b,18.5)
  PhotomeanR21c <- beta.Tmax.lin(ymaxR21c,TminR,ToptR,aR,bR,TempR21c,18.5)
  PhotomeanR21d <- beta.Tmax.lin(ymaxR21d,TminR,ToptR,aR,bR,TempR21d,18.5)
  PhotomeanR21e <- beta.Tmax.lin(ymaxR21e,TminR,ToptR,aR,bR,TempR21e,18.5)
  PhotomeanR21f <- beta.Tmax.lin(ymaxR21f,TminR,ToptR,aR,bR,TempR21f,18.5)
  PhotomeanR21g <- beta.Tmax.lin(ymaxR21g,TminR,ToptR,aR,bR,TempR21g,18.5)
  PhotomeanR21h <- beta.Tmax.lin(ymaxR21h,TminR,ToptR,aR,bR,TempR21h,18.5)
  PhotomeanR26a <- beta.Tmax.lin(ymaxR26a,TminR,ToptR,aR,bR,TempR26a,23.5)
  PhotomeanR26b <- beta.Tmax.lin(ymaxR26b,TminR,ToptR,aR,bR,TempR26b,23.5)
  PhotomeanR26c <- beta.Tmax.lin(ymaxR26c,TminR,ToptR,aR,bR,TempR26c,23.5)
  PhotomeanR26d <- beta.Tmax.lin(ymaxR26d,TminR,ToptR,aR,bR,TempR26d,23.5)
  PhotomeanR26e <- beta.Tmax.lin(ymaxR26e,TminR,ToptR,aR,bR,TempR26e,23.5)
  PhotomeanR26f <- beta.Tmax.lin(ymaxR26f,TminR,ToptR,aR,bR,TempR26f,23.5)
  PhotomeanR31a <- beta.Tmax.lin(ymaxR31a,TminR,ToptR,aR,bR,TempR31a,28.5)
  PhotomeanR31b <- beta.Tmax.lin(ymaxR31b,TminR,ToptR,aR,bR,TempR31b,28.5)
  PhotomeanR31c <- beta.Tmax.lin(ymaxR31c,TminR,ToptR,aR,bR,TempR31c,28.5)
  PhotomeanR31d <- beta.Tmax.lin(ymaxR31d,TminR,ToptR,aR,bR,TempR31d,28.5)
  PhotomeanR31e <- beta.Tmax.lin(ymaxR31e,TminR,ToptR,aR,bR,TempR31e,28.5)
  PhotomeanR31f <- beta.Tmax.lin(ymaxR31f,TminR,ToptR,aR,bR,TempR31f,28.5)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

#NLL function for modified beta without acclimation
Photo_beta_normNLL_all_bin <- function(sdPhoto,ymaxM21a,ymaxM21b,ymaxM21c,ymaxM21d,ymaxM21e,ymaxM21f,
                                       ymaxM26a,ymaxM26b,ymaxM26c,ymaxM26d,ymaxM26e,ymaxM26f,
                                       ymaxM31a,ymaxM31b,ymaxM31c,ymaxM31d,ymaxM31e,ymaxM31f,
                                       ymaxA21a,ymaxA21b,ymaxA21c,ymaxA21d,ymaxA21e,ymaxA21f,
                                       ymaxA26a,ymaxA26b,ymaxA26c,ymaxA26d,ymaxA26e,ymaxA26f,
                                       ymaxA31a,ymaxA31b,ymaxA31c,ymaxA31d,ymaxA31e,ymaxA31f,
                                       ymaxG21a,ymaxG21b,ymaxG21c,ymaxG21d,ymaxG21e,ymaxG21f,
                                       ymaxG26a,ymaxG26b,ymaxG26c,ymaxG26d,ymaxG26e,ymaxG26f,ymaxG26g,
                                       ymaxG31a,ymaxG31b,ymaxG31c,ymaxG31d,ymaxG31e,ymaxG31f,
                                       ymaxR21a,ymaxR21b,ymaxR21c,ymaxR21d,ymaxR21e,ymaxR21f,ymaxR21g,ymaxR21h,
                                       ymaxR26a,ymaxR26b,ymaxR26c,ymaxR26d,ymaxR26e,ymaxR26f,
                                       ymaxR31a,ymaxR31b,ymaxR31c,ymaxR31d,ymaxR31e,ymaxR31f,
                                       TmaxM,
                                       TmaxA,
                                       TmaxG,
                                       TmaxR,
                                       TminM,
                                       TminA,
                                       TminG,
                                       TminR,
                                       ToptM,
                                       ToptA,
                                       ToptG,
                                       ToptR,
                                       TempM21a,TempM21b,TempM21c,TempM21d,TempM21e,TempM21f,
                                       TempM26a,TempM26b,TempM26c,TempM26d,TempM26e,TempM26f,
                                       TempM31a,TempM31b,TempM31c,TempM31d,TempM31e,TempM31f,
                                       TempA21a,TempA21b,TempA21c,TempA21d,TempA21e,TempA21f,
                                       TempA26a,TempA26b,TempA26c,TempA26d,TempA26e,TempA26f,
                                       TempA31a,TempA31b,TempA31c,TempA31d,TempA31e,TempA31f,
                                       TempG21a,TempG21b,TempG21c,TempG21d,TempG21e,TempG21f,
                                       TempG26a,TempG26b,TempG26c,TempG26d,TempG26e,TempG26f,TempG26g,
                                       TempG31a,TempG31b,TempG31c,TempG31d,TempG31e,TempG31f,
                                       TempR21a,TempR21b,TempR21c,TempR21d,TempR21e,TempR21f,TempR21g,TempR21h,
                                       TempR26a,TempR26b,TempR26c,TempR26d,TempR26e,TempR26f,
                                       TempR31a,TempR31b,TempR31c,TempR31d,TempR31e,TempR31f,
                                       PhotodatM21a,PhotodatM21b,PhotodatM21c,PhotodatM21d,PhotodatM21e,PhotodatM21f,
                                       PhotodatM26a,PhotodatM26b,PhotodatM26c,PhotodatM26d,PhotodatM26e,PhotodatM26f,
                                       PhotodatM31a,PhotodatM31b,PhotodatM31c,PhotodatM31d,PhotodatM31e,PhotodatM31f,
                                       PhotodatA21a,PhotodatA21b,PhotodatA21c,PhotodatA21d,PhotodatA21e,PhotodatA21f,
                                       PhotodatA26a,PhotodatA26b,PhotodatA26c,PhotodatA26d,PhotodatA26e,PhotodatA26f,
                                       PhotodatA31a,PhotodatA31b,PhotodatA31c,PhotodatA31d,PhotodatA31e,PhotodatA31f,
                                       PhotodatG21a,PhotodatG21b,PhotodatG21c,PhotodatG21d,PhotodatG21e,PhotodatG21f,
                                       PhotodatG26a,PhotodatG26b,PhotodatG26c,PhotodatG26d,PhotodatG26e,PhotodatG26f,PhotodatG26g,
                                       PhotodatG31a,PhotodatG31b,PhotodatG31c,PhotodatG31d,PhotodatG31e,PhotodatG31f,
                                       PhotodatR21a,PhotodatR21b,PhotodatR21c,PhotodatR21d,PhotodatR21e,PhotodatR21f,PhotodatR21g,PhotodatR21h,
                                       PhotodatR26a,PhotodatR26b,PhotodatR26c,PhotodatR26d,PhotodatR26e,PhotodatR26f,
                                       PhotodatR31a,PhotodatR31b,PhotodatR31c,PhotodatR31d,PhotodatR31e,PhotodatR31f){
  PhotomeanM21a <- beta(ymaxM21a,TminM,ToptM,TmaxM,TempM21a)
  PhotomeanM21b <- beta(ymaxM21b,TminM,ToptM,TmaxM,TempM21b)
  PhotomeanM21c <- beta(ymaxM21c,TminM,ToptM,TmaxM,TempM21c)
  PhotomeanM21d <- beta(ymaxM21d,TminM,ToptM,TmaxM,TempM21d)
  PhotomeanM21e <- beta(ymaxM21e,TminM,ToptM,TmaxM,TempM21e)
  PhotomeanM21f <- beta(ymaxM21f,TminM,ToptM,TmaxM,TempM21f)
  PhotomeanM26a <- beta(ymaxM26a,TminM,ToptM,TmaxM,TempM26a)
  PhotomeanM26b <- beta(ymaxM26b,TminM,ToptM,TmaxM,TempM26b)
  PhotomeanM26c <- beta(ymaxM26c,TminM,ToptM,TmaxM,TempM26c)
  PhotomeanM26d <- beta(ymaxM26d,TminM,ToptM,TmaxM,TempM26d)
  PhotomeanM26e <- beta(ymaxM26e,TminM,ToptM,TmaxM,TempM26e)
  PhotomeanM26f <- beta(ymaxM26f,TminM,ToptM,TmaxM,TempM26f)
  PhotomeanM31a <- beta(ymaxM31a,TminM,ToptM,TmaxM,TempM31a)
  PhotomeanM31b <- beta(ymaxM31b,TminM,ToptM,TmaxM,TempM31b)
  PhotomeanM31c <- beta(ymaxM31c,TminM,ToptM,TmaxM,TempM31c)
  PhotomeanM31d <- beta(ymaxM31d,TminM,ToptM,TmaxM,TempM31d)
  PhotomeanM31e <- beta(ymaxM31e,TminM,ToptM,TmaxM,TempM31e)
  PhotomeanM31f <- beta(ymaxM31f,TminM,ToptM,TmaxM,TempM31f)
  PhotomeanA21a <- beta(ymaxA21a,TminA,ToptA,TmaxA,TempA21a)
  PhotomeanA21b <- beta(ymaxA21b,TminA,ToptA,TmaxA,TempA21b)
  PhotomeanA21c <- beta(ymaxA21c,TminA,ToptA,TmaxA,TempA21c)
  PhotomeanA21d <- beta(ymaxA21d,TminA,ToptA,TmaxA,TempA21d)
  PhotomeanA21e <- beta(ymaxA21e,TminA,ToptA,TmaxA,TempA21e)
  PhotomeanA21f <- beta(ymaxA21f,TminA,ToptA,TmaxA,TempA21f)
  PhotomeanA26a <- beta(ymaxA26a,TminA,ToptA,TmaxA,TempA26a)
  PhotomeanA26b <- beta(ymaxA26b,TminA,ToptA,TmaxA,TempA26b)
  PhotomeanA26c <- beta(ymaxA26c,TminA,ToptA,TmaxA,TempA26c)
  PhotomeanA26d <- beta(ymaxA26d,TminA,ToptA,TmaxA,TempA26d)
  PhotomeanA26e <- beta(ymaxA26e,TminA,ToptA,TmaxA,TempA26e)
  PhotomeanA26f <- beta(ymaxA26f,TminA,ToptA,TmaxA,TempA26f)
  PhotomeanA31a <- beta(ymaxA31a,TminA,ToptA,TmaxA,TempA31a)
  PhotomeanA31b <- beta(ymaxA31b,TminA,ToptA,TmaxA,TempA31b)
  PhotomeanA31c <- beta(ymaxA31c,TminA,ToptA,TmaxA,TempA31c)
  PhotomeanA31d <- beta(ymaxA31d,TminA,ToptA,TmaxA,TempA31d)
  PhotomeanA31e <- beta(ymaxA31e,TminA,ToptA,TmaxA,TempA31e)
  PhotomeanA31f <- beta(ymaxA31f,TminA,ToptA,TmaxA,TempA31f)
  PhotomeanG21a <- beta(ymaxG21a,TminG,ToptG,TmaxG,TempG21a)
  PhotomeanG21b <- beta(ymaxG21b,TminG,ToptG,TmaxG,TempG21b)
  PhotomeanG21c <- beta(ymaxG21c,TminG,ToptG,TmaxG,TempG21c)
  PhotomeanG21d <- beta(ymaxG21d,TminG,ToptG,TmaxG,TempG21d)
  PhotomeanG21e <- beta(ymaxG21e,TminG,ToptG,TmaxG,TempG21e)
  PhotomeanG21f <- beta(ymaxG21f,TminG,ToptG,TmaxG,TempG21f)
  PhotomeanG26a <- beta(ymaxG26a,TminG,ToptG,TmaxG,TempG26a)
  PhotomeanG26b <- beta(ymaxG26b,TminG,ToptG,TmaxG,TempG26b)
  PhotomeanG26c <- beta(ymaxG26c,TminG,ToptG,TmaxG,TempG26c)
  PhotomeanG26d <- beta(ymaxG26d,TminG,ToptG,TmaxG,TempG26d)
  PhotomeanG26e <- beta(ymaxG26e,TminG,ToptG,TmaxG,TempG26e)
  PhotomeanG26f <- beta(ymaxG26f,TminG,ToptG,TmaxG,TempG26f)
  PhotomeanG26g <- beta(ymaxG26g,TminG,ToptG,TmaxG,TempG26g)
  PhotomeanG31a <- beta(ymaxG31a,TminG,ToptG,TmaxG,TempG31a)
  PhotomeanG31b <- beta(ymaxG31b,TminG,ToptG,TmaxG,TempG31b)
  PhotomeanG31c <- beta(ymaxG31c,TminG,ToptG,TmaxG,TempG31c)
  PhotomeanG31d <- beta(ymaxG31d,TminG,ToptG,TmaxG,TempG31d)
  PhotomeanG31e <- beta(ymaxG31e,TminG,ToptG,TmaxG,TempG31e)
  PhotomeanG31f <- beta(ymaxG31f,TminG,ToptG,TmaxG,TempG31f)
  PhotomeanR21a <- beta(ymaxR21a,TminR,ToptR,TmaxR,TempR21a)
  PhotomeanR21b <- beta(ymaxR21b,TminR,ToptR,TmaxR,TempR21b)
  PhotomeanR21c <- beta(ymaxR21c,TminR,ToptR,TmaxR,TempR21c)
  PhotomeanR21d <- beta(ymaxR21d,TminR,ToptR,TmaxR,TempR21d)
  PhotomeanR21e <- beta(ymaxR21e,TminR,ToptR,TmaxR,TempR21e)
  PhotomeanR21f <- beta(ymaxR21f,TminR,ToptR,TmaxR,TempR21f)
  PhotomeanR21g <- beta(ymaxR21g,TminR,ToptR,TmaxR,TempR21g)
  PhotomeanR21h <- beta(ymaxR21h,TminR,ToptR,TmaxR,TempR21h)
  PhotomeanR26a <- beta(ymaxR26a,TminR,ToptR,TmaxR,TempR26a)
  PhotomeanR26b <- beta(ymaxR26b,TminR,ToptR,TmaxR,TempR26b)
  PhotomeanR26c <- beta(ymaxR26c,TminR,ToptR,TmaxR,TempR26c)
  PhotomeanR26d <- beta(ymaxR26d,TminR,ToptR,TmaxR,TempR26d)
  PhotomeanR26e <- beta(ymaxR26e,TminR,ToptR,TmaxR,TempR26e)
  PhotomeanR26f <- beta(ymaxR26f,TminR,ToptR,TmaxR,TempR26f)
  PhotomeanR31a <- beta(ymaxR31a,TminR,ToptR,TmaxR,TempR31a)
  PhotomeanR31b <- beta(ymaxR31b,TminR,ToptR,TmaxR,TempR31b)
  PhotomeanR31c <- beta(ymaxR31c,TminR,ToptR,TmaxR,TempR31c)
  PhotomeanR31d <- beta(ymaxR31d,TminR,ToptR,TmaxR,TempR31d)
  PhotomeanR31e <- beta(ymaxR31e,TminR,ToptR,TmaxR,TempR31e)
  PhotomeanR31f <- beta(ymaxR31f,TminR,ToptR,TmaxR,TempR31f)
  -(sum(dnorm(PhotodatM21a,mean=PhotomeanM21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21b,mean=PhotomeanM21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21c,mean=PhotomeanM21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM21d,mean=PhotomeanM21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21e,mean=PhotomeanM21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM21f,mean=PhotomeanM21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26a,mean=PhotomeanM26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26b,mean=PhotomeanM26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26c,mean=PhotomeanM26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM26d,mean=PhotomeanM26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26e,mean=PhotomeanM26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM26f,mean=PhotomeanM26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31a,mean=PhotomeanM31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31b,mean=PhotomeanM31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31c,mean=PhotomeanM31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatM31d,mean=PhotomeanM31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31e,mean=PhotomeanM31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatM31f,mean=PhotomeanM31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21a,mean=PhotomeanA21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21b,mean=PhotomeanA21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21c,mean=PhotomeanA21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA21d,mean=PhotomeanA21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21e,mean=PhotomeanA21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA21f,mean=PhotomeanA21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26a,mean=PhotomeanA26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26b,mean=PhotomeanA26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26c,mean=PhotomeanA26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA26d,mean=PhotomeanA26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26e,mean=PhotomeanA26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA26f,mean=PhotomeanA26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31a,mean=PhotomeanA31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31b,mean=PhotomeanA31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31c,mean=PhotomeanA31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatA31d,mean=PhotomeanA31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31e,mean=PhotomeanA31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatA31f,mean=PhotomeanA31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21a,mean=PhotomeanG21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21b,mean=PhotomeanG21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21c,mean=PhotomeanG21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG21d,mean=PhotomeanG21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21e,mean=PhotomeanG21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG21f,mean=PhotomeanG21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26a,mean=PhotomeanG26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26b,mean=PhotomeanG26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26c,mean=PhotomeanG26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26d,mean=PhotomeanG26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26e,mean=PhotomeanG26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG26f,mean=PhotomeanG26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG26g,mean=PhotomeanG26g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31a,mean=PhotomeanG31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31b,mean=PhotomeanG31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31c,mean=PhotomeanG31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatG31d,mean=PhotomeanG31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31e,mean=PhotomeanG31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatG31f,mean=PhotomeanG31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21a,mean=PhotomeanR21a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21b,mean=PhotomeanR21b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21c,mean=PhotomeanR21c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21d,mean=PhotomeanR21d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21e,mean=PhotomeanR21e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21f,mean=PhotomeanR21f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR21g,mean=PhotomeanR21g,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR21h,mean=PhotomeanR21h,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26a,mean=PhotomeanR26a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26b,mean=PhotomeanR26b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26c,mean=PhotomeanR26c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR26d,mean=PhotomeanR26d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26e,mean=PhotomeanR26e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR26f,mean=PhotomeanR26f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31a,mean=PhotomeanR31a,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31b,mean=PhotomeanR31b,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31c,mean=PhotomeanR31c,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) +
      sum(dnorm(PhotodatR31d,mean=PhotomeanR31d,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31e,mean=PhotomeanR31e,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(PhotodatR31f,mean=PhotomeanR31f,sd=exp(sdPhoto),log=TRUE),na.rm=TRUE))
}

###
#Maximum likelihood fits
###

###
#A275
###

#Acclimation of all parameters
fit_Photo_beta_all_lin.all_A275 <- mle2(Photo_beta_normNLL_all_lin.all,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                             ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                             ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                             ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                             ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                             ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                             ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                             ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                             ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                             ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                             ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                             ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                             aM=0,
                                                                             aA=0,
                                                                             aG=0,
                                                                             aR=0,
                                                                             bM=0.05,
                                                                             bA=0.05,
                                                                             bG=0.05,
                                                                             bR=0.05,
                                                                             cM=26,
                                                                             cA=18.5,
                                                                             cG=14,
                                                                             cR=24.5,
                                                                             dM=0.2,
                                                                             dA=0.4,
                                                                             dG=0.6,
                                                                             dR=0.2,
                                                                             eM=41,
                                                                             eA=41,
                                                                             eG=41,
                                                                             eR=41,
                                                                             fM=0.05,
                                                                             fA=0.07,
                                                                             fG=0.05,
                                                                             fR=0.05),
                                   data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                             TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:16],TempM26e=M26.A275.dat$Temp[17:20],TempM26f=M26.A275.dat$Temp[21:24],
                                             TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                             TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                             TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:16],TempA26e=A26.A275.dat$Temp[17:20],TempA26f=A26.A275.dat$Temp[21:24],
                                             TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:8],TempA31c=A31.A275.dat$Temp[9:13],TempA31d=A31.A275.dat$Temp[14:16],TempA31e=A31.A275.dat$Temp[17:21],TempA31f=A31.A275.dat$Temp[22:24],
                                             TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                             TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:28],
                                             TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:15],TempG31e=G31.A275.dat$Temp[16:21],TempG31f=G31.A275.dat$Temp[22:24],
                                             TempR21a=R21.A275.dat$Temp[1:3],TempR21b=R21.A275.dat$Temp[4:8],TempR21c=R21.A275.dat$Temp[9:11],TempR21d=R21.A275.dat$Temp[12:16],TempR21e=R21.A275.dat$Temp[17:19],TempR21f=R21.A275.dat$Temp[20:24],TempR21g=R21.A275.dat$Temp[25:27],TempR21h=R21.A275.dat$Temp[28:32],
                                             TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                             TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                             PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                             PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:16],PhotodatM26e=M26.A275.dat$A275[17:20],PhotodatM26f=M26.A275.dat$A275[21:24],
                                             PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                             PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                             PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:16],PhotodatA26e=A26.A275.dat$A275[17:20],PhotodatA26f=A26.A275.dat$A275[21:24],
                                             PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:8],PhotodatA31c=A31.A275.dat$A275[9:13],PhotodatA31d=A31.A275.dat$A275[14:16],PhotodatA31e=A31.A275.dat$A275[17:21],PhotodatA31f=A31.A275.dat$A275[22:24],
                                             PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                             PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:28],
                                             PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:15],PhotodatG31e=G31.A275.dat$A275[16:21],PhotodatG31f=G31.A275.dat$A275[22:24],
                                             PhotodatR21a=R21.A275.dat$A275[1:3],PhotodatR21b=R21.A275.dat$A275[4:8],PhotodatR21c=R21.A275.dat$A275[9:11],PhotodatR21d=R21.A275.dat$A275[12:16],PhotodatR21e=R21.A275.dat$A275[17:19],PhotodatR21f=R21.A275.dat$A275[20:24],PhotodatR21g=R21.A275.dat$A275[25:27],PhotodatR21h=R21.A275.dat$A275[28:32],
                                             PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                             PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                   control=list(maxit=20000))
summary(fit_Photo_beta_all_lin.all_A275)
#Did not converge (unrealistic estimates)

#Acclimation of Tmin and Topt
fit_Photo_beta_all_Tmin.Topt.lin_A275 <- mle2(Photo_beta_normNLL_all_Tmin.Topt.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                              ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                              ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                              ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                              ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                              ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                              ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                              ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                              ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                              ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                              ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                              ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                              TmaxM=42,
                                                                                              TmaxA=43,
                                                                                              TmaxG=43,
                                                                                              TmaxR=42,
                                                                                              aM=21,
                                                                                              aA=24,
                                                                                              aG=-1,
                                                                                              aR=-23,
                                                                                              bM=-0.9,
                                                                                              bA=-1.2,
                                                                                              bG=0.5,
                                                                                              bR=1.1,
                                                                                              cM=22.8,
                                                                                              cA=19.1,
                                                                                              cG=12.9,
                                                                                              cR=22.9,
                                                                                              dM=0.24,
                                                                                              dA=0.34,
                                                                                              dG=0.64,
                                                                                              dR=0.19),
                                              data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                                        TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                                        TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                                        TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                                        TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                                        TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                                        TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                                        TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                                        TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                                        TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                                        TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                                        TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                                        PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                                        PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                                        PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                                        PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                                        PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                                        PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                                        PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                                        PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                                        PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                                        PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                                        PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                                        PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                              control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmin.Topt.lin_A275)
#Did not converge (unrealistic estimates)

#Acclimation of Tmin and Tmax
fit_Photo_beta_all_Tmin.Tmax.lin_A275 <- mle2(Photo_beta_normNLL_all_Tmin.Tmax.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                              ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                              ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                              ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                              ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                              ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                              ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                              ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                              ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                              ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                              ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                              ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                              ToptM=28,
                                                                                              ToptA=28,
                                                                                              ToptG=28,
                                                                                              ToptR=28,
                                                                                              aM=21,
                                                                                              aA=24,
                                                                                              aG=-1,
                                                                                              aR=-23,
                                                                                              bM=-0.9,
                                                                                              bA=-1.2,
                                                                                              bG=0.5,
                                                                                              bR=1.1,
                                                                                              cM=41,
                                                                                              cA=41,
                                                                                              cG=41,
                                                                                              cR=41,
                                                                                              dM=0.05,
                                                                                              dA=0.07,
                                                                                              dG=0.05,
                                                                                              dR=0.05),
                                              data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                                        TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                                        TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                                        TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                                        TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                                        TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                                        TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                                        TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                                        TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                                        TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                                        TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                                        TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                                        PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                                        PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                                        PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                                        PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                                        PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                                        PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                                        PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                                        PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                                        PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                                        PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                                        PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                                        PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                              control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmin.Tmax.lin_A275)
#Did not converge (unrealistic estimates)

#Acclimation of Topt and Tmax
fit_Photo_beta_all_Topt.Tmax.lin_A275 <- mle2(Photo_beta_normNLL_all_Topt.Tmax.lin,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=7,ymaxM21e=7,ymaxM21f=9,
                                                                                              ymaxM26a=8,ymaxM26b=6,ymaxM26c=5.5,ymaxM26d=5,ymaxM26e=11,ymaxM26f=6,
                                                                                              ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=4,ymaxM31e=16,ymaxM31f=12,
                                                                                              ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                              ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                              ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                              ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                              ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                              ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                              ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                              ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                              ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                              TminM=1,
                                                                                              TminA=1,
                                                                                              TminG=8,
                                                                                              TminR=4,
                                                                                              aM=22.7,
                                                                                              aA=16.7,
                                                                                              aG=12.6,
                                                                                              aR=24.5,
                                                                                              bM=0.3,
                                                                                              bA=0.52,
                                                                                              bG=0.61,
                                                                                              bR=0.14,
                                                                                              cM=41,
                                                                                              cA=41,
                                                                                              cG=41,
                                                                                              cR=41,
                                                                                              dM=0.05,
                                                                                              dA=0.07,
                                                                                              dG=0.05,
                                                                                              dR=0.05),
                                              data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                                        TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                                        TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                                        TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                                        TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                                        TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                                        TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                                        TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                                        TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                                        TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                                        TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                                        TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                                        PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                                        PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                                        PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                                        PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                                        PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                                        PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                                        PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                                        PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                                        PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                                        PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                                        PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                                        PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                              control=list(maxit=20000))
summary(fit_Photo_beta_all_Topt.Tmax.lin_A275)

#Acclimation of Tmin
fit_Photo_beta_all_Tmin.lin_A275 <- mle2(Photo_beta_normNLL_all_Tmin.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                    ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                    ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                    ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                    ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                    ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                    ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                    ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                    ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                    ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                    ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                    ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                    TmaxM=42,
                                                                                    TmaxA=43,
                                                                                    TmaxG=43,
                                                                                    TmaxR=42,
                                                                                    aM=5,
                                                                                    aA=24,
                                                                                    aG=-1,
                                                                                    aR=-23,
                                                                                    bM=-0.1,
                                                                                    bA=-1.2,
                                                                                    bG=0.5,
                                                                                    bR=1.1,
                                                                                    ToptM=29,
                                                                                    ToptA=29,
                                                                                    ToptG=28,
                                                                                    ToptR=28),
                                         data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                                   TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                                   TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                                   TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                                   TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                                   TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                                   TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                                   TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                                   TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                                   TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                                   TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                                   TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                                   PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                                   PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                                   PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                                   PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                                   PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                                   PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                                   PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                                   PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                                   PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                                   PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                                   PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                                   PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                         control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmin.lin_A275)
#Did not converge (unrealistic estimates)

#Acclimation of Topt
fit_Photo_beta_all_Topt.lin_A275 <- mle2(Photo_beta_normNLL_all_Topt.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                    ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                    ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                    ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                    ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                    ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                    ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                    ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                    ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                    ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                    ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                    ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                    TmaxM=42,
                                                                                    TmaxA=43,
                                                                                    TmaxG=42,
                                                                                    TmaxR=42,
                                                                                    TminM=1,
                                                                                    TminA=-3,
                                                                                    TminG=8,
                                                                                    TminR=4,
                                                                                    aM=22.7,
                                                                                    aA=16.7,
                                                                                    aG=12.6,
                                                                                    aR=24.5,
                                                                                    bM=0.3,
                                                                                    bA=0.52,
                                                                                    bG=0.61,
                                                                                    bR=0.14),
                                         data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                                   TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                                   TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                                   TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                                   TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                                   TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                                   TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                                   TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                                   TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                                   TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                                   TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                                   TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                                   PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                                   PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                                   PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                                   PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                                   PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                                   PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                                   PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                                   PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                                   PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                                   PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                                   PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                                   PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                         control=list(maxit=20000))
summary(fit_Photo_beta_all_Topt.lin_A275)

#Acclimation of Tmax
fit_Photo_beta_all_Tmax.lin_A275 <- mle2(Photo_beta_normNLL_all_Tmax.lin,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=7,ymaxM21e=7,ymaxM21f=9,
                                                                                    ymaxM26a=8,ymaxM26b=6,ymaxM26c=5.5,ymaxM26d=5,ymaxM26e=11,ymaxM26f=6,
                                                                                    ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=4,ymaxM31e=16,ymaxM31f=12,
                                                                                    ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                    ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                    ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                    ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                    ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                    ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                    ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                    ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                    ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                    aM=39,
                                                                                    aA=40,
                                                                                    aG=41,
                                                                                    aR=44,
                                                                                    bM=0.1,
                                                                                    bA=0.1,
                                                                                    bG=0.05,
                                                                                    bR=-0.1,
                                                                                    TminM=2,
                                                                                    TminA=0.2,
                                                                                    TminG=8,
                                                                                    TminR=4,
                                                                                    ToptM=28,
                                                                                    ToptA=28,
                                                                                    ToptG=28,
                                                                                    ToptR=28),
                                         data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                                   TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                                   TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                                   TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                                   TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                                   TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                                   TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                                   TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                                   TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                                   TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                                   TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                                   TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                                   PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                                   PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                                   PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                                   PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                                   PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                                   PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                                   PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                                   PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                                   PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                                   PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                                   PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                                   PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                         control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmax.lin_A275)

#No acclimation
fit_Photo_beta_all_bin_A275 <- mle2(Photo_beta_normNLL_all_bin,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=7,ymaxM21e=7,ymaxM21f=9,
                                                                          ymaxM26a=8,ymaxM26b=6,ymaxM26c=5.5,ymaxM26d=5,ymaxM26e=11,ymaxM26f=6,
                                                                          ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=4,ymaxM31e=16,ymaxM31f=12,
                                                                          ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                          ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                          ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                          ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                          ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                          ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                          ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                          ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                          ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                          TmaxM=42,
                                                                          TmaxA=42,
                                                                          TmaxG=43,
                                                                          TmaxR=42,
                                                                          TminM=2,
                                                                          TminA=0.2,
                                                                          TminG=8,
                                                                          TminR=4,
                                                                          ToptM=28,
                                                                          ToptA=28,
                                                                          ToptG=28,
                                                                          ToptR=28),
                                    data=list(TempM21a=M21.A275.dat$Temp[1:3],TempM21b=M21.A275.dat$Temp[4:8],TempM21c=M21.A275.dat$Temp[9:11],TempM21d=M21.A275.dat$Temp[12:16],TempM21e=M21.A275.dat$Temp[17:19],TempM21f=M21.A275.dat$Temp[20:24],
                                              TempM26a=M26.A275.dat$Temp[1:4],TempM26b=M26.A275.dat$Temp[5:8],TempM26c=M26.A275.dat$Temp[9:12],TempM26d=M26.A275.dat$Temp[13:15],TempM26e=M26.A275.dat$Temp[16:19],TempM26f=M26.A275.dat$Temp[20:23],
                                              TempM31a=M31.A275.dat$Temp[1:5],TempM31b=M31.A275.dat$Temp[6:8],TempM31c=M31.A275.dat$Temp[9:12],TempM31d=M31.A275.dat$Temp[13:15],TempM31e=M31.A275.dat$Temp[16:20],TempM31f=M31.A275.dat$Temp[21:23],
                                              TempA21a=A21.A275.dat$Temp[1:3],TempA21b=A21.A275.dat$Temp[4:8],TempA21c=A21.A275.dat$Temp[9:11],TempA21d=A21.A275.dat$Temp[12:16],TempA21e=A21.A275.dat$Temp[17:19],TempA21f=A21.A275.dat$Temp[20:24],
                                              TempA26a=A26.A275.dat$Temp[1:4],TempA26b=A26.A275.dat$Temp[5:8],TempA26c=A26.A275.dat$Temp[9:12],TempA26d=A26.A275.dat$Temp[13:15],TempA26e=A26.A275.dat$Temp[16:19],TempA26f=A26.A275.dat$Temp[20:22],
                                              TempA31a=A31.A275.dat$Temp[1:5],TempA31b=A31.A275.dat$Temp[6:7],TempA31c=A31.A275.dat$Temp[8:12],TempA31d=A31.A275.dat$Temp[13:15],TempA31e=A31.A275.dat$Temp[16:20],TempA31f=A31.A275.dat$Temp[21:23],
                                              TempG21a=G21.A275.dat$Temp[1:3],TempG21b=G21.A275.dat$Temp[4:8],TempG21c=G21.A275.dat$Temp[9:11],TempG21d=G21.A275.dat$Temp[12:16],TempG21e=G21.A275.dat$Temp[17:19],TempG21f=G21.A275.dat$Temp[20:24],
                                              TempG26a=G26.A275.dat$Temp[1:4],TempG26b=G26.A275.dat$Temp[5:8],TempG26c=G26.A275.dat$Temp[9:12],TempG26d=G26.A275.dat$Temp[13:16],TempG26e=G26.A275.dat$Temp[17:20],TempG26f=G26.A275.dat$Temp[21:24],TempG26g=G26.A275.dat$Temp[25:27],
                                              TempG31a=G31.A275.dat$Temp[1:4],TempG31b=G31.A275.dat$Temp[5:7],TempG31c=G31.A275.dat$Temp[8:12],TempG31d=G31.A275.dat$Temp[13:14],TempG31e=G31.A275.dat$Temp[15:20],TempG31f=G31.A275.dat$Temp[21:23],
                                              TempR21a=R21.A275.dat$Temp[1:2],TempR21b=R21.A275.dat$Temp[3:6],TempR21c=R21.A275.dat$Temp[7:9],TempR21d=R21.A275.dat$Temp[10:14],TempR21e=R21.A275.dat$Temp[15:17],TempR21f=R21.A275.dat$Temp[18:22],TempR21g=R21.A275.dat$Temp[23:25],TempR21h=R21.A275.dat$Temp[26:30],
                                              TempR26a=R26.A275.dat$Temp[1:4],TempR26b=R26.A275.dat$Temp[5:8],TempR26c=R26.A275.dat$Temp[9:13],TempR26d=R26.A275.dat$Temp[14:17],TempR26e=R26.A275.dat$Temp[18:21],TempR26f=R26.A275.dat$Temp[22:25],
                                              TempR31a=R31.A275.dat$Temp[1:5],TempR31b=R31.A275.dat$Temp[6:8],TempR31c=R31.A275.dat$Temp[9:13],TempR31d=R31.A275.dat$Temp[14:16],TempR31e=R31.A275.dat$Temp[17:21],TempR31f=R31.A275.dat$Temp[22:24],
                                              PhotodatM21a=M21.A275.dat$A275[1:3],PhotodatM21b=M21.A275.dat$A275[4:8],PhotodatM21c=M21.A275.dat$A275[9:11],PhotodatM21d=M21.A275.dat$A275[12:16],PhotodatM21e=M21.A275.dat$A275[17:19],PhotodatM21f=M21.A275.dat$A275[20:24],
                                              PhotodatM26a=M26.A275.dat$A275[1:4],PhotodatM26b=M26.A275.dat$A275[5:8],PhotodatM26c=M26.A275.dat$A275[9:12],PhotodatM26d=M26.A275.dat$A275[13:15],PhotodatM26e=M26.A275.dat$A275[16:19],PhotodatM26f=M26.A275.dat$A275[20:23],
                                              PhotodatM31a=M31.A275.dat$A275[1:5],PhotodatM31b=M31.A275.dat$A275[6:8],PhotodatM31c=M31.A275.dat$A275[9:12],PhotodatM31d=M31.A275.dat$A275[13:15],PhotodatM31e=M31.A275.dat$A275[16:20],PhotodatM31f=M31.A275.dat$A275[21:23],
                                              PhotodatA21a=A21.A275.dat$A275[1:3],PhotodatA21b=A21.A275.dat$A275[4:8],PhotodatA21c=A21.A275.dat$A275[9:11],PhotodatA21d=A21.A275.dat$A275[12:16],PhotodatA21e=A21.A275.dat$A275[17:19],PhotodatA21f=A21.A275.dat$A275[20:24],
                                              PhotodatA26a=A26.A275.dat$A275[1:4],PhotodatA26b=A26.A275.dat$A275[5:8],PhotodatA26c=A26.A275.dat$A275[9:12],PhotodatA26d=A26.A275.dat$A275[13:15],PhotodatA26e=A26.A275.dat$A275[16:19],PhotodatA26f=A26.A275.dat$A275[20:22],
                                              PhotodatA31a=A31.A275.dat$A275[1:5],PhotodatA31b=A31.A275.dat$A275[6:7],PhotodatA31c=A31.A275.dat$A275[8:12],PhotodatA31d=A31.A275.dat$A275[13:15],PhotodatA31e=A31.A275.dat$A275[16:20],PhotodatA31f=A31.A275.dat$A275[21:23],
                                              PhotodatG21a=G21.A275.dat$A275[1:3],PhotodatG21b=G21.A275.dat$A275[4:8],PhotodatG21c=G21.A275.dat$A275[9:11],PhotodatG21d=G21.A275.dat$A275[12:16],PhotodatG21e=G21.A275.dat$A275[17:19],PhotodatG21f=G21.A275.dat$A275[20:24],
                                              PhotodatG26a=G26.A275.dat$A275[1:4],PhotodatG26b=G26.A275.dat$A275[5:8],PhotodatG26c=G26.A275.dat$A275[9:12],PhotodatG26d=G26.A275.dat$A275[13:16],PhotodatG26e=G26.A275.dat$A275[17:20],PhotodatG26f=G26.A275.dat$A275[21:24],PhotodatG26g=G26.A275.dat$A275[25:27],
                                              PhotodatG31a=G31.A275.dat$A275[1:4],PhotodatG31b=G31.A275.dat$A275[5:7],PhotodatG31c=G31.A275.dat$A275[8:12],PhotodatG31d=G31.A275.dat$A275[13:14],PhotodatG31e=G31.A275.dat$A275[15:20],PhotodatG31f=G31.A275.dat$A275[21:23],
                                              PhotodatR21a=R21.A275.dat$A275[1:2],PhotodatR21b=R21.A275.dat$A275[3:6],PhotodatR21c=R21.A275.dat$A275[7:9],PhotodatR21d=R21.A275.dat$A275[10:14],PhotodatR21e=R21.A275.dat$A275[15:17],PhotodatR21f=R21.A275.dat$A275[18:22],PhotodatR21g=R21.A275.dat$A275[23:25],PhotodatR21h=R21.A275.dat$A275[26:30],
                                              PhotodatR26a=R26.A275.dat$A275[1:4],PhotodatR26b=R26.A275.dat$A275[5:8],PhotodatR26c=R26.A275.dat$A275[9:13],PhotodatR26d=R26.A275.dat$A275[14:17],PhotodatR26e=R26.A275.dat$A275[18:21],PhotodatR26f=R26.A275.dat$A275[22:25],
                                              PhotodatR31a=R31.A275.dat$A275[1:5],PhotodatR31b=R31.A275.dat$A275[6:8],PhotodatR31c=R31.A275.dat$A275[9:13],PhotodatR31d=R31.A275.dat$A275[14:16],PhotodatR31e=R31.A275.dat$A275[17:21],PhotodatR31f=R31.A275.dat$A275[22:24]),
                                    control=list(maxit=20000))
summary(fit_Photo_beta_all_bin_A275)

#Delta AIC
AICctab(fit_Photo_beta_all_Topt.Tmax.lin_A275,fit_Photo_beta_all_Topt.lin_A275,
        fit_Photo_beta_all_Tmax.lin_A275,fit_Photo_beta_all_bin_A275,nobs=292)

###
#Asat
###

#Acclimation of all parameters
fit_Photo_beta_all_lin.all_Asat <- mle2(Photo_beta_normNLL_all_lin.all,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                aM=21,
                                                                                aA=24,
                                                                                aG=-1,
                                                                                aR=-23,
                                                                                bM=-0.9,
                                                                                bA=-1.2,
                                                                                bG=0.5,
                                                                                bR=1.1,
                                                                                cM=26,
                                                                                cA=18.5,
                                                                                cG=14,
                                                                                cR=24.5,
                                                                                dM=0.2,
                                                                                dA=0.4,
                                                                                dG=0.6,
                                                                                dR=0.2,
                                                                                eM=41,
                                                                                eA=41,
                                                                                eG=41,
                                                                                eR=41,
                                                                                fM=0.05,
                                                                                fA=0.07,
                                                                                fG=0.05,
                                                                                fR=0.05),
                                       data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                 TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                 TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                 TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                 TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                 TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                 TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                 TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                 TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                 TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                 TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                 TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                 PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                 PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                 PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                 PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                 PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                 PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                 PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                 PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                 PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                 PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                 PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                 PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                       control=list(maxit=20000))
summary(fit_Photo_beta_all_lin.all_Asat)
#Did not converge (unrealistic estimates)

#Acclimation of Tmin and Topt
fit_Photo_beta_all_Tmin.Topt.lin_Asat <- mle2(Photo_beta_normNLL_all_Tmin.Topt.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                            ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                            ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                            ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                            ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                            ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                            ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                            ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                            ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                            ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                            ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                            ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                            TmaxM=42,
                                                                                            TmaxA=43,
                                                                                            TmaxG=43,
                                                                                            TmaxR=42,
                                                                                            aM=0,
                                                                                            aA=0,
                                                                                            aG=0,
                                                                                            aR=0,
                                                                                            bM=0.05,
                                                                                            bA=0.05,
                                                                                            bG=0.05,
                                                                                            bR=0.05,
                                                                                            cM=26,
                                                                                            cA=18.5,
                                                                                            cG=14,
                                                                                            cR=24.5,
                                                                                            dM=0.2,
                                                                                            dA=0.4,
                                                                                            dG=0.6,
                                                                                            dR=0.2),
                                             data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                       TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                       TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                       TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                       TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                       TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                       TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                       TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                       TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                       TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                       TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                       TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                       PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                       PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                       PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                       PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                       PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                       PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                       PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                       PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                       PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                       PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                       PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                       PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                             control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmin.Topt.lin_Asat)

#Acclimation of Tmin and Tmax
fit_Photo_beta_all_Tmin.Tmax.lin_Asat <- mle2(Photo_beta_normNLL_all_Tmin.Tmax.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                            ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                            ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                            ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                            ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                            ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                            ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                            ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                            ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                            ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                            ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                            ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                            aM=21,
                                                                                            aA=24,
                                                                                            aG=-1,
                                                                                            aR=-23,
                                                                                            bM=-0.9,
                                                                                            bA=-1.2,
                                                                                            bG=0.5,
                                                                                            bR=1.1,
                                                                                            ToptM=29,
                                                                                            ToptA=29,
                                                                                            ToptG=28,
                                                                                            ToptR=28,
                                                                                            cM=41,
                                                                                            cA=41,
                                                                                            cG=41,
                                                                                            cR=41,
                                                                                            dM=0.05,
                                                                                            dA=0.07,
                                                                                            dG=0.05,
                                                                                            dR=0.05),
                                             data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                       TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                       TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                       TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                       TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                       TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                       TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                       TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                       TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                       TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                       TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                       TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                       PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                       PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                       PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                       PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                       PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                       PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                       PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                       PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                       PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                       PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                       PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                       PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                             control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmin.Tmax.lin_Asat)

#Acclimation of Topt and Tmax
fit_Photo_beta_all_Topt.Tmax.lin_Asat <- mle2(Photo_beta_normNLL_all_Topt.Tmax.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                            ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                            ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                            ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                            ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                            ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                            ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                            ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                            ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                            ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                            ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                            ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                            TminM=1,
                                                                                            TminA=-3,
                                                                                            TminG=8,
                                                                                            TminR=4,
                                                                                            aM=22.7,
                                                                                            aA=16.7,
                                                                                            aG=12.6,
                                                                                            aR=24.5,
                                                                                            bM=0.3,
                                                                                            bA=0.52,
                                                                                            bG=0.61,
                                                                                            bR=0.14,
                                                                                            cM=41,
                                                                                            cA=41,
                                                                                            cG=41,
                                                                                            cR=41,
                                                                                            dM=0.05,
                                                                                            dA=0.07,
                                                                                            dG=0.05,
                                                                                            dR=0.05),
                                             data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                       TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                       TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                       TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                       TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                       TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                       TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                       TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                       TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                       TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                       TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                       TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                       PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                       PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                       PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                       PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                       PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                       PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                       PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                       PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                       PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                       PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                       PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                       PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                             control=list(maxit=20000))
summary(fit_Photo_beta_all_Topt.Tmax.lin_Asat)

#Acclimation of Tmin
fit_Photo_beta_all_Tmin.lin_Asat <- mle2(Photo_beta_normNLL_all_Tmin.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                  ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                  ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                  ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                  ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                  ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                  ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                  ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                  ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                  ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                  ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                  ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                  TmaxM=42,
                                                                                  TmaxA=43,
                                                                                  TmaxG=43,
                                                                                  TmaxR=42,
                                                                                  aM=21,
                                                                                  aA=24,
                                                                                  aG=-1,
                                                                                  aR=-23,
                                                                                  bM=-0.9,
                                                                                  bA=-1.2,
                                                                                  bG=0.5,
                                                                                  bR=1.1,
                                                                                  ToptM=29,
                                                                                  ToptA=29,
                                                                                  ToptG=28,
                                                                                  ToptR=28),
                                        data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                  TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                  TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                  TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                  TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                  TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                  TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                  TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                  TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                  TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                  TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                  TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                  PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                  PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                  PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                  PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                  PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                  PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                  PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                  PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                  PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                  PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                  PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                  PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                        control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmin.lin_Asat)

#Acclimation of Topt
fit_Photo_beta_all_Topt.lin_Asat <- mle2(Photo_beta_normNLL_all_Topt.lin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                                  ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                                  ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                                  ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                  ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                  ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                  ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                  ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                  ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                  ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                  ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                  ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                  TmaxM=42,
                                                                                  TmaxA=43,
                                                                                  TmaxG=42,
                                                                                  TmaxR=42,
                                                                                  TminM=1,
                                                                                  TminA=-3,
                                                                                  TminG=8,
                                                                                  TminR=4,
                                                                                  aM=22.7,
                                                                                  aA=16.7,
                                                                                  aG=12.6,
                                                                                  aR=24.5,
                                                                                  bM=0.3,
                                                                                  bA=0.52,
                                                                                  bG=0.61,
                                                                                  bR=0.14),
                                        data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                  TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                  TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                  TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                  TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                  TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                  TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                  TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                  TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                  TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                  TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                  TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                  PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                  PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                  PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                  PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                  PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                  PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                  PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                  PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                  PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                  PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                  PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                  PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                        control=list(maxit=20000))
summary(fit_Photo_beta_all_Topt.lin_Asat)

#Acclimation of Tmax
fit_Photo_beta_all_Tmax.lin_Asat <- mle2(Photo_beta_normNLL_all_Tmax.lin,start=list(sdPhoto=-1,ymaxM21a=5,ymaxM21b=5,ymaxM21c=5,ymaxM21d=7,ymaxM21e=7,ymaxM21f=9,
                                                                                  ymaxM26a=8,ymaxM26b=6,ymaxM26c=5.5,ymaxM26d=5,ymaxM26e=11,ymaxM26f=6,
                                                                                  ymaxM31a=8,ymaxM31b=6,ymaxM31c=6,ymaxM31d=4,ymaxM31e=16,ymaxM31f=12,
                                                                                  ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                                  ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                                  ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                                  ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                                  ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                                  ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                                  ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                                  ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                                  ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                                  aM=39,
                                                                                  aA=40,
                                                                                  aG=41,
                                                                                  aR=44,
                                                                                  bM=0.1,
                                                                                  bA=0.1,
                                                                                  bG=0.05,
                                                                                  bR=-0.1,
                                                                                  TminM=2,
                                                                                  TminA=0.2,
                                                                                  TminG=8,
                                                                                  TminR=4,
                                                                                  ToptM=28,
                                                                                  ToptA=28,
                                                                                  ToptG=28,
                                                                                  ToptR=28),
                                        data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                                  TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                                  TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                                  TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                                  TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                                  TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                                  TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                                  TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                                  TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                                  TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                                  TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                                  TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                                  PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                                  PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                                  PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                                  PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                                  PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                                  PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                                  PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                                  PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                                  PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                                  PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                                  PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                                  PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                        control=list(maxit=20000))
summary(fit_Photo_beta_all_Tmax.lin_Asat)

#No acclimation
fit_Photo_beta_all_bin_Asat <- mle2(Photo_beta_normNLL_all_bin,start=list(sdPhoto=-1,ymaxM21a=6,ymaxM21b=6,ymaxM21c=6,ymaxM21d=6,ymaxM21e=8,ymaxM21f=10,
                                                                        ymaxM26a=9,ymaxM26b=7,ymaxM26c=6.5,ymaxM26d=6,ymaxM26e=11,ymaxM26f=7,
                                                                        ymaxM31a=8,ymaxM31b=6,ymaxM31c=7,ymaxM31d=5,ymaxM31e=18,ymaxM31f=14,
                                                                        ymaxA21a=14.5,ymaxA21b=14,ymaxA21c=13,ymaxA21d=13,ymaxA21e=13,ymaxA21f=20,
                                                                        ymaxA26a=14,ymaxA26b=13,ymaxA26c=9,ymaxA26d=7,ymaxA26e=5,ymaxA26f=4,
                                                                        ymaxA31a=8,ymaxA31b=8.5,ymaxA31c=13,ymaxA31d=13,ymaxA31e=15,ymaxA31f=14,
                                                                        ymaxG21a=15,ymaxG21b=8,ymaxG21c=10,ymaxG21d=13,ymaxG21e=6,ymaxG21f=9,
                                                                        ymaxG26a=7.5,ymaxG26b=7,ymaxG26c=5,ymaxG26d=8.5,ymaxG26e=4,ymaxG26f=5,ymaxG26g=5,
                                                                        ymaxG31a=11,ymaxG31b=14,ymaxG31c=6,ymaxG31d=4,ymaxG31e=7,ymaxG31f=7,
                                                                        ymaxR21a=6,ymaxR21b=7,ymaxR21c=8,ymaxR21d=5,ymaxR21e=7,ymaxR21f=6,ymaxR21g=7,ymaxR21h=13,
                                                                        ymaxR26a=8,ymaxR26b=7,ymaxR26c=8,ymaxR26d=11,ymaxR26e=6,ymaxR26f=8,
                                                                        ymaxR31a=14,ymaxR31b=11,ymaxR31c=8,ymaxR31d=7,ymaxR31e=8,ymaxR31f=9,
                                                                        TmaxM=42,
                                                                        TmaxA=43,
                                                                        TmaxG=43,
                                                                        TmaxR=42,
                                                                        TminM=-0.5,
                                                                        TminA=-2,
                                                                        TminG=8,
                                                                        TminR=4,
                                                                        ToptM=29,
                                                                        ToptA=29,
                                                                        ToptG=28,
                                                                        ToptR=28),
                                   data=list(TempM21a=M21.Asat.dat$Temp[1:3],TempM21b=M21.Asat.dat$Temp[4:8],TempM21c=M21.Asat.dat$Temp[9:11],TempM21d=M21.Asat.dat$Temp[12:16],TempM21e=M21.Asat.dat$Temp[17:19],TempM21f=M21.Asat.dat$Temp[20:24],
                                             TempM26a=M26.Asat.dat$Temp[1:4],TempM26b=M26.Asat.dat$Temp[5:8],TempM26c=M26.Asat.dat$Temp[9:12],TempM26d=M26.Asat.dat$Temp[13:16],TempM26e=M26.Asat.dat$Temp[17:20],TempM26f=M26.Asat.dat$Temp[21:24],
                                             TempM31a=M31.Asat.dat$Temp[1:5],TempM31b=M31.Asat.dat$Temp[6:8],TempM31c=M31.Asat.dat$Temp[9:12],TempM31d=M31.Asat.dat$Temp[13:15],TempM31e=M31.Asat.dat$Temp[16:20],TempM31f=M31.Asat.dat$Temp[21:23],
                                             TempA21a=A21.Asat.dat$Temp[1:3],TempA21b=A21.Asat.dat$Temp[4:8],TempA21c=A21.Asat.dat$Temp[9:11],TempA21d=A21.Asat.dat$Temp[12:16],TempA21e=A21.Asat.dat$Temp[17:19],TempA21f=A21.Asat.dat$Temp[20:24],
                                             TempA26a=A26.Asat.dat$Temp[1:4],TempA26b=A26.Asat.dat$Temp[5:8],TempA26c=A26.Asat.dat$Temp[9:12],TempA26d=A26.Asat.dat$Temp[13:16],TempA26e=A26.Asat.dat$Temp[17:20],TempA26f=A26.Asat.dat$Temp[21:24],
                                             TempA31a=A31.Asat.dat$Temp[1:5],TempA31b=A31.Asat.dat$Temp[6:8],TempA31c=A31.Asat.dat$Temp[9:13],TempA31d=A31.Asat.dat$Temp[14:16],TempA31e=A31.Asat.dat$Temp[17:21],TempA31f=A31.Asat.dat$Temp[22:24],
                                             TempG21a=G21.Asat.dat$Temp[1:3],TempG21b=G21.Asat.dat$Temp[4:8],TempG21c=G21.Asat.dat$Temp[9:11],TempG21d=G21.Asat.dat$Temp[12:16],TempG21e=G21.Asat.dat$Temp[17:19],TempG21f=G21.Asat.dat$Temp[20:24],
                                             TempG26a=G26.Asat.dat$Temp[1:4],TempG26b=G26.Asat.dat$Temp[5:8],TempG26c=G26.Asat.dat$Temp[9:12],TempG26d=G26.Asat.dat$Temp[13:16],TempG26e=G26.Asat.dat$Temp[17:20],TempG26f=G26.Asat.dat$Temp[21:24],TempG26g=G26.Asat.dat$Temp[25:28],
                                             TempG31a=G31.Asat.dat$Temp[1:4],TempG31b=G31.Asat.dat$Temp[5:7],TempG31c=G31.Asat.dat$Temp[8:12],TempG31d=G31.Asat.dat$Temp[13:15],TempG31e=G31.Asat.dat$Temp[16:21],TempG31f=G31.Asat.dat$Temp[22:24],
                                             TempR21a=R21.Asat.dat$Temp[1:3],TempR21b=R21.Asat.dat$Temp[4:8],TempR21c=R21.Asat.dat$Temp[9:11],TempR21d=R21.Asat.dat$Temp[12:16],TempR21e=R21.Asat.dat$Temp[17:19],TempR21f=R21.Asat.dat$Temp[20:24],TempR21g=R21.Asat.dat$Temp[25:27],TempR21h=R21.Asat.dat$Temp[28:32],
                                             TempR26a=R26.Asat.dat$Temp[1:4],TempR26b=R26.Asat.dat$Temp[5:8],TempR26c=R26.Asat.dat$Temp[9:13],TempR26d=R26.Asat.dat$Temp[14:17],TempR26e=R26.Asat.dat$Temp[18:21],TempR26f=R26.Asat.dat$Temp[22:25],
                                             TempR31a=R31.Asat.dat$Temp[1:5],TempR31b=R31.Asat.dat$Temp[6:8],TempR31c=R31.Asat.dat$Temp[9:13],TempR31d=R31.Asat.dat$Temp[14:16],TempR31e=R31.Asat.dat$Temp[17:21],TempR31f=R31.Asat.dat$Temp[22:24],
                                             PhotodatM21a=M21.Asat.dat$A400[1:3],PhotodatM21b=M21.Asat.dat$A400[4:8],PhotodatM21c=M21.Asat.dat$A400[9:11],PhotodatM21d=M21.Asat.dat$A400[12:16],PhotodatM21e=M21.Asat.dat$A400[17:19],PhotodatM21f=M21.Asat.dat$A400[20:24],
                                             PhotodatM26a=M26.Asat.dat$A400[1:4],PhotodatM26b=M26.Asat.dat$A400[5:8],PhotodatM26c=M26.Asat.dat$A400[9:12],PhotodatM26d=M26.Asat.dat$A400[13:16],PhotodatM26e=M26.Asat.dat$A400[17:20],PhotodatM26f=M26.Asat.dat$A400[21:24],
                                             PhotodatM31a=M31.Asat.dat$A400[1:5],PhotodatM31b=M31.Asat.dat$A400[6:8],PhotodatM31c=M31.Asat.dat$A400[9:12],PhotodatM31d=M31.Asat.dat$A400[13:15],PhotodatM31e=M31.Asat.dat$A400[16:20],PhotodatM31f=M31.Asat.dat$A400[21:23],
                                             PhotodatA21a=A21.Asat.dat$A400[1:3],PhotodatA21b=A21.Asat.dat$A400[4:8],PhotodatA21c=A21.Asat.dat$A400[9:11],PhotodatA21d=A21.Asat.dat$A400[12:16],PhotodatA21e=A21.Asat.dat$A400[17:19],PhotodatA21f=A21.Asat.dat$A400[20:24],
                                             PhotodatA26a=A26.Asat.dat$A400[1:4],PhotodatA26b=A26.Asat.dat$A400[5:8],PhotodatA26c=A26.Asat.dat$A400[9:12],PhotodatA26d=A26.Asat.dat$A400[13:16],PhotodatA26e=A26.Asat.dat$A400[17:20],PhotodatA26f=A26.Asat.dat$A400[21:24],
                                             PhotodatA31a=A31.Asat.dat$A400[1:5],PhotodatA31b=A31.Asat.dat$A400[6:8],PhotodatA31c=A31.Asat.dat$A400[9:13],PhotodatA31d=A31.Asat.dat$A400[14:16],PhotodatA31e=A31.Asat.dat$A400[17:21],PhotodatA31f=A31.Asat.dat$A400[22:24],
                                             PhotodatG21a=G21.Asat.dat$A400[1:3],PhotodatG21b=G21.Asat.dat$A400[4:8],PhotodatG21c=G21.Asat.dat$A400[9:11],PhotodatG21d=G21.Asat.dat$A400[12:16],PhotodatG21e=G21.Asat.dat$A400[17:19],PhotodatG21f=G21.Asat.dat$A400[20:24],
                                             PhotodatG26a=G26.Asat.dat$A400[1:4],PhotodatG26b=G26.Asat.dat$A400[5:8],PhotodatG26c=G26.Asat.dat$A400[9:12],PhotodatG26d=G26.Asat.dat$A400[13:16],PhotodatG26e=G26.Asat.dat$A400[17:20],PhotodatG26f=G26.Asat.dat$A400[21:24],PhotodatG26g=G26.Asat.dat$A400[25:28],
                                             PhotodatG31a=G31.Asat.dat$A400[1:4],PhotodatG31b=G31.Asat.dat$A400[5:7],PhotodatG31c=G31.Asat.dat$A400[8:12],PhotodatG31d=G31.Asat.dat$A400[13:15],PhotodatG31e=G31.Asat.dat$A400[16:21],PhotodatG31f=G31.Asat.dat$A400[22:24],
                                             PhotodatR21a=R21.Asat.dat$A400[1:3],PhotodatR21b=R21.Asat.dat$A400[4:8],PhotodatR21c=R21.Asat.dat$A400[9:11],PhotodatR21d=R21.Asat.dat$A400[12:16],PhotodatR21e=R21.Asat.dat$A400[17:19],PhotodatR21f=R21.Asat.dat$A400[20:24],PhotodatR21g=R21.Asat.dat$A400[25:27],PhotodatR21h=R21.Asat.dat$A400[28:32],
                                             PhotodatR26a=R26.Asat.dat$A400[1:4],PhotodatR26b=R26.Asat.dat$A400[5:8],PhotodatR26c=R26.Asat.dat$A400[9:13],PhotodatR26d=R26.Asat.dat$A400[14:17],PhotodatR26e=R26.Asat.dat$A400[18:21],PhotodatR26f=R26.Asat.dat$A400[22:25],
                                             PhotodatR31a=R31.Asat.dat$A400[1:5],PhotodatR31b=R31.Asat.dat$A400[6:8],PhotodatR31c=R31.Asat.dat$A400[9:13],PhotodatR31d=R31.Asat.dat$A400[14:16],PhotodatR31e=R31.Asat.dat$A400[17:21],PhotodatR31f=R31.Asat.dat$A400[22:24]),
                                   control=list(maxit=20000))
summary(fit_Photo_beta_all_bin_Asat)

#Delta AIC
AICctab(fit_Photo_beta_all_Tmin.Topt.lin_Asat,fit_Photo_beta_all_Tmin.Tmax.lin_Asat,
        fit_Photo_beta_all_Topt.Tmax.lin_Asat,fit_Photo_beta_all_Tmin.lin_Asat,fit_Photo_beta_all_Topt.lin_Asat,
        fit_Photo_beta_all_Tmax.lin_Asat,fit_Photo_beta_all_bin_Asat,nobs=300)

###############################################################################################################
#Generate parameter estimates and calculate 95% CI for Supplementary Table 2
###############################################################################################################

####
#A275
####

#Parameter estimates
coef(fit_Photo_beta_all_Topt.lin_A275)

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.A275=mvrnorm(1000,mu=coef(fit_Photo_beta_all_Topt.lin_A275),Sigma=vcov(fit_Photo_beta_all_Topt.lin_A275))

###
#Morella
###

#Tmin
dist.m.Tmin.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.m.Tmin.A275[i]=vmat.A275[i,81]
}
quantile(dist.m.Tmin.A275,0.025)
quantile(dist.m.Tmin.A275,0.975)

#Topt intercept
dist.m.a.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.m.a.A275[i]=vmat.A275[i,85]
}
quantile(dist.m.a.A275,0.025)
quantile(dist.m.a.A275,0.975)

#Topt slope
dist.m.b.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.m.b.A275[i]=vmat.A275[i,89]
}
quantile(dist.m.b.A275,0.025)
quantile(dist.m.b.A275,0.975)

#Tmax
dist.m.Tmax.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.m.Tmax.A275[i]=vmat.A275[i,77]
}
quantile(dist.m.Tmax.A275,0.025)
quantile(dist.m.Tmax.A275,0.975)

###
#Alnus
###

#Tmin
dist.a.Tmin.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.a.Tmin.A275[i]=vmat.A275[i,82]
}
quantile(dist.a.Tmin.A275,0.025)
quantile(dist.a.Tmin.A275,0.975)

#Topt intercept
dist.a.a.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.a.a.A275[i]=vmat.A275[i,86]
}
quantile(dist.a.a.A275,0.025)
quantile(dist.a.a.A275,0.975)

#Topt slope
dist.a.b.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.a.b.A275[i]=vmat.A275[i,90]
}
quantile(dist.a.b.A275,0.025)
quantile(dist.a.b.A275,0.975)

#Tmax
dist.a.Tmax.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.a.Tmax.A275[i]=vmat.A275[i,78]
}
quantile(dist.a.Tmax.A275,0.025)
quantile(dist.a.Tmax.A275,0.975)

###
#Gliricidia
###

#Tmin
dist.g.Tmin.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.g.Tmin.A275[i]=vmat.A275[i,83]
}
quantile(dist.g.Tmin.A275,0.025)
quantile(dist.g.Tmin.A275,0.975)

#Topt intercept
dist.g.a.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.g.a.A275[i]=vmat.A275[i,87]
}
quantile(dist.g.a.A275,0.025)
quantile(dist.g.a.A275,0.975)

#Topt slope
dist.g.b.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.g.b.A275[i]=vmat.A275[i,91]
}
quantile(dist.g.b.A275,0.025)
quantile(dist.g.b.A275,0.975)

#Tmax
dist.g.Tmax.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.g.Tmax.A275[i]=vmat.A275[i,79]
}
quantile(dist.g.Tmax.A275,0.025)
quantile(dist.g.Tmax.A275,0.975)

###
#Robinia
###

#Tmin
dist.r.Tmin.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.r.Tmin.A275[i]=vmat.A275[i,84]
}
quantile(dist.r.Tmin.A275,0.025)
quantile(dist.r.Tmin.A275,0.975)

#Topt intercept
dist.r.a.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.r.a.A275[i]=vmat.A275[i,88]
}
quantile(dist.r.a.A275,0.025)
quantile(dist.r.a.A275,0.975)

#Topt slope
dist.r.b.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.r.b.A275[i]=vmat.A275[i,92]
}
quantile(dist.r.b.A275,0.025)
quantile(dist.r.b.A275,0.975)

#Tmax
dist.r.Tmax.A275<-rep(NA,1000)
for(i in 1:1000){
  dist.r.Tmax.A275[i]=vmat.A275[i,80]
}
quantile(dist.r.Tmax.A275,0.025)
quantile(dist.r.Tmax.A275,0.975)

####
#Asat
####

#Parameter estimates
coef(fit_Photo_beta_all_Topt.lin_Asat)

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.Asat=mvrnorm(1000,mu=coef(fit_Photo_beta_all_Topt.lin_Asat),Sigma=vcov(fit_Photo_beta_all_Topt.lin_Asat))

###
#Morella
###

#Tmin
dist.m.Tmin.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.m.Tmin.Asat[i]=vmat.Asat[i,81]
}
quantile(dist.m.Tmin.Asat,0.025)
quantile(dist.m.Tmin.Asat,0.975)

#Topt intercept
dist.m.a.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.m.a.Asat[i]=vmat.Asat[i,85]
}
quantile(dist.m.a.Asat,0.025)
quantile(dist.m.a.Asat,0.975)

#Topt slope
dist.m.b.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.m.b.Asat[i]=vmat.Asat[i,89]
}
quantile(dist.m.b.Asat,0.025)
quantile(dist.m.b.Asat,0.975)

#Tmax
dist.m.Tmax.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.m.Tmax.Asat[i]=vmat.Asat[i,77]
}
quantile(dist.m.Tmax.Asat,0.025)
quantile(dist.m.Tmax.Asat,0.975)

###
#Alnus
###

#Tmin
dist.a.Tmin.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.a.Tmin.Asat[i]=vmat.Asat[i,82]
}
quantile(dist.a.Tmin.Asat,0.025)
quantile(dist.a.Tmin.Asat,0.975)

#Topt intercept
dist.a.a.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.a.a.Asat[i]=vmat.Asat[i,86]
}
quantile(dist.a.a.Asat,0.025)
quantile(dist.a.a.Asat,0.975)

#Topt slope
dist.a.b.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.a.b.Asat[i]=vmat.Asat[i,90]
}
quantile(dist.a.b.Asat,0.025)
quantile(dist.a.b.Asat,0.975)

#Tmax
dist.a.Tmax.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.a.Tmax.Asat[i]=vmat.Asat[i,78]
}
quantile(dist.a.Tmax.Asat,0.025)
quantile(dist.a.Tmax.Asat,0.975)

###
#Gliricidia
###

#Tmin
dist.g.Tmin.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.g.Tmin.Asat[i]=vmat.Asat[i,83]
}
quantile(dist.g.Tmin.Asat,0.025)
quantile(dist.g.Tmin.Asat,0.975)

#Topt intercept
dist.g.a.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.g.a.Asat[i]=vmat.Asat[i,87]
}
quantile(dist.g.a.Asat,0.025)
quantile(dist.g.a.Asat,0.975)

#Topt slope
dist.g.b.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.g.b.Asat[i]=vmat.Asat[i,91]
}
quantile(dist.g.b.Asat,0.025)
quantile(dist.g.b.Asat,0.975)

#Tmax
dist.g.Tmax.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.g.Tmax.Asat[i]=vmat.Asat[i,79]
}
quantile(dist.g.Tmax.Asat,0.025)
quantile(dist.g.Tmax.Asat,0.975)

###
#Robinia
###

#Tmin
dist.r.Tmin.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.r.Tmin.Asat[i]=vmat.Asat[i,84]
}
quantile(dist.r.Tmin.Asat,0.025)
quantile(dist.r.Tmin.Asat,0.975)

#Topt intercept
dist.r.a.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.r.a.Asat[i]=vmat.Asat[i,88]
}
quantile(dist.r.a.Asat,0.025)
quantile(dist.r.a.Asat,0.975)

#Topt slope
dist.r.b.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.r.b.Asat[i]=vmat.Asat[i,92]
}
quantile(dist.r.b.Asat,0.025)
quantile(dist.r.b.Asat,0.975)

#Tmax
dist.r.Tmax.Asat<-rep(NA,1000)
for(i in 1:1000){
  dist.r.Tmax.Asat[i]=vmat.Asat[i,80]
}
quantile(dist.r.Tmax.Asat,0.025)
quantile(dist.r.Tmax.Asat,0.975)