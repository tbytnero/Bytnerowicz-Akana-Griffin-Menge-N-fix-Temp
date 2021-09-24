###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figures 7-10
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
#Read in output data from A-Ci curve statistical fits
#See A-Ci_Calculations.R for A-Ci curve statistical fits
####

A21.dat<-read.csv("ACi.dat.A21.csv")
M21.dat<-read.csv("ACi.dat.M21.csv")
G21.dat<-read.csv("ACi.dat.G21.csv")
R21.dat<-read.csv("ACi.dat.R21.csv")
A26.dat<-read.csv("ACi.dat.A26.csv")
M26.dat<-read.csv("ACi.dat.M26.csv")
G26.dat<-read.csv("ACi.dat.G26.csv")
R26.dat<-read.csv("ACi.dat.R26.csv")
A31.dat<-read.csv("ACi.dat.A31.csv")
M31.dat<-read.csv("ACi.dat.M31.csv")
G31.dat<-read.csv("ACi.dat.G31.csv")
R31.dat<-read.csv("ACi.dat.R31.csv")

#Calculate temperature in Kelvin
A21.dat$TsK<-A21.dat$Temp+273.15
M21.dat$TsK<-M21.dat$Temp+273.15
G21.dat$TsK<-G21.dat$Temp+273.15
R21.dat$TsK<-R21.dat$Temp+273.15
A26.dat$TsK<-A26.dat$Temp+273.15
M26.dat$TsK<-M26.dat$Temp+273.15
G26.dat$TsK<-G26.dat$Temp+273.15
R26.dat$TsK<-R26.dat$Temp+273.15
A31.dat$TsK<-A31.dat$Temp+273.15
M31.dat$TsK<-M31.dat$Temp+273.15
G31.dat$TsK<-G31.dat$Temp+273.15
R31.dat$TsK<-R31.dat$Temp+273.15

####
#Define function
####

#Peaked Arrhenius function (equation 4)
##Here, delta S is reparameterized as Hd/Topt+0.008314*log(Ea/(Hd-Ea)),
##following Medlyn et al. 2002 (reference 57)
##this helped with estimating starting parameter values and made it easier to fit the model
peak.ar <- function(k25,Ea,Topt,Hd,Tk){
  y <- k25*(Tk/298.15)*exp((Ea*(Tk - 298.15))/(298.15*0.008314*Tk)) * 
    (1+exp((298.15*(Hd/Topt+0.008314*log(Ea/(Hd-Ea))) - Hd)/(298.15*0.008314))) / 
    (1+exp((Tk*(Hd/Topt+0.008314*log(Ea/(Hd-Ea)))-Hd)/(Tk*0.008314)))
  y
}

####
#Negative log-likelihood (NLL) functions
####

#NLL function for everything except Gliricidia at 26:20 deg. C and Robinia at 21:15 deg. C
normNLL_f_Hd <- function(sdNase,k25a,k25b,k25c,k25d,k25e,k25f,
                                     Ea,
                                     Topt,
                                     Tka,Tkb,Tkc,Tkd,Tke,Tkf,
                                     Nasedata,Nasedatb,Nasedatc,Nasedatd,Nasedate,Nasedatf){
  Nasemeana <- peak.ar(k25a,Ea,Topt,200,Tka)
  Nasemeanb <- peak.ar(k25b,Ea,Topt,200,Tkb)
  Nasemeanc <- peak.ar(k25c,Ea,Topt,200,Tkc)
  Nasemeand <- peak.ar(k25d,Ea,Topt,200,Tkd)
  Nasemeane <- peak.ar(k25e,Ea,Topt,200,Tke)
  Nasemeanf <- peak.ar(k25f,Ea,Topt,200,Tkf)
  -(sum(dnorm(Nasedata,mean=Nasemeana,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatb,mean=Nasemeanb,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatc,mean=Nasemeanc,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedatd,mean=Nasemeand,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedate,mean=Nasemeane,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatf,mean=Nasemeanf,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for Gliricidia at 26:20 deg. C
normNLL_g_Hd <- function(sdNase,k25a,k25b,k25c,k25d,k25e,k25f,k25g,
                                     Ea,
                                     Topt,
                                     Tka,Tkb,Tkc,Tkd,Tke,Tkf,Tkg,
                                     Nasedata,Nasedatb,Nasedatc,Nasedatd,Nasedate,Nasedatf,Nasedatg){
  Nasemeana <- peak.ar(k25a,Ea,Topt,200,Tka)
  Nasemeanb <- peak.ar(k25b,Ea,Topt,200,Tkb)
  Nasemeanc <- peak.ar(k25c,Ea,Topt,200,Tkc)
  Nasemeand <- peak.ar(k25d,Ea,Topt,200,Tkd)
  Nasemeane <- peak.ar(k25e,Ea,Topt,200,Tke)
  Nasemeanf <- peak.ar(k25f,Ea,Topt,200,Tkf)
  Nasemeang <- peak.ar(k25g,Ea,Topt,200,Tkg)
  -(sum(dnorm(Nasedata,mean=Nasemeana,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatb,mean=Nasemeanb,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatc,mean=Nasemeanc,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedatd,mean=Nasemeand,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedate,mean=Nasemeane,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatf,mean=Nasemeanf,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedatg,mean=Nasemeang,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#NLL function for Robinia at 21:15 deg. C
normNLL_h_Hd <- function(sdNase,k25a,k25b,k25c,k25d,k25e,k25f,k25g,k25h,
                                     Ea,
                                     Topt,
                                     Tka,Tkb,Tkc,Tkd,Tke,Tkf,Tkg,Tkh,
                                     Nasedata,Nasedatb,Nasedatc,Nasedatd,Nasedate,Nasedatf,Nasedatg,Nasedath){
  Nasemeana <- peak.ar(k25a,Ea,Topt,200,Tka)
  Nasemeanb <- peak.ar(k25b,Ea,Topt,200,Tkb)
  Nasemeanc <- peak.ar(k25c,Ea,Topt,200,Tkc)
  Nasemeand <- peak.ar(k25d,Ea,Topt,200,Tkd)
  Nasemeane <- peak.ar(k25e,Ea,Topt,200,Tke)
  Nasemeanf <- peak.ar(k25f,Ea,Topt,200,Tkf)
  Nasemeang <- peak.ar(k25g,Ea,Topt,200,Tkg)
  Nasemeanh <- peak.ar(k25h,Ea,Topt,200,Tkh)
  -(sum(dnorm(Nasedata,mean=Nasemeana,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatb,mean=Nasemeanb,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatc,mean=Nasemeanc,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedatd,mean=Nasemeand,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedate,mean=Nasemeane,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedatf,mean=Nasemeanf,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedatg,mean=Nasemeang,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedath,mean=Nasemeanh,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

####
#Maximum likelihood fits
####

###
#Morella
###

##
#Vcmax
##

#21:15 deg. C growing temperature
fit_M21_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=45*.5,k25b=45*.5,k25c=50*.5,k25d=45*.5,k25e=60*.5,k25f=80*.5,
                                                                  Ea=80,
                                                                  Topt=309),
                              data=list(Tka=M21.dat$TsK[1:3],Tkb=M21.dat$TsK[4:8],Tkc=M21.dat$TsK[9:11],Tkd=M21.dat$TsK[12:16],Tke=M21.dat$TsK[17:19],Tkf=M21.dat$TsK[20:24],
                                        Nasedata=M21.dat$Vcmax[1:3],Nasedatb=M21.dat$Vcmax[4:8],Nasedatc=M21.dat$Vcmax[9:11],Nasedatd=M21.dat$Vcmax[12:16],Nasedate=M21.dat$Vcmax[17:19],Nasedatf=M21.dat$Vcmax[20:24]),
                              control=list(maxit=20000))
summary(fit_M21_VcmaxHd)

#26:20 deg. C growing temperature
fit_M26_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=30,k25b=25,k25c=20,k25d=20,k25e=40,k25f=20,
                                                                  Ea=80,
                                                                  Topt=309),
                              data=list(Tka=M26.dat$TsK[1:4],Tkb=M26.dat$TsK[5:8],Tkc=M26.dat$TsK[9:12],Tkd=M26.dat$TsK[13:15],Tke=M26.dat$TsK[16:19],Tkf=M26.dat$TsK[20:23],
                                        Nasedata=M26.dat$Vcmax[1:4],Nasedatb=M26.dat$Vcmax[5:8],Nasedatc=M26.dat$Vcmax[9:12],Nasedatd=M26.dat$Vcmax[13:15],Nasedate=M26.dat$Vcmax[16:19],Nasedatf=M26.dat$Vcmax[20:23]),
                              control=list(maxit=20000))
summary(fit_M26_VcmaxHd)

#31:25 deg. C growing temperature
fit_M31_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=30,k25b=23,k25c=39,k25d=20,k25e=63,k25f=60,
                                                                  Ea=40,
                                                                  Topt=309),
                              data=list(Tka=M31.dat$TsK[1:5],Tkb=M31.dat$TsK[6:8],Tkc=M31.dat$TsK[9:12],Tkd=M31.dat$TsK[13:15],Tke=M31.dat$TsK[16:20],Tkf=M31.dat$TsK[21:23],
                                        Nasedata=M31.dat$Vcmax[1:5],Nasedatb=M31.dat$Vcmax[6:8],Nasedatc=M31.dat$Vcmax[9:12],Nasedatd=M31.dat$Vcmax[13:15],Nasedate=M31.dat$Vcmax[16:20],Nasedatf=M31.dat$Vcmax[21:23]),
                              control=list(maxit=20000))
summary(fit_M31_VcmaxHd)

##
#Jmax
##

#21:15 deg. C growing temperature
fit_M21_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=70*.65,k25b=65*.65,k25c=75*.65,k25d=65*.65,k25e=75*.65,k25f=100*.65,
                                                                 Ea=70,
                                                                 Topt=306),
                             data=list(Tka=M21.dat$TsK[1:3],Tkb=M21.dat$TsK[4:8],Tkc=M21.dat$TsK[9:11],Tkd=M21.dat$TsK[12:16],Tke=M21.dat$TsK[17:19],Tkf=M21.dat$TsK[20:24],
                                       Nasedata=M21.dat$Jmax[1:3],Nasedatb=M21.dat$Jmax[4:8],Nasedatc=M21.dat$Jmax[9:11],Nasedatd=M21.dat$Jmax[12:16],Nasedate=M21.dat$Jmax[17:19],Nasedatf=M21.dat$Jmax[20:24]),
                             control=list(maxit=20000))
summary(fit_M21_JmaxHd)

#26:20 deg. C growing temperature
fit_M26_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=90*.65,k25b=75*.65,k25c=65*.65,k25d=50*.65,k25e=110*.65,k25f=65*.65,
                                                                 Ea=80,
                                                                 Topt=306),
                             data=list(Tka=M26.dat$TsK[1:4],Tkb=M26.dat$TsK[5:8],Tkc=M26.dat$TsK[9:12],Tkd=M26.dat$TsK[13:15],Tke=M26.dat$TsK[16:19],Tkf=M26.dat$TsK[20:23],
                                       Nasedata=M26.dat$Jmax[1:4],Nasedatb=M26.dat$Jmax[5:8],Nasedatc=M26.dat$Jmax[9:12],Nasedatd=M26.dat$Jmax[13:15],Nasedate=M26.dat$Jmax[16:19],Nasedatf=M26.dat$Jmax[20:23]),
                             control=list(maxit=20000))
summary(fit_M26_JmaxHd)

#31:25 deg. C growing temperature
fit_M31_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=85*.65,k25b=65*.65,k25c=70*.65,k25d=45*.65,k25e=170*.65,k25f=120*.65,
                                                                 Ea=40,
                                                                 Topt=307),
                             data=list(Tka=M31.dat$TsK[1:5],Tkb=M31.dat$TsK[6:8],Tkc=M31.dat$TsK[9:12],Tkd=M31.dat$TsK[13:15],Tke=M31.dat$TsK[16:20],Tkf=M31.dat$TsK[21:23],
                                       Nasedata=M31.dat$Jmax[1:5],Nasedatb=M31.dat$Jmax[6:8],Nasedatc=M31.dat$Jmax[9:12],Nasedatd=M31.dat$Jmax[13:15],Nasedate=M31.dat$Jmax[16:20],Nasedatf=M31.dat$Jmax[21:23]),
                             control=list(maxit=20000))
summary(fit_M31_JmaxHd)

###
#Alnus
###

##
#Vcmax
##

#21:15 deg. C growing temperature
fit_A21_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=60,k25b=66,k25c=55,k25d=61,k25e=55,k25f=81,
                                                                  Ea=80,
                                                                  Topt=307),
                              data=list(Tka=A21.dat$TsK[1:3],Tkb=A21.dat$TsK[4:8],Tkc=A21.dat$TsK[9:11],Tkd=A21.dat$TsK[12:16],Tke=A21.dat$TsK[17:19],Tkf=A21.dat$TsK[20:24],
                                        Nasedata=A21.dat$Vcmax[1:3],Nasedatb=A21.dat$Vcmax[4:8],Nasedatc=A21.dat$Vcmax[9:11],Nasedatd=A21.dat$Vcmax[12:16],Nasedate=A21.dat$Vcmax[17:19],Nasedatf=A21.dat$Vcmax[20:24]),
                              control=list(maxit=20000))
summary(fit_A21_VcmaxHd)

#26:20 deg. C growing temperature
fit_A26_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=53,k25b=52,k25c=37,k25d=35,k25e=17,k25f=19,
                                                                  Ea=120,
                                                                  Topt=308),
                              data=list(Tka=A26.dat$TsK[1:4],Tkb=A26.dat$TsK[5:8],Tkc=A26.dat$TsK[9:12],Tkd=A26.dat$TsK[13:15],Tke=A26.dat$TsK[16:19],Tkf=A26.dat$TsK[20:22],
                                        Nasedata=A26.dat$Vcmax[1:4],Nasedatb=A26.dat$Vcmax[5:8],Nasedatc=A26.dat$Vcmax[9:12],Nasedatd=A26.dat$Vcmax[13:15],Nasedate=A26.dat$Vcmax[16:19],Nasedatf=A26.dat$Vcmax[20:22]),
                              control=list(maxit=20000))
summary(fit_A26_VcmaxHd)

#31:25 deg. C growing temperature
fit_A31_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=31,k25b=35,k25c=48,k25d=55,k25e=55,k25f=56,
                                                                  Ea=30,
                                                                  Topt=307),
                              data=list(Tka=A31.dat$TsK[1:5],Tkb=A31.dat$TsK[6:7],Tkc=A31.dat$TsK[8:12],Tkd=A31.dat$TsK[13:15],Tke=A31.dat$TsK[16:20],Tkf=A31.dat$TsK[21:23],
                                        Nasedata=A31.dat$Vcmax[1:5],Nasedatb=A31.dat$Vcmax[6:7],Nasedatc=A31.dat$Vcmax[8:12],Nasedatd=A31.dat$Vcmax[13:15],Nasedate=A31.dat$Vcmax[16:20],Nasedatf=A31.dat$Vcmax[21:23]),
                              control=list(maxit=20000))
summary(fit_A31_VcmaxHd)

##
#Jmax
##

#21:15 deg. C growing temperature
fit_A21_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=170*.65,k25b=200*.65,k25c=170*.65,k25d=110*.65,k25e=170*.65,k25f=150*.65,
                                                                 Ea=60,
                                                                 Topt=305),
                             data=list(Tka=A21.dat$TsK[1:3],Tkb=A21.dat$TsK[4:8],Tkc=A21.dat$TsK[9:11],Tkd=A21.dat$TsK[12:16],Tke=A21.dat$TsK[17:19],Tkf=A21.dat$TsK[20:24],
                                       Nasedata=A21.dat$Jmax[1:3],Nasedatb=A21.dat$Jmax[4:8],Nasedatc=A21.dat$Jmax[9:11],Nasedatd=A21.dat$Jmax[12:16],Nasedate=A21.dat$Jmax[17:19],Nasedatf=A21.dat$Jmax[20:24]),
                             control=list(maxit=20000))
summary(fit_A21_JmaxHd) 

#26:20 deg. C growing temperature
fit_A26_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=140*.65,k25b=130*.65,k25c=90*.65,k25d=70*.65,k25e=50*.65,k25f=45*.65,
                                                                 Ea=65,
                                                                 Topt=307),
                             data=list(Tka=A26.dat$TsK[1:4],Tkb=A26.dat$TsK[5:8],Tkc=A26.dat$TsK[9:12],Tkd=A26.dat$TsK[13:15],Tke=A26.dat$TsK[16:19],Tkf=A26.dat$TsK[20:22],
                                       Nasedata=A26.dat$Jmax[1:4],Nasedatb=A26.dat$Jmax[5:8],Nasedatc=A26.dat$Jmax[9:12],Nasedatd=A26.dat$Jmax[13:15],Nasedate=A26.dat$Jmax[16:19],Nasedatf=A26.dat$Jmax[20:22]),
                             control=list(maxit=20000))
summary(fit_A26_JmaxHd)

#31:25 deg. C growing temperature
fit_A31_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=33,k25b=40,k25c=50,k25d=60,k25e=56,k25f=60,
                                                                 Ea=40,
                                                                 Topt=308),
                             data=list(Tka=A31.dat$TsK[1:5],Tkb=A31.dat$TsK[6:8],Tkc=A31.dat$TsK[9:12],Tkd=A31.dat$TsK[13:15],Tke=A31.dat$TsK[16:20],Tkf=A31.dat$TsK[21:23],
                                       Nasedata=A31.dat$Jmax[1:5],Nasedatb=A31.dat$Jmax[6:8],Nasedatc=A31.dat$Jmax[9:12],Nasedatd=A31.dat$Jmax[13:15],Nasedate=A31.dat$Jmax[16:20],Nasedatf=A31.dat$Jmax[21:23]),
                             control=list(maxit=20000))
summary(fit_A31_JmaxHd)

###
#Gliricidia
###

##
#Vcmax
##

#21:15 deg. C growing temperature
fit_G21_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=120*.5,k25b=60*.5,k25c=85*.5,k25d=115*.5,k25e=65*.5,k25f=87*.5,
                                                                  Ea=90,
                                                                  Topt=309),
                              data=list(Tka=G21.dat$TsK[1:3],Tkb=G21.dat$TsK[4:8],Tkc=G21.dat$TsK[9:11],Tkd=G21.dat$TsK[12:16],Tke=G21.dat$TsK[17:19],Tkf=G21.dat$TsK[20:24],
                                        Nasedata=G21.dat$Vcmax[1:3],Nasedatb=G21.dat$Vcmax[4:8],Nasedatc=G21.dat$Vcmax[9:11],Nasedatd=G21.dat$Vcmax[12:16],Nasedate=G21.dat$Vcmax[17:19],Nasedatf=G21.dat$Vcmax[20:24]),
                              control=list(maxit=20000))
summary(fit_G21_VcmaxHd)

#26:20 deg. C growing temperature
fit_G26_VcmaxHd <- mle2(normNLL_g_Hd,start=list(sdNase=-1,k25a=70*.5,k25b=55*.5,k25c=50*.5,k25d=70*.5,k25e=50*.5,k25f=50*.5,k25g=45*.5,
                                                                  Ea=50,
                                                                  Topt=306),
                              data=list(Tka=G26.dat$TsK[1:4],Tkb=G26.dat$TsK[5:8],Tkc=G26.dat$TsK[9:12],Tkd=G26.dat$TsK[13:16],Tke=G26.dat$TsK[17:20],Tkf=G26.dat$TsK[21:24],Tkg=G26.dat$TsK[25:27],
                                        Nasedata=G26.dat$Vcmax[1:4],Nasedatb=G26.dat$Vcmax[5:8],Nasedatc=G26.dat$Vcmax[9:12],Nasedatd=G26.dat$Vcmax[13:16],Nasedate=G26.dat$Vcmax[17:20],Nasedatf=G26.dat$Vcmax[21:24],Nasedatg=G26.dat$Vcmax[25:27]),
                              control=list(maxit=20000))
summary(fit_G26_VcmaxHd)

#31:25 deg. C growing temperature
fit_G31_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=1.281581,k25a=35.665629,k25b=57.017459,k25c=25.062748,k25d=16.341698,k25e=26.640753,k25f=21.283619,
                                                                  Ea=191.711238,
                                                                  Topt=313.616944),
                              data=list(Tka=G31.dat$TsK[1:4],Tkb=G31.dat$TsK[5:7],Tkc=G31.dat$TsK[8:12],Tkd=G31.dat$TsK[13:14],Tke=G31.dat$TsK[15:20],Tkf=G31.dat$TsK[21:23],
                                        Nasedata=G31.dat$Vcmax[1:4],Nasedatb=G31.dat$Vcmax[5:7],Nasedatc=G31.dat$Vcmax[8:12],Nasedatd=G31.dat$Vcmax[13:14],Nasedate=G31.dat$Vcmax[15:20],Nasedatf=G31.dat$Vcmax[21:23]),
                              control=list(maxit=20000))
summary(fit_G31_VcmaxHd)

##
#Jmax
##

#21:15 deg. C growing temperature
fit_G21_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=200*.7,k25b=100*.7,k25c=120*.7,k25d=190*.7,k25e=85*.7,k25f=100*.7,
                                                                 Ea=95,
                                                                 Topt=306),
                             data=list(Tka=G21.dat$TsK[1:3],Tkb=G21.dat$TsK[4:8],Tkc=G21.dat$TsK[9:11],Tkd=G21.dat$TsK[12:16],Tke=G21.dat$TsK[17:19],Tkf=G21.dat$TsK[20:24],
                                       Nasedata=G21.dat$Jmax[1:3],Nasedatb=G21.dat$Jmax[4:8],Nasedatc=G21.dat$Jmax[9:11],Nasedatd=G21.dat$Jmax[12:16],Nasedate=G21.dat$Jmax[17:19],Nasedatf=G21.dat$Jmax[20:24]),
                             control=list(maxit=20000))
summary(fit_G21_JmaxHd)

#26:20 deg. C growing temperature
fit_G26_JmaxHd <- mle2(normNLL_g_Hd,start=list(sdNase=-1,k25a=95*.7,k25b=90*.7,k25c=70*.7,k25d=100*.7,k25e=60*.7,k25f=70*.7,k25g=70*.7,
                                                                 Ea=40,
                                                                 Topt=306),
                             data=list(Tka=G26.dat$TsK[1:4],Tkb=G26.dat$TsK[5:8],Tkc=G26.dat$TsK[9:12],Tkd=G26.dat$TsK[13:16],Tke=G26.dat$TsK[17:20],Tkf=G26.dat$TsK[21:24],Tkg=G26.dat$TsK[25:27],
                                       Nasedata=G26.dat$Jmax[1:4],Nasedatb=G26.dat$Jmax[5:8],Nasedatc=G26.dat$Jmax[9:12],Nasedatd=G26.dat$Jmax[13:16],Nasedate=G26.dat$Jmax[17:20],Nasedatf=G26.dat$Jmax[21:24],Nasedatg=G26.dat$Jmax[25:27]),
                             control=list(maxit=20000))
summary(fit_G26_JmaxHd)

#31:25 deg. C growing temperature
fit_G31_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=100*.4,k25b=160*.4,k25c=70*.4,k25d=50*.4,k25e=75*.4,k25f=60*.4,
                                                                 Ea=50,
                                                                 Topt=308),
                             data=list(Tka=G31.dat$TsK[1:4],Tkb=G31.dat$TsK[5:7],Tkc=G31.dat$TsK[8:12],Tkd=G31.dat$TsK[13:14],Tke=G31.dat$TsK[15:20],Tkf=G31.dat$TsK[21:23],
                                       Nasedata=G31.dat$Jmax[1:4],Nasedatb=G31.dat$Jmax[5:7],Nasedatc=G31.dat$Jmax[8:12],Nasedatd=G31.dat$Jmax[13:14],Nasedate=G31.dat$Jmax[15:20],Nasedatf=G31.dat$Jmax[21:23]),
                             control=list(maxit=20000))
summary(fit_G31_JmaxHd)

###
#Robinia
###

##
#Vcmax
##

#21:15 deg. C growing temperature
fit_R21_VcmaxHd <- mle2(normNLL_h_Hd,start=list(sdNase=-1,k25a=89.668284,k25b=42.131186,k25c=66.637371,k25d=21.117929,k25e=59.765637,k25f=28.397229,k25g=62.282305,k25h=61.890062,
                                                                  Ea=30,
                                                                  Topt=308),
                              data=list(Tka=R21.dat$TsK[1:2],Tkb=R21.dat$TsK[3:6],Tkc=R21.dat$TsK[7:9],Tkd=R21.dat$TsK[10:14],Tke=R21.dat$TsK[15:17],Tkf=R21.dat$TsK[18:22],Tkg=R21.dat$TsK[23:25],Tkh=R21.dat$TsK[26:30],
                                        Nasedata=R21.dat$Vcmax[1:2],Nasedatb=R21.dat$Vcmax[3:6],Nasedatc=R21.dat$Vcmax[7:9],Nasedatd=R21.dat$Vcmax[10:14],Nasedate=R21.dat$Vcmax[15:17],Nasedatf=R21.dat$Vcmax[18:22],Nasedatg=R21.dat$Vcmax[23:25],Nasedath=R21.dat$Vcmax[26:30]),
                              control=list(maxit=20000))
summary(fit_R21_VcmaxHd)

#26:20 deg. C growing temperature
fit_R26_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=44.598846,k25b=39.011687,k25c=37.819034,k25d=60.926606,k25e=29.727416,k25f=39.052503,
                                                                  Ea=100,
                                                                  Topt=305.887533),
                              data=list(Tka=R26.dat$TsK[1:4],Tkb=R26.dat$TsK[5:8],Tkc=R26.dat$TsK[9:13],Tkd=R26.dat$TsK[14:17],Tke=R26.dat$TsK[18:21],Tkf=R26.dat$TsK[22:25],
                                        Nasedata=R26.dat$Vcmax[1:4],Nasedatb=R26.dat$Vcmax[5:8],Nasedatc=R26.dat$Vcmax[9:13],Nasedatd=R26.dat$Vcmax[14:17],Nasedate=R26.dat$Vcmax[18:21],Nasedatf=R26.dat$Vcmax[22:25]),
                              control=list(maxit=20000))
summary(fit_R26_VcmaxHd)

#31:25 deg. C growing temperature
fit_R31_VcmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=71.814186,k25b=57.450557,k25c=30.773607,k25d=27.749574,k25e=40.056782,k25f=43.359829,
                                                                  Ea=90,
                                                                  Topt=307),
                              data=list(Tka=R31.dat$TsK[1:5],Tkb=R31.dat$TsK[6:8],Tkc=R31.dat$TsK[9:13],Tkd=R31.dat$TsK[14:16],Tke=R31.dat$TsK[17:21],Tkf=R31.dat$TsK[22:24],
                                        Nasedata=R31.dat$Vcmax[1:5],Nasedatb=R31.dat$Vcmax[6:8],Nasedatc=R31.dat$Vcmax[9:13],Nasedatd=R31.dat$Vcmax[14:16],Nasedate=R31.dat$Vcmax[17:21],Nasedatf=R31.dat$Vcmax[22:24]),
                              control=list(maxit=20000))
summary(fit_R31_VcmaxHd)

##
#Jmax
##

#21:15 deg. C growing temperature
fit_R21_JmaxHd <- mle2(normNLL_h_Hd,start=list(sdNase=-1,k25a=64,k25b=72,k25c=70,k25d=43,k25e=28,k25f=62,k25g=56,k25h=104,
                                                                  Ea=70,
                                                                  Topt=307),
                              data=list(Tka=R21.dat$TsK[1:2],Tkb=R21.dat$TsK[3:6],Tkc=R21.dat$TsK[7:9],Tkd=R21.dat$TsK[10:14],Tke=R21.dat$TsK[15:17],Tkf=R21.dat$TsK[18:22],Tkg=R21.dat$TsK[23:25],Tkh=R21.dat$TsK[26:30],
                                        Nasedata=R21.dat$Jmax[1:2],Nasedatb=R21.dat$Jmax[3:6],Nasedatc=R21.dat$Jmax[7:9],Nasedatd=R21.dat$Jmax[10:14],Nasedate=R21.dat$Jmax[15:17],Nasedatf=R21.dat$Jmax[18:22],Nasedatg=R21.dat$Jmax[23:25],Nasedath=R21.dat$Jmax[26:30]),
                              control=list(maxit=20000))
summary(fit_R21_JmaxHd)

#26:20 deg. C growing temperature
fit_R26_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=41,k25b=73,k25c=72,k25d=95,k25e=60,k25f=71,
                                                                  Ea=89,
                                                                  Topt=304),
                              data=list(Tka=R26.dat$TsK[1:4],Tkb=R26.dat$TsK[5:8],Tkc=R26.dat$TsK[9:13],Tkd=R26.dat$TsK[14:17],Tke=R26.dat$TsK[18:21],Tkf=R26.dat$TsK[22:25],
                                        Nasedata=R26.dat$Jmax[1:4],Nasedatb=R26.dat$Jmax[5:8],Nasedatc=R26.dat$Jmax[9:13],Nasedatd=R26.dat$Jmax[14:17],Nasedate=R26.dat$Jmax[18:21],Nasedatf=R26.dat$Jmax[22:25]),
                              control=list(maxit=20000))
summary(fit_R26_JmaxHd)

#31:25 deg. C growing temperature
fit_R31_JmaxHd <- mle2(normNLL_f_Hd,start=list(sdNase=-1,k25a=115,k25b=91,k25c=51,k25d=56,k25e=48,k25f=64,
                                                                  Ea=90,
                                                                  Topt=305),
                              data=list(Tka=R31.dat$TsK[1:5],Tkb=R31.dat$TsK[6:8],Tkc=R31.dat$TsK[9:13],Tkd=R31.dat$TsK[14:16],Tke=R31.dat$TsK[17:21],Tkf=R31.dat$TsK[22:24],
                                        Nasedata=R31.dat$Jmax[1:5],Nasedatb=R31.dat$Jmax[6:8],Nasedatc=R31.dat$Jmax[9:13],Nasedatd=R31.dat$Jmax[14:16],Nasedate=R31.dat$Jmax[17:21],Nasedatf=R31.dat$Jmax[22:24]),
                              control=list(maxit=20000))
summary(fit_R31_JmaxHd)

####
#Parameter estimates
####

f.m21.v<-coef(fit_M21_VcmaxHd)
f.m26.v<-coef(fit_M26_VcmaxHd)
f.m31.v<-coef(fit_M31_VcmaxHd)
f.m21.j<-coef(fit_M21_JmaxHd)
f.m26.j<-coef(fit_M26_JmaxHd)
f.m31.j<-coef(fit_M31_JmaxHd)
f.a21.v<-coef(fit_A21_VcmaxHd)
f.a26.v<-coef(fit_A26_VcmaxHd)
f.a31.v<-coef(fit_A31_VcmaxHd)
f.a21.j<-coef(fit_A21_JmaxHd)
f.a26.j<-coef(fit_A26_JmaxHd)
f.a31.j<-coef(fit_A31_JmaxHd)
f.g21.v<-coef(fit_G21_VcmaxHd)
f.g26.v<-coef(fit_G26_VcmaxHd)
f.g31.v<-coef(fit_G31_VcmaxHd)
f.g21.j<-coef(fit_G21_JmaxHd)
f.g26.j<-coef(fit_G26_JmaxHd)
f.g31.j<-coef(fit_G31_JmaxHd)
f.r21.v<-coef(fit_R21_VcmaxHd)
f.r26.v<-coef(fit_R26_VcmaxHd)
f.r31.v<-coef(fit_R31_VcmaxHd)
f.r21.j<-coef(fit_R21_JmaxHd)
f.r26.j<-coef(fit_R26_JmaxHd)
f.r31.j<-coef(fit_R31_JmaxHd)

####
#Extract Topt estimate
####


#Vector of measurement temperatures
Tsim<-seq(10,40,0.01)

###
#Vcmax
###

#Morella
m21.vec<-peak.ar(1,f.m21.v[8],f.m21.v[9],200,Tsim+273.15)
Topt.m21.v<-Tsim[which(m21.vec==max(m21.vec))]
m26.vec<-peak.ar(1,f.m26.v[8],f.m26.v[9],200,Tsim+273.15)
Topt.m26.v<-Tsim[which(m26.vec==max(m26.vec))]
m31.vec<-peak.ar(1,f.m31.v[8],f.m31.v[9],200,Tsim+273.15)
Topt.m31.v<-Tsim[which(m31.vec==max(m31.vec))]

#Alnus
a21.vec<-peak.ar(1,f.a21.v[8],f.a21.v[9],200,Tsim+273.15)
Topt.a21.v<-Tsim[which(a21.vec==max(a21.vec))]
a26.vec<-peak.ar(1,f.a26.v[8],f.a26.v[9],200,Tsim+273.15)
Topt.a26.v<-Tsim[which(a26.vec==max(a26.vec))]
a31.vec<-peak.ar(1,f.a31.v[8],f.a31.v[9],200,Tsim+273.15)
Topt.a31.v<-Tsim[which(a31.vec==max(a31.vec))]

#Gliricidia
g21.vec<-peak.ar(1,f.g21.v[8],f.g21.v[9],200,Tsim+273.15)
Topt.g21.v<-Tsim[which(g21.vec==max(g21.vec))]
g26.vec<-peak.ar(1,f.g26.v[9],f.g26.v[10],200,Tsim+273.15)
Topt.g26.v<-Tsim[which(g26.vec==max(g26.vec))]
g31.vec<-peak.ar(1,f.g31.v[8],f.g31.v[9],200,Tsim+273.15)
Topt.g31.v<-Tsim[which(g31.vec==max(g31.vec))]

#Robinia
r21.vec<-peak.ar(1,f.r21.v[10],f.r21.v[11],200,Tsim+273.15)
Topt.r21.v<-Tsim[which(r21.vec==max(r21.vec))]
r26.vec<-peak.ar(1,f.r26.v[8],f.r26.v[9],200,Tsim+273.15)
Topt.r26.v<-Tsim[which(r26.vec==max(r26.vec))]
r31.vec<-peak.ar(1,f.r31.v[8],f.r31.v[9],200,Tsim+273.15)
Topt.r31.v<-Tsim[which(r31.vec==max(r31.vec))]

###
#Jmax
###

#Morella
m21.vec<-peak.ar(1,f.m21.j[8],f.m21.j[9],200,Tsim+273.15)
Topt.m21.j<-Tsim[which(m21.vec==max(m21.vec))]
m26.vec<-peak.ar(1,f.m26.j[8],f.m26.j[9],200,Tsim+273.15)
Topt.m26.j<-Tsim[which(m26.vec==max(m26.vec))]
m31.vec<-peak.ar(1,f.m31.j[8],f.m31.j[9],200,Tsim+273.15)
Topt.m31.j<-Tsim[which(m31.vec==max(m31.vec))]

#Alnus
a21.vec<-peak.ar(1,f.a21.j[8],f.a21.j[9],200,Tsim+273.15)
Topt.a21.j<-Tsim[which(a21.vec==max(a21.vec))]
a26.vec<-peak.ar(1,f.a26.j[8],f.a26.j[9],200,Tsim+273.15)
Topt.a26.j<-Tsim[which(a26.vec==max(a26.vec))]
a31.vec<-peak.ar(1,f.a31.j[8],f.a31.j[9],200,Tsim+273.15)
Topt.a31.j<-Tsim[which(a31.vec==max(a31.vec))]

#Gliricidia
g21.vec<-peak.ar(1,f.g21.j[8],f.g21.j[9],200,Tsim+273.15)
Topt.g21.j<-Tsim[which(g21.vec==max(g21.vec))]
g26.vec<-peak.ar(1,f.g26.j[9],f.g26.j[10],200,Tsim+273.15)
Topt.g26.j<-Tsim[which(g26.vec==max(g26.vec))]
g31.vec<-peak.ar(1,f.g31.j[8],f.g31.j[9],200,Tsim+273.15)
Topt.g31.j<-Tsim[which(g31.vec==max(g31.vec))]

#Robinia
r21.vec<-peak.ar(1,f.r21.j[10],f.r21.j[11],200,Tsim+273.15)
Topt.r21.j<-Tsim[which(r21.vec==max(r21.vec))]
r26.vec<-peak.ar(1,f.r26.j[8],f.r26.j[9],200,Tsim+273.15)
Topt.r26.j<-Tsim[which(r26.vec==max(r26.vec))]
r31.vec<-peak.ar(1,f.r31.j[8],f.r31.j[9],200,Tsim+273.15)
Topt.r31.j<-Tsim[which(r31.vec==max(r31.vec))]

####
#Compile Topt and Ea estimates by species
####

#Topt
Topt.m.v<-c(Topt.m21.v,Topt.m26.v,Topt.m31.v)
Topt.m.j<-c(Topt.m21.j,Topt.m26.j,Topt.m31.j)
Topt.a.v<-c(Topt.a21.v,Topt.a26.v,Topt.a31.v)
Topt.a.j<-c(Topt.a21.j,Topt.a26.j,Topt.a31.j)
Topt.g.v<-c(Topt.g21.v,Topt.g26.v,Topt.g31.v)
Topt.g.j<-c(Topt.g21.j,Topt.g26.j,Topt.g31.j)
Topt.r.v<-c(Topt.r21.v,Topt.r26.v,Topt.r31.v)
Topt.r.j<-c(Topt.r21.j,Topt.r26.j,Topt.r31.j)

#Ea
Ea.m.v<-c(f.m21.v[8],f.m26.v[8],f.m31.v[8])
Ea.m.j<-c(f.m21.j[8],f.m26.j[8],f.m31.j[8])
Ea.a.v<-c(f.a21.v[8],f.a26.v[8],f.a31.v[8])
Ea.a.j<-c(f.a21.j[8],f.a26.j[8],f.a31.j[8])
Ea.g.v<-c(f.g21.v[8],f.g26.v[9],f.g31.v[8])
Ea.g.j<-c(f.g21.j[8],f.g26.j[9],f.g31.j[8])
Ea.r.v<-c(f.r21.v[10],f.r26.v[8],f.r31.v[8])
Ea.r.j<-c(f.r21.j[10],f.r26.j[8],f.r31.j[8])

####
#95% CI
####

###
#Morella
###

##
#Vcmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m21.v=mvrnorm(1000,mu=f.m21.v,Sigma=vcov(fit_M21_VcmaxHd))

#Vcmax ~ temperature
dist.m21.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m21.v<-rep(NA,length(Tsim))
high.m21.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m21.v[i,j]=peak.ar(1,vmat.m21.v[i,8],vmat.m21.v[i,9],200,Tsim[j]+273.15)
  }
  low.m21.v[j]<-quantile(na.omit(dist.m21.v[,j]),0.025)
  high.m21.v[j]<-quantile(na.omit(dist.m21.v[,j]),0.975)
}

#Ea
dist.m21.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.m21.v.Ea[i]=vmat.m21.v[i,8]
}
low.Ea.m21.v<-quantile(na.omit(dist.m21.v.Ea),0.025)
high.Ea.m21.v<-quantile(na.omit(dist.m21.v.Ea),0.975)

#Topt
dist.m21.v.Topt=rep(NA,1000)
for(i in 1:1000){
  m21.vec<-peak.ar(1,vmat.m21.v[i,8],vmat.m21.v[i,9],200,Tsim+273.15)
  dist.m21.v.Topt[i]<-Tsim[which(m21.vec==max(m21.vec))]
}
low.Topt.m21.v<-quantile(na.omit(dist.m21.v.Topt),0.025)
high.Topt.m21.v<-quantile(na.omit(dist.m21.v.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m26.v=mvrnorm(1000,mu=f.m26.v,Sigma=vcov(fit_M26_VcmaxHd))

#Vcmax ~ temperature
dist.m26.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m26.v<-rep(NA,length(Tsim))
high.m26.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m26.v[i,j]=peak.ar(1,vmat.m26.v[i,8],vmat.m26.v[i,9],200,Tsim[j]+273.15)
  }
  low.m26.v[j]<-quantile(na.omit(dist.m26.v[,j]),0.025)
  high.m26.v[j]<-quantile(na.omit(dist.m26.v[,j]),0.975)
}

#Ea
dist.m26.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.m26.v.Ea[i]=vmat.m26.v[i,8]
}
low.Ea.m26.v<-quantile(na.omit(dist.m26.v.Ea),0.025)
high.Ea.m26.v<-quantile(na.omit(dist.m26.v.Ea),0.975)

#Topt
dist.m26.v.Topt=rep(NA,1000)
for(i in 1:1000){
  m26.vec<-peak.ar(1,vmat.m26.v[i,8],vmat.m26.v[i,9],200,Tsim+273.15)
  dist.m26.v.Topt[i]<-Tsim[which(m26.vec==max(m26.vec))]
}
low.Topt.m26.v<-quantile(na.omit(dist.m26.v.Topt),0.025)
high.Topt.m26.v<-quantile(na.omit(dist.m26.v.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m31.v=mvrnorm(1000,mu=f.m31.v,Sigma=vcov(fit_M31_VcmaxHd))

#Vcmax ~ temperature
dist.m31.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m31.v<-rep(NA,length(Tsim))
high.m31.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m31.v[i,j]=peak.ar(1,vmat.m31.v[i,8],vmat.m31.v[i,9],200,Tsim[j]+273.15)
  }
  low.m31.v[j]<-quantile(na.omit(dist.m31.v[,j]),0.025)
  high.m31.v[j]<-quantile(na.omit(dist.m31.v[,j]),0.975)
}

#Ea
dist.m31.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.m31.v.Ea[i]=vmat.m31.v[i,8]
}
low.Ea.m31.v<-quantile(na.omit(dist.m31.v.Ea),0.025)
high.Ea.m31.v<-quantile(na.omit(dist.m31.v.Ea),0.975)

#Topt
dist.m31.v.Topt=rep(NA,1000)
for(i in 1:1000){
  m31.vec<-peak.ar(1,vmat.m31.v[i,8],vmat.m31.v[i,9],200,Tsim+273.15)
  dist.m31.v.Topt[i]<-Tsim[which(m31.vec==max(m31.vec))]
}
low.Topt.m31.v<-quantile(na.omit(dist.m31.v.Topt),0.025)
high.Topt.m31.v<-quantile(na.omit(dist.m31.v.Topt),0.975)

##
#Jmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m21.j=mvrnorm(1000,mu=f.m21.j,Sigma=vcov(fit_M21_JmaxHd))

#Jmax ~ temperature
dist.m21.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m21.j<-rep(NA,length(Tsim))
high.m21.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m21.j[i,j]=peak.ar(1,vmat.m21.j[i,8],vmat.m21.j[i,9],200,Tsim[j]+273.15)
  }
  low.m21.j[j]<-quantile(na.omit(dist.m21.j[,j]),0.025)
  high.m21.j[j]<-quantile(na.omit(dist.m21.j[,j]),0.975)
}

#Ea
dist.m21.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.m21.j.Ea[i]=vmat.m21.j[i,8]
}
low.Ea.m21.j<-quantile(na.omit(dist.m21.j.Ea),0.025)
high.Ea.m21.j<-quantile(na.omit(dist.m21.j.Ea),0.975)

#Topt
dist.m21.j.Topt=rep(NA,1000)
for(i in 1:1000){
  m21.vec<-peak.ar(1,vmat.m21.j[i,8],vmat.m21.j[i,9],200,Tsim+273.15)
  dist.m21.j.Topt[i]<-Tsim[which(m21.vec==max(m21.vec))]
}
low.Topt.m21.j<-quantile(na.omit(dist.m21.j.Topt),0.025)
high.Topt.m21.j<-quantile(na.omit(dist.m21.j.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m26.j=mvrnorm(1000,mu=f.m26.j,Sigma=vcov(fit_M26_JmaxHd))

#Jmax ~ temperature
dist.m26.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m26.j<-rep(NA,length(Tsim))
high.m26.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m26.j[i,j]=peak.ar(1,vmat.m26.j[i,8],vmat.m26.j[i,9],200,Tsim[j]+273.15)
  }
  low.m26.j[j]<-quantile(na.omit(dist.m26.j[,j]),0.025)
  high.m26.j[j]<-quantile(na.omit(dist.m26.j[,j]),0.975)
}

#Ea
dist.m26.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.m26.j.Ea[i]=vmat.m26.j[i,8]
}
low.Ea.m26.j<-quantile(na.omit(dist.m26.j.Ea),0.025)
high.Ea.m26.j<-quantile(na.omit(dist.m26.j.Ea),0.975)

#Topt
dist.m26.j.Topt=rep(NA,1000)
for(i in 1:1000){
  m26.vec<-peak.ar(1,vmat.m26.j[i,8],vmat.m26.j[i,9],200,Tsim+273.15)
  dist.m26.j.Topt[i]<-Tsim[which(m26.vec==max(m26.vec))]
}
low.Topt.m26.j<-quantile(na.omit(dist.m26.j.Topt),0.025)
high.Topt.m26.j<-quantile(na.omit(dist.m26.j.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m31.j=mvrnorm(1000,mu=f.m31.j,Sigma=vcov(fit_M31_JmaxHd))

#Jmax ~ temperature
dist.m31.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m31.j<-rep(NA,length(Tsim))
high.m31.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m31.j[i,j]=peak.ar(1,vmat.m31.j[i,8],vmat.m31.j[i,9],200,Tsim[j]+273.15)
  }
  low.m31.j[j]<-quantile(na.omit(dist.m31.j[,j]),0.025)
  high.m31.j[j]<-quantile(na.omit(dist.m31.j[,j]),0.975)
}

#Ea
dist.m31.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.m31.j.Ea[i]=vmat.m31.j[i,8]
}
low.Ea.m31.j<-quantile(na.omit(dist.m31.j.Ea),0.025)
high.Ea.m31.j<-quantile(na.omit(dist.m31.j.Ea),0.975)

#Topt
dist.m31.j.Topt=rep(NA,1000)
for(i in 1:1000){
  m31.vec<-peak.ar(1,vmat.m31.j[i,8],vmat.m31.j[i,9],200,Tsim+273.15)
  dist.m31.j.Topt[i]<-Tsim[which(m31.vec==max(m31.vec))]
}
low.Topt.m31.j<-quantile(na.omit(dist.m31.j.Topt),0.025)
high.Topt.m31.j<-quantile(na.omit(dist.m31.j.Topt),0.975)

###
#Alnus
###

##
#Vcmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a21.v=mvrnorm(1000,mu=f.a21.v,Sigma=vcov(fit_A21_VcmaxHd))

#Vcmax ~ temperature
dist.a21.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a21.v<-rep(NA,length(Tsim))
high.a21.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a21.v[i,j]=peak.ar(1,vmat.a21.v[i,8],vmat.a21.v[i,9],200,Tsim[j]+273.15)
  }
  low.a21.v[j]<-quantile(na.omit(dist.a21.v[,j]),0.025)
  high.a21.v[j]<-quantile(na.omit(dist.a21.v[,j]),0.975)
}

#Ea
dist.a21.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.a21.v.Ea[i]=vmat.a21.v[i,8]
}
low.Ea.a21.v<-quantile(na.omit(dist.a21.v.Ea),0.025)
high.Ea.a21.v<-quantile(na.omit(dist.a21.v.Ea),0.975)

#Topt
dist.a21.v.Topt=rep(NA,1000)
for(i in 1:1000){
  a21.vec<-peak.ar(1,vmat.a21.v[i,8],vmat.a21.v[i,9],200,Tsim+273.15)
  dist.a21.v.Topt[i]<-Tsim[which(a21.vec==max(a21.vec))]
}
low.Topt.a21.v<-quantile(na.omit(dist.a21.v.Topt),0.025)
high.Topt.a21.v<-quantile(na.omit(dist.a21.v.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a26.v=mvrnorm(1000,mu=f.a26.v,Sigma=vcov(fit_A26_VcmaxHd))

#Vcmax ~ temperature
dist.a26.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a26.v<-rep(NA,length(Tsim))
high.a26.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a26.v[i,j]=peak.ar(1,vmat.a26.v[i,8],vmat.a26.v[i,9],200,Tsim[j]+273.15)
  }
  low.a26.v[j]<-quantile(na.omit(dist.a26.v[,j]),0.025)
  high.a26.v[j]<-quantile(na.omit(dist.a26.v[,j]),0.975)
}

#Ea
dist.a26.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.a26.v.Ea[i]=vmat.a26.v[i,8]
}
low.Ea.a26.v<-quantile(na.omit(dist.a26.v.Ea),0.025)
high.Ea.a26.v<-quantile(na.omit(dist.a26.v.Ea),0.975)

#Topt
dist.a26.v.Topt=rep(NA,1000)
for(i in 1:1000){
  a26.vec<-peak.ar(1,vmat.a26.v[i,8],vmat.a26.v[i,9],200,Tsim+273.15)
  dist.a26.v.Topt[i]<-Tsim[which(a26.vec==max(a26.vec))]
}
low.Topt.a26.v<-quantile(na.omit(dist.a26.v.Topt),0.025)
high.Topt.a26.v<-quantile(na.omit(dist.a26.v.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a31.v=mvrnorm(1000,mu=f.a31.v,Sigma=vcov(fit_A31_VcmaxHd))

#Vcmax ~ temperature
dist.a31.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a31.v<-rep(NA,length(Tsim))
high.a31.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a31.v[i,j]=peak.ar(1,vmat.a31.v[i,8],vmat.a31.v[i,9],200,Tsim[j]+273.15)
  }
  low.a31.v[j]<-quantile(na.omit(dist.a31.v[,j]),0.025)
  high.a31.v[j]<-quantile(na.omit(dist.a31.v[,j]),0.975)
}

#Ea
dist.a31.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.a31.v.Ea[i]=vmat.a31.v[i,8]
}
low.Ea.a31.v<-quantile(na.omit(dist.a31.v.Ea),0.025)
high.Ea.a31.v<-quantile(na.omit(dist.a31.v.Ea),0.975)

#Topt
dist.a31.v.Topt=rep(NA,1000)
for(i in 1:1000){
  a31.vec<-peak.ar(1,vmat.a31.v[i,8],vmat.a31.v[i,9],200,Tsim+273.15)
  dist.a31.v.Topt[i]<-Tsim[which(a31.vec==max(a31.vec))]
}
low.Topt.a31.v<-quantile(na.omit(dist.a31.v.Topt),0.025)
high.Topt.a31.v<-quantile(na.omit(dist.a31.v.Topt),0.975)

##
#Jmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a21.j=mvrnorm(1000,mu=f.a21.j,Sigma=vcov(fit_A21_JmaxHd))

#Jmax ~ temperature
dist.a21.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a21.j<-rep(NA,length(Tsim))
high.a21.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a21.j[i,j]=peak.ar(1,vmat.a21.j[i,8],vmat.a21.j[i,9],200,Tsim[j]+273.15)
  }
  low.a21.j[j]<-quantile(na.omit(dist.a21.j[,j]),0.025)
  high.a21.j[j]<-quantile(na.omit(dist.a21.j[,j]),0.975)
}

#Ea
dist.a21.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.a21.j.Ea[i]=vmat.a21.j[i,8]
}
low.Ea.a21.j<-quantile(na.omit(dist.a21.j.Ea),0.025)
high.Ea.a21.j<-quantile(na.omit(dist.a21.j.Ea),0.975)

#Topt
dist.a21.j.Topt=rep(NA,1000)
for(i in 1:1000){
  a21.vec<-peak.ar(1,vmat.a21.j[i,8],vmat.a21.j[i,9],200,Tsim+273.15)
  dist.a21.j.Topt[i]<-Tsim[which(a21.vec==max(a21.vec))]
}
low.Topt.a21.j<-quantile(na.omit(dist.a21.j.Topt),0.025)
high.Topt.a21.j<-quantile(na.omit(dist.a21.j.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a26.j=mvrnorm(1000,mu=f.a26.j,Sigma=vcov(fit_A26_JmaxHd))

#Jmax ~ temperature
dist.a26.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a26.j<-rep(NA,length(Tsim))
high.a26.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a26.j[i,j]=peak.ar(1,vmat.a26.j[i,8],vmat.a26.j[i,9],200,Tsim[j]+273.15)
  }
  low.a26.j[j]<-quantile(na.omit(dist.a26.j[,j]),0.025)
  high.a26.j[j]<-quantile(na.omit(dist.a26.j[,j]),0.975)
}

#Ea
dist.a26.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.a26.j.Ea[i]=vmat.a26.j[i,8]
}
low.Ea.a26.j<-quantile(na.omit(dist.a26.j.Ea),0.025)
high.Ea.a26.j<-quantile(na.omit(dist.a26.j.Ea),0.975)

#Topt
dist.a26.j.Topt=rep(NA,1000)
for(i in 1:1000){
  a26.vec<-peak.ar(1,vmat.a26.j[i,8],vmat.a26.j[i,9],200,Tsim+273.15)
  dist.a26.j.Topt[i]<-Tsim[which(a26.vec==max(a26.vec))]
}
low.Topt.a26.j<-quantile(na.omit(dist.a26.j.Topt),0.025)
high.Topt.a26.j<-quantile(na.omit(dist.a26.j.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a31.j=mvrnorm(1000,mu=f.a31.j,Sigma=vcov(fit_A31_JmaxHd))

#Jmax ~ temperature
dist.a31.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a31.j<-rep(NA,length(Tsim))
high.a31.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a31.j[i,j]=peak.ar(1,vmat.a31.j[i,8],vmat.a31.j[i,9],200,Tsim[j]+273.15)
  }
  low.a31.j[j]<-quantile(na.omit(dist.a31.j[,j]),0.025)
  high.a31.j[j]<-quantile(na.omit(dist.a31.j[,j]),0.975)
}

#Ea
dist.a31.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.a31.j.Ea[i]=vmat.a31.j[i,8]
}
low.Ea.a31.j<-quantile(na.omit(dist.a31.j.Ea),0.025)
high.Ea.a31.j<-quantile(na.omit(dist.a31.j.Ea),0.975)

#Topt
dist.a31.j.Topt=rep(NA,1000)
for(i in 1:1000){
  a31.vec<-peak.ar(1,vmat.a31.j[i,8],vmat.a31.j[i,9],200,Tsim+273.15)
  dist.a31.j.Topt[i]<-Tsim[which(a31.vec==max(a31.vec))]
}
low.Topt.a31.j<-quantile(na.omit(dist.a31.j.Topt),0.025)
high.Topt.a31.j<-quantile(na.omit(dist.a31.j.Topt),0.975)

###
#Gliricidia
###

##
#Vcmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g21.v=mvrnorm(1000,mu=f.g21.v,Sigma=vcov(fit_G21_VcmaxHd))

#Vcmax ~ temperature
dist.g21.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g21.v<-rep(NA,length(Tsim))
high.g21.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g21.v[i,j]=peak.ar(1,vmat.g21.v[i,8],vmat.g21.v[i,9],200,Tsim[j]+273.15)
  }
  low.g21.v[j]<-quantile(na.omit(dist.g21.v[,j]),0.025)
  high.g21.v[j]<-quantile(na.omit(dist.g21.v[,j]),0.975)
}

#Ea
dist.g21.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.g21.v.Ea[i]=vmat.g21.v[i,8]
}
low.Ea.g21.v<-quantile(na.omit(dist.g21.v.Ea),0.025)
high.Ea.g21.v<-quantile(na.omit(dist.g21.v.Ea),0.975)

#Topt
dist.g21.v.Topt=rep(NA,1000)
for(i in 1:1000){
  g21.vec<-peak.ar(1,vmat.g21.v[i,8],vmat.g21.v[i,9],200,Tsim+273.15)
  dist.g21.v.Topt[i]<-Tsim[which(g21.vec==max(g21.vec))]
}
low.Topt.g21.v<-quantile(na.omit(dist.g21.v.Topt),0.025)
high.Topt.g21.v<-quantile(na.omit(dist.g21.v.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g26.v=mvrnorm(1000,mu=f.g26.v,Sigma=vcov(fit_G26_VcmaxHd))

#Vcmax ~ temperature
dist.g26.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g26.v<-rep(NA,length(Tsim))
high.g26.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g26.v[i,j]=peak.ar(1,vmat.g26.v[i,9],vmat.g26.v[i,10],200,Tsim[j]+273.15)
  }
  low.g26.v[j]<-quantile(na.omit(dist.g26.v[,j]),0.025)
  high.g26.v[j]<-quantile(na.omit(dist.g26.v[,j]),0.975)
}

#Ea
dist.g26.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.g26.v.Ea[i]=vmat.g26.v[i,9]
}
low.Ea.g26.v<-quantile(na.omit(dist.g26.v.Ea),0.025)
high.Ea.g26.v<-quantile(na.omit(dist.g26.v.Ea),0.975)

#Topt
dist.g26.v.Topt=rep(NA,1000)
for(i in 1:1000){
  g26.vec<-peak.ar(1,vmat.g26.v[i,9],vmat.g26.v[i,10],200,Tsim+273.15)
  dist.g26.v.Topt[i]<-Tsim[which(g26.vec==max(g26.vec))]
}
low.Topt.g26.v<-quantile(na.omit(dist.g26.v.Topt),0.025)
high.Topt.g26.v<-quantile(na.omit(dist.g26.v.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g31.v=mvrnorm(1000,mu=f.g31.v,Sigma=vcov(fit_G31_VcmaxHd))

#Vcmax ~ temperature
dist.g31.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g31.v<-rep(NA,length(Tsim))
high.g31.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g31.v[i,j]=peak.ar(1,vmat.g31.v[i,8],vmat.g31.v[i,9],200,Tsim[j]+273.15)
  }
  low.g31.v[j]<-quantile(na.omit(dist.g31.v[,j]),0.025)
  high.g31.v[j]<-quantile(na.omit(dist.g31.v[,j]),0.975)
}

#Ea
dist.g31.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.g31.v.Ea[i]=vmat.g31.v[i,8]
}
low.Ea.g31.v<-quantile(na.omit(dist.g31.v.Ea),0.025)
high.Ea.g31.v<-quantile(na.omit(dist.g31.v.Ea),0.975)

#Topt
dist.g31.v.Topt=rep(NA,1000)
for(i in 1:1000){
  g31.vec<-peak.ar(1,vmat.g31.v[i,8],vmat.g31.v[i,9],200,Tsim+273.15)
  dist.g31.v.Topt[i]<-Tsim[which(g31.vec==max(g31.vec))]
}
low.Topt.g31.v<-quantile(na.omit(dist.g31.v.Topt),0.025)
high.Topt.g31.v<-quantile(na.omit(dist.g31.v.Topt),0.975)

##
#Jmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g21.j=mvrnorm(1000,mu=f.g21.j,Sigma=vcov(fit_G21_JmaxHd))

#Jmax ~ temperature
dist.g21.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g21.j<-rep(NA,length(Tsim))
high.g21.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g21.j[i,j]=peak.ar(1,vmat.g21.j[i,8],vmat.g21.j[i,9],200,Tsim[j]+273.15)
  }
  low.g21.j[j]<-quantile(na.omit(dist.g21.j[,j]),0.025)
  high.g21.j[j]<-quantile(na.omit(dist.g21.j[,j]),0.975)
}

#Ea
dist.g21.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.g21.j.Ea[i]=vmat.g21.j[i,8]
}
low.Ea.g21.j<-quantile(na.omit(dist.g21.j.Ea),0.025)
high.Ea.g21.j<-quantile(na.omit(dist.g21.j.Ea),0.975)

#Topt
dist.g21.j.Topt=rep(NA,1000)
for(i in 1:1000){
  g21.vec<-peak.ar(1,vmat.g21.j[i,8],vmat.g21.j[i,9],200,Tsim+273.15)
  dist.g21.j.Topt[i]<-Tsim[which(g21.vec==max(g21.vec))]
}
low.Topt.g21.j<-quantile(na.omit(dist.g21.j.Topt),0.025)
high.Topt.g21.j<-quantile(na.omit(dist.g21.j.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g26.j=mvrnorm(1000,mu=f.g26.j,Sigma=vcov(fit_G26_JmaxHd))

#Jmax ~ temperature
dist.g26.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g26.j<-rep(NA,length(Tsim))
high.g26.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g26.j[i,j]=peak.ar(1,vmat.g26.j[i,9],vmat.g26.j[i,10],200,Tsim[j]+273.15)
  }
  low.g26.j[j]<-quantile(na.omit(dist.g26.j[,j]),0.025)
  high.g26.j[j]<-quantile(na.omit(dist.g26.j[,j]),0.975)
}

#Ea
dist.g26.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.g26.j.Ea[i]=vmat.g26.j[i,9]
}
low.Ea.g26.j<-quantile(na.omit(dist.g26.j.Ea),0.025)
high.Ea.g26.j<-quantile(na.omit(dist.g26.j.Ea),0.975)

#Topt
dist.g26.j.Topt=rep(NA,1000)
for(i in 1:1000){
  g26.vec<-peak.ar(1,vmat.g26.j[i,9],vmat.g26.j[i,10],200,Tsim+273.15)
  dist.g26.j.Topt[i]<-Tsim[which(g26.vec==max(g26.vec))]
}
low.Topt.g26.j<-quantile(na.omit(dist.g26.j.Topt),0.025)
high.Topt.g26.j<-quantile(na.omit(dist.g26.j.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g31.j=mvrnorm(1000,mu=f.g31.j,Sigma=vcov(fit_G31_JmaxHd))

#Jmax ~ temperature
dist.g31.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g31.j<-rep(NA,length(Tsim))
high.g31.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g31.j[i,j]=peak.ar(1,vmat.g31.j[i,8],vmat.g31.j[i,9],200,Tsim[j]+273.15)
  }
  low.g31.j[j]<-quantile(na.omit(dist.g31.j[,j]),0.025)
  high.g31.j[j]<-quantile(na.omit(dist.g31.j[,j]),0.975)
}

#Ea
dist.g31.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.g31.j.Ea[i]=vmat.g31.j[i,8]
}
low.Ea.g31.j<-quantile(na.omit(dist.g31.j.Ea),0.025)
high.Ea.g31.j<-quantile(na.omit(dist.g31.j.Ea),0.975)

#Topt
dist.g31.j.Topt=rep(NA,1000)
for(i in 1:1000){
  g31.vec<-peak.ar(1,vmat.g31.j[i,8],vmat.g31.j[i,9],200,Tsim+273.15)
  dist.g31.j.Topt[i]<-Tsim[which(g31.vec==max(g31.vec))]
}
low.Topt.g31.j<-quantile(na.omit(dist.g31.j.Topt),0.025)
high.Topt.g31.j<-quantile(na.omit(dist.g31.j.Topt),0.975)

###
#Robinia
###

##
#Vcmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r21.v=mvrnorm(1000,mu=f.r21.v,Sigma=vcov(fit_R21_VcmaxHd))

#Vcmax ~ temperature
dist.r21.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r21.v<-rep(NA,length(Tsim))
high.r21.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r21.v[i,j]=peak.ar(1,vmat.r21.v[i,10],vmat.r21.v[i,11],200,Tsim[j]+273.15)
  }
  low.r21.v[j]<-quantile(na.omit(dist.r21.v[,j]),0.025)
  high.r21.v[j]<-quantile(na.omit(dist.r21.v[,j]),0.975)
}

#Ea
dist.r21.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.r21.v.Ea[i]=vmat.r21.v[i,10]
}
low.Ea.r21.v<-quantile(na.omit(dist.r21.v.Ea),0.025)
high.Ea.r21.v<-quantile(na.omit(dist.r21.v.Ea),0.975)

#Topt
dist.r21.v.Topt=rep(NA,1000)
for(i in 1:1000){
  r21.vec<-peak.ar(1,vmat.r21.v[i,10],vmat.r21.v[i,11],200,Tsim+273.15)
  dist.r21.v.Topt[i]<-Tsim[which(r21.vec==max(r21.vec))]
}
low.Topt.r21.v<-quantile(na.omit(dist.r21.v.Topt),0.025)
high.Topt.r21.v<-quantile(na.omit(dist.r21.v.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r26.v=mvrnorm(1000,mu=f.r26.v,Sigma=vcov(fit_R26_VcmaxHd))

#Vcmax ~ temperature
dist.r26.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r26.v<-rep(NA,length(Tsim))
high.r26.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r26.v[i,j]=peak.ar(1,vmat.r26.v[i,8],vmat.r26.v[i,9],200,Tsim[j]+273.15)
  }
  low.r26.v[j]<-quantile(na.omit(dist.r26.v[,j]),0.025)
  high.r26.v[j]<-quantile(na.omit(dist.r26.v[,j]),0.975)
}

#Ea
dist.r26.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.r26.v.Ea[i]=vmat.r26.v[i,8]
}
low.Ea.r26.v<-quantile(na.omit(dist.r26.v.Ea),0.025)
high.Ea.r26.v<-quantile(na.omit(dist.r26.v.Ea),0.975)

#Topt
dist.r26.v.Topt=rep(NA,1000)
for(i in 1:1000){
  r26.vec<-peak.ar(1,vmat.r26.v[i,8],vmat.r26.v[i,9],200,Tsim+273.15)
  dist.r26.v.Topt[i]<-Tsim[which(na.omit(r26.vec)==max(na.omit(r26.vec)))]
}
low.Topt.r26.v<-quantile(na.omit(dist.r26.v.Topt),0.025)
high.Topt.r26.v<-quantile(na.omit(dist.r26.v.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r31.v=mvrnorm(1000,mu=f.r31.v,Sigma=vcov(fit_R31_VcmaxHd))

#Vcmax ~ temperature
dist.r31.v=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r31.v<-rep(NA,length(Tsim))
high.r31.v<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r31.v[i,j]=peak.ar(1,vmat.r31.v[i,8],vmat.r31.v[i,9],200,Tsim[j]+273.15)
  }
  low.r31.v[j]<-quantile(na.omit(dist.r31.v[,j]),0.025)
  high.r31.v[j]<-quantile(na.omit(dist.r31.v[,j]),0.975)
}

#Ea
dist.r31.v.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.r31.v.Ea[i]=vmat.r31.v[i,8]
}
low.Ea.r31.v<-quantile(na.omit(dist.r31.v.Ea),0.025)
high.Ea.r31.v<-quantile(na.omit(dist.r31.v.Ea),0.975)

#Topt
dist.r31.v.Topt=rep(NA,1000)
for(i in 1:1000){
  r31.vec<-peak.ar(1,vmat.r31.v[i,8],vmat.r31.v[i,9],200,Tsim+273.15)
  dist.r31.v.Topt[i]<-Tsim[which(r31.vec==max(r31.vec))]
}
low.Topt.r31.v<-quantile(na.omit(dist.r31.v.Topt),0.025)
high.Topt.r31.v<-quantile(na.omit(dist.r31.v.Topt),0.975)

##
#Jmax
##

#
#21:15 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r21.j=mvrnorm(1000,mu=f.r21.j,Sigma=vcov(fit_R21_JmaxHd))

#Jmax ~ temperature
dist.r21.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r21.j<-rep(NA,length(Tsim))
high.r21.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r21.j[i,j]=peak.ar(1,vmat.r21.j[i,10],vmat.r21.j[i,11],200,Tsim[j]+273.15)
  }
  low.r21.j[j]<-quantile(na.omit(dist.r21.j[,j]),0.025)
  high.r21.j[j]<-quantile(na.omit(dist.r21.j[,j]),0.975)
}

#Ea
dist.r21.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.r21.j.Ea[i]=vmat.r21.j[i,10]
}
low.Ea.r21.j<-quantile(na.omit(dist.r21.j.Ea),0.025)
high.Ea.r21.j<-quantile(na.omit(dist.r21.j.Ea),0.975)

#Topt
dist.r21.j.Topt=rep(NA,1000)
for(i in 1:1000){
  r21.vec<-peak.ar(1,vmat.r21.j[i,10],vmat.r21.j[i,11],200,Tsim+273.15)
  dist.r21.j.Topt[i]<-Tsim[which(r21.vec==max(r21.vec))]
}
low.Topt.r21.j<-quantile(na.omit(dist.r21.j.Topt),0.025)
high.Topt.r21.j<-quantile(na.omit(dist.r21.j.Topt),0.975)

#
#26:20 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r26.j=mvrnorm(1000,mu=f.r26.j,Sigma=vcov(fit_R26_JmaxHd))

#Jmax ~ temperature
dist.r26.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r26.j<-rep(NA,length(Tsim))
high.r26.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r26.j[i,j]=peak.ar(1,vmat.r26.j[i,8],vmat.r26.j[i,9],200,Tsim[j]+273.15)
  }
  low.r26.j[j]<-quantile(na.omit(dist.r26.j[,j]),0.025)
  high.r26.j[j]<-quantile(na.omit(dist.r26.j[,j]),0.975)
}

#Ea
dist.r26.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.r26.j.Ea[i]=vmat.r26.j[i,8]
}
low.Ea.r26.j<-quantile(na.omit(dist.r26.j.Ea),0.025)
high.Ea.r26.j<-quantile(na.omit(dist.r26.j.Ea),0.975)

#Topt
dist.r26.j.Topt=rep(NA,1000)
for(i in 1:1000){
  r26.vec<-peak.ar(1,vmat.r26.j[i,8],vmat.r26.j[i,9],200,Tsim+273.15)
  dist.r26.j.Topt[i]<-Tsim[which(r26.vec==max(r26.vec))]
}
low.Topt.r26.j<-quantile(na.omit(dist.r26.j.Topt),0.025)
high.Topt.r26.j<-quantile(na.omit(dist.r26.j.Topt),0.975)

#
#31:25 deg. C growing temperature
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r31.j=mvrnorm(1000,mu=f.r31.j,Sigma=vcov(fit_R31_JmaxHd))

#Jmax ~ temperature
dist.r31.j=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r31.j<-rep(NA,length(Tsim))
high.r31.j<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r31.j[i,j]=peak.ar(1,vmat.r31.j[i,8],vmat.r31.j[i,9],200,Tsim[j]+273.15)
  }
  low.r31.j[j]<-quantile(na.omit(dist.r31.j[,j]),0.025)
  high.r31.j[j]<-quantile(na.omit(dist.r31.j[,j]),0.975)
}

#Ea
dist.r31.j.Ea=rep(NA,1000)
for(i in 1:1000){
  dist.r31.j.Ea[i]=vmat.r31.j[i,8]
}
low.Ea.r31.j<-quantile(na.omit(dist.r31.j.Ea),0.025)
high.Ea.r31.j<-quantile(na.omit(dist.r31.j.Ea),0.975)

#Topt
dist.r31.j.Topt=rep(NA,1000)
for(i in 1:1000){
  r31.vec<-peak.ar(1,vmat.r31.j[i,8],vmat.r31.j[i,9],200,Tsim+273.15)
  dist.r31.j.Topt[i]<-Tsim[which(r31.vec==max(r31.vec))]
}
low.Topt.r31.j<-quantile(na.omit(dist.r31.j.Topt),0.025)
high.Topt.r31.j<-quantile(na.omit(dist.r31.j.Topt),0.975)

####
#Compile Topt and Ea 95% CI by species
####

###
#Vcmax
###

#Topt
l.Topt.m.v<-c(low.Topt.m21.v,low.Topt.m26.v,low.Topt.m31.v)
h.Topt.m.v<-c(high.Topt.m21.v,high.Topt.m26.v,high.Topt.m31.v)
l.Topt.a.v<-c(low.Topt.a21.v,low.Topt.a26.v,low.Topt.a31.v)
h.Topt.a.v<-c(high.Topt.a21.v,high.Topt.a26.v,high.Topt.a31.v)
l.Topt.g.v<-c(low.Topt.g21.v,low.Topt.g26.v,NA)
h.Topt.g.v<-c(high.Topt.g21.v,high.Topt.g26.v,high.Topt.g31.v)
l.Topt.r.v<-c(low.Topt.r21.v,low.Topt.r26.v,low.Topt.r31.v)
h.Topt.r.v<-c(high.Topt.r21.v,high.Topt.r26.v,high.Topt.r31.v)

#Ea
l.Ea.m.v<-c(low.Ea.m21.v,low.Ea.m26.v,low.Ea.m31.v)
h.Ea.m.v<-c(high.Ea.m21.v,high.Ea.m26.v,high.Ea.m31.v)
l.Ea.a.v<-c(low.Ea.a21.v,low.Ea.a26.v,low.Ea.a31.v)
h.Ea.a.v<-c(high.Ea.a21.v,high.Ea.a26.v,high.Ea.a31.v)
l.Ea.g.v<-c(low.Ea.g21.v,low.Ea.g26.v,low.Ea.g31.v)
h.Ea.g.v<-c(high.Ea.g21.v,high.Ea.g26.v,NA)
l.Ea.r.v<-c(low.Ea.r21.v,low.Ea.r26.v,low.Ea.r31.v)
h.Ea.r.v<-c(high.Ea.r21.v,high.Ea.r26.v,high.Ea.r31.v)

###
#Jmax
###

#Topt
l.Topt.m.j<-c(low.Topt.m21.j,low.Topt.m26.j,low.Topt.m31.j)
h.Topt.m.j<-c(high.Topt.m21.j,high.Topt.m26.j,high.Topt.m31.j)
l.Topt.a.j<-c(low.Topt.a21.j,low.Topt.a26.j,low.Topt.a31.j)
h.Topt.a.j<-c(high.Topt.a21.j,high.Topt.a26.j,high.Topt.a31.j)
l.Topt.g.j<-c(low.Topt.g21.j,low.Topt.g26.j,low.Topt.g31.j)
h.Topt.g.j<-c(high.Topt.g21.j,high.Topt.g26.j,high.Topt.g31.j)
l.Topt.r.j<-c(low.Topt.r21.j,low.Topt.r26.j,low.Topt.r31.j)
h.Topt.r.j<-c(high.Topt.r21.j,high.Topt.r26.j,high.Topt.r31.j)

#Ea
l.Ea.m.j<-c(low.Ea.m21.j,low.Ea.m26.j,low.Ea.m31.j)
h.Ea.m.j<-c(high.Ea.m21.j,high.Ea.m26.j,high.Ea.m31.j)
l.Ea.a.j<-c(low.Ea.a21.j,low.Ea.a26.j,low.Ea.a31.j)
h.Ea.a.j<-c(high.Ea.a21.j,high.Ea.a26.j,high.Ea.a31.j)
l.Ea.g.j<-c(low.Ea.g21.j,low.Ea.g26.j,low.Ea.g31.j)
h.Ea.g.j<-c(high.Ea.g21.j,high.Ea.g26.j,high.Ea.g31.j)
l.Ea.r.j<-c(low.Ea.r21.j,low.Ea.r26.j,low.Ea.r31.j)
h.Ea.r.j<-c(high.Ea.r21.j,high.Ea.r26.j,high.Ea.r31.j)

###############################################################################################################
#Supplementary Figure 7
###############################################################################################################

#PDF dimension is 7x13 inches  

#Plotting Settings
par(mfrow=c(2,1))
par(mar=c(0.3,0.8,0.8,0.5))
par(oma=c(5,5.5,1,0.5))
par(pty="s")

#S7a
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(283.15,313.15),ylim=c(0,3.3),cex.lab=1.5,cex.axis=1.2,xaxt="n")
axis(1,at=c(283.15,288.15,293.15,298.15,303.15,308.15,313.15),labels=F,cex.axis=1.2)
curve(peak.ar(1,f.m21.v[8],f.m21.v[9],200,x),from=283.15,to=313.15,lty=1,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.a21.v[8],f.a21.v[9],200,x),from=283.15,to=313.15,lty=2,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.g21.v[8],f.g21.v[9],200,x),from=283.15,to=313.15,lty=3,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.r21.v[10],f.r21.v[11],200,x),from=283.15,to=313.15,lty=4,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.m26.v[8],f.m26.v[9],200,x),from=283.15,to=313.15,lty=1,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.a26.v[8],f.a26.v[9],200,x),from=283.15,to=313.15,lty=2,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.g26.v[9],f.g26.v[10],200,x),from=283.15,to=313.15,lty=3,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.r26.v[8],f.r26.v[9],200,x),from=283.15,to=313.15,lty=4,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.m31.v[8],f.m31.v[9],200,x),from=283.15,to=313.15,lty=1,col="orangered3",add=T,lwd=3)
curve(peak.ar(1,f.a31.v[8],f.a31.v[9],200,x),from=283.15,to=313.15,lty=2,col="orangered3",add=T,lwd=3)
curve(peak.ar(1,f.g31.v[8],f.g31.v[9],200,x),from=283.15,to=313.15,lty=3,col="orangered3",add=T,lwd=3)
curve(peak.ar(1,f.r31.v[8],f.r31.v[9],200,x),from=283.15,to=313.15,lty=4,col="orangered3",add=T,lwd=3)
legend(273.15+8,3.5,c(expression(underline(bolditalic("Morella"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                      expression(underline(bolditalic("Alnus"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                      expression(underline(bolditalic("Gliricidia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                      expression(underline(bolditalic("Robinia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C')),
       col=c(NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3"),
       lty=c(NA,1,1,1,NA,2,2,2,NA,3,3,3,NA,4,4,4),bty="n",lwd=3,y.intersp = 0.6,cex=1,seg.len=2,x.intersp = 0.5)
mtext(expression(italic('V')[cmax]),side=2,line=4.5,cex=1.5)
mtext(expression('(normalized to 1 at 25 '*degree*'C)'),side=2,line=3,cex=1.5)
mtext(text="a",side=3,cex=1.4,adj=0)

#S7b
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(283.15,313.15),ylim=c(0,3.3),cex.lab=1.5,cex.axis=1.2,xaxt="n")
axis(1,at=c(283.15,288.15,293.15,298.15,303.15,308.15,313.15),labels=c("10","15","20","25","30","35","40"),cex.axis=1.2)
curve(peak.ar(1,f.m21.j[8],f.m21.j[9],200,x),from=283.15,to=313.15,lty=1,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.a21.j[8],f.a21.j[9],200,x),from=283.15,to=313.15,lty=2,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.g21.j[8],f.g21.j[9],200,x),from=283.15,to=313.15,lty=3,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.r21.j[10],f.r21.j[11],200,x),from=283.15,to=313.15,lty=4,col="dodgerblue1",add=T,lwd=3)
curve(peak.ar(1,f.m26.j[8],f.m26.j[9],200,x),from=283.15,to=313.15,lty=1,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.a26.j[8],f.a26.j[9],200,x),from=283.15,to=313.15,lty=2,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.g26.j[9],f.g26.j[10],200,x),from=283.15,to=313.15,lty=3,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.r26.j[8],f.r26.j[9],200,x),from=283.15,to=313.15,lty=4,col="gold1",add=T,lwd=3)
curve(peak.ar(1,f.m31.j[8],f.m31.j[9],200,x),from=283.15,to=313.15,lty=1,col="orangered3",add=T,lwd=3)
curve(peak.ar(1,f.a31.j[8],f.a31.j[9],200,x),from=283.15,to=313.15,lty=2,col="orangered3",add=T,lwd=3)
curve(peak.ar(1,f.g31.j[8],f.g31.j[9],200,x),from=283.15,to=313.15,lty=3,col="orangered3",add=T,lwd=3)
curve(peak.ar(1,f.r31.j[8],f.r31.j[9],200,x),from=283.15,to=313.15,lty=4,col="orangered3",add=T,lwd=3)
mtext(expression(italic('J')[max]),side=2,line=4.5,cex=1.5)
mtext(expression('(normalized to 1 at 25 '*degree*'C)'),side=2,line=3,cex=1.5)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=F)
mtext(text="b",side=3,cex=1.4,adj=0)

###############################################################################################################
#Supplementary Figure 8
###############################################################################################################

#PDF dimension is 9x7 inches  

#Plotting Settings
par(pty="s")
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

#S8a
plot(M21.dat$Temp[1:3],M21.dat$Vcmax[1:3]/f.m21.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m21.v,rev(high.m21.v)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(M21.dat$Temp[1:3],M21.dat$Vcmax[1:3]/f.m21.v[2],pch=1,cex=1.5,col="black")
points(M21.dat$Temp[4:8],M21.dat$Vcmax[4:8]/f.m21.v[3],pch=2,cex=1.5,col="black")
points(M21.dat$Temp[9:11],M21.dat$Vcmax[9:11]/f.m21.v[4],pch=3,cex=1.5,col="black")
points(M21.dat$Temp[12:16],M21.dat$Vcmax[12:16]/f.m21.v[5],pch=4,cex=1.5,col="black")
points(M21.dat$Temp[17:19],M21.dat$Vcmax[17:19]/f.m21.v[6],pch=5,cex=1.5,col="black")
points(M21.dat$Temp[20:24],M21.dat$Vcmax[20:24]/f.m21.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.m21.v[8],f.m21.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)

#S8b
plot(A21.dat$Temp[1:3],A21.dat$Vcmax[1:3]/f.a21.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a21.v,rev(high.a21.v)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(A21.dat$Temp[1:3],A21.dat$Vcmax[1:3]/f.a21.v[2],pch=1,cex=1.5,col="black")
points(A21.dat$Temp[4:8],A21.dat$Vcmax[4:8]/f.a21.v[3],pch=2,cex=1.5,col="black")
points(A21.dat$Temp[9:11],A21.dat$Vcmax[9:11]/f.a21.v[4],pch=3,cex=1.5,col="black")
points(A21.dat$Temp[12:16],A21.dat$Vcmax[12:16]/f.a21.v[5],pch=4,cex=1.5,col="black")
points(A21.dat$Temp[17:19],A21.dat$Vcmax[17:19]/f.a21.v[6],pch=5,cex=1.5,col="black")
points(A21.dat$Temp[20:24],A21.dat$Vcmax[20:24]/f.a21.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.a21.v[8],f.a21.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)

#S8c
plot(G21.dat$Temp[1:3],G21.dat$Vcmax[1:3]/f.g21.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g21.v,rev(high.g21.v)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(G21.dat$Temp[1:3],G21.dat$Vcmax[1:3]/f.g21.v[2],pch=1,cex=1.5,col="black")
points(G21.dat$Temp[4:8],G21.dat$Vcmax[4:8]/f.g21.v[3],pch=2,cex=1.5,col="black")
points(G21.dat$Temp[9:11],G21.dat$Vcmax[9:11]/f.g21.v[4],pch=3,cex=1.5,col="black")
points(G21.dat$Temp[12:16],G21.dat$Vcmax[12:16]/f.g21.v[5],pch=4,cex=1.5,col="black")
points(G21.dat$Temp[17:19],G21.dat$Vcmax[17:19]/f.g21.v[6],pch=5,cex=1.5,col="black")
points(G21.dat$Temp[20:24],G21.dat$Vcmax[20:24]/f.g21.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.g21.v[8],f.g21.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)

#S8d
plot(R21.dat$Temp[1:2],R21.dat$Vcmax[1:2]/f.r21.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r21.v,rev(high.r21.v)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(R21.dat$Temp[1:2],R21.dat$Vcmax[1:2]/f.r21.v[2],pch=1,cex=1.5,col="black")
points(R21.dat$Temp[3:6],R21.dat$Vcmax[3:6]/f.r21.v[3],pch=2,cex=1.5,col="black")
points(R21.dat$Temp[7:9],R21.dat$Vcmax[7:9]/f.r21.v[4],pch=3,cex=1.5,col="black")
points(R21.dat$Temp[10:14],R21.dat$Vcmax[10:14]/f.r21.v[5],pch=4,cex=1.5,col="black")
points(R21.dat$Temp[15:17],R21.dat$Vcmax[15:17]/f.r21.v[6],pch=5,cex=1.5,col="black")
points(R21.dat$Temp[18:22],R21.dat$Vcmax[18:22]/f.r21.v[7],pch=6,cex=1.5,col="black")
points(R21.dat$Temp[23:25],R21.dat$Vcmax[23:25]/f.r21.v[8],pch=7,cex=1.5,col="black")
points(R21.dat$Temp[26:30],R21.dat$Vcmax[26:30]/f.r21.v[9],pch=8,cex=1.5,col="black")
curve(peak.ar(1,f.r21.v[10],f.r21.v[11],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)

#S8e
plot(M26.dat$Temp[1:4],M26.dat$Vcmax[1:4]/f.m26.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m26.v,rev(high.m26.v)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(M26.dat$Temp[1:4],M26.dat$Vcmax[1:4]/f.m26.v[2],pch=1,cex=1.5,col="black")
points(M26.dat$Temp[5:8],M26.dat$Vcmax[5:8]/f.m26.v[3],pch=2,cex=1.5,col="black")
points(M26.dat$Temp[9:12],M26.dat$Vcmax[9:12]/f.m26.v[4],pch=3,cex=1.5,col="black")
points(M26.dat$Temp[13:15],M26.dat$Vcmax[13:15]/f.m26.v[5],pch=4,cex=1.5,col="black")
points(M26.dat$Temp[16:19],M26.dat$Vcmax[16:19]/f.m26.v[6],pch=5,cex=1.5,col="black")
points(M26.dat$Temp[20:23],M26.dat$Vcmax[20:23]/f.m26.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.m26.v[8],f.m26.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S8f
plot(A26.dat$Temp[1:4],A26.dat$Vcmax[1:4]/f.a26.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a26.v,rev(high.a26.v)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(A26.dat$Temp[1:4],A26.dat$Vcmax[1:4]/f.a26.v[2],pch=1,cex=1.5,col="black")
points(A26.dat$Temp[5:8],A26.dat$Vcmax[5:8]/f.a26.v[3],pch=2,cex=1.5,col="black")
points(A26.dat$Temp[9:12],A26.dat$Vcmax[9:12]/f.a26.v[4],pch=3,cex=1.5,col="black")
points(A26.dat$Temp[13:15],A26.dat$Vcmax[13:15]/f.a26.v[5],pch=4,cex=1.5,col="black")
points(A26.dat$Temp[16:19],A26.dat$Vcmax[16:19]/f.a26.v[6],pch=5,cex=1.5,col="black")
points(A26.dat$Temp[20:22],A26.dat$Vcmax[20:22]/f.a26.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.a26.v[8],f.a26.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S8g
plot(G26.dat$Temp[1:4],G26.dat$Vcmax[1:4]/f.g26.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g26.v,rev(high.g26.v)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(G26.dat$Temp[1:4],G26.dat$Vcmax[1:4]/f.g26.v[2],pch=1,cex=1.5,col="black")
points(G26.dat$Temp[5:8],G26.dat$Vcmax[5:8]/f.g26.v[3],pch=2,cex=1.5,col="black")
points(G26.dat$Temp[9:12],G26.dat$Vcmax[9:12]/f.g26.v[4],pch=3,cex=1.5,col="black")
points(G26.dat$Temp[13:16],G26.dat$Vcmax[13:16]/f.g26.v[5],pch=4,cex=1.5,col="black")
points(G26.dat$Temp[17:20],G26.dat$Vcmax[17:20]/f.g26.v[6],pch=5,cex=1.5,col="black")
points(G26.dat$Temp[21:24],G26.dat$Vcmax[21:24]/f.g26.v[7],pch=6,cex=1.5,col="black")
points(G26.dat$Temp[25:27],G26.dat$Vcmax[25:27]/f.g26.v[8],pch=7,cex=1.5,col="black")
curve(peak.ar(1,f.g26.v[9],f.g26.v[10],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S8h
plot(R26.dat$Temp[1:4],R26.dat$Vcmax[1:4]/f.r26.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r26.v,rev(high.r26.v)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(R26.dat$Temp[1:4],R26.dat$Vcmax[1:4]/f.r26.v[2],pch=1,cex=1.5,col="black")
points(R26.dat$Temp[5:8],R26.dat$Vcmax[5:8]/f.r26.v[3],pch=2,cex=1.5,col="black")
points(R26.dat$Temp[9:13],R26.dat$Vcmax[9:13]/f.r26.v[4],pch=3,cex=1.5,col="black")
points(R26.dat$Temp[14:17],R26.dat$Vcmax[14:17]/f.r26.v[5],pch=4,cex=1.5,col="black")
points(R26.dat$Temp[18:21],R26.dat$Vcmax[18:21]/f.r26.v[6],pch=5,cex=1.5,col="black")
points(R26.dat$Temp[22:25],R26.dat$Vcmax[22:25]/f.r26.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.r26.v[8],f.r26.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

#S8i
plot(M31.dat$Temp[1:5],M31.dat$Vcmax[1:5]/f.m31.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m31.v,rev(high.m31.v)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(M31.dat$Temp[1:5],M31.dat$Vcmax[1:5]/f.m31.v[2],pch=1,cex=1.5,col="black")
points(M31.dat$Temp[6:8],M31.dat$Vcmax[6:8]/f.m31.v[3],pch=2,cex=1.5,col="black")
points(M31.dat$Temp[9:12],M31.dat$Vcmax[9:12]/f.m31.v[4],pch=3,cex=1.5,col="black")
points(M31.dat$Temp[13:15],M31.dat$Vcmax[13:15]/f.m31.v[5],pch=4,cex=1.5,col="black")
points(M31.dat$Temp[16:20],M31.dat$Vcmax[16:20]/f.m31.v[6],pch=5,cex=1.5,col="black")
points(M31.dat$Temp[21:23],M31.dat$Vcmax[21:23]/f.m31.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.m31.v[8],f.m31.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S8j
plot(A31.dat$Temp[1:5],A31.dat$Vcmax[1:5]/f.a31.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a31.v,rev(high.a31.v)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(A31.dat$Temp[1:5],A31.dat$Vcmax[1:5]/f.a31.v[2],pch=1,cex=1.5,col="black")
points(A31.dat$Temp[6:7],A31.dat$Vcmax[6:7]/f.a31.v[3],pch=2,cex=1.5,col="black")
points(A31.dat$Temp[8:12],A31.dat$Vcmax[8:12]/f.a31.v[4],pch=3,cex=1.5,col="black")
points(A31.dat$Temp[13:15],A31.dat$Vcmax[13:15]/f.a31.v[5],pch=4,cex=1.5,col="black")
points(A31.dat$Temp[16:20],A31.dat$Vcmax[16:20]/f.a31.v[6],pch=5,cex=1.5,col="black")
points(A31.dat$Temp[21:23],A31.dat$Vcmax[21:23]/f.a31.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.a31.v[8],f.a31.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S8k
plot(G31.dat$Temp[1:4],G31.dat$Vcmax[1:4]/f.g31.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g31.v,rev(high.g31.v)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(G31.dat$Temp[1:4],G31.dat$Vcmax[1:4]/f.g31.v[2],pch=1,cex=1.5,col="black")
points(G31.dat$Temp[5:7],G31.dat$Vcmax[5:7]/f.g31.v[3],pch=2,cex=1.5,col="black")
points(G31.dat$Temp[8:12],G31.dat$Vcmax[8:12]/f.g31.v[4],pch=3,cex=1.5,col="black")
points(G31.dat$Temp[13:14],G31.dat$Vcmax[13:14]/f.g31.v[5],pch=4,cex=1.5,col="black")
points(G31.dat$Temp[15:20],G31.dat$Vcmax[15:20]/f.g31.v[6],pch=5,cex=1.5,col="black")
points(G31.dat$Temp[21:23],G31.dat$Vcmax[21:23]/f.g31.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.g31.v[8],f.g31.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S8l
plot(R31.dat$Temp[1:5],R31.dat$Vcmax[1:5]/f.r31.v[2],pch=1,cex=1.5,col="black",ylim=c(0,3),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r31.v,rev(high.r31.v)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(R31.dat$Temp[1:5],R31.dat$Vcmax[1:5]/f.r31.v[2],pch=1,cex=1.5,col="black")
points(R31.dat$Temp[6:8],R31.dat$Vcmax[6:8]/f.r31.v[3],pch=2,cex=1.5,col="black")
points(R31.dat$Temp[9:13],R31.dat$Vcmax[9:13]/f.r31.v[4],pch=3,cex=1.5,col="black")
points(R31.dat$Temp[14:16],R31.dat$Vcmax[14:16]/f.r31.v[5],pch=4,cex=1.5,col="black")
points(R31.dat$Temp[17:21],R31.dat$Vcmax[17:21]/f.r31.v[6],pch=5,cex=1.5,col="black")
points(R31.dat$Temp[22:24],R31.dat$Vcmax[22:24]/f.r31.v[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.r31.v[8],f.r31.v[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression(italic('V')[cmax]*' (normalized to 1 at 25 '*degree*'C)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)

###############################################################################################################
#Supplementary Figure 9
###############################################################################################################

#PDF dimension is 9x7 inches  

#Plotting Settings
par(pty="s")
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

#S9a
plot(M21.dat$Temp[1:3],M21.dat$Jmax[1:3]/f.m21.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m21.j,rev(high.m21.j)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(M21.dat$Temp[1:3],M21.dat$Jmax[1:3]/f.m21.j[2],pch=1,cex=1.5,col="black")
points(M21.dat$Temp[4:8],M21.dat$Jmax[4:8]/f.m21.j[3],pch=2,cex=1.5,col="black")
points(M21.dat$Temp[9:11],M21.dat$Jmax[9:11]/f.m21.j[4],pch=3,cex=1.5,col="black")
points(M21.dat$Temp[12:16],M21.dat$Jmax[12:16]/f.m21.j[5],pch=4,cex=1.5,col="black")
points(M21.dat$Temp[17:19],M21.dat$Jmax[17:19]/f.m21.j[6],pch=5,cex=1.5,col="black")
points(M21.dat$Temp[20:24],M21.dat$Jmax[20:24]/f.m21.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.m21.j[8],f.m21.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)

#S9b
plot(A21.dat$Temp[1:3],A21.dat$Jmax[1:3]/f.a21.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a21.j,rev(high.a21.j)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(A21.dat$Temp[1:3],A21.dat$Jmax[1:3]/f.a21.j[2],pch=1,cex=1.5,col="black")
points(A21.dat$Temp[4:8],A21.dat$Jmax[4:8]/f.a21.j[3],pch=2,cex=1.5,col="black")
points(A21.dat$Temp[9:11],A21.dat$Jmax[9:11]/f.a21.j[4],pch=3,cex=1.5,col="black")
points(A21.dat$Temp[12:16],A21.dat$Jmax[12:16]/f.a21.j[5],pch=4,cex=1.5,col="black")
points(A21.dat$Temp[17:19],A21.dat$Jmax[17:19]/f.a21.j[6],pch=5,cex=1.5,col="black")
points(A21.dat$Temp[20:24],A21.dat$Jmax[20:24]/f.a21.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.a21.j[8],f.a21.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)

#S9c
plot(G21.dat$Temp[1:3],G21.dat$Jmax[1:3]/f.g21.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g21.j,rev(high.g21.j)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(G21.dat$Temp[1:3],G21.dat$Jmax[1:3]/f.g21.j[2],pch=1,cex=1.5,col="black")
points(G21.dat$Temp[4:8],G21.dat$Jmax[4:8]/f.g21.j[3],pch=2,cex=1.5,col="black")
points(G21.dat$Temp[9:11],G21.dat$Jmax[9:11]/f.g21.j[4],pch=3,cex=1.5,col="black")
points(G21.dat$Temp[12:16],G21.dat$Jmax[12:16]/f.g21.j[5],pch=4,cex=1.5,col="black")
points(G21.dat$Temp[17:19],G21.dat$Jmax[17:19]/f.g21.j[6],pch=5,cex=1.5,col="black")
points(G21.dat$Temp[20:24],G21.dat$Jmax[20:24]/f.g21.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.g21.j[8],f.g21.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)

#S9d
plot(R21.dat$Temp[1:2],R21.dat$Jmax[1:2]/f.r21.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r21.j,rev(high.r21.j)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(R21.dat$Temp[1:2],R21.dat$Jmax[1:2]/f.r21.j[2],pch=1,cex=1.5,col="black")
points(R21.dat$Temp[3:6],R21.dat$Jmax[3:6]/f.r21.j[3],pch=2,cex=1.5,col="black")
points(R21.dat$Temp[7:9],R21.dat$Jmax[7:9]/f.r21.j[4],pch=3,cex=1.5,col="black")
points(R21.dat$Temp[10:14],R21.dat$Jmax[10:14]/f.r21.j[5],pch=4,cex=1.5,col="black")
points(R21.dat$Temp[15:17],R21.dat$Jmax[15:17]/f.r21.j[6],pch=5,cex=1.5,col="black")
points(R21.dat$Temp[18:22],R21.dat$Jmax[18:22]/f.r21.j[7],pch=6,cex=1.5,col="black")
points(R21.dat$Temp[23:25],R21.dat$Jmax[23:25]/f.r21.j[8],pch=7,cex=1.5,col="black")
points(R21.dat$Temp[26:30],R21.dat$Jmax[26:30]/f.r21.j[9],pch=8,cex=1.5,col="black")
curve(peak.ar(1,f.r21.j[10],f.r21.j[11],200,x+273.15),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)

#S9e
plot(M26.dat$Temp[1:4],M26.dat$Jmax[1:4]/f.m26.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m26.j,rev(high.m26.j)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(M26.dat$Temp[1:4],M26.dat$Jmax[1:4]/f.m26.j[2],pch=1,cex=1.5,col="black")
points(M26.dat$Temp[5:8],M26.dat$Jmax[5:8]/f.m26.j[3],pch=2,cex=1.5,col="black")
points(M26.dat$Temp[9:12],M26.dat$Jmax[9:12]/f.m26.j[4],pch=3,cex=1.5,col="black")
points(M26.dat$Temp[13:15],M26.dat$Jmax[13:15]/f.m26.j[5],pch=4,cex=1.5,col="black")
points(M26.dat$Temp[16:19],M26.dat$Jmax[16:19]/f.m26.j[6],pch=5,cex=1.5,col="black")
points(M26.dat$Temp[20:23],M26.dat$Jmax[20:23]/f.m26.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.m26.j[8],f.m26.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S9f
plot(A26.dat$Temp[1:4],A26.dat$Jmax[1:4]/f.a26.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a26.j,rev(high.a26.j)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(A26.dat$Temp[1:4],A26.dat$Jmax[1:4]/f.a26.j[2],pch=1,cex=1.5,col="black")
points(A26.dat$Temp[5:8],A26.dat$Jmax[5:8]/f.a26.j[3],pch=2,cex=1.5,col="black")
points(A26.dat$Temp[9:12],A26.dat$Jmax[9:12]/f.a26.j[4],pch=3,cex=1.5,col="black")
points(A26.dat$Temp[13:15],A26.dat$Jmax[13:15]/f.a26.j[5],pch=4,cex=1.5,col="black")
points(A26.dat$Temp[16:19],A26.dat$Jmax[16:19]/f.a26.j[6],pch=5,cex=1.5,col="black")
points(A26.dat$Temp[20:22],A26.dat$Jmax[20:22]/f.a26.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.a26.j[8],f.a26.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S9g
plot(G26.dat$Temp[1:4],G26.dat$Jmax[1:4]/f.g26.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g26.j,rev(high.g26.j)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(G26.dat$Temp[1:4],G26.dat$Jmax[1:4]/f.g26.j[2],pch=1,cex=1.5,col="black")
points(G26.dat$Temp[5:8],G26.dat$Jmax[5:8]/f.g26.j[3],pch=2,cex=1.5,col="black")
points(G26.dat$Temp[9:12],G26.dat$Jmax[9:12]/f.g26.j[4],pch=3,cex=1.5,col="black")
points(G26.dat$Temp[13:16],G26.dat$Jmax[13:16]/f.g26.j[5],pch=4,cex=1.5,col="black")
points(G26.dat$Temp[17:20],G26.dat$Jmax[17:20]/f.g26.j[6],pch=5,cex=1.5,col="black")
points(G26.dat$Temp[21:24],G26.dat$Jmax[21:24]/f.g26.j[7],pch=6,cex=1.5,col="black")
points(G26.dat$Temp[25:27],G26.dat$Jmax[25:27]/f.g26.j[8],pch=7,cex=1.5,col="black")
curve(peak.ar(1,f.g26.j[9],f.g26.j[10],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S9h
plot(R26.dat$Temp[1:4],R26.dat$Jmax[1:4]/f.r26.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r26.j,rev(high.r26.j)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(R26.dat$Temp[1:4],R26.dat$Jmax[1:4]/f.r26.j[2],pch=1,cex=1.5,col="black")
points(R26.dat$Temp[5:8],R26.dat$Jmax[5:8]/f.r26.j[3],pch=2,cex=1.5,col="black")
points(R26.dat$Temp[9:13],R26.dat$Jmax[9:13]/f.r26.j[4],pch=3,cex=1.5,col="black")
points(R26.dat$Temp[14:17],R26.dat$Jmax[14:17]/f.r26.j[5],pch=4,cex=1.5,col="black")
points(R26.dat$Temp[18:21],R26.dat$Jmax[18:21]/f.r26.j[6],pch=5,cex=1.5,col="black")
points(R26.dat$Temp[22:25],R26.dat$Jmax[22:25]/f.r26.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.r26.j[8],f.r26.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

#S9i
plot(M31.dat$Temp[1:5],M31.dat$Jmax[1:5]/f.m31.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",cex.axis=1.5)
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m31.j,rev(high.m31.j)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(M31.dat$Temp[1:5],M31.dat$Jmax[1:5]/f.m31.j[2],pch=1,cex=1.5,col="black")
points(M31.dat$Temp[6:8],M31.dat$Jmax[6:8]/f.m31.j[3],pch=2,cex=1.5,col="black")
points(M31.dat$Temp[9:12],M31.dat$Jmax[9:12]/f.m31.j[4],pch=3,cex=1.5,col="black")
points(M31.dat$Temp[13:15],M31.dat$Jmax[13:15]/f.m31.j[5],pch=4,cex=1.5,col="black")
points(M31.dat$Temp[16:20],M31.dat$Jmax[16:20]/f.m31.j[6],pch=5,cex=1.5,col="black")
points(M31.dat$Temp[21:23],M31.dat$Jmax[21:23]/f.m31.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.m31.j[8],f.m31.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S9j
plot(A31.dat$Temp[1:5],A31.dat$Jmax[1:5]/f.a31.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a31.j,rev(high.a31.j)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(A31.dat$Temp[1:5],A31.dat$Jmax[1:5]/f.a31.j[2],pch=1,cex=1.5,col="black")
points(A31.dat$Temp[6:7],A31.dat$Jmax[6:7]/f.a31.j[3],pch=2,cex=1.5,col="black")
points(A31.dat$Temp[8:12],A31.dat$Jmax[8:12]/f.a31.j[4],pch=3,cex=1.5,col="black")
points(A31.dat$Temp[13:15],A31.dat$Jmax[13:15]/f.a31.j[5],pch=4,cex=1.5,col="black")
points(A31.dat$Temp[16:20],A31.dat$Jmax[16:20]/f.a31.j[6],pch=5,cex=1.5,col="black")
points(A31.dat$Temp[21:23],A31.dat$Jmax[21:23]/f.a31.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.a31.j[8],f.a31.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S9k
plot(G31.dat$Temp[1:4],G31.dat$Jmax[1:4]/f.g31.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g31.j,rev(high.g31.j)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(G31.dat$Temp[1:4],G31.dat$Jmax[1:4]/f.g31.j[2],pch=1,cex=1.5,col="black")
points(G31.dat$Temp[5:7],G31.dat$Jmax[5:7]/f.g31.j[3],pch=2,cex=1.5,col="black")
points(G31.dat$Temp[8:12],G31.dat$Jmax[8:12]/f.g31.j[4],pch=3,cex=1.5,col="black")
points(G31.dat$Temp[13:14],G31.dat$Jmax[13:14]/f.g31.j[5],pch=4,cex=1.5,col="black")
points(G31.dat$Temp[15:20],G31.dat$Jmax[15:20]/f.g31.j[6],pch=5,cex=1.5,col="black")
points(G31.dat$Temp[21:23],G31.dat$Jmax[21:23]/f.g31.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.g31.j[8],f.g31.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S9l
plot(R31.dat$Temp[1:5],R31.dat$Jmax[1:5]/f.r31.j[2],pch=1,cex=1.5,col="black",ylim=c(0,2.5),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.5,1,1.5,2,2.5,3),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r31.j,rev(high.r31.j)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(R31.dat$Temp[1:5],R31.dat$Jmax[1:5]/f.r31.j[2],pch=1,cex=1.5,col="black")
points(R31.dat$Temp[6:8],R31.dat$Jmax[6:8]/f.r31.j[3],pch=2,cex=1.5,col="black")
points(R31.dat$Temp[9:13],R31.dat$Jmax[9:13]/f.r31.j[4],pch=3,cex=1.5,col="black")
points(R31.dat$Temp[14:16],R31.dat$Jmax[14:16]/f.r31.j[5],pch=4,cex=1.5,col="black")
points(R31.dat$Temp[17:21],R31.dat$Jmax[17:21]/f.r31.j[6],pch=5,cex=1.5,col="black")
points(R31.dat$Temp[22:24],R31.dat$Jmax[22:24]/f.r31.j[7],pch=6,cex=1.5,col="black")
curve(peak.ar(1,f.r31.j[8],f.r31.j[9],200,x+273.15),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression(italic('J')[max]*' (normalized to 1 at 25 '*degree*'C)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)

###############################################################################################################
#Supplementary Figure 10
###############################################################################################################

#PDF dimension is 8x6 inches  

#Add random horizontal noise to predicted values
Tg.jitter<-jitter(rep(c(18.5,23.5,28.5),4))

#Plotting Settings
par(pty="s")
nf<-layout(matrix(c(1,2,5,3,4,5),2,3,byrow=T),c(3,3,1,3,3,1),c(3,3,6,3,3,6),T)
layout.show(nf)
par(mar=c(4,6,2,1))
par(oma=c(0,0,0,0))

#S10a
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(30,45),cex.lab=1.5,cex.axis=1.2)
mtext(text="a",side=3,cex=1.2,adj=0)
mtext(expression(italic('T')['opt,'*italic('V')[cmax]]*' ('*degree*'C)'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
abline(v=18.5,col="gray")
abline(v=23.5,col="gray")
abline(v=28.5,col="gray")
points(Topt.m.v~Tg.jitter[1:3],pch=16,col="darkorange1",cex=1.5,lwd=1)
points(Topt.a.v~Tg.jitter[4:6],pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(Topt.g.v~Tg.jitter[7:9],pch=17,col="orangered2",cex=1.5,lwd=1)
points(Topt.r.v~Tg.jitter[10:12],pch=17,col="dodgerblue3",cex=1.5,lwd=1)
arrows(Tg.jitter[1:3], l.Topt.m.v, Tg.jitter[1:3],h.Topt.m.v, length=0.05, angle=0, code=3,lwd=1.2,col="darkorange1")
arrows(Tg.jitter[4:6], l.Topt.a.v, Tg.jitter[4:6],h.Topt.a.v, length=0.05, angle=0, code=3,lwd=1.2,col="darkturquoise")
arrows(Tg.jitter[7:9], l.Topt.g.v, Tg.jitter[7:9],h.Topt.g.v, length=0.05, angle=0, code=3,lwd=1.2,col="orangered2")
arrows(Tg.jitter[10:12], l.Topt.r.v, Tg.jitter[10:12],h.Topt.r.v, length=0.05, angle=0, code=3,lwd=1.2,col="dodgerblue3")

#S10b
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(30,45),cex.lab=1.5,cex.axis=1.2)
mtext(text="b",side=3,cex=1.2,adj=0)
mtext(expression(italic('T')['opt,'*italic('J')[max]]*' ('*degree*'C)'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
abline(v=18.5,col="gray")
abline(v=23.5,col="gray")
abline(v=28.5,col="gray")
points(Topt.m.j~Tg.jitter[1:3],pch=16,col="darkorange1",cex=1.5,lwd=1)
points(Topt.a.j~Tg.jitter[4:6],pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(Topt.g.j~Tg.jitter[7:9],pch=17,col="orangered2",cex=1.5,lwd=1)
points(Topt.r.j~Tg.jitter[10:12],pch=17,col="dodgerblue3",cex=1.5,lwd=1)
arrows(Tg.jitter[1:3], l.Topt.m.j, Tg.jitter[1:3],h.Topt.m.j, length=0.05, angle=0, code=3,lwd=1.2,col="darkorange1")
arrows(Tg.jitter[4:6], l.Topt.a.j, Tg.jitter[4:6],h.Topt.a.j, length=0.05, angle=0, code=3,lwd=1.2,col="darkturquoise")
arrows(Tg.jitter[7:9], l.Topt.g.j, Tg.jitter[7:9],h.Topt.g.j, length=0.05, angle=0, code=3,lwd=1.2,col="orangered2")
arrows(Tg.jitter[10:12], l.Topt.r.j, Tg.jitter[10:12],h.Topt.r.j, length=0.05, angle=0, code=3,lwd=1.2,col="dodgerblue3")

#S10c
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(0,200),cex.lab=1.5,cex.axis=1.2)
mtext(text="c",side=3,cex=1.2,adj=0)
mtext(expression(italic('E')['a,'*italic('V')[cmax]]*' (kJ mol'^-1*')'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
abline(v=18.5,col="gray")
abline(v=23.5,col="gray")
abline(v=28.5,col="gray")
points(Ea.m.v~Tg.jitter[1:3],pch=16,col="darkorange1",cex=1.5,lwd=1)
points(Ea.a.v~Tg.jitter[4:6],pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(Ea.g.v~Tg.jitter[7:9],pch=17,col="orangered2",cex=1.5,lwd=1)
points(Ea.r.v~Tg.jitter[10:12],pch=17,col="dodgerblue3",cex=1.5,lwd=1)
arrows(Tg.jitter[1:3], l.Ea.m.v, Tg.jitter[1:3],h.Ea.m.v, length=0.05, angle=0, code=3,lwd=1.2,col="darkorange1")
arrows(Tg.jitter[4:6], l.Ea.a.v, Tg.jitter[4:6],h.Ea.a.v, length=0.05, angle=0, code=3,lwd=1.2,col="darkturquoise")
arrows(Tg.jitter[7:9], l.Ea.g.v, Tg.jitter[7:9],h.Ea.g.v, length=0.05, angle=0, code=3,lwd=1.2,col="orangered2")
arrows(Tg.jitter[10:12], l.Ea.r.v, Tg.jitter[10:12],h.Ea.r.v, length=0.05, angle=0, code=3,lwd=1.2,col="dodgerblue3")

#S10d
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(0,200),cex.lab=1.5,cex.axis=1.2)
mtext(text="d",side=3,cex=1.2,adj=0)
mtext(expression(italic('E')['a,'*italic('J')[cmax]]*' (kJ mol'^-1*')'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
abline(v=18.5,col="gray")
abline(v=23.5,col="gray")
abline(v=28.5,col="gray")
points(Ea.m.j~Tg.jitter[1:3],pch=16,col="darkorange1",cex=1.5,lwd=1)
points(Ea.a.j~Tg.jitter[4:6],pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(Ea.g.j~Tg.jitter[7:9],pch=17,col="orangered2",cex=1.5,lwd=1)
points(Ea.r.j~Tg.jitter[10:12],pch=17,col="dodgerblue3",cex=1.5,lwd=1)
arrows(Tg.jitter[1:3], l.Ea.m.j, Tg.jitter[1:3],h.Ea.m.j, length=0.05, angle=0, code=3,lwd=1.2,col="darkorange1")
arrows(Tg.jitter[4:6], l.Ea.a.j, Tg.jitter[4:6],h.Ea.a.j, length=0.05, angle=0, code=3,lwd=1.2,col="darkturquoise")
arrows(Tg.jitter[7:9], l.Ea.g.j, Tg.jitter[7:9],h.Ea.g.j, length=0.05, angle=0, code=3,lwd=1.2,col="orangered2")
arrows(Tg.jitter[10:12], l.Ea.r.j, Tg.jitter[10:12],h.Ea.r.j, length=0.05, angle=0, code=3,lwd=1.2,col="dodgerblue3")

#Change plotting Settings for legend
par(mar=c(0,0,0,0))
par(xpd=TRUE)

#Legend
plot(0:10, 0:10, type='n', bty='n', xaxt='n', yaxt='n',xlab=NA,ylab=NA)
legend("left",legend=c("Morella","Alnus","Gliricidia","Robinia"),text.font=c(3,3,3,3),
       col=c("darkorange1","darkturquoise","orangered2","dodgerblue3"),pch=c(16,16,17,17),bty="n",pt.cex=1.5,cex=1.2)
par(xpd=FALSE)