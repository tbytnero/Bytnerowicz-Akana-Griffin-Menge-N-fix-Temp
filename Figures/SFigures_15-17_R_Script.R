###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figures 15-17
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
#Read data for Alnus in from "SNF_Temp" folder (there was no linear acclimation of Km for Alnus)
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

####
#Read data for other species in from "SNF_Temp_linKm" folder
####

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

####
#Read photosynthesis data (A275) in from "Photo_Temp" folder 
####

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

####
#Define functions
####

#No acclimation (used for calculating N-fixation rates at 25 deg. C [Supplementary Figure 16])
beta <- function(ymax,Tmin,Topt,Tmax,T){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of all parameters (used for N-fixation with Morella and Alnus)
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

#Acclimation of Tmin and Tmax (used for N-fixation with Robinia)
beta.lin.Tmin.Tmax <- function(ymax,Topt,a,b,c,d,T,Tgrow){
  y <- pmax(0,ymax*((c+d*Tgrow)-T)/((c+d*Tgrow)-Topt)*(((T-(a+b*Tgrow))/(Topt-(a+b*Tgrow)))^((Topt-(a+b*Tgrow))/((c+d*Tgrow)-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of Topt (used for photosynthesis [A275])
beta.Topt.lin <- function(ymax,Tmin,a,b,Tmax,T,Tgrow){
  y <- ymax*(Tmax-T)/(Tmax-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/(Tmax-(a+b*Tgrow))))
  y
}

#Houlton et al. 2008 function (Supplementary Figure 16)
Houlton_func <- function(Ts){
  s <- sqrt(-2*(-3.62))/0.27
  Topt <- 25.15
  F <- exp(-(Ts - Topt)^2/(2*s*s))
  F
}

#Normalized Houlton et al. 2008 function (Supplementary Figure 16)
Houlton_func_norm <- function(Ts){
  s <- sqrt(-2*(-3.62))/0.27
  Topt <- 25.15
  F <- exp(-(Ts - Topt)^2/(2*s*s))
  F <- F/max(F)
}

####
#Negative log-likelihood (NLL) functions
####

#NLL function for acclimation of all parameters (used for N-fixation with Morella and Alnus)
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

#NLL function for acclimation of Tmin and Tmax (used for N-fixation with Robinia)
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

#NLL function for acclimation of Topt (used for photosynthesis)
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

####
#Maximum likelihood fits
####

###
#N-fixation
###

#Morella
fit_Nase_beta_linall_MOCE_linKm <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
                                                                            ymax26a=1,ymax26b=1,ymax26c=1,
                                                                            ymax31a=1,ymax31b=1,ymax31c=1,
                                                                            a=-20,b=1,c=19,d=0.67,e=44,f=0.03),
                                        data=list(T21a=MOCE21a$Temperature,T21b=MOCE21b$Temperature,T21c=MOCE21c$Temperature,
                                                  T26a=MOCE26a$Temperature,T26b=MOCE26b$Temperature,T26c=MOCE26c$Temperature,
                                                  T31a=MOCE31a$Temperature,T31b=MOCE31b$Temperature,T31c=MOCE31c$Temperature,
                                                  Nasedat21a=MOCE21a$Vmax/9960.998949,Nasedat21b=MOCE21b$Vmax/3833.121977,Nasedat21c=MOCE21c$Vmax/3095.456985,
                                                  Nasedat26a=MOCE26a$Vmax/6492.677546,Nasedat26b=MOCE26b$Vmax/5974.076428,Nasedat26c=MOCE26c$Vmax/7890.763063,
                                                  Nasedat31a=MOCE31a$Vmax/8826.601393,Nasedat31b=MOCE31b$Vmax/2371.923752,Nasedat31c=MOCE31c$Vmax/629.7544222),
                                        control=list(maxit=20000))
summary(fit_Nase_beta_linall_MOCE_linKm)

#Alnus
fit_Nase_beta_linall_ALRU <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-2.39845357,ymax21a=0.99686505,ymax21b=1.04547194,ymax21c=0.97011417,
                                                                      ymax26a=1.00547291,ymax26b=1.00533434,ymax26c=0.91749066,
                                                                      ymax31a=1.17651704,ymax31b=1.08128675,ymax31c=1.08072962,
                                                                      a=-19.92168066,b=0.48755757,c=31.23940575,d=0.06415693,e=41.92691426,f=0.02707867),
                                  data=list(T21a=ALRU21a$Temperature,T21b=ALRU21b$Temperature,T21c=ALRU21c$Temperature,
                                            T26a=ALRU26a$Temperature,T26b=ALRU26b$Temperature,T26c=ALRU26c$Temperature,
                                            T31a=ALRU31a$Temperature,T31b=ALRU31b$Temperature,T31c=ALRU31c$Temperature,
                                            Nasedat21a=ALRU21a$Vmax/10035.24751,Nasedat21b=ALRU21b$Vmax/7990.765838,Nasedat21c=ALRU21c$Vmax/4046.082544,
                                            Nasedat26a=ALRU26a$Vmax/2236.979133,Nasedat26b=ALRU26b$Vmax/2744.252462,Nasedat26c=ALRU26c$Vmax/1831.763172,
                                            Nasedat31a=ALRU31a$Vmax/406.3524389,Nasedat31b=ALRU31b$Vmax/1832.723413,Nasedat31c=ALRU31c$Vmax/1350.729381),
                                  control=list(maxit=20000))
summary(fit_Nase_beta_linall_ALRU)

#Gliricidia
fit_Nase_beta_linTminTopt_GLSE_linKm <- mle2(Nase_beta_linTminTopt_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,
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
summary(fit_Nase_beta_linTminTopt_GLSE_linKm)

#Robinia
fit_Nase_beta_linTminTmax_ROPS_linKm <- mle2(Nase_beta_linTminTmax_normNLL,start=list(sdNase=-1,ymax21a=0.96320475,ymax21b=0.82996072,ymax21c=1.08051641,
                                                                                      ymax26a=0.89390710,ymax26b=0.95916073,ymax26c=0.88695838,
                                                                                      ymax31a=0.89351822,ymax31b=0.69375718,ymax31c=0.88715637,
                                                                                      Topt=33,a=-16.33352131,b=0.81903181,c=44.75447753,d=0.03978480),
                                             data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                                       T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                                       T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                                       Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                                       Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                                       Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                             control=list(maxit=20000))
summary(fit_Nase_beta_linTminTmax_ROPS_linKm)

###
#Photosynthesis (A275)
###

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

####
#Save N-fixation parameter estimates
####

#Morella
ft.m_linKm <- coef(fit_Nase_beta_linall_MOCE_linKm)

#Alnus
ft.a <- coef(fit_Nase_beta_linall_ALRU)

#Gliricidia
ft.g_linKm <- coef(fit_Nase_beta_linTminTopt_GLSE_linKm)

#Robinia
ft.r_linKm <- coef(fit_Nase_beta_linTminTmax_ROPS_linKm)

####
#Calculate N-fixation rates at 15 and 40 deg. C as a function of growing temperature
####

#Vector of mean growing and measurement temperatures
Tg.seq<-seq(15,35,0.1) #Growing temperature
Tsim<-seq(0,50,0.1) #Measurement temperature

###
#Rates at 15 deg. C (Figure 3a)
###

#Morella
M15_linKm<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  M15_linKm[i]<-beta(1,ft.m_linKm[12]*Tg.seq[i]+ft.m_linKm[11],ft.m_linKm[14]*Tg.seq[i]+ft.m_linKm[13],ft.m_linKm[16]*Tg.seq[i]+ft.m_linKm[15],15)
}
#At high growing temperatures Tmin is greater than 15 deg. C causing NA values to appear
#This makes those NA values equal to zero
M15_linKm[is.nan(M15_linKm)]<-0 

#Alnus
A15<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  A15[i]<-beta(1,ft.a[12]*Tg.seq[i]+ft.a[11],ft.a[14]*Tg.seq[i]+ft.a[13],ft.a[16]*Tg.seq[i]+ft.a[15],15)
}

#Gliricidia
G15_linKm<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  G15_linKm[i]<-beta(1,ft.g_linKm[13]*Tg.seq[i]+ft.g_linKm[12],ft.g_linKm[15]*Tg.seq[i]+ft.g_linKm[14],ft.g_linKm[11],15)
}

#Robinia
R15_linKm<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  R15_linKm[i]<-beta(1,ft.r_linKm[13]*Tg.seq[i]+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*Tg.seq[i]+ft.r_linKm[14],15)
}

###
#Rates at 40 deg. C (Figure 3c)
###

#Morella
M40_linKm<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  M40_linKm[i]<-beta(1,ft.m_linKm[12]*Tg.seq[i]+ft.m_linKm[11],ft.m_linKm[14]*Tg.seq[i]+ft.m_linKm[13],ft.m_linKm[16]*Tg.seq[i]+ft.m_linKm[15],40)
}

#Alnus
A40<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  A40[i]<-beta(1,ft.a[12]*Tg.seq[i]+ft.a[11],ft.a[14]*Tg.seq[i]+ft.a[13],ft.a[16]*Tg.seq[i]+ft.a[15],40)
}

#Gliricidia
G40_linKm<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  G40_linKm[i]<-beta(1,ft.g_linKm[13]*Tg.seq[i]+ft.g_linKm[12],ft.g_linKm[15]*Tg.seq[i]+ft.g_linKm[14],ft.g_linKm[11],40)
}

#Robinia
R40_linKm<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  R40_linKm[i]<-beta(1,ft.r_linKm[13]*Tg.seq[i]+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*Tg.seq[i]+ft.r_linKm[14],40)
}

####
#Calculate N-fixation 95% CI
####

##
#Morella
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m_linKm=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_MOCE_linKm),Sigma=vcov(fit_Nase_beta_linall_MOCE_linKm))

#N-fixation ~ measurement temperature for 21:15 deg. C growing temperature
dist.m21_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m21_linKm<-rep(NA,length(Tsim))
high.m21_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m21_linKm[i,j]=beta.lin.all(1,vmat.m_linKm[i,11],vmat.m_linKm[i,12],vmat.m_linKm[i,13],vmat.m_linKm[i,14],vmat.m_linKm[i,15],vmat.m_linKm[i,16],Tsim[j],18.5)
  }
  low.m21_linKm[j]<-quantile(dist.m21_linKm[,j],0.025)
  high.m21_linKm[j]<-quantile(dist.m21_linKm[,j],0.975)
}

#N-fixation ~ measurement temperature for 26:20 deg. C growing temperature
dist.m26_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m26_linKm<-rep(NA,length(Tsim))
high.m26_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m26_linKm[i,j]=beta.lin.all(1,vmat.m_linKm[i,11],vmat.m_linKm[i,12],vmat.m_linKm[i,13],vmat.m_linKm[i,14],vmat.m_linKm[i,15],vmat.m_linKm[i,16],Tsim[j],23.5)
  }
  low.m26_linKm[j]<-quantile(dist.m26_linKm[,j],0.025)
  high.m26_linKm[j]<-quantile(dist.m26_linKm[,j],0.975)
}

#N-fixation ~ measurement temperature for 31:25 deg. C growing temperature
dist.m31_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m31_linKm<-rep(NA,length(Tsim))
high.m31_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m31_linKm[i,j]=beta.lin.all(1,vmat.m_linKm[i,11],vmat.m_linKm[i,12],vmat.m_linKm[i,13],vmat.m_linKm[i,14],vmat.m_linKm[i,15],vmat.m_linKm[i,16],Tsim[j],28.5)
  }
  low.m31_linKm[j]<-quantile(dist.m31_linKm[,j],0.025)
  high.m31_linKm[j]<-quantile(dist.m31_linKm[,j],0.975)
}

#Rates at 15 deg. C 95% CI
dist.m15_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m15_linKm<-rep(NA,length(Tg.seq))
high.m15_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m15_linKm[i,j]=beta.lin.all(1,vmat.m_linKm[i,11],vmat.m_linKm[i,12],vmat.m_linKm[i,13],vmat.m_linKm[i,14],vmat.m_linKm[i,15],vmat.m_linKm[i,16],15,Tg.seq[j])
  }
  low.m15_linKm[j]<-quantile(dist.m15_linKm[,j],0.025)
  high.m15_linKm[j]<-quantile(dist.m15_linKm[,j],0.975)
}

#Topt 95% CI
dist.mTopt_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.mTopt_linKm<-rep(NA,length(Tg.seq))
high.mTopt_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.mTopt_linKm[i,j]=vmat.m_linKm[i,13]+vmat.m_linKm[i,14]*Tg.seq[j]
  }
  low.mTopt_linKm[j]<-quantile(dist.mTopt_linKm[,j],0.025)
  high.mTopt_linKm[j]<-quantile(dist.mTopt_linKm[,j],0.975)
}

#Rates at 40 deg. C 95% CI
dist.m40_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m40_linKm<-rep(NA,length(Tg.seq))
high.m40_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m40_linKm[i,j]=beta.lin.all(1,vmat.m_linKm[i,11],vmat.m_linKm[i,12],vmat.m_linKm[i,13],vmat.m_linKm[i,14],vmat.m_linKm[i,15],vmat.m_linKm[i,16],40,Tg.seq[j])
  }
  low.m40_linKm[j]<-quantile(dist.m40_linKm[,j],0.025)
  high.m40_linKm[j]<-quantile(dist.m40_linKm[,j],0.975)
}

##
#Alnus
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_ALRU),Sigma=vcov(fit_Nase_beta_linall_ALRU))

#N-fixation ~ measurement temperature for 21:15 deg. C growing temperature
dist.a21=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a21<-rep(NA,length(Tsim))
high.a21<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a21[i,j]=beta.lin.all(1,vmat.a[i,11],vmat.a[i,12],vmat.a[i,13],vmat.a[i,14],vmat.a[i,15],vmat.a[i,16],Tsim[j],18.5)
  }
  low.a21[j]<-quantile(dist.a21[,j],0.025)
  high.a21[j]<-quantile(dist.a21[,j],0.975)
}

#N-fixation ~ measurement temperature for 26:20 deg. C growing temperature
dist.a26=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a26<-rep(NA,length(Tsim))
high.a26<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a26[i,j]=beta.lin.all(1,vmat.a[i,11],vmat.a[i,12],vmat.a[i,13],vmat.a[i,14],vmat.a[i,15],vmat.a[i,16],Tsim[j],23.5)
  }
  low.a26[j]<-quantile(dist.a26[,j],0.025)
  high.a26[j]<-quantile(dist.a26[,j],0.975)
}

#N-fixation ~ measurement temperature for 31:25 deg. C growing temperature
dist.a31=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a31<-rep(NA,length(Tsim))
high.a31<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a31[i,j]=beta.lin.all(1,vmat.a[i,11],vmat.a[i,12],vmat.a[i,13],vmat.a[i,14],vmat.a[i,15],vmat.a[i,16],Tsim[j],28.5)
  }
  low.a31[j]<-quantile(dist.a31[,j],0.025)
  high.a31[j]<-quantile(dist.a31[,j],0.975)
}

#Rates at 15 deg. C 95% CI
dist.a15=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.a15<-rep(NA,length(Tg.seq))
high.a15<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.a15[i,j]=beta.lin.all(1,vmat.a[i,11],vmat.a[i,12],vmat.a[i,13],vmat.a[i,14],vmat.a[i,15],vmat.a[i,16],15,Tg.seq[j])
  }
  low.a15[j]<-quantile(dist.a15[,j],0.025)
  high.a15[j]<-quantile(dist.a15[,j],0.975)
}

#Topt 95% CI
dist.aTopt=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.aTopt<-rep(NA,length(Tg.seq))
high.aTopt<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.aTopt[i,j]=vmat.a[i,13]+vmat.a[i,14]*Tg.seq[j]
  }
  low.aTopt[j]<-quantile(dist.aTopt[,j],0.025)
  high.aTopt[j]<-quantile(dist.aTopt[,j],0.975)
}

#Rates at 40 deg. C 95% CI
dist.a40=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.a40<-rep(NA,length(Tg.seq))
high.a40<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.a40[i,j]=beta.lin.all(1,vmat.a[i,11],vmat.a[i,12],vmat.a[i,13],vmat.a[i,14],vmat.a[i,15],vmat.a[i,16],40,Tg.seq[j])
  }
  low.a40[j]<-quantile(dist.a40[,j],0.025)
  high.a40[j]<-quantile(dist.a40[,j],0.975)
}

##
#Gliricidia
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g_linKm=mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_GLSE_linKm),Sigma=vcov(fit_Nase_beta_linTminTopt_GLSE_linKm))

#N-fixation ~ measurement temperature for 21:15 deg. C growing temperature
dist.g21_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g21_linKm<-rep(NA,length(Tsim))
high.g21_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g21_linKm[i,j]=beta.lin.Tmin.Topt(1,vmat.g_linKm[i,11],vmat.g_linKm[i,12],vmat.g_linKm[i,13],vmat.g_linKm[i,14],vmat.g_linKm[i,15],Tsim[j],18.5)
  }
  low.g21_linKm[j]<-quantile(dist.g21_linKm[,j],0.025)
  high.g21_linKm[j]<-quantile(dist.g21_linKm[,j],0.975)
}

#N-fixation ~ measurement temperature for 26:20 deg. C growing temperature
dist.g26_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g26_linKm<-rep(NA,length(Tsim))
high.g26_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g26_linKm[i,j]=beta.lin.Tmin.Topt(1,vmat.g_linKm[i,11],vmat.g_linKm[i,12],vmat.g_linKm[i,13],vmat.g_linKm[i,14],vmat.g_linKm[i,15],Tsim[j],23.5)
  }
  low.g26_linKm[j]<-quantile(dist.g26_linKm[,j],0.025)
  high.g26_linKm[j]<-quantile(dist.g26_linKm[,j],0.975)
}

#N-fixation ~ measurement temperature for 31:25 deg. C growing temperature
dist.g31_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g31_linKm<-rep(NA,length(Tsim))
high.g31_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g31_linKm[i,j]=beta.lin.Tmin.Topt(1,vmat.g_linKm[i,11],vmat.g_linKm[i,12],vmat.g_linKm[i,13],vmat.g_linKm[i,14],vmat.g_linKm[i,15],Tsim[j],28.5)
  }
  low.g31_linKm[j]<-quantile(dist.g31_linKm[,j],0.025)
  high.g31_linKm[j]<-quantile(dist.g31_linKm[,j],0.975)
}

#Rates at 15 deg. C 95% CI
dist.g15_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g15_linKm<-rep(NA,length(Tg.seq))
high.g15_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g15_linKm[i,j]=beta.lin.Tmin.Topt(1,vmat.g_linKm[i,11],vmat.g_linKm[i,12],vmat.g_linKm[i,13],vmat.g_linKm[i,14],vmat.g_linKm[i,15],15,Tg.seq[j])
  }
  low.g15_linKm[j]<-quantile(dist.g15_linKm[,j],0.025)
  high.g15_linKm[j]<-quantile(dist.g15_linKm[,j],0.975)
}

#Topt 95% CI
dist.gTopt_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.gTopt_linKm<-rep(NA,length(Tg.seq))
high.gTopt_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.gTopt_linKm[i,j]=vmat.g_linKm[i,14]+vmat.g_linKm[i,15]*Tg.seq[j]
  }
  low.gTopt_linKm[j]<-quantile(dist.gTopt_linKm[,j],0.025)
  high.gTopt_linKm[j]<-quantile(dist.gTopt_linKm[,j],0.975)
}

#Rates at 40 deg. C 95% CI
dist.g40_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g40_linKm<-rep(NA,length(Tg.seq))
high.g40_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g40_linKm[i,j]=beta.lin.Tmin.Topt(1,vmat.g_linKm[i,11],vmat.g_linKm[i,12],vmat.g_linKm[i,13],vmat.g_linKm[i,14],vmat.g_linKm[i,15],40,Tg.seq[j])
  }
  low.g40_linKm[j]<-quantile(dist.g40_linKm[,j],0.025)
  high.g40_linKm[j]<-quantile(dist.g40_linKm[,j],0.975)
}

##
#Robinia
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r_linKm=mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTmax_ROPS_linKm),Sigma=vcov(fit_Nase_beta_linTminTmax_ROPS_linKm))

#N-fixation ~ measurement temperature for 21:15 deg. C growing temperature
dist.r21_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r21_linKm<-rep(NA,length(Tsim))
high.r21_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r21_linKm[i,j]=beta.lin.Tmin.Tmax(1,vmat.r_linKm[i,11],vmat.r_linKm[i,12],vmat.r_linKm[i,13],vmat.r_linKm[i,14],vmat.r_linKm[i,15],Tsim[j],18.5)
  }
  low.r21_linKm[j]<-quantile(dist.r21_linKm[,j],0.025)
  high.r21_linKm[j]<-quantile(dist.r21_linKm[,j],0.975)
}

#N-fixation ~ measurement temperature for 26:20 deg. C growing temperature
dist.r26_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r26_linKm<-rep(NA,length(Tsim))
high.r26_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r26_linKm[i,j]=beta.lin.Tmin.Tmax(1,vmat.r_linKm[i,11],vmat.r_linKm[i,12],vmat.r_linKm[i,13],vmat.r_linKm[i,14],vmat.r_linKm[i,15],Tsim[j],23.5)
  }
  low.r26_linKm[j]<-quantile(dist.r26_linKm[,j],0.025)
  high.r26_linKm[j]<-quantile(dist.r26_linKm[,j],0.975)
}

#N-fixation ~ measurement temperature for 31:25 deg. C growing temperature
dist.r31_linKm=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r31_linKm<-rep(NA,length(Tsim))
high.r31_linKm<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r31_linKm[i,j]=beta.lin.Tmin.Tmax(1,vmat.r_linKm[i,11],vmat.r_linKm[i,12],vmat.r_linKm[i,13],vmat.r_linKm[i,14],vmat.r_linKm[i,15],Tsim[j],28.5)
  }
  low.r31_linKm[j]<-quantile(dist.r31_linKm[,j],0.025)
  high.r31_linKm[j]<-quantile(dist.r31_linKm[,j],0.975)
}

#Rates at 15 deg. C 95% CI
dist.r15_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r15_linKm<-rep(NA,length(Tg.seq))
high.r15_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r15_linKm[i,j]=beta.lin.Tmin.Tmax(1,vmat.r_linKm[i,11],vmat.r_linKm[i,12],vmat.r_linKm[i,13],vmat.r_linKm[i,14],vmat.r_linKm[i,15],15,Tg.seq[j])
  }
  low.r15_linKm[j]<-quantile(dist.r15_linKm[,j],0.025)
  high.r15_linKm[j]<-quantile(dist.r15_linKm[,j],0.975)
}

#Rates at 40 deg. C 95% CI
dist.r40_linKm=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r40_linKm<-rep(NA,length(Tg.seq))
high.r40_linKm<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r40_linKm[i,j]=beta.lin.Tmin.Tmax(1,vmat.r_linKm[i,11],vmat.r_linKm[i,12],vmat.r_linKm[i,13],vmat.r_linKm[i,14],vmat.r_linKm[i,15],40,Tg.seq[j])
  }
  low.r40_linKm[j]<-quantile(dist.r40_linKm[,j],0.025)
  high.r40_linKm[j]<-quantile(dist.r40_linKm[,j],0.975)
}

####
#Calculate photosynthesis rates at 15 and 40 deg. C as a function of growing temperature
####

###
#Rates at 15 deg. C (Figure 3a)
###

#Morella
m15.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  m15.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[85],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[89],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[77],15,Tg.seq[i])
}

#Alnus
a15.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  a15.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[86],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[90],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[78],15,Tg.seq[i])
}

#Gliricidia
g15.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  g15.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[87],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[91],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[79],15,Tg.seq[i])
}

#Robinia
r15.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  r15.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[88],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[92],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[80],15,Tg.seq[i])
}

###
#Rates at 40 deg. C (Figure 3c)
###

#Morella
m40.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  m40.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[85],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[89],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[77],40,Tg.seq[i])
}

#Alnus
a40.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  a40.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[86],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[90],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[78],40,Tg.seq[i])
}

#Gliricidia
g40.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  g40.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[87],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[91],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[79],40,Tg.seq[i])
}

#Robinia
r40.A275=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  r40.A275[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[88],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[92],
                            coef(fit_Photo_beta_all_Topt.lin_A275)[80],40,Tg.seq[i])
}


####
#Calculate photosynthesis 95% CI
####

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.A275=mvrnorm(1000,mu=coef(fit_Photo_beta_all_Topt.lin_A275),Sigma=vcov(fit_Photo_beta_all_Topt.lin_A275))

###
#Morella
###

#Rates at 15 deg. C 95% CI (Figure 3a)
dist.m15.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m15.A275<-rep(NA,length(Tg.seq))
high.m15.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m15.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,81],vmat.A275[i,85],vmat.A275[i,89],vmat.A275[i,77],15,Tg.seq[j])
  }
  low.m15.A275[j]<-quantile(dist.m15.A275[,j],0.025)
  high.m15.A275[j]<-quantile(dist.m15.A275[,j],0.975)
}

#Topt 95% CI (Figure 3b)
dist.mTopt.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.mTopt.A275<-rep(NA,length(Tg.seq))
high.mTopt.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.mTopt.A275[i,j]=vmat.A275[i,85]+vmat.A275[i,89]*Tg.seq[j]
  }
  low.mTopt.A275[j]<-quantile(dist.mTopt.A275[,j],0.025)
  high.mTopt.A275[j]<-quantile(dist.mTopt.A275[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Figure 3c)
dist.m40.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m40.A275<-rep(NA,length(Tg.seq))
high.m40.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m40.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,81],vmat.A275[i,85],vmat.A275[i,89],vmat.A275[i,77],40,Tg.seq[j])
  }
  low.m40.A275[j]<-quantile(dist.m40.A275[,j],0.025)
  high.m40.A275[j]<-quantile(dist.m40.A275[,j],0.975)
}

###
#Alnus
###

#Rates at 15 deg. C 95% CI (Figure 3a)
dist.a15.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.a15.A275<-rep(NA,length(Tg.seq))
high.a15.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.a15.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,82],vmat.A275[i,86],vmat.A275[i,90],vmat.A275[i,78],15,Tg.seq[j])
  }
  low.a15.A275[j]<-quantile(dist.a15.A275[,j],0.025)
  high.a15.A275[j]<-quantile(dist.a15.A275[,j],0.975)
}

#Topt 95% CI (Figure 3b)
dist.aTopt.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.aTopt.A275<-rep(NA,length(Tg.seq))
high.aTopt.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.aTopt.A275[i,j]=vmat.A275[i,86]+vmat.A275[i,90]*Tg.seq[j]
  }
  low.aTopt.A275[j]<-quantile(dist.aTopt.A275[,j],0.025)
  high.aTopt.A275[j]<-quantile(dist.aTopt.A275[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Figure 3c)
dist.a40.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.a40.A275<-rep(NA,length(Tg.seq))
high.a40.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.a40.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,82],vmat.A275[i,86],vmat.A275[i,90],vmat.A275[i,78],40,Tg.seq[j])
  }
  low.a40.A275[j]<-quantile(dist.a40.A275[,j],0.025)
  high.a40.A275[j]<-quantile(dist.a40.A275[,j],0.975)
}

###
#Gliricidia
###

#Rates at 15 deg. C 95% CI (Figure 3a)
dist.g15.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g15.A275<-rep(NA,length(Tg.seq))
high.g15.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g15.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,83],vmat.A275[i,87],vmat.A275[i,91],vmat.A275[i,79],15,Tg.seq[j])
  }
  low.g15.A275[j]<-quantile(dist.g15.A275[,j],0.025)
  high.g15.A275[j]<-quantile(dist.g15.A275[,j],0.975)
}

#Topt 95% CI (Figure 3b)
dist.gTopt.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.gTopt.A275<-rep(NA,length(Tg.seq))
high.gTopt.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.gTopt.A275[i,j]=vmat.A275[i,87]+vmat.A275[i,91]*Tg.seq[j]
  }
  low.gTopt.A275[j]<-quantile(dist.gTopt.A275[,j],0.025)
  high.gTopt.A275[j]<-quantile(dist.gTopt.A275[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Figure 3c)
dist.g40.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g40.A275<-rep(NA,length(Tg.seq))
high.g40.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g40.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,83],vmat.A275[i,87],vmat.A275[i,91],vmat.A275[i,79],40,Tg.seq[j])
  }
  low.g40.A275[j]<-quantile(dist.g40.A275[,j],0.025)
  high.g40.A275[j]<-quantile(dist.g40.A275[,j],0.975)
}

###
#Robinia
###

#Rates at 15 deg. C 95% CI (Figure 3a)
dist.r15.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r15.A275<-rep(NA,length(Tg.seq))
high.r15.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r15.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,84],vmat.A275[i,88],vmat.A275[i,92],vmat.A275[i,80],15,Tg.seq[j])
  }
  low.r15.A275[j]<-quantile(dist.r15.A275[,j],0.025)
  high.r15.A275[j]<-quantile(dist.r15.A275[,j],0.975)
}

#Topt 95% CI (Figure 3b)
dist.rTopt.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.rTopt.A275<-rep(NA,length(Tg.seq))
high.rTopt.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.rTopt.A275[i,j]=vmat.A275[i,88]+vmat.A275[i,92]*Tg.seq[j]
  }
  low.rTopt.A275[j]<-quantile(dist.rTopt.A275[,j],0.025)
  high.rTopt.A275[j]<-quantile(dist.rTopt.A275[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Figure 3c)
dist.r40.A275=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r40.A275<-rep(NA,length(Tg.seq))
high.r40.A275<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r40.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,84],vmat.A275[i,88],vmat.A275[i,92],vmat.A275[i,80],40,Tg.seq[j])
  }
  low.r40.A275[j]<-quantile(dist.r40.A275[,j],0.025)
  high.r40.A275[j]<-quantile(dist.r40.A275[,j],0.975)
}

####
#Calculate range of N-fixation data
####

#Morella
M.min<-min(c(MOCE21a$Temperature,MOCE21b$Temperature,MOCE21c$Temperature,
             MOCE26a$Temperature,MOCE26b$Temperature,MOCE26c$Temperature,
             MOCE31a$Temperature,MOCE31b$Temperature,MOCE31c$Temperature))
M.max<-max(c(MOCE21a$Temperature,MOCE21b$Temperature,MOCE21c$Temperature,
             MOCE26a$Temperature,MOCE26b$Temperature,MOCE26c$Temperature,
             MOCE31a$Temperature,MOCE31b$Temperature,MOCE31c$Temperature))

#Alnus
A.min<-min(c(ALRU21a$Temperature,ALRU21b$Temperature,ALRU21c$Temperature,
             ALRU26a$Temperature,ALRU26b$Temperature,ALRU26c$Temperature,
             ALRU31a$Temperature,ALRU31b$Temperature,ALRU31c$Temperature))
A.max<-max(c(ALRU21a$Temperature,ALRU21b$Temperature,ALRU21c$Temperature,
             ALRU26a$Temperature,ALRU26b$Temperature,ALRU26c$Temperature,
             ALRU31a$Temperature,ALRU31b$Temperature,ALRU31c$Temperature))

#Gliricidia
G.min<-min(c(GLSE21a$Temperature,GLSE21b$Temperature,GLSE21c$Temperature,
             GLSE26a$Temperature,GLSE26b$Temperature,GLSE26c$Temperature,
             GLSE31a$Temperature,GLSE31b$Temperature,GLSE31c$Temperature))
G.max<-max(c(GLSE21a$Temperature,GLSE21b$Temperature,GLSE21c$Temperature,
             GLSE26a$Temperature,GLSE26b$Temperature,GLSE26c$Temperature,
             GLSE31a$Temperature,GLSE31b$Temperature,GLSE31c$Temperature))

#Robinia
R.min<-min(c(ROPS21a$Temperature,ROPS21b$Temperature,ROPS21c$Temperature,
             ROPS26a$Temperature,ROPS26b$Temperature,ROPS26c$Temperature,
             ROPS31a$Temperature,ROPS31b$Temperature,ROPS31c$Temperature))
R.max<-max(c(ROPS21a$Temperature,ROPS21b$Temperature,ROPS21c$Temperature,
             ROPS26a$Temperature,ROPS26b$Temperature,ROPS26c$Temperature,
             ROPS31a$Temperature,ROPS31b$Temperature,ROPS31c$Temperature))

###############################################################################################################
#Supplementary Figure 15
###############################################################################################################

#PDF dimension is 9x7 inches  

#Plotting Settings
par(pty="s")
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

#S15a
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
points(MOCE21a$Temperature,MOCE21a$Vmax/(9960.998949*ft.m_linKm[2]),pch=16,cex=0.001,col="gray")
points(MOCE21b$Temperature,MOCE21b$Vmax/(3833.121977*ft.m_linKm[3]),pch=16,cex=0.001,col="gray")
points(MOCE21c$Temperature,MOCE21c$Vmax/(3095.456985*ft.m_linKm[4]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.m_linKm[12]*18.5+ft.m_linKm[11],ft.m_linKm[14]*18.5+ft.m_linKm[13],ft.m_linKm[16]*18.5+ft.m_linKm[15],x),
      from=max(M.min,ft.m_linKm[12]*18.5+ft.m_linKm[11]),to=min(M.max,ft.m_linKm[16]*18.5+ft.m_linKm[15]),lwd=2,col="dodgerblue1",lty=1,add=TRUE)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)

#S15b
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(ALRU21a$Temperature,ALRU21a$Vmax/(10035.24751*ft.a[2]),pch=16,cex=0.001,col="gray")
points(ALRU21b$Temperature,ALRU21b$Vmax/(7990.765838*ft.a[3]),pch=16,cex=0.001,col="gray")
points(ALRU21c$Temperature,ALRU21c$Vmax/(4046.082544*ft.a[4]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.a[12]*18.5+ft.a[11],ft.a[14]*18.5+ft.a[13],ft.a[16]*18.5+ft.a[15],x),
      from=max(A.min,ft.a[12]*18.5+ft.a[11]),to=min(A.max,ft.a[16]*18.5+ft.a[15]),lwd=2,col="dodgerblue1",lty=1,add=TRUE)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)

#S15c
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(GLSE21a$Temperature,GLSE21a$Vmax/(308.3421323*ft.g_linKm[2]),pch=16,cex=0.001,col="gray")
points(GLSE21b$Temperature,GLSE21b$Vmax/(2028.57371*ft.g_linKm[3]),pch=16,cex=0.001,col="gray")
points(GLSE21c$Temperature,GLSE21c$Vmax/(879.0395611*ft.g_linKm[4]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.g_linKm[13]*18.5+ft.g_linKm[12],ft.g_linKm[15]*18.5+ft.g_linKm[14],ft.g_linKm[11],x),
      from=max(G.min,ft.g_linKm[13]*18.5+ft.g_linKm[12]),to=min(G.max,ft.g_linKm[11]),lwd=2,col="dodgerblue1",lty=1,add=TRUE)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)

#S15d
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(ROPS21a$Temperature,ROPS21a$Vmax/(1922.245351*ft.r_linKm[2]),pch=16,cex=0.001,col="gray")
points(ROPS21b$Temperature,ROPS21b$Vmax/(1690.867422*ft.r_linKm[3]),pch=16,cex=0.001,col="gray")
points(ROPS21c$Temperature,ROPS21c$Vmax/(450.8292229*ft.r_linKm[4]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.r_linKm[13]*18.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*18.5+ft.r_linKm[14],x),
      from=max(R.min,ft.r_linKm[13]*18.5+ft.r_linKm[12]),to=min(R.max,ft.r_linKm[15]*18.5+ft.r_linKm[14]),lwd=2,col="dodgerblue1",lty=1,add=TRUE)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)

#S15e
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
points(MOCE26a$Temperature,MOCE26a$Vmax/(6492.677546*ft.m_linKm[5]),pch=16,cex=0.001,col="gray")
points(MOCE26b$Temperature,MOCE26b$Vmax/(5974.076428*ft.m_linKm[6]),pch=16,cex=0.001,col="gray")
points(MOCE26c$Temperature,MOCE26c$Vmax/(7890.763063*ft.m_linKm[7]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.m_linKm[12]*23.5+ft.m_linKm[11],ft.m_linKm[14]*23.5+ft.m_linKm[13],ft.m_linKm[16]*23.5+ft.m_linKm[15],x),
      from=max(M.min,ft.m_linKm[12]*23.5+ft.m_linKm[11]),to=min(M.max,ft.m_linKm[16]*23.5+ft.m_linKm[15]),lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S15f
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(ALRU26a$Temperature,ALRU26a$Vmax/(2236.979133*ft.a[5]),pch=16,cex=0.001,col="gray")
points(ALRU26b$Temperature,ALRU26b$Vmax/(2744.252462*ft.a[6]),pch=16,cex=0.001,col="gray")
points(ALRU26c$Temperature,ALRU26c$Vmax/(1831.763172*ft.a[7]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.a[12]*23.5+ft.a[11],ft.a[14]*23.5+ft.a[13],ft.a[16]*23.5+ft.a[15],x),
      from=max(A.min,ft.a[12]*23.5+ft.a[11]),to=min(A.max,ft.a[16]*23.5+ft.a[15]),lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S15g
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(GLSE26a$Temperature,GLSE26a$Vmax/(909.4624356*ft.g_linKm[5]),pch=16,cex=0.001,col="gray")
points(GLSE26b$Temperature,GLSE26b$Vmax/(1373.906634*ft.g_linKm[6]),pch=16,cex=0.001,col="gray")
points(GLSE26c$Temperature,GLSE26c$Vmax/(705.9735103*ft.g_linKm[7]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.g_linKm[13]*23.5+ft.g_linKm[12],ft.g_linKm[15]*23.5+ft.g_linKm[14],ft.g_linKm[11],x),
      from=max(G.min,ft.g_linKm[13]*23.5+ft.g_linKm[12]),to=min(G.max,ft.g_linKm[11]),lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S15h
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(ROPS26a$Temperature,ROPS26a$Vmax/(2785.452785*ft.r_linKm[5]),pch=16,cex=0.001,col="gray")
points(ROPS26b$Temperature,ROPS26b$Vmax/(1572.439251*ft.r_linKm[6]),pch=16,cex=0.001,col="gray")
points(ROPS26c$Temperature,ROPS26c$Vmax/(1686.352*ft.r_linKm[7]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.r_linKm[13]*23.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*23.5+ft.r_linKm[14],x),
      from=max(R.min,ft.r_linKm[13]*23.5+ft.r_linKm[12]),to=min(R.max,ft.r_linKm[15]*23.5+ft.r_linKm[14]),lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

#S15i
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
points(MOCE31a$Temperature,MOCE31a$Vmax/(8826.601393*ft.m_linKm[8]),pch=16,cex=0.001,col="gray")
points(MOCE31b$Temperature,MOCE31b$Vmax/(2371.923752*ft.m_linKm[9]),pch=16,cex=0.001,col="gray")
points(MOCE31c$Temperature,MOCE31c$Vmax/(629.7544222*ft.m_linKm[10]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.m_linKm[12]*28.5+ft.m_linKm[11],ft.m_linKm[14]*28.5+ft.m_linKm[13],ft.m_linKm[16]*28.5+ft.m_linKm[15],x),
      from=max(M.min,ft.m_linKm[12]*28.5+ft.m_linKm[11]),to=min(M.max,ft.m_linKm[16]*28.5+ft.m_linKm[15]),lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S15j
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(ALRU31a$Temperature,ALRU31a$Vmax/(406.3524389*ft.a[8]),pch=16,cex=0.001,col="gray")
points(ALRU31b$Temperature,ALRU31b$Vmax/(1832.723413*ft.a[9]),pch=16,cex=0.001,col="gray")
points(ALRU31c$Temperature,ALRU31c$Vmax/(1350.729381*ft.a[10]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.a[12]*28.5+ft.a[11],ft.a[14]*28.5+ft.a[13],ft.a[16]*28.5+ft.a[15],x),
      from=max(A.min,ft.a[12]*28.5+ft.a[11]),to=min(A.max,ft.a[16]*28.5+ft.a[15]),lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S15k
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(GLSE31a$Temperature,GLSE31a$Vmax/(1196.3736273*ft.g_linKm[8]),pch=16,cex=0.001,col="gray")
points(GLSE31b$Temperature,GLSE31b$Vmax/(312.4607561*ft.g_linKm[9]),pch=16,cex=0.001,col="gray")
points(GLSE31c$Temperature,GLSE31c$Vmax/(209.4577018*ft.g_linKm[10]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.g_linKm[13]*28.5+ft.g_linKm[12],ft.g_linKm[15]*28.5+ft.g_linKm[14],ft.g_linKm[11],x),
      from=max(G.min,ft.g_linKm[13]*28.5+ft.g_linKm[12]),to=min(G.max,ft.g_linKm[11]),lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S15l
plot(0:10,0:10,t="n",ylim=c(0,1.6),xlim=c(5,50),xlab=NA,ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
points(ROPS31a$Temperature,ROPS31a$Vmax/(1646.769302*ft.r_linKm[8]),pch=16,cex=0.001,col="gray")
points(ROPS31b$Temperature,ROPS31b$Vmax/(920.7416901*ft.r_linKm[9]),pch=16,cex=0.001,col="gray")
points(ROPS31c$Temperature,ROPS31c$Vmax/(1978.187071*ft.r_linKm[10]),pch=16,cex=0.001,col="gray")
curve(beta(1,ft.r_linKm[13]*28.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*28.5+ft.r_linKm[14],x),
      from=max(R.min,ft.r_linKm[13]*28.5+ft.r_linKm[12]),to=min(R.max,ft.r_linKm[15]*28.5+ft.r_linKm[14]),lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression('Nitrogen fixation (% of max)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)

###############################################################################################################
#Supplementary Figure 16
###############################################################################################################

####
#Calculate relative N-fixation rate at 25C
####

#Houlton et al. 2008 function
Houl.25<-Houlton_func(25)

#Morella
SNF.M21.25_linKm<-beta(1,ft.m_linKm[12]*18.5+ft.m_linKm[11],ft.m_linKm[14]*18.5+ft.m_linKm[13],ft.m_linKm[16]*18.5+ft.m_linKm[15],25)
SNF.M26.25_linKm<-beta(1,ft.m_linKm[12]*23.5+ft.m_linKm[11],ft.m_linKm[14]*23.5+ft.m_linKm[13],ft.m_linKm[16]*23.5+ft.m_linKm[15],25)
SNF.M31.25_linKm<-beta(1,ft.m_linKm[12]*28.5+ft.m_linKm[11],ft.m_linKm[14]*28.5+ft.m_linKm[13],ft.m_linKm[16]*28.5+ft.m_linKm[15],25)

#Alnus
SNF.A21.25<-beta(1,ft.a[12]*18.5+ft.a[11],ft.a[14]*18.5+ft.a[13],ft.a[16]*18.5+ft.a[15],25)
SNF.A26.25<-beta(1,ft.a[12]*23.5+ft.a[11],ft.a[14]*23.5+ft.a[13],ft.a[16]*23.5+ft.a[15],25)
SNF.A31.25<-beta(1,ft.a[12]*28.5+ft.a[11],ft.a[14]*28.5+ft.a[13],ft.a[16]*28.5+ft.a[15],25)

#Gliricidia
SNF.G21.25_linKm<-beta(1,ft.g_linKm[13]*18.5+ft.g_linKm[12],ft.g_linKm[15]*18.5+ft.g_linKm[14],ft.g_linKm[11],25)
SNF.G26.25_linKm<-beta(1,ft.g_linKm[13]*23.5+ft.g_linKm[12],ft.g_linKm[15]*23.5+ft.g_linKm[14],ft.g_linKm[11],25)
SNF.G31.25_linKm<-beta(1,ft.g_linKm[13]*28.5+ft.g_linKm[12],ft.g_linKm[15]*28.5+ft.g_linKm[14],ft.g_linKm[11],25)

#Robinia
SNF.R21.25_linKm<-beta(1,ft.r_linKm[13]*18.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*18.5+ft.r_linKm[14],25)
SNF.R26.25_linKm<-beta(1,ft.r_linKm[13]*23.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*23.5+ft.r_linKm[14],25)
SNF.R31.25_linKm<-beta(1,ft.r_linKm[13]*28.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*28.5+ft.r_linKm[14],25)

####
#Calculate relative photosynthesis rate at 25C
####

#Morella
A275.M21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],25,18.5)
A275.M26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],25,23.5)
A275.M31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],25,28.5)

#Alnus
A275.A21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],25,18.5)
A275.A26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],25,23.5)
A275.A31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],25,28.5)

#Gliricidia
A275.G21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],25,18.5)
A275.G26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],25,23.5)
A275.G31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],25,28.5)

#Robinia
A275.R21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],25,18.5)
A275.R26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],25,23.5)
A275.R31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],25,28.5)

#PDF dimension is 7x13 inches  

#Plotting Settings
par(pty="s")
par(mfrow=c(2,1))
par(mar=c(0.3,0.8,0.8,0.5))
par(oma=c(5,5.5,1,0.5))

#S16a
plot(0:50,(0:50)/50,col="white",xlab=NA,ylab=NA,las=1,ylim=c(0,3.5),cex.axis=1.2,xaxt="n")
axis(1,at=c(0,10,20,30,40,50),labels=F)
mtext(expression('Nitrogen fixation'),side=2,line=4.5,cex=1.5)
mtext(expression('(normalized to 1 at 25 '*degree*'C)'),side=2,line=3,cex=1.5)
curve(Houlton_func_norm(x)/Houl.25,from=0,to=50,lwd=3,col="black",add=TRUE)
curve(beta(1/SNF.M21.25_linKm,ft.m_linKm[12]*18.5+ft.m_linKm[11],ft.m_linKm[14]*18.5+ft.m_linKm[13],ft.m_linKm[16]*18.5+ft.m_linKm[15],x),from=max(0,ft.m_linKm[12]*18.5+ft.m_linKm[11]),to=ft.m_linKm[16]*18.5+ft.m_linKm[15],lwd=3,col="dodgerblue1",lty=1,add=TRUE)
curve(beta(1/SNF.M26.25_linKm,ft.m_linKm[12]*23.5+ft.m_linKm[11],ft.m_linKm[14]*23.5+ft.m_linKm[13],ft.m_linKm[16]*23.5+ft.m_linKm[15],x),from=max(0,ft.m_linKm[12]*23.5+ft.m_linKm[11]),to=ft.m_linKm[16]*23.5+ft.m_linKm[15],lwd=3,col="gold1",lty=1,add=TRUE)
curve(beta(1/SNF.M31.25_linKm,ft.m_linKm[12]*28.5+ft.m_linKm[11],ft.m_linKm[14]*28.5+ft.m_linKm[13],ft.m_linKm[16]*28.5+ft.m_linKm[15],x),from=max(0,ft.m_linKm[12]*28.5+ft.m_linKm[11]),to=ft.m_linKm[16]*28.5+ft.m_linKm[15],lwd=3,col="orangered3",lty=1,add=TRUE)
curve(beta(1/SNF.A21.25,ft.a[12]*18.5+ft.a[11],ft.a[14]*18.5+ft.a[13],ft.a[16]*18.5+ft.a[15],x),from=max(0,ft.a[12]*18.5+ft.a[11]),to=ft.a[16]*18.5+ft.a[15],lwd=3,col="dodgerblue1",lty=2,add=TRUE)
curve(beta(1/SNF.A26.25,ft.a[12]*23.5+ft.a[11],ft.a[14]*23.5+ft.a[13],ft.a[16]*23.5+ft.a[15],x),from=max(0,ft.a[12]*23.5+ft.a[11]),to=ft.a[16]*23.5+ft.a[15],lwd=3,col="gold1",lty=2,add=TRUE)
curve(beta(1/SNF.A31.25,ft.a[12]*28.5+ft.a[11],ft.a[14]*28.5+ft.a[13],ft.a[16]*28.5+ft.a[15],x),from=max(0,ft.a[12]*28.5+ft.a[11]),to=ft.a[16]*28.5+ft.a[15],lwd=3,col="orangered3",lty=2,add=TRUE)
curve(beta(1/SNF.G21.25_linKm,ft.g_linKm[13]*18.5+ft.g_linKm[12],ft.g_linKm[15]*18.5+ft.g_linKm[14],ft.g_linKm[11],x),from=max(0,ft.g_linKm[13]*18.5+ft.g_linKm[12]),to=ft.g_linKm[11],lwd=3,col="dodgerblue1",lty=3,add=TRUE)
curve(beta(1/SNF.G26.25_linKm,ft.g_linKm[13]*23.5+ft.g_linKm[12],ft.g_linKm[15]*23.5+ft.g_linKm[14],ft.g_linKm[11],x),from=max(0,ft.g_linKm[13]*23.5+ft.g_linKm[12]),to=ft.g_linKm[11],lwd=3,col="gold1",lty=3,add=TRUE)
curve(beta(1/SNF.G31.25_linKm,ft.g_linKm[13]*28.5+ft.g_linKm[12],ft.g_linKm[15]*28.5+ft.g_linKm[14],ft.g_linKm[11],x),from=max(0,ft.g_linKm[13]*28.5+ft.g_linKm[12]),to=ft.g_linKm[11],lwd=3,col="orangered3",lty=3,add=TRUE)
curve(beta(1/SNF.R21.25_linKm,ft.r_linKm[13]*18.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*18.5+ft.r_linKm[14],x),from=max(0,ft.r_linKm[12]*18.5+ft.r_linKm[11]),to=ft.r_linKm[15]*18.5+ft.r_linKm[14],lwd=3,col="dodgerblue1",lty=4,add=TRUE)
curve(beta(1/SNF.R26.25_linKm,ft.r_linKm[13]*23.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*23.5+ft.r_linKm[14],x),from=max(0,ft.r_linKm[12]*23.5+ft.r_linKm[11]),to=ft.r_linKm[15]*23.5+ft.r_linKm[14],lwd=3,col="gold1",lty=4,add=TRUE)
curve(beta(1/SNF.R31.25_linKm,ft.r_linKm[13]*28.5+ft.r_linKm[12],ft.r_linKm[11],ft.r_linKm[15]*28.5+ft.r_linKm[14],x),from=max(0,ft.r_linKm[12]*28.5+ft.r_linKm[11]),to=ft.r_linKm[15]*28.5+ft.r_linKm[14],lwd=3,col="orangered3",lty=4,add=TRUE)
legend(-3,3.7,c(expression('Houlton '*italic('et al.')*' (2008)'),expression(underline(bolditalic("Morella"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Alnus"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Gliricidia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Robinia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C')),
       col=c("black",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3"),
       lty=c(1,NA,1,1,1,NA,2,2,2,NA,3,3,3,NA,4,4,4),bty="n",lwd=3,y.intersp = 0.5,cex=1,seg.len=1.2,x.intersp = 0.5)
mtext(text="a",side=3,cex=1.4,adj=0)

#S16b
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(0,50),ylim=c(0,1.3),cex.lab=1.5,cex.axis=1.2) #,xaxt="n"
curve(beta.Topt.lin(1/A275.M21.25,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],x,18.5),from=10,to=40,lty=1,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.A21.25,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],x,18.5),from=10,to=40,lty=2,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.G21.25,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],x,18.5),from=10,to=40,lty=3,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.R21.25,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],x,18.5),from=10,to=40,lty=4,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.M26.25,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],x,23.5),from=10,to=40,lty=1,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.A26.25,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],x,23.5),from=10,to=40,lty=2,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.G26.25,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],x,23.5),from=10,to=40,lty=3,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.R26.25,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],x,23.5),from=10,to=40,lty=4,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.M31.25,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],x,28.5),from=10,to=40,lty=1,col="orangered3",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.A31.25,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],x,28.5),from=10,to=40,lty=2,col="orangered3",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.G31.25,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],x,28.5),from=10,to=40,lty=3,col="orangered3",add=T,lwd=3)
curve(beta.Topt.lin(1/A275.R31.25,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],x,28.5),from=10,to=40,lty=4,col="orangered3",add=T,lwd=3)
mtext(expression('Photosynthesis'),side=2,line=4.5,cex=1.5)
mtext(expression('(normalized to 1 at 25 '*degree*'C)'),side=2,line=3,cex=1.5)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=F)
mtext(text="b",side=3,cex=1.4,adj=0)

###############################################################################################################
#Figure 17
###############################################################################################################

#PDF dimension is 10.5x3 inches  

#Plotting Settings
nf<-layout(matrix(c(1,2,3,4),1,4,byrow=T),c(3,3,3,1.3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,4,0))
par(pty="s")

#S17a
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(0,100),cex.lab=1.5,cex.axis=1.2)
mtext(text="a",side=3,cex=1.2,adj=0)
mtext(expression('Rate at 15 '*degree*'C'),side=2,cex=1.2,line=4.25)
mtext(expression('(% of max)'),side=2,cex=1.2,line=2.5)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m15_linKm,rev(high.m15_linKm))*100,
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a15,rev(high.a15))*100,
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g15_linKm,rev(high.g15_linKm))*100,
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r15_linKm,rev(high.r15_linKm))*100,
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m15.A275*100,rev(high.m15.A275*100)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a15.A275*100,rev(high.a15.A275*100)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g15.A275*100,rev(high.g15.A275*100)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r15.A275*100,rev(high.r15.A275*100)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
points(Tg.seq,M15_linKm*100,type="l",lty=1,col="darkorange1",lwd=1)
points(Tg.seq,A15*100,type="l",lty=1,col="darkturquoise",lwd=1)
points(Tg.seq,G15_linKm*100,type="l",lty=1,col="orangered2",lwd=1)
points(Tg.seq,R15_linKm*100,type="l",lty=1,col="dodgerblue3",lwd=1)
points(Tg.seq,m15.A275*100,type="l",lty=2,col="darkorange1",lwd=1)
points(Tg.seq,a15.A275*100,type="l",lty=2,col="darkturquoise",lwd=1)
points(Tg.seq,g15.A275*100,type="l",lty=2,col="orangered2",lwd=1)
points(Tg.seq,r15.A275*100,type="l",lty=2,col="dodgerblue3",lwd=1)
points(c(M15_linKm[which(Tg.seq==18.5)],M15_linKm[which(Tg.seq==23.5)],M15_linKm[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkorange1",cex=1.5,lwd=1)
points(c(A15[which(Tg.seq==18.5)],A15[which(Tg.seq==23.5)],A15[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(c(G15_linKm[which(Tg.seq==18.5)],G15_linKm[which(Tg.seq==23.5)],G15_linKm[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="orangered2",cex=1.5,lwd=1)
points(c(R15_linKm[which(Tg.seq==18.5)],R15_linKm[which(Tg.seq==23.5)],R15_linKm[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="dodgerblue3",cex=1.5,lwd=1)
points(c(m15.A275[which(Tg.seq==18.5)],m15.A275[which(Tg.seq==23.5)],m15.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkorange1",cex=1.5,lwd=1)
points(c(a15.A275[which(Tg.seq==18.5)],a15.A275[which(Tg.seq==23.5)],a15.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkturquoise",cex=1.5,lwd=1)
points(c(r15.A275[which(Tg.seq==18.5)],r15.A275[which(Tg.seq==23.5)],r15.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="dodgerblue3",cex=1.5,lwd=1)
points(c(g15.A275[which(Tg.seq==18.5)],g15.A275[which(Tg.seq==23.5)],g15.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="orangered2",cex=1.5,lwd=1)

#S17b
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(20,40),cex.lab=1.5,cex.axis=1.2)
mtext(text="b",side=3,cex=1.2,adj=0)
mtext(expression(italic('T')[opt]*' ('*degree*'C)'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.mTopt_linKm,rev(high.mTopt_linKm)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.aTopt,rev(high.aTopt)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.gTopt_linKm,rev(high.gTopt_linKm)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(rep(32.99,201),rev(rep(33.03,201))),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.mTopt.A275,rev(high.mTopt.A275)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.aTopt.A275,rev(high.aTopt.A275)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.gTopt.A275,rev(high.gTopt.A275)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.rTopt.A275,rev(high.rTopt.A275)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
abline(a=19.74,b=0.633,col="darkorange1",lwd=1)
abline(a=31.25,b=0.061,col="darkturquoise",lwd=1)
abline(a=26.64,b=0.323,col="orangered2",lwd=1)
abline(h=33.01,col="dodgerblue3",lwd=1)
abline(a=22.83,b=0.238,col="darkorange1",lwd=1,lty=2)
abline(a=19.10,b=0.339,col="darkturquoise",lwd=1,lty=2)
abline(a=22.87,b=0.194,col="dodgerblue3",lwd=1,lty=2)
abline(a=12.90,b=0.637,col="orangered2",lwd=1,lty=2)
points(19.74+0.633*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=16,col="darkorange1",cex=1.5,lwd=1)
points(31.25+0.061*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(26.64+0.323*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=17,col="orangered2",cex=1.5,lwd=1)
points(33.01+0*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=17,col="dodgerblue3",cex=1.5,lwd=1)
points(22.83+0.238*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=1,col="darkorange1",cex=1.5,lwd=1)
points(19.10+0.339*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=1,col="darkturquoise",cex=1.5,lwd=1)
points(22.87+0.194*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=2,col="dodgerblue3",cex=1.5,lwd=1)
points(12.90+0.637*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=2,col="orangered2",cex=1.5,lwd=1)

#S17c
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(0,100),cex.lab=1.5,cex.axis=1.2)
mtext(text="c",side=3,cex=1.2,adj=0)
mtext(expression('Rate at 40 '*degree*'C'),side=2,cex=1.2,line=4.25)
mtext(expression('(% of max)'),side=2,cex=1.2,line=2.5)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m40_linKm,rev(high.m40_linKm))*100,
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a40,rev(high.a40))*100,
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g40_linKm,rev(high.g40_linKm))*100,
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r40_linKm,rev(high.r40_linKm))*100,
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m40.A275*100,rev(high.m40.A275*100)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a40.A275*100,rev(high.a40.A275*100)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g40.A275*100,rev(high.g40.A275*100)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r40.A275*100,rev(high.r40.A275*100)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
points(Tg.seq,M40_linKm*100,type="l",lty=1,col="darkorange1",lwd=1)
points(Tg.seq,A40*100,type="l",lty=1,col="darkturquoise",lwd=1)
points(Tg.seq,G40_linKm*100,type="l",lty=1,col="orangered2",lwd=1)
points(Tg.seq,R40_linKm*100,type="l",lty=1,col="dodgerblue3",lwd=1)
points(Tg.seq,m40.A275*100,type="l",lty=2,col="darkorange1",lwd=1)
points(Tg.seq,a40.A275*100,type="l",lty=2,col="darkturquoise",lwd=1)
points(Tg.seq,g40.A275*100,type="l",lty=2,col="orangered2",lwd=1)
points(Tg.seq,r40.A275*100,type="l",lty=2,col="dodgerblue3",lwd=1)
points(c(M40_linKm[which(Tg.seq==18.5)],M40_linKm[which(Tg.seq==23.5)],M40_linKm[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkorange1",cex=1.5,lwd=1)
points(c(A40[which(Tg.seq==18.5)],A40[which(Tg.seq==23.5)],A40[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(c(G40_linKm[which(Tg.seq==18.5)],G40_linKm[which(Tg.seq==23.5)],G40_linKm[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="orangered2",cex=1.5,lwd=1)
points(c(R40_linKm[which(Tg.seq==18.5)],R40_linKm[which(Tg.seq==23.5)],R40_linKm[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="dodgerblue3",cex=1.5,lwd=1)
points(c(m40.A275[which(Tg.seq==18.5)],m40.A275[which(Tg.seq==23.5)],m40.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkorange1",cex=1.5,lwd=1)
points(c(a40.A275[which(Tg.seq==18.5)],a40.A275[which(Tg.seq==23.5)],a40.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkturquoise",cex=1.5,lwd=1)
points(c(r40.A275[which(Tg.seq==18.5)],r40.A275[which(Tg.seq==23.5)],r40.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="dodgerblue3",cex=1.5,lwd=1)
points(c(g40.A275[which(Tg.seq==18.5)],g40.A275[which(Tg.seq==23.5)],g40.A275[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="orangered2",cex=1.5,lwd=1)

#Change plotting settings for legend
par(mar=c(0,0,0,0))
par(xpd=TRUE)

#Legend
plot(0:10, 0:10, type='n', bty='n', xaxt='n', yaxt='n',xlab=NA,ylab=NA)
legend("left",legend=c(expression(underline(bold('N-fixation'))),expression(italic(Morella)),expression(italic(Alnus)),expression(italic(Gliricidia)),expression(italic(Robinia)),
                       expression(underline(bold('Photosynthesis'))),expression(italic(Morella)),expression(italic(Alnus)),expression(italic(Gliricidia)),expression(italic(Robinia))),
       col=c(NA,"darkorange1","darkturquoise","orangered2","dodgerblue3",NA,"darkorange1","darkturquoise","orangered2","dodgerblue3"),
       lty=c(NA,1,1,1,1,NA,2,2,2,2),pch=c(NA,16,16,17,17,NA,1,1,2,2),bty="n",pt.cex=1.5,lwd=1,x.intersp = 0.2,y.intersp = 1,
       cex=1.1)
par(xpd=FALSE)