###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figures 1,3, and 4
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
#Read N-fixation data in from "SNF_Temp" folder
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

####
#Read photosynthesis data (Asat) in from "Photo_Temp" folder 
####

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

####
#Define functions
####

#No acclimation (used for calculating N-fixation rates at 15 and 40 deg. C [Supplementary Figure 4])
beta <- function(ymax,Tmin,Topt,Tmax,T){
  y <- pmax(0,ymax*(Tmax-T)/(Tmax-Topt)*(((T-Tmin)/(Topt-Tmin))^((Topt-Tmin)/(Tmax-Topt))))
  for(i in 1:length(y)){
    if(is.nan(y[i])){
      y[i] <- 0
    }
  }
  y
}

#Acclimation of all parameters (used for N-fixation with Morella, Alnus, and Robinia)
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

#Acclimation of Topt (used for photosynthesis [Asat])
beta.Topt.lin <- function(ymax,Tmin,a,b,Tmax,T,Tgrow){
  y <- ymax*(Tmax-T)/(Tmax-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/(Tmax-(a+b*Tgrow))))
  y
}

####
#Negative log-likelihood (NLL) functions
####

#NLL function for acclimation of all parameters (used for N-fixation with Morella, Alnus, and Robinia)
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

#Robinia
fit_Nase_beta_linall_ROPS <- mle2(Nase_beta_linall_normNLL,start=list(sdNase=-1,ymax21a=0.96320475,ymax21b=0.82996072,ymax21c=1.08051641,
                                                                      ymax26a=0.89390710,ymax26b=0.95916073,ymax26c=0.88695838,
                                                                      ymax31a=0.89351822,ymax31b=0.69375718,ymax31c=0.88715637,
                                                                      a=-16.33352131,b=0.81903181,c=31.04814934,d=0.05409619,e=44.75447753,f=0.03978480),
                                  data=list(T21a=ROPS21a$Temperature,T21b=ROPS21b$Temperature,T21c=ROPS21c$Temperature,
                                            T26a=ROPS26a$Temperature,T26b=ROPS26b$Temperature,T26c=ROPS26c$Temperature,
                                            T31a=ROPS31a$Temperature,T31b=ROPS31b$Temperature,T31c=ROPS31c$Temperature,
                                            Nasedat21a=ROPS21a$Vmax/1922.245351,Nasedat21b=ROPS21b$Vmax/1690.867422,Nasedat21c=ROPS21c$Vmax/450.8292229,
                                            Nasedat26a=ROPS26a$Vmax/2785.452785,Nasedat26b=ROPS26b$Vmax/1572.439251,Nasedat26c=ROPS26c$Vmax/1686.352,
                                            Nasedat31a=ROPS31a$Vmax/1646.769302,Nasedat31b=ROPS31b$Vmax/920.7416901,Nasedat31c=ROPS31c$Vmax/1978.187071),
                                  control=list(maxit=50000))
summary(fit_Nase_beta_linall_ROPS)

###
#Photosynthesis (Asat)
###

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

####
#Save N-fixation parameter estimates
####

#Morella
ft.m <- coef(fit_Nase_beta_linall_MOCE)

#Alnus
ft.a <- coef(fit_Nase_beta_linall_ALRU)

#Gliricidia
ft.g <- coef(fit_Nase_beta_linTminTopt_GLSE)

#Robinia
ft.r <- coef(fit_Nase_beta_linall_ROPS)

####
#Calculate N-fixation rates at 15 and 40 deg. C as a function of growing temperature
####

#Vector of mean growing and measurement temperatures
Tg.seq<-seq(15,35,0.1) #Growing temperature
Tsim<-seq(10,40,0.1) #Measurement temperature

###
#Rates at 15 deg. C (Supplementary Figure 4a)
###

#Morella
M15<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  M15[i]<-beta(1,ft.m[12]*Tg.seq[i]+ft.m[11],ft.m[14]*Tg.seq[i]+ft.m[13],ft.m[16]*Tg.seq[i]+ft.m[15],15)
}
#At high growing temperatures Tmin is greater than 15 deg. C causing NA values to appear
#This makes those NA values equal to zero
M15[is.nan(M15)]<-0 

#Alnus
A15<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  A15[i]<-beta(1,ft.a[12]*Tg.seq[i]+ft.a[11],ft.a[14]*Tg.seq[i]+ft.a[13],ft.a[16]*Tg.seq[i]+ft.a[15],15)
}

#Gliricidia
G15<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  G15[i]<-beta(1,ft.g[13]*Tg.seq[i]+ft.g[12],ft.g[15]*Tg.seq[i]+ft.g[14],ft.g[11],15)
}

#Robinia
R15<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  R15[i]<-beta(1,ft.r[12]*Tg.seq[i]+ft.r[11],ft.r[14]*Tg.seq[i]+ft.r[13],ft.r[16]*Tg.seq[i]+ft.r[15],15)
}

###
#Rates at 40 deg. C (Supplementary Figure 4c)
###

#Morella
M40<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  M40[i]<-beta(1,ft.m[12]*Tg.seq[i]+ft.m[11],ft.m[14]*Tg.seq[i]+ft.m[13],ft.m[16]*Tg.seq[i]+ft.m[15],40)
}

#Alnus
A40<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  A40[i]<-beta(1,ft.a[12]*Tg.seq[i]+ft.a[11],ft.a[14]*Tg.seq[i]+ft.a[13],ft.a[16]*Tg.seq[i]+ft.a[15],40)
}

#Gliricidia
G40<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  G40[i]<-beta(1,ft.g[13]*Tg.seq[i]+ft.g[12],ft.g[15]*Tg.seq[i]+ft.g[14],ft.g[11],40)
}

#Robinia
R40<-rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  R40[i]<-beta(1,ft.r[12]*Tg.seq[i]+ft.r[11],ft.r[14]*Tg.seq[i]+ft.r[13],ft.r[16]*Tg.seq[i]+ft.r[15],40)
}

####
#Calculate N-fixation 95% CI
####

##
#Morella
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_MOCE),Sigma=vcov(fit_Nase_beta_linall_MOCE))

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.m15=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m15<-rep(NA,length(Tg.seq))
high.m15<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m15[i,j]=beta.lin.all(1,vmat.m[i,11],vmat.m[i,12],vmat.m[i,13],vmat.m[i,14],vmat.m[i,15],vmat.m[i,16],15,Tg.seq[j])
  }
  low.m15[j]<-quantile(dist.m15[,j],0.025)
  high.m15[j]<-quantile(dist.m15[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.mTopt=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.mTopt<-rep(NA,length(Tg.seq))
high.mTopt<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.mTopt[i,j]=vmat.m[i,13]+vmat.m[i,14]*Tg.seq[j]
  }
  low.mTopt[j]<-quantile(dist.mTopt[,j],0.025)
  high.mTopt[j]<-quantile(dist.mTopt[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.m40=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m40<-rep(NA,length(Tg.seq))
high.m40<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m40[i,j]=beta.lin.all(1,vmat.m[i,11],vmat.m[i,12],vmat.m[i,13],vmat.m[i,14],vmat.m[i,15],vmat.m[i,16],40,Tg.seq[j])
  }
  low.m40[j]<-quantile(dist.m40[,j],0.025)
  high.m40[j]<-quantile(dist.m40[,j],0.975)
}

##
#Alnus
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.a=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_ALRU),Sigma=vcov(fit_Nase_beta_linall_ALRU))

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
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

#Topt 95% CI (Supplementary Figure 4b)
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

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
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
vmat.g=mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_GLSE),Sigma=vcov(fit_Nase_beta_linTminTopt_GLSE))

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.g15=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g15<-rep(NA,length(Tg.seq))
high.g15<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g15[i,j]=beta.lin.Tmin.Topt(1,vmat.g[i,11],vmat.g[i,12],vmat.g[i,13],vmat.g[i,14],vmat.g[i,15],15,Tg.seq[j])
  }
  low.g15[j]<-quantile(dist.g15[,j],0.025)
  high.g15[j]<-quantile(dist.g15[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.gTopt=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.gTopt<-rep(NA,length(Tg.seq))
high.gTopt<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.gTopt[i,j]=vmat.g[i,14]+vmat.g[i,15]*Tg.seq[j]
  }
  low.gTopt[j]<-quantile(dist.gTopt[,j],0.025)
  high.gTopt[j]<-quantile(dist.gTopt[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.g40=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g40<-rep(NA,length(Tg.seq))
high.g40<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g40[i,j]=beta.lin.Tmin.Topt(1,vmat.g[i,11],vmat.g[i,12],vmat.g[i,13],vmat.g[i,14],vmat.g[i,15],40,Tg.seq[j])
  }
  low.g40[j]<-quantile(dist.g40[,j],0.025)
  high.g40[j]<-quantile(dist.g40[,j],0.975)
}

##
#Robinia
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_ROPS),Sigma=vcov(fit_Nase_beta_linall_ROPS))

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.r15=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r15<-rep(NA,length(Tg.seq))
high.r15<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r15[i,j]=beta.lin.all(1,vmat.r[i,11],vmat.r[i,12],vmat.r[i,13],vmat.r[i,14],vmat.r[i,15],vmat.r[i,16],15,Tg.seq[j])
  }
  low.r15[j]<-quantile(dist.r15[,j],0.025)
  high.r15[j]<-quantile(dist.r15[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.rTopt=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.rTopt<-rep(NA,length(Tg.seq))
high.rTopt<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.rTopt[i,j]=vmat.r[i,13]+vmat.r[i,14]*Tg.seq[j]
  }
  low.rTopt[j]<-quantile(dist.rTopt[,j],0.025)
  high.rTopt[j]<-quantile(dist.rTopt[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.r40=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r40<-rep(NA,length(Tg.seq))
high.r40<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r40[i,j]=beta.lin.all(1,vmat.r[i,11],vmat.r[i,12],vmat.r[i,13],vmat.r[i,14],vmat.r[i,15],vmat.r[i,16],40,Tg.seq[j])
  }
  low.r40[j]<-quantile(dist.r40[,j],0.025)
  high.r40[j]<-quantile(dist.r40[,j],0.975)
}

####
#Calculate photosynthesis rates at 15 and 40 deg. C as a function of growing temperature
####

###
#Rates at 15 deg. C (Supplementary Figure 4a)
###

#Morella
m15.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  m15.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[85],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[89],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[77],15,Tg.seq[i])
}

#Alnus
a15.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  a15.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[86],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[90],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[78],15,Tg.seq[i])
}

#Gliricidia
g15.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  g15.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[87],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[91],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[79],15,Tg.seq[i])
}

#Robinia
r15.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  r15.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[88],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[92],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[80],15,Tg.seq[i])
}

###
#Rates at 40 deg. C (Supplementary Figure 4c)
###

#Morella
m40.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  m40.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[85],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[89],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[77],40,Tg.seq[i])
}

#Alnus
a40.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  a40.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[86],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[90],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[78],40,Tg.seq[i])
}

#Gliricidia
g40.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  g40.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[87],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[91],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[79],40,Tg.seq[i])
}

#Robinia
r40.Asat=rep(NA,length(Tg.seq))
for(i in 1:length(Tg.seq)){
  r40.Asat[i]=beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[88],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[92],
                            coef(fit_Photo_beta_all_Topt.lin_Asat)[80],40,Tg.seq[i])
}

####
#Calculate photosynthesis 95% CI
####

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.Asat=mvrnorm(1000,mu=coef(fit_Photo_beta_all_Topt.lin_Asat),Sigma=vcov(fit_Photo_beta_all_Topt.lin_Asat))

###
#Morella
###

#Asat ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.m21.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m21.Asat<-rep(NA,length(Tsim))
high.m21.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m21.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,81],vmat.Asat[i,85],vmat.Asat[i,89],vmat.Asat[i,77],Tsim[j],18.5)
  }
  low.m21.Asat[j]<-quantile(na.omit(dist.m21.Asat[,j]),0.025)
  high.m21.Asat[j]<-quantile(na.omit(dist.m21.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.m26.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m26.Asat<-rep(NA,length(Tsim))
high.m26.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m26.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,81],vmat.Asat[i,85],vmat.Asat[i,89],vmat.Asat[i,77],Tsim[j],23.5)
  }
  low.m26.Asat[j]<-quantile(na.omit(dist.m26.Asat[,j]),0.025)
  high.m26.Asat[j]<-quantile(na.omit(dist.m26.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.m31.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m31.Asat<-rep(NA,length(Tsim))
high.m31.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m31.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,81],vmat.Asat[i,85],vmat.Asat[i,89],vmat.Asat[i,77],Tsim[j],28.5)
  }
  low.m31.Asat[j]<-quantile(na.omit(dist.m31.Asat[,j]),0.025)
  high.m31.Asat[j]<-quantile(na.omit(dist.m31.Asat[,j]),0.975)
}

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.m15.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m15.Asat<-rep(NA,length(Tg.seq))
high.m15.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m15.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,81],vmat.Asat[i,85],vmat.Asat[i,89],vmat.Asat[i,77],15,Tg.seq[j])
  }
  low.m15.Asat[j]<-quantile(dist.m15.Asat[,j],0.025)
  high.m15.Asat[j]<-quantile(dist.m15.Asat[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.mTopt.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.mTopt.Asat<-rep(NA,length(Tg.seq))
high.mTopt.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.mTopt.Asat[i,j]=vmat.Asat[i,85]+vmat.Asat[i,89]*Tg.seq[j]
  }
  low.mTopt.Asat[j]<-quantile(dist.mTopt.Asat[,j],0.025)
  high.mTopt.Asat[j]<-quantile(dist.mTopt.Asat[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.m40.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.m40.Asat<-rep(NA,length(Tg.seq))
high.m40.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.m40.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,81],vmat.Asat[i,85],vmat.Asat[i,89],vmat.Asat[i,77],40,Tg.seq[j])
  }
  low.m40.Asat[j]<-quantile(dist.m40.Asat[,j],0.025)
  high.m40.Asat[j]<-quantile(dist.m40.Asat[,j],0.975)
}

###
#Alnus
###

#Asat ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.a21.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a21.Asat<-rep(NA,length(Tsim))
high.a21.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a21.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,82],vmat.Asat[i,86],vmat.Asat[i,90],vmat.Asat[i,78],Tsim[j],18.5)
  }
  low.a21.Asat[j]<-quantile(na.omit(dist.a21.Asat[,j]),0.025)
  high.a21.Asat[j]<-quantile(na.omit(dist.a21.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.a26.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a26.Asat<-rep(NA,length(Tsim))
high.a26.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a26.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,82],vmat.Asat[i,86],vmat.Asat[i,90],vmat.Asat[i,78],Tsim[j],23.5)
  }
  low.a26.Asat[j]<-quantile(na.omit(dist.a26.Asat[,j]),0.025)
  high.a26.Asat[j]<-quantile(na.omit(dist.a26.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.a31.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a31.Asat<-rep(NA,length(Tsim))
high.a31.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a31.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,82],vmat.Asat[i,86],vmat.Asat[i,90],vmat.Asat[i,78],Tsim[j],28.5)
  }
  low.a31.Asat[j]<-quantile(na.omit(dist.a31.Asat[,j]),0.025)
  high.a31.Asat[j]<-quantile(na.omit(dist.a31.Asat[,j]),0.975)
}

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.a15.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.a15.Asat<-rep(NA,length(Tg.seq))
high.a15.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.a15.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,82],vmat.Asat[i,86],vmat.Asat[i,90],vmat.Asat[i,78],15,Tg.seq[j])
  }
  low.a15.Asat[j]<-quantile(dist.a15.Asat[,j],0.025)
  high.a15.Asat[j]<-quantile(dist.a15.Asat[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.aTopt.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.aTopt.Asat<-rep(NA,length(Tg.seq))
high.aTopt.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.aTopt.Asat[i,j]=vmat.Asat[i,86]+vmat.Asat[i,90]*Tg.seq[j]
  }
  low.aTopt.Asat[j]<-quantile(dist.aTopt.Asat[,j],0.025)
  high.aTopt.Asat[j]<-quantile(dist.aTopt.Asat[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.a40.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.a40.Asat<-rep(NA,length(Tg.seq))
high.a40.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.a40.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,82],vmat.Asat[i,86],vmat.Asat[i,90],vmat.Asat[i,78],40,Tg.seq[j])
  }
  low.a40.Asat[j]<-quantile(dist.a40.Asat[,j],0.025)
  high.a40.Asat[j]<-quantile(dist.a40.Asat[,j],0.975)
}

###
#Gliricidia
###

#Asat ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.g21.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g21.Asat<-rep(NA,length(Tsim))
high.g21.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g21.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,83],vmat.Asat[i,87],vmat.Asat[i,91],vmat.Asat[i,79],Tsim[j],18.5)
  }
  low.g21.Asat[j]<-quantile(na.omit(dist.g21.Asat[,j]),0.025)
  high.g21.Asat[j]<-quantile(na.omit(dist.g21.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.g26.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g26.Asat<-rep(NA,length(Tsim))
high.g26.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g26.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,83],vmat.Asat[i,87],vmat.Asat[i,91],vmat.Asat[i,79],Tsim[j],23.5)
  }
  low.g26.Asat[j]<-quantile(na.omit(dist.g26.Asat[,j]),0.025)
  high.g26.Asat[j]<-quantile(na.omit(dist.g26.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.g31.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g31.Asat<-rep(NA,length(Tsim))
high.g31.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g31.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,83],vmat.Asat[i,87],vmat.Asat[i,91],vmat.Asat[i,79],Tsim[j],28.5)
  }
  low.g31.Asat[j]<-quantile(na.omit(dist.g31.Asat[,j]),0.025)
  high.g31.Asat[j]<-quantile(na.omit(dist.g31.Asat[,j]),0.975)
}

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.g15.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g15.Asat<-rep(NA,length(Tg.seq))
high.g15.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g15.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,83],vmat.Asat[i,87],vmat.Asat[i,91],vmat.Asat[i,79],15,Tg.seq[j])
  }
  low.g15.Asat[j]<-quantile(dist.g15.Asat[,j],0.025)
  high.g15.Asat[j]<-quantile(dist.g15.Asat[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.gTopt.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.gTopt.Asat<-rep(NA,length(Tg.seq))
high.gTopt.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.gTopt.Asat[i,j]=vmat.Asat[i,87]+vmat.Asat[i,91]*Tg.seq[j]
  }
  low.gTopt.Asat[j]<-quantile(dist.gTopt.Asat[,j],0.025)
  high.gTopt.Asat[j]<-quantile(dist.gTopt.Asat[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.g40.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.g40.Asat<-rep(NA,length(Tg.seq))
high.g40.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.g40.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,83],vmat.Asat[i,87],vmat.Asat[i,91],vmat.Asat[i,79],40,Tg.seq[j])
  }
  low.g40.Asat[j]<-quantile(dist.g40.Asat[,j],0.025)
  high.g40.Asat[j]<-quantile(dist.g40.Asat[,j],0.975)
}

###
#Robinia
###

#Asat ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.r21.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r21.Asat<-rep(NA,length(Tsim))
high.r21.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r21.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,84],vmat.Asat[i,88],vmat.Asat[i,92],vmat.Asat[i,80],Tsim[j],18.5)
  }
  low.r21.Asat[j]<-quantile(na.omit(dist.r21.Asat[,j]),0.025)
  high.r21.Asat[j]<-quantile(na.omit(dist.r21.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.r26.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r26.Asat<-rep(NA,length(Tsim))
high.r26.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r26.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,84],vmat.Asat[i,88],vmat.Asat[i,92],vmat.Asat[i,80],Tsim[j],23.5)
  }
  low.r26.Asat[j]<-quantile(na.omit(dist.r26.Asat[,j]),0.025)
  high.r26.Asat[j]<-quantile(na.omit(dist.r26.Asat[,j]),0.975)
}

#Asat ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.r31.Asat=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r31.Asat<-rep(NA,length(Tsim))
high.r31.Asat<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r31.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,84],vmat.Asat[i,88],vmat.Asat[i,92],vmat.Asat[i,80],Tsim[j],28.5)
  }
  low.r31.Asat[j]<-quantile(na.omit(dist.r31.Asat[,j]),0.025)
  high.r31.Asat[j]<-quantile(na.omit(dist.r31.Asat[,j]),0.975)
}

#Rates at 15 deg. C 95% CI (Supplementary Figure 4a)
dist.r15.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r15.Asat<-rep(NA,length(Tg.seq))
high.r15.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r15.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,84],vmat.Asat[i,88],vmat.Asat[i,92],vmat.Asat[i,80],15,Tg.seq[j])
  }
  low.r15.Asat[j]<-quantile(dist.r15.Asat[,j],0.025)
  high.r15.Asat[j]<-quantile(dist.r15.Asat[,j],0.975)
}

#Topt 95% CI (Supplementary Figure 4b)
dist.rTopt.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.rTopt.Asat<-rep(NA,length(Tg.seq))
high.rTopt.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.rTopt.Asat[i,j]=vmat.Asat[i,88]+vmat.Asat[i,92]*Tg.seq[j]
  }
  low.rTopt.Asat[j]<-quantile(dist.rTopt.Asat[,j],0.025)
  high.rTopt.Asat[j]<-quantile(dist.rTopt.Asat[,j],0.975)
}

#Rates at 40 deg. C 95% CI (Supplementary Figure 4c)
dist.r40.Asat=matrix(NA,nrow=1000,ncol=length(Tg.seq))
low.r40.Asat<-rep(NA,length(Tg.seq))
high.r40.Asat<-rep(NA,length(Tg.seq))
for(j in 1:length(Tg.seq)){
  for(i in 1:1000){
    dist.r40.Asat[i,j]=beta.Topt.lin(1,vmat.Asat[i,84],vmat.Asat[i,88],vmat.Asat[i,92],vmat.Asat[i,80],40,Tg.seq[j])
  }
  low.r40.Asat[j]<-quantile(dist.r40.Asat[,j],0.025)
  high.r40.Asat[j]<-quantile(dist.r40.Asat[,j],0.975)
}

####
#Calculate relative photosynthesis rate at 25 deg. C
####

Asat.A21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],25,18.5)
Asat.M21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],25,18.5)
Asat.G21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],25,18.5)
Asat.R21.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],25,18.5)
Asat.A26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],25,23.5)
Asat.M26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],25,23.5)
Asat.G26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],25,23.5)
Asat.R26.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],25,23.5)
Asat.A31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],25,28.5)
Asat.M31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],25,28.5)
Asat.G31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],25,28.5)
Asat.R31.25<-beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],25,28.5)

###############################################################################################################
#Supplementary Figure 1
###############################################################################################################

#PDF dimension is 8x8 inches  

#Plotting Settings
dev.off() #Use default margins
par(pty="s")

plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(0,50),ylim=c(0,1.4),cex.lab=1.5,cex.axis=1.2)
curve(beta.Topt.lin(1/Asat.M21.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],x,18.5),from=10,to=40,lty=1,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.A21.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],x,18.5),from=10,to=40,lty=2,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.G21.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],x,18.5),from=10,to=40,lty=3,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.R21.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],x,18.5),from=10,to=40,lty=4,col="dodgerblue1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.M26.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],x,23.5),from=10,to=40,lty=1,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.A26.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],x,23.5),from=10,to=40,lty=2,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.G26.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],x,23.5),from=10,to=40,lty=3,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.R26.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],x,23.5),from=10,to=40,lty=4,col="gold1",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.M31.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],x,28.5),from=10,to=40,lty=1,col="orangered3",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.A31.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],x,28.5),from=10,to=40,lty=2,col="orangered3",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.G31.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],x,28.5),from=10,to=40,lty=3,col="orangered3",add=T,lwd=3)
curve(beta.Topt.lin(1/Asat.R31.25,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],x,28.5),from=10,to=40,lty=4,col="orangered3",add=T,lwd=3)
mtext(expression(italic('A')[sat]*' (normalized to 1 at 25 '*degree*'C)'),side=2,line=3,cex=1.5)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=F)
legend(-2,1.48,c(expression(underline(bolditalic("Morella"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                 expression(underline(bolditalic("Alnus"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                 expression(underline(bolditalic("Gliricidia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                 expression(underline(bolditalic("Robinia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C')),
       col=c(NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3"),
       lty=c(NA,1,1,1,NA,2,2,2,NA,3,3,3,NA,4,4,4),bty="n",lwd=3,y.intersp = 0.8,cex=1,seg.len=2,x.intersp = 0.5)

###############################################################################################################
#Supplementary Figure 3
###############################################################################################################

#PDF dimension is 9x7 inches  

#Plotting Settings
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))
par(pty="s")

#S3a
plot(M21.Asat.dat$Temp[1:3],M21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[2],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m21.Asat,rev(high.m21.Asat)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(M21.Asat.dat$Temp[1:3],M21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[2],pch=1,cex=1.5,col="black")
points(M21.Asat.dat$Temp[4:8],M21.Asat.dat$A400[4:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[3],pch=2,cex=1.5,col="black")
points(M21.Asat.dat$Temp[9:11],M21.Asat.dat$A400[9:11]/coef(fit_Photo_beta_all_Topt.lin_Asat)[4],pch=3,cex=1.5,col="black")
points(M21.Asat.dat$Temp[12:16],M21.Asat.dat$A400[12:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[5],pch=4,cex=1.5,col="black")
points(M21.Asat.dat$Temp[17:19],M21.Asat.dat$A400[17:19]/coef(fit_Photo_beta_all_Topt.lin_Asat)[6],pch=5,cex=1.5,col="black")
points(M21.Asat.dat$Temp[20:24],M21.Asat.dat$A400[20:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[7],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)

#S3b
plot(A21.Asat.dat$Temp[1:3],A21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[20],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a21.Asat,rev(high.a21.Asat)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(A21.Asat.dat$Temp[1:3],A21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[20],pch=1,cex=1.5,col="black")
points(A21.Asat.dat$Temp[4:8],A21.Asat.dat$A400[4:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[21],pch=2,cex=1.5,col="black")
points(A21.Asat.dat$Temp[9:11],A21.Asat.dat$A400[9:11]/coef(fit_Photo_beta_all_Topt.lin_Asat)[22],pch=3,cex=1.5,col="black")
points(A21.Asat.dat$Temp[12:16],A21.Asat.dat$A400[12:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[23],pch=4,cex=1.5,col="black")
points(A21.Asat.dat$Temp[17:19],A21.Asat.dat$A400[17:19]/coef(fit_Photo_beta_all_Topt.lin_Asat)[24],pch=5,cex=1.5,col="black")
points(A21.Asat.dat$Temp[20:24],A21.Asat.dat$A400[20:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[25],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)

#S3c
plot(G21.Asat.dat$Temp[1:3],G21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[38],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g21.Asat,rev(high.g21.Asat)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(G21.Asat.dat$Temp[1:3],G21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[38],pch=1,cex=1.5,col="black")
points(G21.Asat.dat$Temp[4:8],G21.Asat.dat$A400[4:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[39],pch=2,cex=1.5,col="black")
points(G21.Asat.dat$Temp[9:11],G21.Asat.dat$A400[9:11]/coef(fit_Photo_beta_all_Topt.lin_Asat)[40],pch=3,cex=1.5,col="black")
points(G21.Asat.dat$Temp[12:16],G21.Asat.dat$A400[12:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[41],pch=4,cex=1.5,col="black")
points(G21.Asat.dat$Temp[17:19],G21.Asat.dat$A400[17:19]/coef(fit_Photo_beta_all_Topt.lin_Asat)[42],pch=5,cex=1.5,col="black")
points(G21.Asat.dat$Temp[20:24],G21.Asat.dat$A400[20:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[43],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)

#S3d
plot(R21.Asat.dat$Temp[1:3],R21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[57],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r21.Asat,rev(high.r21.Asat)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(R21.Asat.dat$Temp[1:3],R21.Asat.dat$A400[1:3]/coef(fit_Photo_beta_all_Topt.lin_Asat)[57],pch=1,cex=1.5,col="black")
points(R21.Asat.dat$Temp[4:8],R21.Asat.dat$A400[4:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[58],pch=2,cex=1.5,col="black")
points(R21.Asat.dat$Temp[9:11],R21.Asat.dat$A400[9:11]/coef(fit_Photo_beta_all_Topt.lin_Asat)[59],pch=3,cex=1.5,col="black")
points(R21.Asat.dat$Temp[12:16],R21.Asat.dat$A400[12:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[60],pch=4,cex=1.5,col="black")
points(R21.Asat.dat$Temp[17:19],R21.Asat.dat$A400[17:19]/coef(fit_Photo_beta_all_Topt.lin_Asat)[61],pch=5,cex=1.5,col="black")
points(R21.Asat.dat$Temp[20:24],R21.Asat.dat$A400[20:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[62],pch=6,cex=1.5,col="black")
points(R21.Asat.dat$Temp[25:27],R21.Asat.dat$A400[25:27]/coef(fit_Photo_beta_all_Topt.lin_Asat)[63],pch=7,cex=1.5,col="black")
points(R21.Asat.dat$Temp[28:32],R21.Asat.dat$A400[28:32]/coef(fit_Photo_beta_all_Topt.lin_Asat)[64],pch=8,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)

#S3e
plot(M26.Asat.dat$Temp[1:4],M26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[8],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m26.Asat,rev(high.m26.Asat)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(M26.Asat.dat$Temp[1:4],M26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[8],pch=1,cex=1.5,col="black")
points(M26.Asat.dat$Temp[5:8],M26.Asat.dat$A400[5:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[9],pch=2,cex=1.5,col="black")
points(M26.Asat.dat$Temp[9:12],M26.Asat.dat$A400[9:12]/coef(fit_Photo_beta_all_Topt.lin_Asat)[10],pch=3,cex=1.5,col="black")
points(M26.Asat.dat$Temp[13:16],M26.Asat.dat$A400[13:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[11],pch=4,cex=1.5,col="black")
points(M26.Asat.dat$Temp[17:20],M26.Asat.dat$A400[17:20]/coef(fit_Photo_beta_all_Topt.lin_Asat)[12],pch=5,cex=1.5,col="black")
points(M26.Asat.dat$Temp[21:24],M26.Asat.dat$A400[21:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[13],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S3f
plot(A26.Asat.dat$Temp[1:4],A26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[26],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a26.Asat,rev(high.a26.Asat)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(A26.Asat.dat$Temp[1:4],A26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[26],pch=1,cex=1.5,col="black")
points(A26.Asat.dat$Temp[5:8],A26.Asat.dat$A400[5:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[27],pch=2,cex=1.5,col="black")
points(A26.Asat.dat$Temp[9:12],A26.Asat.dat$A400[9:12]/coef(fit_Photo_beta_all_Topt.lin_Asat)[28],pch=3,cex=1.5,col="black")
points(A26.Asat.dat$Temp[13:16],A26.Asat.dat$A400[13:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[29],pch=4,cex=1.5,col="black")
points(A26.Asat.dat$Temp[17:20],A26.Asat.dat$A400[17:20]/coef(fit_Photo_beta_all_Topt.lin_Asat)[30],pch=5,cex=1.5,col="black")
points(A26.Asat.dat$Temp[21:24],A26.Asat.dat$A400[21:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[31],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S3g
plot(G26.Asat.dat$Temp[1:4],G26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[44],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g26.Asat,rev(high.g26.Asat)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(G26.Asat.dat$Temp[1:4],G26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[44],pch=1,cex=1.5,col="black")
points(G26.Asat.dat$Temp[5:8],G26.Asat.dat$A400[5:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[45],pch=2,cex=1.5,col="black")
points(G26.Asat.dat$Temp[9:12],G26.Asat.dat$A400[9:12]/coef(fit_Photo_beta_all_Topt.lin_Asat)[46],pch=3,cex=1.5,col="black")
points(G26.Asat.dat$Temp[13:16],G26.Asat.dat$A400[13:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[47],pch=4,cex=1.5,col="black")
points(G26.Asat.dat$Temp[17:20],G26.Asat.dat$A400[17:20]/coef(fit_Photo_beta_all_Topt.lin_Asat)[48],pch=5,cex=1.5,col="black")
points(G26.Asat.dat$Temp[21:24],G26.Asat.dat$A400[21:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[49],pch=6,cex=1.5,col="black")
points(G26.Asat.dat$Temp[25:28],G26.Asat.dat$A400[25:28]/coef(fit_Photo_beta_all_Topt.lin_Asat)[50],pch=7,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S3h
plot(R26.Asat.dat$Temp[1:4],R26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[65],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r26.Asat,rev(high.r26.Asat)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(R26.Asat.dat$Temp[1:4],R26.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[65],pch=1,cex=1.5,col="black")
points(R26.Asat.dat$Temp[5:8],R26.Asat.dat$A400[5:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[66],pch=2,cex=1.5,col="black")
points(R26.Asat.dat$Temp[9:13],R26.Asat.dat$A400[9:13]/coef(fit_Photo_beta_all_Topt.lin_Asat)[67],pch=3,cex=1.5,col="black")
points(R26.Asat.dat$Temp[14:17],R26.Asat.dat$A400[14:17]/coef(fit_Photo_beta_all_Topt.lin_Asat)[68],pch=4,cex=1.5,col="black")
points(R26.Asat.dat$Temp[18:21],R26.Asat.dat$A400[18:21]/coef(fit_Photo_beta_all_Topt.lin_Asat)[69],pch=5,cex=1.5,col="black")
points(R26.Asat.dat$Temp[22:25],R26.Asat.dat$A400[22:25]/coef(fit_Photo_beta_all_Topt.lin_Asat)[70],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

#S3i
plot(M31.Asat.dat$Temp[1:5],M31.Asat.dat$A400[1:5]/coef(fit_Photo_beta_all_Topt.lin_Asat)[14],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m31.Asat,rev(high.m31.Asat)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(M31.Asat.dat$Temp[1:5],M31.Asat.dat$A400[1:5]/coef(fit_Photo_beta_all_Topt.lin_Asat)[14],pch=1,cex=1.5,col="black")
points(M31.Asat.dat$Temp[6:8],M31.Asat.dat$A400[6:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[15],pch=2,cex=1.5,col="black")
points(M31.Asat.dat$Temp[9:12],M31.Asat.dat$A400[9:12]/coef(fit_Photo_beta_all_Topt.lin_Asat)[16],pch=3,cex=1.5,col="black")
points(M31.Asat.dat$Temp[13:15],M31.Asat.dat$A400[13:15]/coef(fit_Photo_beta_all_Topt.lin_Asat)[17],pch=4,cex=1.5,col="black")
points(M31.Asat.dat$Temp[16:20],M31.Asat.dat$A400[16:20]/coef(fit_Photo_beta_all_Topt.lin_Asat)[18],pch=5,cex=1.5,col="black")
points(M31.Asat.dat$Temp[21:23],M31.Asat.dat$A400[21:23]/coef(fit_Photo_beta_all_Topt.lin_Asat)[19],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[81],coef(fit_Photo_beta_all_Topt.lin_Asat)[85],coef(fit_Photo_beta_all_Topt.lin_Asat)[89],coef(fit_Photo_beta_all_Topt.lin_Asat)[77],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S3j
plot(A31.Asat.dat$Temp[1:5],A31.Asat.dat$A400[1:5]/coef(fit_Photo_beta_all_Topt.lin_Asat)[32],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a31.Asat,rev(high.a31.Asat)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(A31.Asat.dat$Temp[1:5],A31.Asat.dat$A400[1:5]/coef(fit_Photo_beta_all_Topt.lin_Asat)[32],pch=1,cex=1.5,col="black")
points(A31.Asat.dat$Temp[6:8],A31.Asat.dat$A400[6:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[33],pch=2,cex=1.5,col="black")
points(A31.Asat.dat$Temp[9:13],A31.Asat.dat$A400[9:13]/coef(fit_Photo_beta_all_Topt.lin_Asat)[34],pch=3,cex=1.5,col="black")
points(A31.Asat.dat$Temp[14:16],A31.Asat.dat$A400[14:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[35],pch=4,cex=1.5,col="black")
points(A31.Asat.dat$Temp[17:21],A31.Asat.dat$A400[17:21]/coef(fit_Photo_beta_all_Topt.lin_Asat)[36],pch=5,cex=1.5,col="black")
points(A31.Asat.dat$Temp[22:24],A31.Asat.dat$A400[22:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[37],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[82],coef(fit_Photo_beta_all_Topt.lin_Asat)[86],coef(fit_Photo_beta_all_Topt.lin_Asat)[90],coef(fit_Photo_beta_all_Topt.lin_Asat)[78],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S3k
plot(G31.Asat.dat$Temp[1:4],G31.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[51],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g31.Asat,rev(high.g31.Asat)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(G31.Asat.dat$Temp[1:4],G31.Asat.dat$A400[1:4]/coef(fit_Photo_beta_all_Topt.lin_Asat)[51],pch=1,cex=1.5,col="black")
points(G31.Asat.dat$Temp[5:7],G31.Asat.dat$A400[5:7]/coef(fit_Photo_beta_all_Topt.lin_Asat)[52],pch=2,cex=1.5,col="black")
points(G31.Asat.dat$Temp[8:12],G31.Asat.dat$A400[8:12]/coef(fit_Photo_beta_all_Topt.lin_Asat)[53],pch=3,cex=1.5,col="black")
points(G31.Asat.dat$Temp[13:15],G31.Asat.dat$A400[13:15]/coef(fit_Photo_beta_all_Topt.lin_Asat)[54],pch=4,cex=1.5,col="black")
points(G31.Asat.dat$Temp[16:21],G31.Asat.dat$A400[16:21]/coef(fit_Photo_beta_all_Topt.lin_Asat)[55],pch=5,cex=1.5,col="black")
points(G31.Asat.dat$Temp[22:24],G31.Asat.dat$A400[22:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[56],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[83],coef(fit_Photo_beta_all_Topt.lin_Asat)[87],coef(fit_Photo_beta_all_Topt.lin_Asat)[91],coef(fit_Photo_beta_all_Topt.lin_Asat)[79],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S3l
plot(R31.Asat.dat$Temp[1:5],R31.Asat.dat$A400[1:5]/coef(fit_Photo_beta_all_Topt.lin_Asat)[71],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r31.Asat,rev(high.r31.Asat)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(R31.Asat.dat$Temp[1:5],R31.Asat.dat$A400[1:5]/coef(fit_Photo_beta_all_Topt.lin_Asat)[71],pch=1,cex=1.5,col="black")
points(R31.Asat.dat$Temp[6:8],R31.Asat.dat$A400[6:8]/coef(fit_Photo_beta_all_Topt.lin_Asat)[72],pch=2,cex=1.5,col="black")
points(R31.Asat.dat$Temp[9:13],R31.Asat.dat$A400[9:13]/coef(fit_Photo_beta_all_Topt.lin_Asat)[73],pch=3,cex=1.5,col="black")
points(R31.Asat.dat$Temp[14:16],R31.Asat.dat$A400[14:16]/coef(fit_Photo_beta_all_Topt.lin_Asat)[74],pch=4,cex=1.5,col="black")
points(R31.Asat.dat$Temp[17:21],R31.Asat.dat$A400[17:21]/coef(fit_Photo_beta_all_Topt.lin_Asat)[75],pch=5,cex=1.5,col="black")
points(R31.Asat.dat$Temp[22:24],R31.Asat.dat$A400[22:24]/coef(fit_Photo_beta_all_Topt.lin_Asat)[76],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_Asat)[84],coef(fit_Photo_beta_all_Topt.lin_Asat)[88],coef(fit_Photo_beta_all_Topt.lin_Asat)[92],coef(fit_Photo_beta_all_Topt.lin_Asat)[80],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression('Photosynthesis ('*italic('A')[sat]*'; % of max)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)

###############################################################################################################
#Supplementary Figure 4
###############################################################################################################

#PDF dimension is 10.5x3 inches  

#Plotting Settings
par(pty="m")
nf<-layout(matrix(c(1,2,3,4),1,4,byrow=T),c(3,3,3,1.3),c(3,3,3,3),T)
layout.show(nf)
par(oma=c(0,0,0,0))
par(mar=c(4,5,4,0))
par(pty="s")

#S4a
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(0,100),cex.lab=1.5,cex.axis=1.2)
mtext(text="a",side=3,cex=1.2,adj=0)
mtext(expression('Rate at 15 '*degree*'C'),side=2,cex=1.2,line=4.25)
mtext(expression('(% of max)'),side=2,cex=1.2,line=2.5)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m15,rev(high.m15))*100,
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a15,rev(high.a15))*100,
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g15,rev(high.g15))*100,
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r15,rev(high.r15))*100,
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m15.Asat*100,rev(high.m15.Asat*100)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a15.Asat*100,rev(high.a15.Asat*100)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g15.Asat*100,rev(high.g15.Asat*100)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r15.Asat*100,rev(high.r15.Asat*100)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
points(Tg.seq,M15*100,type="l",lty=1,col="darkorange1",lwd=1)
points(Tg.seq,A15*100,type="l",lty=1,col="darkturquoise",lwd=1)
points(Tg.seq,G15*100,type="l",lty=1,col="orangered2",lwd=1)
points(Tg.seq,R15*100,type="l",lty=1,col="dodgerblue3",lwd=1)
points(Tg.seq,m15.Asat*100,type="l",lty=2,col="darkorange1",lwd=1)
points(Tg.seq,a15.Asat*100,type="l",lty=2,col="darkturquoise",lwd=1)
points(Tg.seq,g15.Asat*100,type="l",lty=2,col="orangered2",lwd=1)
points(Tg.seq,r15.Asat*100,type="l",lty=2,col="dodgerblue3",lwd=1)
points(c(M15[which(Tg.seq==18.5)],M15[which(Tg.seq==23.5)],M15[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkorange1",cex=1.5,lwd=1)
points(c(A15[which(Tg.seq==18.5)],A15[which(Tg.seq==23.5)],A15[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(c(G15[which(Tg.seq==18.5)],G15[which(Tg.seq==23.5)],G15[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="orangered2",cex=1.5,lwd=1)
points(c(R15[which(Tg.seq==18.5)],R15[which(Tg.seq==23.5)],R15[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="dodgerblue3",cex=1.5,lwd=1)
points(c(m15.Asat[which(Tg.seq==18.5)],m15.Asat[which(Tg.seq==23.5)],m15.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkorange1",cex=1.5,lwd=1)
points(c(a15.Asat[which(Tg.seq==18.5)],a15.Asat[which(Tg.seq==23.5)],a15.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkturquoise",cex=1.5,lwd=1)
points(c(r15.Asat[which(Tg.seq==18.5)],r15.Asat[which(Tg.seq==23.5)],r15.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="dodgerblue3",cex=1.5,lwd=1)
points(c(g15.Asat[which(Tg.seq==18.5)],g15.Asat[which(Tg.seq==23.5)],g15.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="orangered2",cex=1.5,lwd=1)

#S4b
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(20,40),cex.lab=1.5,cex.axis=1.2)
mtext(text="b",side=3,cex=1.2,adj=0)
mtext(expression(italic('T')[opt]*' ('*degree*'C)'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.mTopt,rev(high.mTopt)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.aTopt,rev(high.aTopt)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.gTopt,rev(high.gTopt)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.rTopt,rev(high.rTopt)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.mTopt.Asat,rev(high.mTopt.Asat)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.aTopt.Asat,rev(high.aTopt.Asat)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.gTopt.Asat,rev(high.gTopt.Asat)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.rTopt.Asat,rev(high.rTopt.Asat)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
abline(a=14.56,b=0.783,col="darkorange1",lwd=1)
abline(a=31.20,b=0.066,col="darkturquoise",lwd=1)
abline(a=23.88,b=0.416,col="orangered2",lwd=1)
abline(a=31.44,b=0.025,col="dodgerblue3",lwd=1)
abline(a=25.99,b=0.163,col="darkorange1",lwd=1,lty=2)
abline(a=18.45,b=0.417,col="darkturquoise",lwd=1,lty=2)
abline(a=24.62,b=0.160,col="dodgerblue3",lwd=1,lty=2)
abline(a=14.26,b=0.553,col="orangered2",lwd=1,lty=2)
points(14.56+0.783*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=16,col="darkorange1",cex=1.5,lwd=1)
points(31.24+0.064*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(23.88+0.416*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=17,col="orangered2",cex=1.5,lwd=1)
points(31.44+0.025*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=17,col="dodgerblue3",cex=1.5,lwd=1)
points(25.99+0.163*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=1,col="darkorange1",cex=1.5,lwd=1)
points(18.45+0.417*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=1,col="darkturquoise",cex=1.5,lwd=1)
points(24.62+0.160*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=2,col="dodgerblue3",cex=1.5,lwd=1)
points(14.26+0.553*c(18.5,23.5,28.5)~c(18.5,23.5,28.5),pch=2,col="orangered2",cex=1.5,lwd=1)

#S4c
plot(15:40,15:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(17,30),ylim=c(0,100),cex.lab=1.5,cex.axis=1.2)
mtext(text="c",side=3,cex=1.2,adj=0)
mtext(expression('Rate at 40 '*degree*'C'),side=2,cex=1.2,line=4.25)
mtext(expression('(% of max)'),side=2,cex=1.2,line=2.5)
mtext(expression(italic('T')[growth]*' ('*degree*'C)'),side=1,line=3,cex=1.2)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m40,rev(high.m40))*100,
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a40,rev(high.a40))*100,
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g40,rev(high.g40))*100,
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r40,rev(high.r40))*100,
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.m40.Asat*100,rev(high.m40.Asat*100)),
        col=adjustcolor("darkorange1",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.a40.Asat*100,rev(high.a40.Asat*100)),
        col=adjustcolor("darkturquoise",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.g40.Asat*100,rev(high.g40.Asat*100)),
        col=adjustcolor("orangered2",alpha.f = 0.2),border=NA)
polygon(x=c(Tg.seq,rev(Tg.seq)),y=c(low.r40.Asat*100,rev(high.r40.Asat*100)),
        col=adjustcolor("dodgerblue3",alpha.f = 0.2),border=NA)
points(Tg.seq,M40*100,type="l",lty=1,col="darkorange1",lwd=1)
points(Tg.seq,A40*100,type="l",lty=1,col="darkturquoise",lwd=1)
points(Tg.seq,G40*100,type="l",lty=1,col="orangered2",lwd=1)
points(Tg.seq,R40*100,type="l",lty=1,col="dodgerblue3",lwd=1)
points(Tg.seq,m40.Asat*100,type="l",lty=2,col="darkorange1",lwd=1)
points(Tg.seq,a40.Asat*100,type="l",lty=2,col="darkturquoise",lwd=1)
points(Tg.seq,g40.Asat*100,type="l",lty=2,col="orangered2",lwd=1)
points(Tg.seq,r40.Asat*100,type="l",lty=2,col="dodgerblue3",lwd=1)
points(c(M40[which(Tg.seq==18.5)],M40[which(Tg.seq==23.5)],M40[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkorange1",cex=1.5,lwd=1)
points(c(A40[which(Tg.seq==18.5)],A40[which(Tg.seq==23.5)],A40[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=16,col="darkturquoise",cex=1.5,lwd=1)
points(c(G40[which(Tg.seq==18.5)],G40[which(Tg.seq==23.5)],G40[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="orangered2",cex=1.5,lwd=1)
points(c(R40[which(Tg.seq==18.5)],R40[which(Tg.seq==23.5)],R40[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=17,col="dodgerblue3",cex=1.5,lwd=1)
points(c(m40.Asat[which(Tg.seq==18.5)],m40.Asat[which(Tg.seq==23.5)],m40.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkorange1",cex=1.5,lwd=1)
points(c(a40.Asat[which(Tg.seq==18.5)],a40.Asat[which(Tg.seq==23.5)],a40.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=1,col="darkturquoise",cex=1.5,lwd=1)
points(c(r40.Asat[which(Tg.seq==18.5)],r40.Asat[which(Tg.seq==23.5)],r40.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="dodgerblue3",cex=1.5,lwd=1)
points(c(g40.Asat[which(Tg.seq==18.5)],g40.Asat[which(Tg.seq==23.5)],g40.Asat[which(Tg.seq==28.5)])*100~c(18.5,23.5,28.5),pch=2,col="orangered2",cex=1.5,lwd=1)

#Change plotting Settings for legend
par(mar=c(0,0,0,0))
par(xpd=TRUE)

#Legend
plot(0:10, 0:10, type='n', bty='n', xaxt='n', yaxt='n',xlab=NA,ylab=NA)
legend("left",legend=c(expression(underline(bold('N-fixation'))),expression(italic(Morella)),expression(italic(Alnus)),expression(italic(Gliricidia)),expression(italic(Robinia)),
                       expression(underline(bold('Photosynthesis'))),expression(italic(Morella)),expression(italic(Alnus)),expression(italic(Gliricidia)),expression(italic(Robinia))),
       col=c(NA,"darkorange1","darkturquoise","orangered2","dodgerblue3",NA,"darkorange1","darkturquoise","orangered2","dodgerblue3"),
       lty=c(NA,1,1,1,1,NA,2,2,2,2),pch=c(NA,16,16,17,17,NA,1,1,2,2),bty="n",pt.cex=1.5,lwd=1,x.intersp = 0.2,y.intersp = 1,
       cex=1.1)
par(xpd=F)