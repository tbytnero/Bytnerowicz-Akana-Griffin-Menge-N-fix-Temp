###############################################################################################################
###############################################################################################################
#This script generates parameter estimates and 95% CIs for N-fixation 
#with the possibility of linear acclimation of Km (Supplementary Tables 11 and 12) 
###############################################################################################################
###############################################################################################################

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

#NLL function for acclimation of Tmin and Topt (2 species)
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

#NLL function for overall fit (all data)
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

####
#Maximum likelihood fits, parameter estimates, and 95% CI
####

###
#Supplementary Table 11 first
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
coef(fit_Nase_beta_linall_MOCE_linKm)
fit_Nase_beta_linall_MOCE_linKm.pr<-profile(fit_Nase_beta_linall_MOCE_linKm)
confint(fit_Nase_beta_linall_MOCE_linKm.pr)

#See "Nfix_Model_Comparison.R" script for Alnus fit

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
fit_Nase_beta_linTminTopt_GLSE_linKm.pr<-profile(fit_Nase_beta_linTminTopt_GLSE_linKm)
confint(fit_Nase_beta_linTminTopt_GLSE_linKm.pr)
coef(fit_Nase_beta_linTminTopt_GLSE_linKm)

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
coef(fit_Nase_beta_linTminTmax_ROPS_linKm)
fit_Nase_beta_linTminTmax_ROPS_linKm.pr<-profile(fit_Nase_beta_linTminTmax_ROPS_linKm)
confint(fit_Nase_beta_linTminTmax_ROPS_linKm.pr)

#Tropical species (Morella and Gliricidia)
fit_Nase_beta_linTminTopt_Trop_linKm <- mle2(Nase_beta_lin.Tmin.Topt_2sp_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,ymax21d=1,ymax21e=1,ymax21f=1,
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
summary(fit_Nase_beta_linTminTopt_Trop_linKm)
coef(fit_Nase_beta_linTminTopt_Trop_linKm)
fit_Nase_beta_linTminTopt_Trop_linKm.pr<-profile(fit_Nase_beta_linTminTopt_Trop_linKm)
confint(fit_Nase_beta_linTminTopt_Trop_linKm.pr)

#Temperate species (Alnus and Robinia)
fit_Nase_beta_linTminTopt_Temp_linKm <- mle2(Nase_beta_lin.Tmin.Topt_2sp_normNLL,start=list(sdNase=-1,ymax21a=1,ymax21b=1,ymax21c=1,ymax21d=1,ymax21e=1,ymax21f=1,
                                                                                            ymax26a=1,ymax26b=1,ymax26c=1,ymax26d=1,ymax26e=1,ymax26f=1,
                                                                                            ymax31a=1,ymax31b=1,ymax31c=1,ymax31d=1,ymax31e=1,ymax31f=1,
                                                                                            a=-14.4,b=0.68,c=32.1,d=0.03,Tmax=44.3),
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
summary(fit_Nase_beta_linTminTopt_Temp_linKm)
fit_Nase_beta_linTminTopt_Temp_linKm.pr<-profile(fit_Nase_beta_linTminTopt_Temp_linKm)
confint(fit_Nase_beta_linTminTopt_Temp_linKm.pr)

#Overall (all species, no acclimation)
fit_Nase_beta_4sp_bin_linKm <- mle2(Nase_beta_bin_4sp_normNLL,start=list(sdNase=-1.7620180,ymax21a=1.3089718,ymax21b=1.1895483,ymax21c=1.3177746,ymax21d=0.9484360,ymax21e=0.9302366,ymax21f=0.9302773,
                                                                         ymax21g=1.1016426,ymax21h=1.1193006,ymax21i=1.3132580,ymax21j=1.0395706,ymax21k=0.9256191,ymax21l=1.2180454,
                                                                         ymax26a=0.9672031,ymax26b=0.9183402,ymax26c=0.9542005,ymax26d=0.9524491,ymax26e=0.9530968,ymax26f=0.8678270,
                                                                         ymax26g=1.1495894,ymax26h=1.0470307,ymax26i=1.2153890,ymax26j=1.0028608,ymax26k=1.0710558,ymax26l=0.9899612,
                                                                         ymax31a=0.8597254,ymax31b=0.7769438,ymax31c=0.5401770,ymax31d=1.1062547,ymax31e=1.0000478,ymax31f=1.0135986,
                                                                         ymax31g=0.9019845,ymax31h=0.9365425,ymax31i=0.7224562,ymax31j=0.9662776,ymax31k=0.7513595,ymax31l=0.9650542,
                                                                         Tmax=45.0505816,Tmin=1.8373446,Topt=33.4400312),
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
summary(fit_Nase_beta_4sp_bin_linKm)
coef(fit_Nase_beta_4sp_bin_linKm)
fit_Nase_beta_4sp_bin_linKm.pr<-profile(fit_Nase_beta_4sp_bin_linKm)
confint(fit_Nase_beta_4sp_bin_linKm.pr)

###
#Supplementary Table 12
###

##
#Morella
##

#
#Parameter estimates
#

#Tmin
coef(fit_Nase_beta_linall_MOCE_linKm)[11]+coef(fit_Nase_beta_linall_MOCE_linKm)[12]*23.5

#Topt
coef(fit_Nase_beta_linall_MOCE_linKm)[13]+coef(fit_Nase_beta_linall_MOCE_linKm)[14]*23.5

#Tmax
coef(fit_Nase_beta_linall_MOCE_linKm)[15]+coef(fit_Nase_beta_linall_MOCE_linKm)[16]*23.5

#
#95% CI
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.m_linKm=mvrnorm(1000,mu=coef(fit_Nase_beta_linall_MOCE_linKm),Sigma=vcov(fit_Nase_beta_linall_MOCE_linKm))

#Tmin
dist.mTmin26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.mTmin26_linKm[i]=vmat.m_linKm[i,11]+vmat.m_linKm[i,12]*23.5
}
quantile(dist.mTmin26_linKm,0.025)
quantile(dist.mTmin26_linKm,0.975)

#Topt
dist.mTopt26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.mTopt26_linKm[i]=vmat.m_linKm[i,13]+vmat.m_linKm[i,14]*23.5
}
quantile(dist.mTopt26_linKm,0.025)
quantile(dist.mTopt26_linKm,0.975)

#Tmax
dist.mTmax26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.mTmax26_linKm[i]=vmat.m_linKm[i,15]+vmat.m_linKm[i,16]*23.5
}
quantile(dist.mTmax26_linKm,0.025)
quantile(dist.mTmax26_linKm,0.975)

##
#See "Nfix_Model_Comparison.R" script for Alnus parameter estimates and 95% CI
##

##
#Gliricidia
##

#
#Parameter estimates
#

#Tmin
coef(fit_Nase_beta_linTminTopt_GLSE_linKm)[12]+coef(fit_Nase_beta_linTminTopt_GLSE_linKm)[13]*23.5

#Topt
coef(fit_Nase_beta_linTminTopt_GLSE_linKm)[14]+coef(fit_Nase_beta_linTminTopt_GLSE_linKm)[15]*23.5

#Tmax estimate is the same as in Supplementary Table 11

#
#95% CI
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.g_linKm=mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_GLSE_linKm),Sigma=vcov(fit_Nase_beta_linTminTopt_GLSE_linKm))

#Tmin
dist.gTmin26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.gTmin26_linKm[i]=vmat.g_linKm[i,12]+vmat.g_linKm[i,13]*23.5
}
quantile(dist.gTmin26_linKm,0.025)
quantile(dist.gTmin26_linKm,0.975)

#Topt
dist.gTopt26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.gTopt26_linKm[i]=vmat.g_linKm[i,14]+vmat.g_linKm[i,15]*23.5
}
quantile(dist.gTopt26_linKm,0.025)
quantile(dist.gTopt26_linKm,0.975)

#Tmax 95% CI is the same as in Supplementary Table 11

##
#Robinia
##

#
#Parameter estimates
#

#Tmin
coef(fit_Nase_beta_linTminTmax_ROPS_linKm)[12]+coef(fit_Nase_beta_linTminTmax_ROPS_linKm)[13]*18.5

#Topt estimate is the same as in Supplementary Table 11

#Tmax
coef(fit_Nase_beta_linTminTmax_ROPS_linKm)[14]+coef(fit_Nase_beta_linTminTmax_ROPS_linKm)[15]*18.5

##
#95% CI
##

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.r_linKm=mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTmax_ROPS_linKm),Sigma=vcov(fit_Nase_beta_linTminTmax_ROPS_linKm))

#Tmin
dist.rTmin21_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.rTmin21_linKm[i]=vmat.r_linKm[i,12]+vmat.r_linKm[i,13]*18.5
}
quantile(dist.rTmin21_linKm,0.025)
quantile(dist.rTmin21_linKm,0.975)

#Topt 95% CI is the same as in Supplementary Table 11

#Tmax
dist.rTmax21_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.rTmax21_linKm[i]=vmat.r_linKm[i,14]+vmat.r_linKm[i,15]*18.5
}
quantile(dist.rTmax21_linKm,0.025)
quantile(dist.rTmax21_linKm,0.975)

##
#Tropical
##

#
#Parameter estimates
#

#Tmin
coef(fit_Nase_beta_linTminTopt_Trop_linKm)[21]+coef(fit_Nase_beta_linTminTopt_Trop_linKm)[22]*23.5

#Topt
coef(fit_Nase_beta_linTminTopt_Trop_linKm)[23]+coef(fit_Nase_beta_linTminTopt_Trop_linKm)[24]*23.5

#Tmax estimate is the same as in Supplementary Table 11

#
#95% CI
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.Trop_linKm<-mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_Trop_linKm),Sigma=vcov(fit_Nase_beta_linTminTopt_Trop_linKm))

#Tmin
dist.TropTmin26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.TropTmin26_linKm[i]=vmat.Trop_linKm[i,21]+vmat.Trop_linKm[i,22]*23.5
}
quantile(dist.TropTmin26_linKm,0.025)
quantile(dist.TropTmin26_linKm,0.975)

#Topt
dist.TropTopt26_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.TropTopt26_linKm[i]=vmat.Trop_linKm[i,23]+vmat.Trop_linKm[i,24]*23.5
}
quantile(dist.TropTopt26_linKm,0.025)
quantile(dist.TropTopt26_linKm,0.975)

#Tmax 95% CI is the same as in Supplementary Table 11

##
#Temperate
##

##
#Parameter estimates
##
#Tmin
coef(fit_Nase_beta_linall_Temp_linKm)[21]+coef(fit_Nase_beta_linall_Temp_linKm)[22]*18.5

#Topt
coef(fit_Nase_beta_linall_Temp_linKm)[23]+coef(fit_Nase_beta_linall_Temp_linKm)[24]*18.5

#Tmax estimate is the same as in Supplementary Table 11

#
#95% CI
#

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.Temp_linKm<-mvrnorm(1000,mu=coef(fit_Nase_beta_linTminTopt_Temp_linKm),Sigma=vcov(fit_Nase_beta_linTminTopt_Temp_linKm))

#Tmin
dist.TempTmin21_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.TempTmin21_linKm[i]=vmat.Temp_linKm[i,21]+vmat.Temp_linKm[i,22]*18.5
}
quantile(dist.TempTmin21_linKm,0.025)
quantile(dist.TempTmin21_linKm,0.975)

#Topt
dist.TempTopt21_linKm<-rep(NA,1000)
for(i in 1:1000){
  dist.TempTopt21_linKm[i]=vmat.Temp_linKm[i,23]+vmat.Temp_linKm[i,24]*18.5
}
quantile(dist.TempTopt21_linKm,0.025)
quantile(dist.TempTopt21_linKm,0.975)

#Tmax 95% CI is the same as in Supplementary Table 11

##
#See above (lines 493-495) for overall parameter estimates and 95% CI
##