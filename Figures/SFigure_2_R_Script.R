###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figure 2
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
#Read in photosynthesis data (A275)
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
#Define function
####

#Acclimation of Topt
beta.Topt.lin <- function(ymax,Tmin,a,b,Tmax,T,Tgrow){
  y <- ymax*(Tmax-T)/(Tmax-(a+b*Tgrow))*(((T-Tmin)/((a+b*Tgrow)-Tmin))^(((a+b*Tgrow)-Tmin)/(Tmax-(a+b*Tgrow))))
  y
}

####
#Negative log-likelihood (NLL) function
####

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
#Maximum likelihood fit
####

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

#Vector of measurement temperatures
Tsim<-seq(10,40,0.1) #Measurement temperature

####
#Calculate photosynthesis 95% CI
####

#Generate 1000 multivariate normal random deviates for calculating CI
vmat.A275=mvrnorm(1000,mu=coef(fit_Photo_beta_all_Topt.lin_A275),Sigma=vcov(fit_Photo_beta_all_Topt.lin_A275))

###
#Morella
###

#A275 ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.m21.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m21.A275<-rep(NA,length(Tsim))
high.m21.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m21.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,81],vmat.A275[i,85],vmat.A275[i,89],vmat.A275[i,77],Tsim[j],18.5)
  }
  low.m21.A275[j]<-quantile(na.omit(dist.m21.A275[,j]),0.025)
  high.m21.A275[j]<-quantile(na.omit(dist.m21.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.m26.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m26.A275<-rep(NA,length(Tsim))
high.m26.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m26.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,81],vmat.A275[i,85],vmat.A275[i,89],vmat.A275[i,77],Tsim[j],23.5)
  }
  low.m26.A275[j]<-quantile(na.omit(dist.m26.A275[,j]),0.025)
  high.m26.A275[j]<-quantile(na.omit(dist.m26.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.m31.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.m31.A275<-rep(NA,length(Tsim))
high.m31.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.m31.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,81],vmat.A275[i,85],vmat.A275[i,89],vmat.A275[i,77],Tsim[j],28.5)
  }
  low.m31.A275[j]<-quantile(na.omit(dist.m31.A275[,j]),0.025)
  high.m31.A275[j]<-quantile(na.omit(dist.m31.A275[,j]),0.975)
}

###
#Alnus
###

#A275 ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.a21.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a21.A275<-rep(NA,length(Tsim))
high.a21.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a21.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,82],vmat.A275[i,86],vmat.A275[i,90],vmat.A275[i,78],Tsim[j],18.5)
  }
  low.a21.A275[j]<-quantile(na.omit(dist.a21.A275[,j]),0.025)
  high.a21.A275[j]<-quantile(na.omit(dist.a21.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.a26.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a26.A275<-rep(NA,length(Tsim))
high.a26.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a26.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,82],vmat.A275[i,86],vmat.A275[i,90],vmat.A275[i,78],Tsim[j],23.5)
  }
  low.a26.A275[j]<-quantile(na.omit(dist.a26.A275[,j]),0.025)
  high.a26.A275[j]<-quantile(na.omit(dist.a26.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.a31.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.a31.A275<-rep(NA,length(Tsim))
high.a31.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.a31.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,82],vmat.A275[i,86],vmat.A275[i,90],vmat.A275[i,78],Tsim[j],28.5)
  }
  low.a31.A275[j]<-quantile(na.omit(dist.a31.A275[,j]),0.025)
  high.a31.A275[j]<-quantile(na.omit(dist.a31.A275[,j]),0.975)
}

###
#Gliricidia
###

#A275 ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.g21.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g21.A275<-rep(NA,length(Tsim))
high.g21.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g21.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,83],vmat.A275[i,87],vmat.A275[i,91],vmat.A275[i,79],Tsim[j],18.5)
  }
  low.g21.A275[j]<-quantile(na.omit(dist.g21.A275[,j]),0.025)
  high.g21.A275[j]<-quantile(na.omit(dist.g21.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.g26.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g26.A275<-rep(NA,length(Tsim))
high.g26.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g26.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,83],vmat.A275[i,87],vmat.A275[i,91],vmat.A275[i,79],Tsim[j],23.5)
  }
  low.g26.A275[j]<-quantile(na.omit(dist.g26.A275[,j]),0.025)
  high.g26.A275[j]<-quantile(na.omit(dist.g26.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.g31.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.g31.A275<-rep(NA,length(Tsim))
high.g31.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.g31.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,83],vmat.A275[i,87],vmat.A275[i,91],vmat.A275[i,79],Tsim[j],28.5)
  }
  low.g31.A275[j]<-quantile(na.omit(dist.g31.A275[,j]),0.025)
  high.g31.A275[j]<-quantile(na.omit(dist.g31.A275[,j]),0.975)
}

###
#Robinia
###

#A275 ~ measurement temperature for 21:15 deg. C growing temperature (Supplementary Figure 3)
dist.r21.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r21.A275<-rep(NA,length(Tsim))
high.r21.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r21.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,84],vmat.A275[i,88],vmat.A275[i,92],vmat.A275[i,80],Tsim[j],18.5)
  }
  low.r21.A275[j]<-quantile(na.omit(dist.r21.A275[,j]),0.025)
  high.r21.A275[j]<-quantile(na.omit(dist.r21.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 26:20 deg. C growing temperature (Supplementary Figure 3)
dist.r26.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r26.A275<-rep(NA,length(Tsim))
high.r26.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r26.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,84],vmat.A275[i,88],vmat.A275[i,92],vmat.A275[i,80],Tsim[j],23.5)
  }
  low.r26.A275[j]<-quantile(na.omit(dist.r26.A275[,j]),0.025)
  high.r26.A275[j]<-quantile(na.omit(dist.r26.A275[,j]),0.975)
}

#A275 ~ measurement temperature for 31:25 deg. C growing temperature (Supplementary Figure 3)
dist.r31.A275=matrix(NA,nrow=1000,ncol=length(Tsim))
low.r31.A275<-rep(NA,length(Tsim))
high.r31.A275<-rep(NA,length(Tsim))
for(j in 1:length(Tsim)){
  for(i in 1:1000){
    dist.r31.A275[i,j]=beta.Topt.lin(1,vmat.A275[i,84],vmat.A275[i,88],vmat.A275[i,92],vmat.A275[i,80],Tsim[j],28.5)
  }
  low.r31.A275[j]<-quantile(na.omit(dist.r31.A275[,j]),0.025)
  high.r31.A275[j]<-quantile(na.omit(dist.r31.A275[,j]),0.975)
}

####
#Plot Supplementary Figure 2
####

#PDF dimension is 9x7 inches  

#Plotting Settings
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))
par(pty="s")

#S2a
plot(M21.A275.dat$Temp[1:3],M21.A275.dat$A275[1:3]/coef(fit_Photo_beta_all_Topt.lin_A275)[2],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m21.A275,rev(high.m21.A275)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(M21.A275.dat$Temp[1:3],M21.A275.dat$A275[1:3]/coef(fit_Photo_beta_all_Topt.lin_A275)[2],pch=1,cex=1.5,col="black")
points(M21.A275.dat$Temp[4:8],M21.A275.dat$A275[4:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[3],pch=2,cex=1.5,col="black")
points(M21.A275.dat$Temp[9:11],M21.A275.dat$A275[9:11]/coef(fit_Photo_beta_all_Topt.lin_A275)[4],pch=3,cex=1.5,col="black")
points(M21.A275.dat$Temp[12:16],M21.A275.dat$A275[12:16]/coef(fit_Photo_beta_all_Topt.lin_A275)[5],pch=4,cex=1.5,col="black")
points(M21.A275.dat$Temp[17:19],M21.A275.dat$A275[17:19]/coef(fit_Photo_beta_all_Topt.lin_A275)[6],pch=5,cex=1.5,col="black")
points(M21.A275.dat$Temp[20:24],M21.A275.dat$A275[20:24]/coef(fit_Photo_beta_all_Topt.lin_A275)[7],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)

#S2b
plot(A21.A275.dat$Temp[1:3],A21.A275.dat$A275[1:3]/coef(fit_Photo_beta_all_Topt.lin_A275)[20],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a21.A275,rev(high.a21.A275)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(A21.A275.dat$Temp[1:3],A21.A275.dat$A275[1:3]/coef(fit_Photo_beta_all_Topt.lin_A275)[20],pch=1,cex=1.5,col="black")
points(A21.A275.dat$Temp[4:8],A21.A275.dat$A275[4:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[21],pch=2,cex=1.5,col="black")
points(A21.A275.dat$Temp[9:11],A21.A275.dat$A275[9:11]/coef(fit_Photo_beta_all_Topt.lin_A275)[22],pch=3,cex=1.5,col="black")
points(A21.A275.dat$Temp[12:16],A21.A275.dat$A275[12:16]/coef(fit_Photo_beta_all_Topt.lin_A275)[23],pch=4,cex=1.5,col="black")
points(A21.A275.dat$Temp[17:19],A21.A275.dat$A275[17:19]/coef(fit_Photo_beta_all_Topt.lin_A275)[24],pch=5,cex=1.5,col="black")
points(A21.A275.dat$Temp[20:24],A21.A275.dat$A275[20:24]/coef(fit_Photo_beta_all_Topt.lin_A275)[25],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)

#S2c
plot(G21.A275.dat$Temp[1:3],G21.A275.dat$A275[1:3]/coef(fit_Photo_beta_all_Topt.lin_A275)[38],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g21.A275,rev(high.g21.A275)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(G21.A275.dat$Temp[1:3],G21.A275.dat$A275[1:3]/coef(fit_Photo_beta_all_Topt.lin_A275)[38],pch=1,cex=1.5,col="black")
points(G21.A275.dat$Temp[4:8],G21.A275.dat$A275[4:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[39],pch=2,cex=1.5,col="black")
points(G21.A275.dat$Temp[9:11],G21.A275.dat$A275[9:11]/coef(fit_Photo_beta_all_Topt.lin_A275)[40],pch=3,cex=1.5,col="black")
points(G21.A275.dat$Temp[12:16],G21.A275.dat$A275[12:16]/coef(fit_Photo_beta_all_Topt.lin_A275)[41],pch=4,cex=1.5,col="black")
points(G21.A275.dat$Temp[17:19],G21.A275.dat$A275[17:19]/coef(fit_Photo_beta_all_Topt.lin_A275)[42],pch=5,cex=1.5,col="black")
points(G21.A275.dat$Temp[20:24],G21.A275.dat$A275[20:24]/coef(fit_Photo_beta_all_Topt.lin_A275)[43],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)

#S2d
plot(R21.A275.dat$Temp[1:2],R21.A275.dat$A275[1:2]/coef(fit_Photo_beta_all_Topt.lin_A275)[57],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r21.A275,rev(high.r21.A275)),
        col=adjustcolor("dodgerblue1",alpha.f = 0.2),border=NA)
points(R21.A275.dat$Temp[1:2],R21.A275.dat$A275[1:2]/coef(fit_Photo_beta_all_Topt.lin_A275)[57],pch=1,cex=1.5,col="black")
points(R21.A275.dat$Temp[3:6],R21.A275.dat$A275[3:6]/coef(fit_Photo_beta_all_Topt.lin_A275)[58],pch=2,cex=1.5,col="black")
points(R21.A275.dat$Temp[7:9],R21.A275.dat$A275[7:9]/coef(fit_Photo_beta_all_Topt.lin_A275)[59],pch=3,cex=1.5,col="black")
points(R21.A275.dat$Temp[10:14],R21.A275.dat$A275[10:14]/coef(fit_Photo_beta_all_Topt.lin_A275)[60],pch=4,cex=1.5,col="black")
points(R21.A275.dat$Temp[15:17],R21.A275.dat$A275[15:17]/coef(fit_Photo_beta_all_Topt.lin_A275)[61],pch=5,cex=1.5,col="black")
points(R21.A275.dat$Temp[18:22],R21.A275.dat$A275[18:22]/coef(fit_Photo_beta_all_Topt.lin_A275)[62],pch=6,cex=1.5,col="black")
points(R21.A275.dat$Temp[23:25],R21.A275.dat$A275[23:25]/coef(fit_Photo_beta_all_Topt.lin_A275)[63],pch=7,cex=1.5,col="black")
points(R21.A275.dat$Temp[26:30],R21.A275.dat$A275[26:30]/coef(fit_Photo_beta_all_Topt.lin_A275)[64],pch=8,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],x,18.5),
      from=10,to=40,lwd=2,col="dodgerblue1",lty=1,add=TRUE)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)

#S2e
plot(M26.A275.dat$Temp[1:4],M26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[8],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m26.A275,rev(high.m26.A275)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(M26.A275.dat$Temp[1:4],M26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[8],pch=1,cex=1.5,col="black")
points(M26.A275.dat$Temp[5:8],M26.A275.dat$A275[5:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[9],pch=2,cex=1.5,col="black")
points(M26.A275.dat$Temp[9:12],M26.A275.dat$A275[9:12]/coef(fit_Photo_beta_all_Topt.lin_A275)[10],pch=3,cex=1.5,col="black")
points(M26.A275.dat$Temp[13:15],M26.A275.dat$A275[13:15]/coef(fit_Photo_beta_all_Topt.lin_A275)[11],pch=4,cex=1.5,col="black")
points(M26.A275.dat$Temp[16:19],M26.A275.dat$A275[16:19]/coef(fit_Photo_beta_all_Topt.lin_A275)[12],pch=5,cex=1.5,col="black")
points(M26.A275.dat$Temp[20:23],M26.A275.dat$A275[20:23]/coef(fit_Photo_beta_all_Topt.lin_A275)[13],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S2f
plot(A26.A275.dat$Temp[1:4],A26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[26],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a26.A275,rev(high.a26.A275)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(A26.A275.dat$Temp[1:4],A26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[26],pch=1,cex=1.5,col="black")
points(A26.A275.dat$Temp[5:8],A26.A275.dat$A275[5:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[27],pch=2,cex=1.5,col="black")
points(A26.A275.dat$Temp[9:12],A26.A275.dat$A275[9:12]/coef(fit_Photo_beta_all_Topt.lin_A275)[28],pch=3,cex=1.5,col="black")
points(A26.A275.dat$Temp[13:15],A26.A275.dat$A275[13:15]/coef(fit_Photo_beta_all_Topt.lin_A275)[29],pch=4,cex=1.5,col="black")
points(A26.A275.dat$Temp[16:19],A26.A275.dat$A275[16:19]/coef(fit_Photo_beta_all_Topt.lin_A275)[30],pch=5,cex=1.5,col="black")
points(A26.A275.dat$Temp[20:22],A26.A275.dat$A275[20:22]/coef(fit_Photo_beta_all_Topt.lin_A275)[31],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S2g
plot(G26.A275.dat$Temp[1:4],G26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[44],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g26.A275,rev(high.g26.A275)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(G26.A275.dat$Temp[1:4],G26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[44],pch=1,cex=1.5,col="black")
points(G26.A275.dat$Temp[5:8],G26.A275.dat$A275[5:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[45],pch=2,cex=1.5,col="black")
points(G26.A275.dat$Temp[9:12],G26.A275.dat$A275[9:12]/coef(fit_Photo_beta_all_Topt.lin_A275)[46],pch=3,cex=1.5,col="black")
points(G26.A275.dat$Temp[13:16],G26.A275.dat$A275[13:16]/coef(fit_Photo_beta_all_Topt.lin_A275)[47],pch=4,cex=1.5,col="black")
points(G26.A275.dat$Temp[17:20],G26.A275.dat$A275[17:20]/coef(fit_Photo_beta_all_Topt.lin_A275)[48],pch=5,cex=1.5,col="black")
points(G26.A275.dat$Temp[21:24],G26.A275.dat$A275[21:24]/coef(fit_Photo_beta_all_Topt.lin_A275)[49],pch=6,cex=1.5,col="black")
points(G26.A275.dat$Temp[25:27],G26.A275.dat$A275[25:27]/coef(fit_Photo_beta_all_Topt.lin_A275)[50],pch=7,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S2h
plot(R26.A275.dat$Temp[1:4],R26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[65],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=F)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r26.A275,rev(high.r26.A275)),
        col=adjustcolor("gold1",alpha.f = 0.2),border=NA)
points(R26.A275.dat$Temp[1:4],R26.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[65],pch=1,cex=1.5,col="black")
points(R26.A275.dat$Temp[5:8],R26.A275.dat$A275[5:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[66],pch=2,cex=1.5,col="black")
points(R26.A275.dat$Temp[9:13],R26.A275.dat$A275[9:13]/coef(fit_Photo_beta_all_Topt.lin_A275)[67],pch=3,cex=1.5,col="black")
points(R26.A275.dat$Temp[14:17],R26.A275.dat$A275[14:17]/coef(fit_Photo_beta_all_Topt.lin_A275)[68],pch=4,cex=1.5,col="black")
points(R26.A275.dat$Temp[18:21],R26.A275.dat$A275[18:21]/coef(fit_Photo_beta_all_Topt.lin_A275)[69],pch=5,cex=1.5,col="black")
points(R26.A275.dat$Temp[22:25],R26.A275.dat$A275[22:25]/coef(fit_Photo_beta_all_Topt.lin_A275)[70],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],x,23.5),
      from=10,to=40,lwd=2,col="gold1",lty=1,add=TRUE)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)

#S2i
plot(M31.A275.dat$Temp[1:5],M31.A275.dat$A275[1:5]/coef(fit_Photo_beta_all_Topt.lin_A275)[14],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=c("0","25","50","75","100","125","150"),las=1,cex.axis=1.5)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.m31.A275,rev(high.m31.A275)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(M31.A275.dat$Temp[1:5],M31.A275.dat$A275[1:5]/coef(fit_Photo_beta_all_Topt.lin_A275)[14],pch=1,cex=1.5,col="black")
points(M31.A275.dat$Temp[6:8],M31.A275.dat$A275[6:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[15],pch=2,cex=1.5,col="black")
points(M31.A275.dat$Temp[9:12],M31.A275.dat$A275[9:12]/coef(fit_Photo_beta_all_Topt.lin_A275)[16],pch=3,cex=1.5,col="black")
points(M31.A275.dat$Temp[13:15],M31.A275.dat$A275[13:15]/coef(fit_Photo_beta_all_Topt.lin_A275)[17],pch=4,cex=1.5,col="black")
points(M31.A275.dat$Temp[16:20],M31.A275.dat$A275[16:20]/coef(fit_Photo_beta_all_Topt.lin_A275)[18],pch=5,cex=1.5,col="black")
points(M31.A275.dat$Temp[21:23],M31.A275.dat$A275[21:23]/coef(fit_Photo_beta_all_Topt.lin_A275)[19],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[81],coef(fit_Photo_beta_all_Topt.lin_A275)[85],coef(fit_Photo_beta_all_Topt.lin_A275)[89],coef(fit_Photo_beta_all_Topt.lin_A275)[77],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S2j
plot(A31.A275.dat$Temp[1:5],A31.A275.dat$A275[1:5]/coef(fit_Photo_beta_all_Topt.lin_A275)[32],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.a31.A275,rev(high.a31.A275)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(A31.A275.dat$Temp[1:5],A31.A275.dat$A275[1:5]/coef(fit_Photo_beta_all_Topt.lin_A275)[32],pch=1,cex=1.5,col="black")
points(A31.A275.dat$Temp[6:7],A31.A275.dat$A275[6:7]/coef(fit_Photo_beta_all_Topt.lin_A275)[33],pch=2,cex=1.5,col="black")
points(A31.A275.dat$Temp[8:12],A31.A275.dat$A275[8:12]/coef(fit_Photo_beta_all_Topt.lin_A275)[34],pch=3,cex=1.5,col="black")
points(A31.A275.dat$Temp[13:15],A31.A275.dat$A275[13:15]/coef(fit_Photo_beta_all_Topt.lin_A275)[35],pch=4,cex=1.5,col="black")
points(A31.A275.dat$Temp[16:20],A31.A275.dat$A275[16:20]/coef(fit_Photo_beta_all_Topt.lin_A275)[36],pch=5,cex=1.5,col="black")
points(A31.A275.dat$Temp[21:23],A31.A275.dat$A275[21:23]/coef(fit_Photo_beta_all_Topt.lin_A275)[37],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[82],coef(fit_Photo_beta_all_Topt.lin_A275)[86],coef(fit_Photo_beta_all_Topt.lin_A275)[90],coef(fit_Photo_beta_all_Topt.lin_A275)[78],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S2k
plot(G31.A275.dat$Temp[1:4],G31.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[51],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.g31.A275,rev(high.g31.A275)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(G31.A275.dat$Temp[1:4],G31.A275.dat$A275[1:4]/coef(fit_Photo_beta_all_Topt.lin_A275)[51],pch=1,cex=1.5,col="black")
points(G31.A275.dat$Temp[5:7],G31.A275.dat$A275[5:7]/coef(fit_Photo_beta_all_Topt.lin_A275)[52],pch=2,cex=1.5,col="black")
points(G31.A275.dat$Temp[8:12],G31.A275.dat$A275[8:12]/coef(fit_Photo_beta_all_Topt.lin_A275)[53],pch=3,cex=1.5,col="black")
points(G31.A275.dat$Temp[13:14],G31.A275.dat$A275[13:14]/coef(fit_Photo_beta_all_Topt.lin_A275)[54],pch=4,cex=1.5,col="black")
points(G31.A275.dat$Temp[15:20],G31.A275.dat$A275[15:20]/coef(fit_Photo_beta_all_Topt.lin_A275)[55],pch=5,cex=1.5,col="black")
points(G31.A275.dat$Temp[21:23],G31.A275.dat$A275[21:23]/coef(fit_Photo_beta_all_Topt.lin_A275)[56],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[83],coef(fit_Photo_beta_all_Topt.lin_A275)[87],coef(fit_Photo_beta_all_Topt.lin_A275)[91],coef(fit_Photo_beta_all_Topt.lin_A275)[79],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S2l
plot(R31.A275.dat$Temp[1:5],R31.A275.dat$A275[1:5]/coef(fit_Photo_beta_all_Topt.lin_A275)[71],pch=1,cex=1.5,col="black",ylim=c(0,1.6),xlim=c(5,45),xlab=NA,
     ylab=NA,las=1,xaxt="n",yaxt="n")
axis(1,at=c(10,20,30,40,50),labels=T,cex.axis=1.5)
axis(2,at=c(0,.25,.5,.75,1,1.25,1.5),labels=F)
polygon(x=c(Tsim,rev(Tsim)),y=c(low.r31.A275,rev(high.r31.A275)),
        col=adjustcolor("orangered3",alpha.f = 0.2),border=NA)
points(R31.A275.dat$Temp[1:5],R31.A275.dat$A275[1:5]/coef(fit_Photo_beta_all_Topt.lin_A275)[71],pch=1,cex=1.5,col="black")
points(R31.A275.dat$Temp[6:8],R31.A275.dat$A275[6:8]/coef(fit_Photo_beta_all_Topt.lin_A275)[72],pch=2,cex=1.5,col="black")
points(R31.A275.dat$Temp[9:13],R31.A275.dat$A275[9:13]/coef(fit_Photo_beta_all_Topt.lin_A275)[73],pch=3,cex=1.5,col="black")
points(R31.A275.dat$Temp[14:16],R31.A275.dat$A275[14:16]/coef(fit_Photo_beta_all_Topt.lin_A275)[74],pch=4,cex=1.5,col="black")
points(R31.A275.dat$Temp[17:21],R31.A275.dat$A275[17:21]/coef(fit_Photo_beta_all_Topt.lin_A275)[75],pch=5,cex=1.5,col="black")
points(R31.A275.dat$Temp[22:24],R31.A275.dat$A275[22:24]/coef(fit_Photo_beta_all_Topt.lin_A275)[76],pch=6,cex=1.5,col="black")
curve(beta.Topt.lin(1,coef(fit_Photo_beta_all_Topt.lin_A275)[84],coef(fit_Photo_beta_all_Topt.lin_A275)[88],coef(fit_Photo_beta_all_Topt.lin_A275)[92],coef(fit_Photo_beta_all_Topt.lin_A275)[80],x,28.5),
      from=10,to=40,lwd=2,col="orangered3",lty=1,add=TRUE)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)

mtext(expression('Photosynthesis ('*italic('A')[275]*'; % of max)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Temperature ('*degree*'C)'),side=1,line=3,cex=1.5,outer=T)