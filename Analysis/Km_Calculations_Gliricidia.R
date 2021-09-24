###############################################################################################################
###############################################################################################################
#This script calculates the Michaelis-Menten half-saturation constant (Km) and 95% CI for Gliricidia for
#measuring nitrogenase activity with ARACAS with (Supplementary Table 4) and without a possible change of Km 
#with measurement temperature (Supplementary Table 5)
###############################################################################################################
###############################################################################################################

####
#Load Necessary Package
####

library(bbmle)

####
#Read the data and subset for each individual plant within a growing temperature
####

###
#21:15 deg. C growing temperature
###

data21<-read.csv("GLSE21_alltemps_MM.csv")

G21_p1<-data21[data21$Plant=="1",]
G21_p2<-data21[data21$Plant=="2",]
G21_p3<-data21[data21$Plant=="3",]
G21_p4<-data21[data21$Plant=="4",]
G21_p5<-data21[data21$Plant=="5",]
G21_p6<-data21[data21$Plant=="6",]

###
#26:20 deg. C growing temperature
###

data26<-read.csv("GLSE26_alltemps_MM.csv")

G26_p1<-data26[data26$Plant=="1",]
G26_p2<-data26[data26$Plant=="2",]
G26_p3<-data26[data26$Plant=="3",]
G26_p4<-data26[data26$Plant=="4",]

###
#31:25 deg. C growing temperature
###

data31<-read.csv("GLSE31_alltemps_MM.csv")

G31_p1<-data31[data31$Plant=="4",]
G31_p2<-data31[data31$Plant=="5",]
G31_p3<-data31[data31$Plant=="1",]
G31_p4<-data31[data31$Plant=="2",]
G31_p5<-data31[data31$Plant=="3",]
G31_p6<-data31[data31$Plant=="6",]

####
#Define functions
####

#Michaelis-Menten function
MM.func <- function(Vmax,Km,acet){
  y <- Vmax*acet/(Km+acet)
  y
}

#Michaelis-Menten function with linear effect of temperature on Km
MM.lin.func <- function(Vmax,acet,a,b,temp){
  y <- Vmax*acet/(a + b*temp+acet)
  y
}

####
#Negative log-likelihood (NLL) functions
####

#No effect of temperature on Km
MM_normNLL <- function(sdNase,Vmax21a,Vmax21b,Vmax21c,Vmax21d,Vmax21e,Vmax21f,
                       Vmax26a,Vmax26b,Vmax26c,Vmax26d,
                       Vmax31a,Vmax31b,Vmax31c,Vmax31d,Vmax31e,Vmax31f,
                       Km,
                       acet21a,acet21b,acet21c,acet21d,acet21e,acet21f,
                       acet26a,acet26b,acet26c,acet26d,
                       acet31a,acet31b,acet31c,acet31d,acet31e,acet31f,
                       Nasedat21a,Nasedat21b,Nasedat21c,Nasedat21d,Nasedat21e,Nasedat21f,
                       Nasedat26a,Nasedat26b,Nasedat26c,
                       Nasedat31a,Nasedat31b,Nasedat31c,Nasedat31d,Nasedat31e,Nasedat31f){
  Nasemean21a <- MM.func(Vmax21a,Km,acet21a)
  Nasemean21b <- MM.func(Vmax21b,Km,acet21b)
  Nasemean21c <- MM.func(Vmax21c,Km,acet21c)
  Nasemean21d <- MM.func(Vmax21d,Km,acet21d)
  Nasemean21e <- MM.func(Vmax21e,Km,acet21e)
  Nasemean21f <- MM.func(Vmax21f,Km,acet21f)
  Nasemean26a <- MM.func(Vmax26a,Km,acet26a)
  Nasemean26b <- MM.func(Vmax26b,Km,acet26b)
  Nasemean26c <- MM.func(Vmax26c,Km,acet26c)
  Nasemean26d <- MM.func(Vmax26d,Km,acet26d)
  Nasemean31a <- MM.func(Vmax31a,Km,acet31a)
  Nasemean31b <- MM.func(Vmax31b,Km,acet31b)
  Nasemean31c <- MM.func(Vmax31c,Km,acet31c)
  Nasemean31d <- MM.func(Vmax31d,Km,acet31d)
  Nasemean31e <- MM.func(Vmax31e,Km,acet31e)
  Nasemean31f <- MM.func(Vmax31f,Km,acet31f)
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
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31d,mean=Nasemean31d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31e,mean=Nasemean31e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31f,mean=Nasemean31f,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#With temperature effect on Km
MM_lin_NLL <- function(sdNase,Vmax21a,Vmax21b,Vmax21c,Vmax21d,Vmax21e,Vmax21f,
                       Vmax26a,Vmax26b,Vmax26c,Vmax26d,
                       Vmax31a,Vmax31b,Vmax31c,Vmax31d,Vmax31e,Vmax31f,
                       a21,a26,a31,b,
                       acet21a,acet21b,acet21c,acet21d,acet21e,acet21f,
                       acet26a,acet26b,acet26c,acet26d,
                       acet31a,acet31b,acet31c,acet31d,acet31e,acet31f,
                       temp21a,temp21b,temp21c,temp21d,temp21e,temp21f,
                       temp26a,temp26b,temp26c,temp26d,
                       temp31a,temp31b,temp31c,temp31d,temp31e,temp31f,
                       Nasedat21a,Nasedat21b,Nasedat21c,Nasedat21d,Nasedat21e,Nasedat21f,
                       Nasedat26a,Nasedat26b,Nasedat26c,
                       Nasedat31a,Nasedat31b,Nasedat31c,Nasedat31d,Nasedat31e,Nasedat31f){
  Nasemean21a <- MM.lin.func(Vmax21a,acet21a,a21,b,temp21a)
  Nasemean21b <- MM.lin.func(Vmax21b,acet21b,a21,b,temp21b)
  Nasemean21c <- MM.lin.func(Vmax21c,acet21c,a21,b,temp21c)
  Nasemean21d <- MM.lin.func(Vmax21d,acet21d,a21,b,temp21d)
  Nasemean21e <- MM.lin.func(Vmax21e,acet21e,a21,b,temp21e)
  Nasemean21f <- MM.lin.func(Vmax21f,acet21f,a21,b,temp21f)
  Nasemean26a <- MM.lin.func(Vmax26a,acet26a,a26,b,temp26a)
  Nasemean26b <- MM.lin.func(Vmax26b,acet26b,a26,b,temp26b)
  Nasemean26c <- MM.lin.func(Vmax26c,acet26c,a26,b,temp26c)
  Nasemean26d <- MM.lin.func(Vmax26d,acet26d,a26,b,temp26d)
  Nasemean31a <- MM.lin.func(Vmax31a,acet31a,a31,b,temp31a)
  Nasemean31b <- MM.lin.func(Vmax31b,acet31b,a31,b,temp31b)
  Nasemean31c <- MM.lin.func(Vmax31c,acet31c,a31,b,temp31c)
  Nasemean31d <- MM.lin.func(Vmax31d,acet31d,a31,b,temp31d)
  Nasemean31e <- MM.lin.func(Vmax31e,acet31e,a31,b,temp31e)
  Nasemean31f <- MM.lin.func(Vmax31f,acet31f,a31,b,temp31f)
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
      sum(dnorm(Nasedat31a,mean=Nasemean31a,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31b,mean=Nasemean31b,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31c,mean=Nasemean31c,sd=exp(sdNase),log=TRUE),na.rm=TRUE) +
      sum(dnorm(Nasedat31d,mean=Nasemean31d,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31e,mean=Nasemean31e,sd=exp(sdNase),log=TRUE),na.rm=TRUE) + 
      sum(dnorm(Nasedat31f,mean=Nasemean31f,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

####
#Maximum likelihood fits
####

###
#No effect of temperature on Km
###

fit_MM_GLSE <- mle2(MM_normNLL,start=list(sdNase=-1,Vmax21a=160,Vmax21b=500,Vmax21c=50,Vmax21d=140,Vmax21e=350,Vmax21f=100,
                                          Vmax26a=1000,Vmax26b=800,Vmax26c=1400,Vmax26d=2500,
                                          Vmax31a=450,Vmax31b=250,Vmax31c=1300,Vmax31d=600,Vmax31e=1000,Vmax31f=400,
                                          Km=1100),
                    data=list(acet21a=G21_p1$mL_acet_cor,acet21b=G21_p2$mL_acet_cor,acet21c=G21_p3$mL_acet_cor,
                              acet21d=G21_p4$mL_acet_cor,acet21e=G21_p5$mL_acet_cor,acet21f=G21_p6$mL_acet_cor,
                              acet26a=G26_p1$mL_acet_cor,acet26b=G26_p2$mL_acet_cor,acet26c=G26_p3$mL_acet_cor,
                              acet26d=G26_p4$mL_acet_cor,
                              acet31a=G31_p1$mL_acet_cor,acet31b=G31_p2$mL_acet_cor,acet31c=G31_p3$mL_acet_cor,
                              acet31d=G31_p4$mL_acet_cor,acet31e=G31_p5$mL_acet_cor,acet31f=G31_p6$mL_acet_cor,
                              Nasedat21a=G21_p1$dE.dt_SNF,Nasedat21b=G21_p2$dE.dt_SNF,Nasedat21c=G21_p3$dE.dt_SNF,
                              Nasedat21d=G21_p4$dE.dt_SNF,Nasedat21e=G21_p5$dE.dt_SNF,Nasedat21f=G21_p6$dE.dt_SNF,
                              Nasedat26a=G26_p1$dE.dt_SNF,Nasedat26b=G26_p2$dE.dt_SNF,Nasedat26c=G26_p3$dE.dt_SNF,
                              Nasedat26d=G26_p4$dE.dt_SNF,
                              Nasedat31a=G31_p1$dE.dt_SNF,Nasedat31b=G31_p2$dE.dt_SNF,Nasedat31c=G31_p3$dE.dt_SNF,
                              Nasedat31d=G31_p4$dE.dt_SNF,Nasedat31e=G31_p5$dE.dt_SNF,Nasedat31f=G31_p6$dE.dt_SNF),
                    control=list(maxit=20000))
summary(fit_MM_GLSE)

##
#Mean and 95% CI
##

coef(fit_MM_GLSE)[[18]]/405 #divide by 405 to convert from mL to %
confint(fit_MM_GLSE)/405 #divide by 405 to convert from mL to %

###
#With temperature effect on Km
###

fit_MM_lin_GLSE <- mle2(MM_lin_NLL,start=list(sdNase=-1,Vmax21a=200,Vmax21b=500,Vmax21c=60,Vmax21d=150,Vmax21e=350,Vmax21f=120,
                                              Vmax26a=1100,Vmax26b=800,Vmax26c=1400,Vmax26d=2500,
                                              Vmax31a=1200,Vmax31b=600,Vmax31c=1400,Vmax31d=600,Vmax31e=1000,Vmax31f=400,
                                              a21=500,a26=500,a31=500,b=13),
                        data=list(acet21a=G21_p1$mL_acet_cor,acet21b=G21_p2$mL_acet_cor,acet21c=G21_p3$mL_acet_cor,
                                  acet21d=G21_p4$mL_acet_cor,acet21e=G21_p5$mL_acet_cor,acet21f=G21_p6$mL_acet_cor,
                                  acet26a=G26_p1$mL_acet_cor,acet26b=G26_p2$mL_acet_cor,acet26c=G26_p3$mL_acet_cor,
                                  acet26d=G26_p4$mL_acet_cor,
                                  acet31a=G31_p1$mL_acet_cor,acet31b=G31_p2$mL_acet_cor,acet31c=G31_p3$mL_acet_cor,
                                  acet31d=G31_p4$mL_acet_cor,acet31e=G31_p5$mL_acet_cor,acet31f=G31_p6$mL_acet_cor,
                                  temp21a=19.9,temp21b=20.4,temp21c=20.4,temp21d=25.5,temp21e=30,temp21f=35.2,
                                  temp26a=20.6,temp26b=25.7,temp26c=29.7,temp26d=34.7,
                                  temp31a=20.4,temp31b=24.9,temp31c=31,temp31d=31,temp31e=31,temp31f=34.7,
                                  Nasedat21a=G21_p1$dE.dt_SNF,Nasedat21b=G21_p2$dE.dt_SNF,Nasedat21c=G21_p3$dE.dt_SNF,
                                  Nasedat21d=G21_p4$dE.dt_SNF,Nasedat21e=G21_p5$dE.dt_SNF,Nasedat21f=G21_p6$dE.dt_SNF,
                                  Nasedat26a=G26_p1$dE.dt_SNF,Nasedat26b=G26_p2$dE.dt_SNF,Nasedat26c=G26_p3$dE.dt_SNF,
                                  Nasedat26d=G26_p4$dE.dt_SNF,
                                  Nasedat31a=G31_p1$dE.dt_SNF,Nasedat31b=G31_p2$dE.dt_SNF,Nasedat31c=G31_p3$dE.dt_SNF,
                                  Nasedat31d=G31_p4$dE.dt_SNF,Nasedat31e=G31_p5$dE.dt_SNF,Nasedat31f=G31_p6$dE.dt_SNF),
                        control=list(maxit=20000))
summary(fit_MM_lin_GLSE)

##
#Mean and 95% CI of slope and intercept
##

#Extract slope and growing temperature specific intercepts
slope.MMlin<-coef(fit_MM_lin_GLSE)[[21]]
a21.int.MMlin<-coef(fit_MM_lin_GLSE)[[18]]
a26.int.MMlin<-coef(fit_MM_lin_GLSE)[[19]]
a31.int.MMlin<-coef(fit_MM_lin_GLSE)[[20]]

#Mean
mean(c(a21.int.MMlin,a26.int.MMlin,a31.int.MMlin))/405

#95% CI of slope
confint(fit_MM_lin_GLSE)/405

#Use parametric bootstrapping to calculate 95% CI of intercept
n=10000 #number of random draws
a21.sim<-rnorm(n,343.94,195.31)
a26.sim<-rnorm(n,781.97,221.50)
a31.sim<-rnorm(n,439.92,215.42)
a.sim<-(a21.sim+a26.sim+a31.sim)/3

quantile(a.sim,c(0.025,0.975))/405