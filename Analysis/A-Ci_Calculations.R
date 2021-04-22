###############################################################################################################
###############################################################################################################
#This script fits the A-Ci curves in order to calculate A275, Vcmax, and Jmax
###############################################################################################################
###############################################################################################################

#Run Resp_Calculations.R first

#Load necessary package
library(plantecophys)

#####
#Read in data from "Photo_Temp" folder
#####

M21.071418<-read.csv("MOCE21_ACi_071418.csv")
M21.120918<-read.csv("MOCE21_ACi_120918.csv")
M21.012119<-read.csv("MOCE21_ACi_012119.csv")
M26.083019<-read.csv("MOCE26_ACi_083019.csv")
M26.091119<-read.csv("MOCE26_ACi_091119.csv")
M26.092119<-read.csv("MOCE26_ACi_092119.csv")
M31.080218<-read.csv("MOCE31_ACi_080218.csv")
M31.120418<-read.csv("MOCE31_ACi_120418.csv")
M31.040919<-read.csv("MOCE31_ACi_040919.csv")
A21.071218<-read.csv("ALRU21_ACi_071218.csv")
A21.071818<-read.csv("ALRU21_ACi_071818.csv")
A21.101718<-read.csv("ALRU21_ACi_101718.csv")
A26.112419<-read.csv("ALRU26_ACi_112419.csv")
A26.020420<-read.csv("ALRU26_ACi_020420.csv")
A26.020820<-read.csv("ALRU26_ACi_020820.csv")
A31.071618<-read.csv("ALRU31_ACi_071618.csv")
A31.092918<-read.csv("ALRU31_ACi_092918.csv")
A31.112518<-read.csv("ALRU31_ACi_112518.csv")
G21.080318<-read.csv("GLSE21_ACi_080318.csv")
G21.112718<-read.csv("GLSE21_ACi_112718.csv")
G21.082019<-read.csv("GLSE21_ACi_082019.csv")
G21.082819<-read.csv("GLSE21_ACi_082819.csv")
G26.082319<-read.csv("GLSE26_ACi_082319.csv")
G26.090619<-read.csv("GLSE26_ACi_090619.csv")
G26.100319<-read.csv("GLSE26_ACi_100319.csv")
G26.110119<-read.csv("GLSE26_ACi_110119.csv")
G31.050418<-read.csv("GLSE31_ACi_050418.csv")
G31.071118<-read.csv("GLSE31_ACi_071118.csv")
G31.072218<-read.csv("GLSE31_ACi_072218.csv")
R21.073118<-read.csv("ROPS21_ACi_073118.csv")
R21.100518<-read.csv("ROPS21_ACi_100518.csv")
R21.101018<-read.csv("ROPS21_ACi_101018.csv")
R21.091719<-read.csv("ROPS21_ACi_091719.csv")
R26.102219<-read.csv("ROPS26_ACi_102219.csv")
R26.112219<-read.csv("ROPS26_ACi_112219.csv")
R26.120319<-read.csv("ROPS26_ACi_120319.csv")
R31.071318<-read.csv("ROPS31_ACi_071318.csv")
R31.072318<-read.csv("ROPS31_ACi_072318.csv")
R31.101418<-read.csv("ROPS31_ACi_101418.csv")

####
#Morella
####

###
#21:15 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate respiration in the light (RL) for each measurement temperature
#Here and below, RL is coded as Rd in order to work with the fitacis function
p.M21.Rd<-rep(NA,length(M21.071418$Curve))
for (i in 1:length(M21.071418$Curve)){
  p.M21.Rd[i]<-norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,M21.071418$Curve[i])
}
M21.071418$Rd<-p.M21.Rd

#Subset data depending on if they were collected in an ascending or descending order
M21.071418.d<-M21.071418[M21.071418$Direction == "down",]
M21.071418.u<-M21.071418[M21.071418$Direction == "up",]

#
#Fit curves
#

#Descending
m.M21.071418.d<-fitacis(M21.071418.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21.071418.d)
coef(m.M21.071418.d)

#Ascending
m.M21.071418.u<-fitacis(M21.071418.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21.071418.u)
coef(m.M21.071418.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M21.071418.d$`10`$df[,2]~m.M21.071418.d$`10`$df[,3])
summary(lm(m.M21.071418.d$`10`$df[,2]~m.M21.071418.d$`10`$df[,3]))
plot(m.M21.071418.d$`15`$df[,2]~m.M21.071418.d$`15`$df[,3])
summary(lm(m.M21.071418.d$`15`$df[,2]~m.M21.071418.d$`15`$df[,3]))
plot(m.M21.071418.d$`20`$df[,2]~m.M21.071418.d$`20`$df[,3])
summary(lm(m.M21.071418.d$`20`$df[,2]~m.M21.071418.d$`20`$df[,3]))

#Ascending
plot(m.M21.071418.u$`20`$df[,2]~m.M21.071418.u$`20`$df[,3])
summary(lm(m.M21.071418.u$`20`$df[,2]~m.M21.071418.u$`20`$df[,3]))
plot(m.M21.071418.u$`25`$df[,2]~m.M21.071418.u$`25`$df[,3])
summary(lm(m.M21.071418.u$`25`$df[,2]~m.M21.071418.u$`25`$df[,3]))
plot(m.M21.071418.u$`30`$df[,2]~m.M21.071418.u$`30`$df[,3])
summary(lm(m.M21.071418.u$`30`$df[,2]~m.M21.071418.u$`30`$df[,3]))
plot(m.M21.071418.u$`35`$df[,2]~m.M21.071418.u$`35`$df[,3])
summary(lm(m.M21.071418.u$`35`$df[,2]~m.M21.071418.u$`35`$df[,3]))
plot(m.M21.071418.u$`40`$df[,2]~m.M21.071418.u$`40`$df[,3])
summary(lm(m.M21.071418.u$`40`$df[,2]~m.M21.071418.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.M21.Rd<-rep(NA,length(M21.120918$Curve))
for (i in 1:length(M21.120918$Curve)){
  p.M21.Rd[i]<-norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,M21.120918$Curve[i])
}
M21.120918$Rd<-p.M21.Rd

#Subset data depending on if they were collected in an ascending or descending order
M21.120918.d<-M21.120918[M21.120918$Direction == "down",]
M21.120918.u<-M21.120918[M21.120918$Direction == "up",]

#
#Fit curves
#

#Descending
m.M21.120918.d<-fitacis(M21.120918.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21.120918.d)
coef(m.M21.120918.d)

#Ascending
m.M21.120918.u<-fitacis(M21.120918.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21.120918.u)
coef(m.M21.120918.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M21.120918.d$`10`$df[,2]~m.M21.120918.d$`10`$df[,3])
summary(lm(m.M21.120918.d$`10`$df[,2]~m.M21.120918.d$`10`$df[,3]))
plot(m.M21.120918.d$`15`$df[,2]~m.M21.120918.d$`15`$df[,3])
summary(lm(m.M21.120918.d$`15`$df[,2]~m.M21.120918.d$`15`$df[,3]))
plot(m.M21.120918.d$`20`$df[,2]~m.M21.120918.d$`20`$df[,3])
summary(lm(m.M21.120918.d$`20`$df[,2]~m.M21.120918.d$`20`$df[,3]))

#Ascending
plot(m.M21.120918.u$`20`$df[,2]~m.M21.120918.u$`20`$df[,3])
summary(lm(m.M21.120918.u$`20`$df[,2]~m.M21.120918.u$`20`$df[,3]))
plot(m.M21.120918.u$`25`$df[,2]~m.M21.120918.u$`25`$df[,3])
summary(lm(m.M21.120918.u$`25`$df[,2]~m.M21.120918.u$`25`$df[,3]))
plot(m.M21.120918.u$`30`$df[,2]~m.M21.120918.u$`30`$df[,3])
summary(lm(m.M21.120918.u$`30`$df[,2]~m.M21.120918.u$`30`$df[,3]))
plot(m.M21.120918.u$`35`$df[,2]~m.M21.120918.u$`35`$df[,3])
summary(lm(m.M21.120918.u$`35`$df[,2]~m.M21.120918.u$`35`$df[,3]))
plot(m.M21.120918.u$`40`$df[,2]~m.M21.120918.u$`40`$df[,3])
summary(lm(m.M21.120918.u$`40`$df[,2]~m.M21.120918.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.M21.Rd<-rep(NA,length(M21.012119$Curve))
for (i in 1:length(M21.012119$Curve)){
  p.M21.Rd[i]<-norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,M21.012119$Curve[i])
}
M21.012119$Rd<-p.M21.Rd

#Subset data depending on if they were collected in an ascending or descending order
M21.012119.d<-M21.012119[M21.012119$Direction == "down",]
M21.012119.u<-M21.012119[M21.012119$Direction == "up",]

#
#Fit curves
#

#Descending
m.M21.012119.d<-fitacis(M21.012119.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21.012119.d)
coef(m.M21.012119.d)

#Ascending
m.M21.012119.u<-fitacis(M21.012119.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21.012119.u)
coef(m.M21.012119.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M21.012119.d$`10`$df[,2]~m.M21.012119.d$`10`$df[,3])
summary(lm(m.M21.012119.d$`10`$df[,2]~m.M21.012119.d$`10`$df[,3]))
plot(m.M21.012119.d$`15`$df[,2]~m.M21.012119.d$`15`$df[,3])
summary(lm(m.M21.012119.d$`15`$df[,2]~m.M21.012119.d$`15`$df[,3]))
plot(m.M21.012119.d$`20`$df[,2]~m.M21.012119.d$`20`$df[,3])
summary(lm(m.M21.012119.d$`20`$df[,2]~m.M21.012119.d$`20`$df[,3]))

#Ascending
plot(m.M21.012119.u$`20`$df[,2]~m.M21.012119.u$`20`$df[,3])
summary(lm(m.M21.012119.u$`20`$df[,2]~m.M21.012119.u$`20`$df[,3]))
plot(m.M21.012119.u$`25`$df[,2]~m.M21.012119.u$`25`$df[,3])
summary(lm(m.M21.012119.u$`25`$df[,2]~m.M21.012119.u$`25`$df[,3]))
plot(m.M21.012119.u$`30`$df[,2]~m.M21.012119.u$`30`$df[,3])
summary(lm(m.M21.012119.u$`30`$df[,2]~m.M21.012119.u$`30`$df[,3]))
plot(m.M21.012119.u$`35`$df[,2]~m.M21.012119.u$`35`$df[,3])
summary(lm(m.M21.012119.u$`35`$df[,2]~m.M21.012119.u$`35`$df[,3]))
plot(m.M21.012119.u$`40`$df[,2]~m.M21.012119.u$`40`$df[,3])
summary(lm(m.M21.012119.u$`40`$df[,2]~m.M21.012119.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.M26.Rd<-rep(NA,length(M26.083019$Curve))
for (i in 1:length(M26.083019$Curve)){
  p.M26.Rd[i]<-norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,M26.083019$Curve[i])
}
M26.083019$Rd<-p.M26.Rd

#Subset data depending on if they were collected in an ascending or descending order
M26.083019.d<-M26.083019[M26.083019$Direction == "down",]
M26.083019.u<-M26.083019[M26.083019$Direction == "up",]

#
#Fit curves
#

#Descending
m.M26.083019.d<-fitacis(M26.083019.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26.083019.d)
coef(m.M26.083019.d)

#Ascending
m.M26.083019.u<-fitacis(M26.083019.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26.083019.u)
coef(m.M26.083019.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M26.083019.d$`10`$df[,2]~m.M26.083019.d$`10`$df[,3])
summary(lm(m.M26.083019.d$`10`$df[,2]~m.M26.083019.d$`10`$df[,3]))
plot(m.M26.083019.d$`15`$df[,2]~m.M26.083019.d$`15`$df[,3])
summary(lm(m.M26.083019.d$`15`$df[,2]~m.M26.083019.d$`15`$df[,3]))
plot(m.M26.083019.d$`20`$df[,2]~m.M26.083019.d$`20`$df[,3])
summary(lm(m.M26.083019.d$`20`$df[,2]~m.M26.083019.d$`20`$df[,3]))
plot(m.M26.083019.d$`25`$df[,2]~m.M26.083019.d$`25`$df[,3])
summary(lm(m.M26.083019.d$`25`$df[,2]~m.M26.083019.d$`25`$df[,3]))

#Ascending
plot(m.M26.083019.u$`25`$df[,2]~m.M26.083019.u$`25`$df[,3])
summary(lm(m.M26.083019.u$`25`$df[,2]~m.M26.083019.u$`25`$df[,3]))
plot(m.M26.083019.u$`30`$df[,2]~m.M26.083019.u$`30`$df[,3])
summary(lm(m.M26.083019.u$`30`$df[,2]~m.M26.083019.u$`30`$df[,3]))
plot(m.M26.083019.u$`35`$df[,2]~m.M26.083019.u$`35`$df[,3])
summary(lm(m.M26.083019.u$`35`$df[,2]~m.M26.083019.u$`35`$df[,3]))
plot(m.M26.083019.u$`40`$df[,2]~m.M26.083019.u$`40`$df[,3])
summary(lm(m.M26.083019.u$`40`$df[,2]~m.M26.083019.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.M26.Rd<-rep(NA,length(M26.091119$Curve))
for (i in 1:length(M26.091119$Curve)){
  p.M26.Rd[i]<-norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,M26.091119$Curve[i])
}
M26.091119$Rd<-p.M26.Rd

#Subset data depending on if they were collected in an ascending or descending order
M26.091119.d<-M26.091119[M26.091119$Direction == "down",]
M26.091119.u<-M26.091119[M26.091119$Direction == "up",]

#
#Fit curves
#

#Descending
m.M26.091119.d<-fitacis(M26.091119.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26.091119.d)
coef(m.M26.091119.d)

#Ascending
m.M26.091119.u<-fitacis(M26.091119.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26.091119.u)
coef(m.M26.091119.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M26.091119.d$`10`$df[,2]~m.M26.091119.d$`10`$df[,3])
summary(lm(m.M26.091119.d$`10`$df[,2]~m.M26.091119.d$`10`$df[,3]))
plot(m.M26.091119.d$`15`$df[,2]~m.M26.091119.d$`15`$df[,3])
summary(lm(m.M26.091119.d$`15`$df[,2]~m.M26.091119.d$`15`$df[,3]))
plot(m.M26.091119.d$`20`$df[,2]~m.M26.091119.d$`20`$df[,3])
summary(lm(m.M26.091119.d$`20`$df[,2]~m.M26.091119.d$`20`$df[,3]))
plot(m.M26.091119.d$`25`$df[,2]~m.M26.091119.d$`25`$df[,3])
summary(lm(m.M26.091119.d$`25`$df[,2]~m.M26.091119.d$`25`$df[,3]))

#Ascending
plot(m.M26.091119.u$`25`$df[,2]~m.M26.091119.u$`25`$df[,3])
summary(lm(m.M26.091119.u$`25`$df[,2]~m.M26.091119.u$`25`$df[,3]))
plot(m.M26.091119.u$`30`$df[,2]~m.M26.091119.u$`30`$df[,3])
summary(lm(m.M26.091119.u$`30`$df[,2]~m.M26.091119.u$`30`$df[,3]))
plot(m.M26.091119.u$`35`$df[,2]~m.M26.091119.u$`35`$df[,3])
summary(lm(m.M26.091119.u$`35`$df[,2]~m.M26.091119.u$`35`$df[,3]))
plot(m.M26.091119.u$`40`$df[,2]~m.M26.091119.u$`40`$df[,3])
summary(lm(m.M26.091119.u$`40`$df[,2]~m.M26.091119.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.M26.Rd<-rep(NA,length(M26.092119$Curve))
for (i in 1:length(M26.092119$Curve)){
  p.M26.Rd[i]<-norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,M26.092119$Curve[i])
}
M26.092119$Rd<-p.M26.Rd

#Subset data depending on if they were collected in an ascending or descending order
M26.092119.d<-M26.092119[M26.092119$Direction == "down",]
M26.092119.u<-M26.092119[M26.092119$Direction == "up",]

#
#Fit curves
#

#Descending
m.M26.092119.d<-fitacis(M26.092119.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26.092119.d)
coef(m.M26.092119.d)

#Ascending
m.M26.092119.u<-fitacis(M26.092119.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26.092119.u)
coef(m.M26.092119.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M26.092119.d$`10`$df[,2]~m.M26.092119.d$`10`$df[,3])
summary(lm(m.M26.092119.d$`10`$df[,2]~m.M26.092119.d$`10`$df[,3]))
plot(m.M26.092119.d$`15`$df[,2]~m.M26.092119.d$`15`$df[,3])
summary(lm(m.M26.092119.d$`15`$df[,2]~m.M26.092119.d$`15`$df[,3]))
plot(m.M26.092119.d$`20`$df[,2]~m.M26.092119.d$`20`$df[,3])
summary(lm(m.M26.092119.d$`20`$df[,2]~m.M26.092119.d$`20`$df[,3]))
plot(m.M26.092119.d$`25`$df[,2]~m.M26.092119.d$`25`$df[,3])
summary(lm(m.M26.092119.d$`25`$df[,2]~m.M26.092119.d$`25`$df[,3]))

#Ascending
plot(m.M26.092119.u$`25`$df[,2]~m.M26.092119.u$`25`$df[,3])
summary(lm(m.M26.092119.u$`25`$df[,2]~m.M26.092119.u$`25`$df[,3]))
plot(m.M26.092119.u$`30`$df[,2]~m.M26.092119.u$`30`$df[,3])
summary(lm(m.M26.092119.u$`30`$df[,2]~m.M26.092119.u$`30`$df[,3]))
plot(m.M26.092119.u$`35`$df[,2]~m.M26.092119.u$`35`$df[,3])
summary(lm(m.M26.092119.u$`35`$df[,2]~m.M26.092119.u$`35`$df[,3]))
plot(m.M26.092119.u$`40`$df[,2]~m.M26.092119.u$`40`$df[,3])
summary(lm(m.M26.092119.u$`40`$df[,2]~m.M26.092119.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.M31.Rd<-rep(NA,length(M31.080218$Curve))
for (i in 1:length(M31.080218$Curve)){
  p.M31.Rd[i]<-norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,M31.080218$Curve[i])
}
M31.080218$Rd<-p.M31.Rd

#Subset data depending on if they were collected in an ascending or descending order
M31.080218.d<-M31.080218[M31.080218$Direction == "down",]
M31.080218.u<-M31.080218[M31.080218$Direction == "up",]

#
#Fit curves
#

#Descending
m.M31.080218.d<-fitacis(M31.080218.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31.080218.d)
coef(m.M31.080218.d)

#Ascending
m.M31.080218.u<-fitacis(M31.080218.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31.080218.u)
coef(m.M31.080218.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M31.080218.d$`10`$df[,2]~m.M31.080218.d$`10`$df[,3])
summary(lm(m.M31.080218.d$`10`$df[,2]~m.M31.080218.d$`10`$df[,3]))
plot(m.M31.080218.d$`15`$df[,2]~m.M31.080218.d$`15`$df[,3])
summary(lm(m.M31.080218.d$`15`$df[,2]~m.M31.080218.d$`15`$df[,3]))
plot(m.M31.080218.d$`20`$df[,2]~m.M31.080218.d$`20`$df[,3])
summary(lm(m.M31.080218.d$`20`$df[,2]~m.M31.080218.d$`20`$df[,3]))
plot(m.M31.080218.d$`25`$df[,2]~m.M31.080218.d$`25`$df[,3])
summary(lm(m.M31.080218.d$`25`$df[,2]~m.M31.080218.d$`25`$df[,3]))
plot(m.M31.080218.d$`30`$df[,2]~m.M31.080218.d$`30`$df[,3])
summary(lm(m.M31.080218.d$`30`$df[,2]~m.M31.080218.d$`30`$df[,3]))

#Ascending
plot(m.M31.080218.u$`30`$df[,2]~m.M31.080218.u$`30`$df[,3])
summary(lm(m.M31.080218.u$`30`$df[,2]~m.M31.080218.u$`30`$df[,3]))
plot(m.M31.080218.u$`35`$df[,2]~m.M31.080218.u$`35`$df[,3])
summary(lm(m.M31.080218.u$`35`$df[,2]~m.M31.080218.u$`35`$df[,3]))
plot(m.M31.080218.u$`40`$df[,2]~m.M31.080218.u$`40`$df[,3])
summary(lm(m.M31.080218.u$`40`$df[,2]~m.M31.080218.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.M31.Rd<-rep(NA,length(M31.120418$Curve))
for (i in 1:length(M31.120418$Curve)){
  p.M31.Rd[i]<-norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,M31.120418$Curve[i])
}
M31.120418$Rd<-p.M31.Rd

#Subset data depending on if they were collected in an ascending or descending order
M31.120418.d<-M31.120418[M31.120418$Direction == "down",]
M31.120418.u<-M31.120418[M31.120418$Direction == "up",]

#
#Fit curves
#

#Descending
m.M31.120418.d<-fitacis(M31.120418.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31.120418.d)
coef(m.M31.120418.d)

#Ascending
m.M31.120418.u<-fitacis(M31.120418.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31.120418.u)
coef(m.M31.120418.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M31.120418.d$`15`$df[,2]~m.M31.120418.d$`15`$df[,3])
summary(lm(m.M31.120418.d$`15`$df[,2]~m.M31.120418.d$`15`$df[,3]))
plot(m.M31.120418.d$`20`$df[,2]~m.M31.120418.d$`20`$df[,3])
summary(lm(m.M31.120418.d$`20`$df[,2]~m.M31.120418.d$`20`$df[,3]))
plot(m.M31.120418.d$`25`$df[,2]~m.M31.120418.d$`25`$df[,3])
summary(lm(m.M31.120418.d$`25`$df[,2]~m.M31.120418.d$`25`$df[,3]))
plot(m.M31.120418.d$`30`$df[,2]~m.M31.120418.d$`30`$df[,3])
summary(lm(m.M31.120418.d$`30`$df[,2]~m.M31.120418.d$`30`$df[,3]))

#Ascending
plot(m.M31.120418.u$`30`$df[,2]~m.M31.120418.u$`30`$df[,3])
summary(lm(m.M31.120418.u$`30`$df[,2]~m.M31.120418.u$`30`$df[,3]))
plot(m.M31.120418.u$`35`$df[,2]~m.M31.120418.u$`35`$df[,3])
summary(lm(m.M31.120418.u$`35`$df[,2]~m.M31.120418.u$`35`$df[,3]))
plot(m.M31.120418.u$`40`$df[,2]~m.M31.120418.u$`40`$df[,3])
summary(lm(m.M31.120418.u$`40`$df[,2]~m.M31.120418.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.M31.Rd<-rep(NA,length(M31.040919$Curve))
for (i in 1:length(M31.040919$Curve)){
  p.M31.Rd[i]<-norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,M31.040919$Curve[i])
}
M31.040919$Rd<-p.M31.Rd

#Subset data depending on if they were collected in an ascending or descending order
M31.040919.d<-M31.040919[M31.040919$Direction == "down",]
M31.040919.u<-M31.040919[M31.040919$Direction == "up",]

#
#Fit curves
#

#Descending
m.M31.040919.d<-fitacis(M31.040919.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31.040919.d)
coef(m.M31.040919.d)

#Ascending
m.M31.040919.u<-fitacis(M31.040919.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31.040919.u)
coef(m.M31.040919.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M31.040919.d$`10`$df[,2]~m.M31.040919.d$`10`$df[,3])
summary(lm(m.M31.040919.d$`10`$df[,2]~m.M31.040919.d$`10`$df[,3]))
plot(m.M31.040919.d$`15`$df[,2]~m.M31.040919.d$`15`$df[,3])
summary(lm(m.M31.040919.d$`15`$df[,2]~m.M31.040919.d$`15`$df[,3]))
plot(m.M31.040919.d$`20`$df[,2]~m.M31.040919.d$`20`$df[,3])
summary(lm(m.M31.040919.d$`20`$df[,2]~m.M31.040919.d$`20`$df[,3]))
plot(m.M31.040919.d$`25`$df[,2]~m.M31.040919.d$`25`$df[,3])
summary(lm(m.M31.040919.d$`25`$df[,2]~m.M31.040919.d$`25`$df[,3]))
plot(m.M31.040919.d$`30`$df[,2]~m.M31.040919.d$`30`$df[,3])
summary(lm(m.M31.040919.d$`30`$df[,2]~m.M31.040919.d$`30`$df[,3]))

#Ascending
plot(m.M31.040919.u$`30`$df[,2]~m.M31.040919.u$`30`$df[,3])
summary(lm(m.M31.040919.u$`30`$df[,2]~m.M31.040919.u$`30`$df[,3]))
plot(m.M31.040919.u$`35`$df[,2]~m.M31.040919.u$`35`$df[,3])
summary(lm(m.M31.040919.u$`35`$df[,2]~m.M31.040919.u$`35`$df[,3]))
plot(m.M31.040919.u$`40`$df[,2]~m.M31.040919.u$`40`$df[,3])
summary(lm(m.M31.040919.u$`40`$df[,2]~m.M31.040919.u$`40`$df[,3]))

####
#Alnus
####

###
#21:15 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.A21.Rd<-rep(NA,length(A21.071218$Curve))
for (i in 1:length(A21.071218$Curve)){
  p.A21.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,A21.071218$Curve[i])
}
A21.071218$Rd<-p.A21.Rd

#Subset data depending on if they were collected in an ascending or descending order
A21.071218.d<-A21.071218[A21.071218$Direction == "down",]
A21.071218.u<-A21.071218[A21.071218$Direction == "up",]

#
#Fit curves
#

#Descending
m.A21.071218.d<-fitacis(A21.071218.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21.071218.d)
coef(m.A21.071218.d)

#Ascending
m.A21.071218.u<-fitacis(A21.071218.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21.071218.u)
coef(m.A21.071218.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A21.071218.d$`10`$df[,2]~m.A21.071218.d$`10`$df[,3])
summary(lm(m.A21.071218.d$`10`$df[,2]~m.A21.071218.d$`10`$df[,3]))
plot(m.A21.071218.d$`15`$df[,2]~m.A21.071218.d$`15`$df[,3])
summary(lm(m.A21.071218.d$`15`$df[,2]~m.A21.071218.d$`15`$df[,3]))
plot(m.A21.071218.d$`20`$df[,2]~m.A21.071218.d$`20`$df[,3])
summary(lm(m.A21.071218.d$`20`$df[,2]~m.A21.071218.d$`20`$df[,3]))

#Ascending
plot(m.A21.071218.u$`20`$df[,2]~m.A21.071218.u$`20`$df[,3])
summary(lm(m.A21.071218.u$`20`$df[,2]~m.A21.071218.u$`20`$df[,3]))
plot(m.A21.071218.u$`25`$df[,2]~m.A21.071218.u$`25`$df[,3])
summary(lm(m.A21.071218.u$`25`$df[,2]~m.A21.071218.u$`25`$df[,3]))
plot(m.A21.071218.u$`30`$df[,2]~m.A21.071218.u$`30`$df[,3])
summary(lm(m.A21.071218.u$`30`$df[,2]~m.A21.071218.u$`30`$df[,3]))
plot(m.A21.071218.u$`35`$df[,2]~m.A21.071218.u$`35`$df[,3])
summary(lm(m.A21.071218.u$`35`$df[,2]~m.A21.071218.u$`35`$df[,3]))
plot(m.A21.071218.u$`40`$df[,2]~m.A21.071218.u$`40`$df[,3])
summary(lm(m.A21.071218.u$`40`$df[,2]~m.A21.071218.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.A21.Rd<-rep(NA,length(A21.071818$Curve))
for (i in 1:length(A21.071818$Curve)){
  p.A21.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,A21.071818$Curve[i])
}
A21.071818$Rd<-p.A21.Rd

#Subset data depending on if they were collected in an ascending or descending order
A21.071818.d<-A21.071818[A21.071818$Direction == "down",]
A21.071818.u<-A21.071818[A21.071818$Direction == "up",]

#
#Fit curves
#

#Descending
m.A21.071818.d<-fitacis(A21.071818.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21.071818.d)
coef(m.A21.071818.d)

#Ascending
m.A21.071818.u<-fitacis(A21.071818.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21.071818.u)
coef(m.A21.071818.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A21.071818.d$`10`$df[,2]~m.A21.071818.d$`10`$df[,3])
summary(lm(m.A21.071818.d$`10`$df[,2]~m.A21.071818.d$`10`$df[,3]))
plot(m.A21.071818.d$`15`$df[,2]~m.A21.071818.d$`15`$df[,3])
summary(lm(m.A21.071818.d$`15`$df[,2]~m.A21.071818.d$`15`$df[,3]))
plot(m.A21.071818.d$`20`$df[,2]~m.A21.071818.d$`20`$df[,3])
summary(lm(m.A21.071818.d$`20`$df[,2]~m.A21.071818.d$`20`$df[,3]))

#Ascending
plot(m.A21.071818.u$`20`$df[,2]~m.A21.071818.u$`20`$df[,3])
summary(lm(m.A21.071818.u$`20`$df[,2]~m.A21.071818.u$`20`$df[,3]))
plot(m.A21.071818.u$`25`$df[,2]~m.A21.071818.u$`25`$df[,3])
summary(lm(m.A21.071818.u$`25`$df[,2]~m.A21.071818.u$`25`$df[,3]))
plot(m.A21.071818.u$`30`$df[,2]~m.A21.071818.u$`30`$df[,3])
summary(lm(m.A21.071818.u$`30`$df[,2]~m.A21.071818.u$`30`$df[,3]))
plot(m.A21.071818.u$`35`$df[,2]~m.A21.071818.u$`35`$df[,3])
summary(lm(m.A21.071818.u$`35`$df[,2]~m.A21.071818.u$`35`$df[,3]))
plot(m.A21.071818.u$`40`$df[,2]~m.A21.071818.u$`40`$df[,3])
summary(lm(m.A21.071818.u$`40`$df[,2]~m.A21.071818.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.A21.Rd<-rep(NA,length(A21.101718$Curve))
for (i in 1:length(A21.101718$Curve)){
  p.A21.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,A21.101718$Curve[i])
}
A21.101718$Rd<-p.A21.Rd

#Subset data depending on if they were collected in an ascending or descending order
A21.101718.d<-A21.101718[A21.101718$Direction == "down",]
A21.101718.u<-A21.101718[A21.101718$Direction == "up",]

#
#Fit curves
#

#Descending
m.A21.101718.d<-fitacis(A21.101718.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21.101718.d)
coef(m.A21.101718.d)

#Ascending
m.A21.101718.u<-fitacis(A21.101718.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21.101718.u)
coef(m.A21.101718.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A21.101718.d$`10`$df[,2]~m.A21.101718.d$`10`$df[,3])
summary(lm(m.A21.101718.d$`10`$df[,2]~m.A21.101718.d$`10`$df[,3]))
plot(m.A21.101718.d$`15`$df[,2]~m.A21.101718.d$`15`$df[,3])
summary(lm(m.A21.101718.d$`15`$df[,2]~m.A21.101718.d$`15`$df[,3]))
plot(m.A21.101718.d$`20`$df[,2]~m.A21.101718.d$`20`$df[,3])
summary(lm(m.A21.101718.d$`20`$df[,2]~m.A21.101718.d$`20`$df[,3]))

#Ascending
plot(m.A21.101718.u$`20`$df[,2]~m.A21.101718.u$`20`$df[,3])
summary(lm(m.A21.101718.u$`20`$df[,2]~m.A21.101718.u$`20`$df[,3]))
plot(m.A21.101718.u$`25`$df[,2]~m.A21.101718.u$`25`$df[,3])
summary(lm(m.A21.101718.u$`25`$df[,2]~m.A21.101718.u$`25`$df[,3]))
plot(m.A21.101718.u$`30`$df[,2]~m.A21.101718.u$`30`$df[,3])
summary(lm(m.A21.101718.u$`30`$df[,2]~m.A21.101718.u$`30`$df[,3]))
plot(m.A21.101718.u$`35`$df[,2]~m.A21.101718.u$`35`$df[,3])
summary(lm(m.A21.101718.u$`35`$df[,2]~m.A21.101718.u$`35`$df[,3]))
plot(m.A21.101718.u$`40`$df[,2]~m.A21.101718.u$`40`$df[,3])
summary(lm(m.A21.101718.u$`40`$df[,2]~m.A21.101718.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.A26.Rd<-rep(NA,length(A26.112419$Curve))
for (i in 1:length(A26.112419$Curve)){
  p.A26.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,A26.112419$Curve[i])
}
A26.112419$Rd<-p.A26.Rd

#Subset data depending on if they were collected in an ascending or descending order
A26.112419.d<-A26.112419[A26.112419$Direction == "down",]
A26.112419.u<-A26.112419[A26.112419$Direction == "up",]

#
#Fit curves
#

#Descending
m.A26.112419.d<-fitacis(A26.112419.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26.112419.d)
coef(m.A26.112419.d)

#Ascending
m.A26.112419.u<-fitacis(A26.112419.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26.112419.u)
coef(m.A26.112419.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A26.112419.d$`10`$df[,2]~m.A26.112419.d$`10`$df[,3])
summary(lm(m.A26.112419.d$`10`$df[,2]~m.A26.112419.d$`10`$df[,3]))
plot(m.A26.112419.d$`15`$df[,2]~m.A26.112419.d$`15`$df[,3])
summary(lm(m.A26.112419.d$`15`$df[,2]~m.A26.112419.d$`15`$df[,3]))
plot(m.A26.112419.d$`20`$df[,2]~m.A26.112419.d$`20`$df[,3])
summary(lm(m.A26.112419.d$`20`$df[,2]~m.A26.112419.d$`20`$df[,3]))
plot(m.A26.112419.d$`25`$df[,2]~m.A26.112419.d$`25`$df[,3])
summary(lm(m.A26.112419.d$`25`$df[,2]~m.A26.112419.d$`25`$df[,3]))

#Ascending
plot(m.A26.112419.u$`25`$df[,2]~m.A26.112419.u$`25`$df[,3])
summary(lm(m.A26.112419.u$`25`$df[,2]~m.A26.112419.u$`25`$df[,3]))
plot(m.A26.112419.u$`30`$df[,2]~m.A26.112419.u$`30`$df[,3])
summary(lm(m.A26.112419.u$`30`$df[,2]~m.A26.112419.u$`30`$df[,3]))
plot(m.A26.112419.u$`35`$df[,2]~m.A26.112419.u$`35`$df[,3])
summary(lm(m.A26.112419.u$`35`$df[,2]~m.A26.112419.u$`35`$df[,3]))
plot(m.A26.112419.u$`40`$df[,2]~m.A26.112419.u$`40`$df[,3])
summary(lm(m.A26.112419.u$`40`$df[,2]~m.A26.112419.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.A26.Rd<-rep(NA,length(A26.020420$Curve))
for (i in 1:length(A26.020420$Curve)){
  p.A26.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,A26.020420$Curve[i])
}
A26.020420$Rd<-p.A26.Rd

#Subset data depending on if they were collected in an ascending or descending order
A26.020420.d<-A26.020420[A26.020420$Direction == "down",]
A26.020420.u<-A26.020420[A26.020420$Direction == "up",]

#
#Fit curves
#

#Descending
m.A26.020420.d<-fitacis(A26.020420.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26.020420.d)
coef(m.A26.020420.d) 

#Ascending
m.A26.020420.u<-fitacis(A26.020420.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26.020420.u)
coef(m.A26.020420.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A26.020420.d$`10`$df[,2]~m.A26.020420.d$`10`$df[,3])
summary(lm(m.A26.020420.d$`10`$df[,2]~m.A26.020420.d$`10`$df[,3]))
plot(m.A26.020420.d$`15`$df[,2]~m.A26.020420.d$`15`$df[,3])
summary(lm(m.A26.020420.d$`15`$df[,2]~m.A26.020420.d$`15`$df[,3]))
plot(m.A26.020420.d$`20`$df[,2]~m.A26.020420.d$`20`$df[,3])
summary(lm(m.A26.020420.d$`20`$df[,2]~m.A26.020420.d$`20`$df[,3]))
plot(m.A26.020420.d$`25`$df[,2]~m.A26.020420.d$`25`$df[,3])
summary(lm(m.A26.020420.d$`25`$df[,2]~m.A26.020420.d$`25`$df[,3]))

#Ascending
plot(m.A26.020420.u$`25`$df[,2]~m.A26.020420.u$`25`$df[,3])
summary(lm(m.A26.020420.u$`25`$df[,2]~m.A26.020420.u$`25`$df[,3]))
plot(m.A26.020420.u$`30`$df[,2]~m.A26.020420.u$`30`$df[,3])
summary(lm(m.A26.020420.u$`30`$df[,2]~m.A26.020420.u$`30`$df[,3]))
plot(m.A26.020420.u$`35`$df[,2]~m.A26.020420.u$`35`$df[,3])
summary(lm(m.A26.020420.u$`35`$df[,2]~m.A26.020420.u$`35`$df[,3]))
plot(m.A26.020420.u$`40`$df[,2]~m.A26.020420.u$`40`$df[,3])
summary(lm(m.A26.020420.u$`40`$df[,2]~m.A26.020420.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.A26.Rd<-rep(NA,length(A26.020820$Curve))
for (i in 1:length(A26.020820$Curve)){
  p.A26.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,A26.020820$Curve[i])
}
A26.020820$Rd<-p.A26.Rd

#Subset data depending on if they were collected in an ascending or descending order
A26.020820.d<-A26.020820[A26.020820$Direction == "down",]
A26.020820.u<-A26.020820[A26.020820$Direction == "up",]

#
#Fit curves
#

#Descending
m.A26.020820.d<-fitacis(A26.020820.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26.020820.d)
coef(m.A26.020820.d) 

#Ascending
m.A26.020820.u<-fitacis(A26.020820.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26.020820.u)
coef(m.A26.020820.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A26.020820.d$`10`$df[,2]~m.A26.020820.d$`10`$df[,3])
summary(lm(m.A26.020820.d$`10`$df[,2]~m.A26.020820.d$`10`$df[,3]))
plot(m.A26.020820.d$`15`$df[,2]~m.A26.020820.d$`15`$df[,3])
summary(lm(m.A26.020820.d$`15`$df[,2]~m.A26.020820.d$`15`$df[,3]))
plot(m.A26.020820.d$`20`$df[,2]~m.A26.020820.d$`20`$df[,3])
summary(lm(m.A26.020820.d$`20`$df[,2]~m.A26.020820.d$`20`$df[,3]))
plot(m.A26.020820.d$`25`$df[,2]~m.A26.020820.d$`25`$df[,3])
summary(lm(m.A26.020820.d$`25`$df[,2]~m.A26.020820.d$`25`$df[,3]))

#Ascending
plot(m.A26.020820.u$`25`$df[,2]~m.A26.020820.u$`25`$df[,3])
summary(lm(m.A26.020820.u$`25`$df[,2]~m.A26.020820.u$`25`$df[,3]))
plot(m.A26.020820.u$`30`$df[,2]~m.A26.020820.u$`30`$df[,3])
summary(lm(m.A26.020820.u$`30`$df[,2]~m.A26.020820.u$`30`$df[,3]))
plot(m.A26.020820.u$`35`$df[,2]~m.A26.020820.u$`35`$df[,3])
summary(lm(m.A26.020820.u$`35`$df[,2]~m.A26.020820.u$`35`$df[,3]))
plot(m.A26.020820.u$`40`$df[,2]~m.A26.020820.u$`40`$df[,3])
summary(lm(m.A26.020820.u$`40`$df[,2]~m.A26.020820.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.A31.Rd<-rep(NA,length(A31.071618$Curve))
for (i in 1:length(A31.071618$Curve)){
  p.A31.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,A31.071618$Curve[i])
}
A31.071618$Rd<-p.A31.Rd

#Subset data depending on if they were collected in an ascending or descending order
A31.071618.d<-A31.071618[A31.071618$Direction == "down",]
A31.071618.u<-A31.071618[A31.071618$Direction == "up",]

#
#Fit curves
#

#Descending
m.A31.071618.d<-fitacis(A31.071618.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31.071618.d)
coef(m.A31.071618.d)

#Ascending
m.A31.071618.u<-fitacis(A31.071618.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31.071618.u)
coef(m.A31.071618.u) 

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A31.071618.d$`10`$df[,2]~m.A31.071618.d$`10`$df[,3])
summary(lm(m.A31.071618.d$`10`$df[,2]~m.A31.071618.d$`10`$df[,3]))
plot(m.A31.071618.d$`15`$df[,2]~m.A31.071618.d$`15`$df[,3])
summary(lm(m.A31.071618.d$`15`$df[,2]~m.A31.071618.d$`15`$df[,3]))
plot(m.A31.071618.d$`20`$df[,2]~m.A31.071618.d$`20`$df[,3])
summary(lm(m.A31.071618.d$`20`$df[,2]~m.A31.071618.d$`20`$df[,3]))
plot(m.A31.071618.d$`25`$df[,2]~m.A31.071618.d$`25`$df[,3])
summary(lm(m.A31.071618.d$`25`$df[,2]~m.A31.071618.d$`25`$df[,3]))
plot(m.A31.071618.d$`30`$df[,2]~m.A31.071618.d$`30`$df[,3])
summary(lm(m.A31.071618.d$`30`$df[,2]~m.A31.071618.d$`30`$df[,3]))

#Ascending
plot(m.A31.071618.u$`30`$df[,2]~m.A31.071618.u$`30`$df[,3])
summary(lm(m.A31.071618.u$`30`$df[,2]~m.A31.071618.u$`30`$df[,3]))
plot(m.A31.071618.u$`35`$df[,2]~m.A31.071618.u$`35`$df[,3])
summary(lm(m.A31.071618.u$`35`$df[,2]~m.A31.071618.u$`35`$df[,3]))
plot(m.A31.071618.u$`40`$df[,2]~m.A31.071618.u$`40`$df[,3])
summary(lm(m.A31.071618.u$`40`$df[,2]~m.A31.071618.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.A31.Rd<-rep(NA,length(A31.092918$Curve))
for (i in 1:length(A31.092918$Curve)){
  p.A31.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,A31.092918$Curve[i])
}
A31.092918$Rd<-p.A31.Rd

#Subset data depending on if they were collected in an ascending or descending order
A31.092918.d<-A31.092918[A31.092918$Direction == "down",]
A31.092918.u<-A31.092918[A31.092918$Direction == "up",]

#
#Fit curves
#

#Descending
m.A31.092918.d<-fitacis(A31.092918.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31.092918.d)
coef(m.A31.092918.d)

#Ascending
m.A31.092918.u<-fitacis(A31.092918.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31.092918.u)
coef(m.A31.092918.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A31.092918.d$`10`$df[,2]~m.A31.092918.d$`10`$df[,3])
summary(lm(m.A31.092918.d$`10`$df[,2]~m.A31.092918.d$`10`$df[,3]))
plot(m.A31.092918.d$`15`$df[,2]~m.A31.092918.d$`15`$df[,3])
summary(lm(m.A31.092918.d$`15`$df[,2]~m.A31.092918.d$`15`$df[,3]))
plot(m.A31.092918.d$`20`$df[,2]~m.A31.092918.d$`20`$df[,3])
summary(lm(m.A31.092918.d$`20`$df[,2]~m.A31.092918.d$`20`$df[,3]))
plot(m.A31.092918.d$`25`$df[,2]~m.A31.092918.d$`25`$df[,3])
summary(lm(m.A31.092918.d$`25`$df[,2]~m.A31.092918.d$`25`$df[,3]))
plot(m.A31.092918.d$`30`$df[,2]~m.A31.092918.d$`30`$df[,3])
summary(lm(m.A31.092918.d$`30`$df[,2]~m.A31.092918.d$`30`$df[,3]))

#Ascending
plot(m.A31.092918.u$`30`$df[,2]~m.A31.092918.u$`30`$df[,3])
summary(lm(m.A31.092918.u$`30`$df[,2]~m.A31.092918.u$`30`$df[,3]))
plot(m.A31.092918.u$`35`$df[,2]~m.A31.092918.u$`35`$df[,3])
summary(lm(m.A31.092918.u$`35`$df[,2]~m.A31.092918.u$`35`$df[,3]))
plot(m.A31.092918.u$`40`$df[,2]~m.A31.092918.u$`40`$df[,3])
summary(lm(m.A31.092918.u$`40`$df[,2]~m.A31.092918.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.A31.Rd<-rep(NA,length(A31.112518$Curve))
for (i in 1:length(A31.112518$Curve)){
  p.A31.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,A31.112518$Curve[i])
}
A31.112518$Rd<-p.A31.Rd

#Subset data depending on if they were collected in an ascending or descending order
A31.112518.d<-A31.112518[A31.112518$Direction == "down",]
A31.112518.u<-A31.112518[A31.112518$Direction == "up",]

#
#Fit curves
#

#Descending
m.A31.112518.d<-fitacis(A31.112518.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31.112518.d)
coef(m.A31.112518.d)

#Ascending
m.A31.112518.u<-fitacis(A31.112518.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31.112518.u)
coef(m.A31.112518.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A31.112518.d$`10`$df[,2]~m.A31.112518.d$`10`$df[,3])
summary(lm(m.A31.112518.d$`10`$df[,2]~m.A31.112518.d$`10`$df[,3]))
plot(m.A31.112518.d$`15`$df[,2]~m.A31.112518.d$`15`$df[,3])
summary(lm(m.A31.112518.d$`15`$df[,2]~m.A31.112518.d$`15`$df[,3]))
plot(m.A31.112518.d$`20`$df[,2]~m.A31.112518.d$`20`$df[,3])
summary(lm(m.A31.112518.d$`20`$df[,2]~m.A31.112518.d$`20`$df[,3]))
plot(m.A31.112518.d$`25`$df[,2]~m.A31.112518.d$`25`$df[,3])
summary(lm(m.A31.112518.d$`25`$df[,2]~m.A31.112518.d$`25`$df[,3]))
plot(m.A31.112518.d$`30`$df[,2]~m.A31.112518.d$`30`$df[,3])
summary(lm(m.A31.112518.d$`30`$df[,2]~m.A31.112518.d$`30`$df[,3]))

#Ascending
plot(m.A31.112518.u$`30`$df[,2]~m.A31.112518.u$`30`$df[,3])
summary(lm(m.A31.112518.u$`30`$df[,2]~m.A31.112518.u$`30`$df[,3]))
plot(m.A31.112518.u$`35`$df[,2]~m.A31.112518.u$`35`$df[,3])
summary(lm(m.A31.112518.u$`35`$df[,2]~m.A31.112518.u$`35`$df[,3]))
plot(m.A31.112518.u$`40`$df[,2]~m.A31.112518.u$`40`$df[,3])
summary(lm(m.A31.112518.u$`40`$df[,2]~m.A31.112518.u$`40`$df[,3]))

####
#Gliricidia
####

###
#21:15 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21.080318$Curve))
for (i in 1:length(G21.080318$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21.080318$Curve[i])
}
G21.080318$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21.080318.d<-G21.080318[G21.080318$Direction == "down",]
G21.080318.u<-G21.080318[G21.080318$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21.080318.d<-fitacis(G21.080318.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.080318.d)
coef(m.G21.080318.d)

#Ascending
m.G21.080318.u<-fitacis(G21.080318.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.080318.u)
coef(m.G21.080318.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21.080318.d$`10`$df[,2]~m.G21.080318.d$`10`$df[,3])
summary(lm(m.G21.080318.d$`10`$df[,2]~m.G21.080318.d$`10`$df[,3]))
plot(m.G21.080318.d$`15`$df[,2]~m.G21.080318.d$`15`$df[,3])
summary(lm(m.G21.080318.d$`15`$df[,2]~m.G21.080318.d$`15`$df[,3]))
plot(m.G21.080318.d$`20`$df[,2]~m.G21.080318.d$`20`$df[,3])
summary(lm(m.G21.080318.d$`20`$df[,2]~m.G21.080318.d$`20`$df[,3]))

#Ascending
plot(m.G21.080318.u$`20`$df[,2]~m.G21.080318.u$`20`$df[,3])
summary(lm(m.G21.080318.u$`20`$df[,2]~m.G21.080318.u$`20`$df[,3]))
plot(m.G21.080318.u$`25`$df[,2]~m.G21.080318.u$`25`$df[,3])
summary(lm(m.G21.080318.u$`25`$df[,2]~m.G21.080318.u$`25`$df[,3]))
plot(m.G21.080318.u$`30`$df[,2]~m.G21.080318.u$`30`$df[,3])
summary(lm(m.G21.080318.u$`30`$df[,2]~m.G21.080318.u$`30`$df[,3]))
plot(m.G21.080318.u$`35`$df[,2]~m.G21.080318.u$`35`$df[,3])
summary(lm(m.G21.080318.u$`35`$df[,2]~m.G21.080318.u$`35`$df[,3]))
plot(m.G21.080318.u$`40`$df[,2]~m.G21.080318.u$`40`$df[,3])
summary(lm(m.G21.080318.u$`40`$df[,2]~m.G21.080318.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21.112718$Curve))
for (i in 1:length(G21.112718$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21.112718$Curve[i])
}
G21.112718$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21.112718.d<-G21.112718[G21.112718$Direction == "down",]
G21.112718.u<-G21.112718[G21.112718$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21.112718.d<-fitacis(G21.112718.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.112718.d)
coef(m.G21.112718.d)

#Ascending
m.G21.112718.u<-fitacis(G21.112718.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.112718.u)
coef(m.G21.112718.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21.112718.d$`10`$df[,2]~m.G21.112718.d$`10`$df[,3])
summary(lm(m.G21.112718.d$`10`$df[,2]~m.G21.112718.d$`10`$df[,3]))
plot(m.G21.112718.d$`15`$df[,2]~m.G21.112718.d$`15`$df[,3])
summary(lm(m.G21.112718.d$`15`$df[,2]~m.G21.112718.d$`15`$df[,3]))
plot(m.G21.112718.d$`20`$df[,2]~m.G21.112718.d$`20`$df[,3])
summary(lm(m.G21.112718.d$`20`$df[,2]~m.G21.112718.d$`20`$df[,3]))

#Ascending
plot(m.G21.112718.u$`20`$df[,2]~m.G21.112718.u$`20`$df[,3])
summary(lm(m.G21.112718.u$`20`$df[,2]~m.G21.112718.u$`20`$df[,3]))
plot(m.G21.112718.u$`25`$df[,2]~m.G21.112718.u$`25`$df[,3])
summary(lm(m.G21.112718.u$`25`$df[,2]~m.G21.112718.u$`25`$df[,3]))
plot(m.G21.112718.u$`30`$df[,2]~m.G21.112718.u$`30`$df[,3])
summary(lm(m.G21.112718.u$`30`$df[,2]~m.G21.112718.u$`30`$df[,3]))
plot(m.G21.112718.u$`35`$df[,2]~m.G21.112718.u$`35`$df[,3])
summary(lm(m.G21.112718.u$`35`$df[,2]~m.G21.112718.u$`35`$df[,3]))
plot(m.G21.112718.u$`40`$df[,2]~m.G21.112718.u$`40`$df[,3])
summary(lm(m.G21.112718.u$`40`$df[,2]~m.G21.112718.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21.082019$Curve))
for (i in 1:length(G21.082019$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21.082019$Curve[i])
}
G21.082019$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21.082019.d<-G21.082019[G21.082019$Direction == "down",]
G21.082019.u<-G21.082019[G21.082019$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21.082019.d<-fitacis(G21.082019.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.082019.d)
coef(m.G21.082019.d)

#Ascending
m.G21.082019.u<-fitacis(G21.082019.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.082019.u)
coef(m.G21.082019.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21.082019.d$`10`$df[,2]~m.G21.082019.d$`10`$df[,3])
summary(lm(m.G21.082019.d$`10`$df[,2]~m.G21.082019.d$`10`$df[,3]))
plot(m.G21.082019.d$`15`$df[,2]~m.G21.082019.d$`15`$df[,3])
summary(lm(m.G21.082019.d$`15`$df[,2]~m.G21.082019.d$`15`$df[,3]))
plot(m.G21.082019.d$`20`$df[,2]~m.G21.082019.d$`20`$df[,3])
summary(lm(m.G21.082019.d$`20`$df[,2]~m.G21.082019.d$`20`$df[,3]))

#Ascending
plot(m.G21.082019.u$`20`$df[,2]~m.G21.082019.u$`20`$df[,3])
summary(lm(m.G21.082019.u$`20`$df[,2]~m.G21.082019.u$`20`$df[,3]))
plot(m.G21.082019.u$`25`$df[,2]~m.G21.082019.u$`25`$df[,3])
summary(lm(m.G21.082019.u$`25`$df[,2]~m.G21.082019.u$`25`$df[,3]))
plot(m.G21.082019.u$`30`$df[,2]~m.G21.082019.u$`30`$df[,3])
summary(lm(m.G21.082019.u$`30`$df[,2]~m.G21.082019.u$`30`$df[,3]))
plot(m.G21.082019.u$`35`$df[,2]~m.G21.082019.u$`35`$df[,3])
summary(lm(m.G21.082019.u$`35`$df[,2]~m.G21.082019.u$`35`$df[,3]))
plot(m.G21.082019.u$`40`$df[,2]~m.G21.082019.u$`40`$df[,3])
summary(lm(m.G21.082019.u$`40`$df[,2]~m.G21.082019.u$`40`$df[,3]))

##
#Replicate 4
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21.082819$Curve))
for (i in 1:length(G21.082819$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21.082819$Curve[i])
}
G21.082819$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21.082819.d<-G21.082819[G21.082819$Direction == "down",]
G21.082819.u<-G21.082819[G21.082819$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21.082819.d<-fitacis(G21.082819.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.082819.d)
coef(m.G21.082819.d)

#Ascending
m.G21.082819.u<-fitacis(G21.082819.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21.082819.u)
coef(m.G21.082819.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21.082819.d$`10`$df[,2]~m.G21.082819.d$`10`$df[,3])
summary(lm(m.G21.082819.d$`10`$df[,2]~m.G21.082819.d$`10`$df[,3]))
plot(m.G21.082819.d$`15`$df[,2]~m.G21.082819.d$`15`$df[,3])
summary(lm(m.G21.082819.d$`15`$df[,2]~m.G21.082819.d$`15`$df[,3]))
plot(m.G21.082819.d$`20`$df[,2]~m.G21.082819.d$`20`$df[,3])
summary(lm(m.G21.082819.d$`20`$df[,2]~m.G21.082819.d$`20`$df[,3]))

#Ascending
plot(m.G21.082819.u$`20`$df[,2]~m.G21.082819.u$`20`$df[,3])
summary(lm(m.G21.082819.u$`20`$df[,2]~m.G21.082819.u$`20`$df[,3]))
plot(m.G21.082819.u$`25`$df[,2]~m.G21.082819.u$`25`$df[,3])
summary(lm(m.G21.082819.u$`25`$df[,2]~m.G21.082819.u$`25`$df[,3]))
plot(m.G21.082819.u$`30`$df[,2]~m.G21.082819.u$`30`$df[,3])
summary(lm(m.G21.082819.u$`30`$df[,2]~m.G21.082819.u$`30`$df[,3]))
plot(m.G21.082819.u$`35`$df[,2]~m.G21.082819.u$`35`$df[,3])
summary(lm(m.G21.082819.u$`35`$df[,2]~m.G21.082819.u$`35`$df[,3]))
plot(m.G21.082819.u$`40`$df[,2]~m.G21.082819.u$`40`$df[,3])
summary(lm(m.G21.082819.u$`40`$df[,2]~m.G21.082819.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26.082319$Curve))
for (i in 1:length(G26.082319$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26.082319$Curve[i])
}
G26.082319$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26.082319.d<-G26.082319[G26.082319$Direction == "down",]
G26.082319.u<-G26.082319[G26.082319$Direction == "up",]

#
#Fit curves
#

#Descending
m.G26.082319.d<-fitacis(G26.082319.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.082319.d)
coef(m.G26.082319.d)

#Ascending
m.G26.082319.u<-fitacis(G26.082319.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.082319.u)
coef(m.G26.082319.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G26.082319.d$`10`$df[,2]~m.G26.082319.d$`10`$df[,3])
summary(lm(m.G26.082319.d$`10`$df[,2]~m.G26.082319.d$`10`$df[,3]))
plot(m.G26.082319.d$`15`$df[,2]~m.G26.082319.d$`15`$df[,3])
summary(lm(m.G26.082319.d$`15`$df[,2]~m.G26.082319.d$`15`$df[,3]))
plot(m.G26.082319.d$`20`$df[,2]~m.G26.082319.d$`20`$df[,3])
summary(lm(m.G26.082319.d$`20`$df[,2]~m.G26.082319.d$`20`$df[,3]))
plot(m.G26.082319.d$`25`$df[,2]~m.G26.082319.d$`25`$df[,3])
summary(lm(m.G26.082319.d$`25`$df[,2]~m.G26.082319.d$`25`$df[,3]))

#Ascending
plot(m.G26.082319.u$`25`$df[,2]~m.G26.082319.u$`25`$df[,3])
summary(lm(m.G26.082319.u$`25`$df[,2]~m.G26.082319.u$`25`$df[,3]))
plot(m.G26.082319.u$`30`$df[,2]~m.G26.082319.u$`30`$df[,3])
summary(lm(m.G26.082319.u$`30`$df[,2]~m.G26.082319.u$`30`$df[,3]))
plot(m.G26.082319.u$`35`$df[,2]~m.G26.082319.u$`35`$df[,3])
summary(lm(m.G26.082319.u$`35`$df[,2]~m.G26.082319.u$`35`$df[,3]))
plot(m.G26.082319.u$`40`$df[,2]~m.G26.082319.u$`40`$df[,3])
summary(lm(m.G26.082319.u$`40`$df[,2]~m.G26.082319.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26.090619$Curve))
for (i in 1:length(G26.090619$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26.090619$Curve[i])
}
G26.090619$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26.090619.d<-G26.090619[G26.090619$Direction == "down",]
G26.090619.u<-G26.090619[G26.090619$Direction == "up",]

#
#Fit curves
#

#Descending
m.G26.090619.d<-fitacis(G26.090619.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.090619.d)
coef(m.G26.090619.d)

#Ascending
m.G26.090619.u<-fitacis(G26.090619.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.090619.u)
coef(m.G26.090619.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G26.090619.d$`10`$df[,2]~m.G26.090619.d$`10`$df[,3])
summary(lm(m.G26.090619.d$`10`$df[,2]~m.G26.090619.d$`10`$df[,3]))
plot(m.G26.090619.d$`15`$df[,2]~m.G26.090619.d$`15`$df[,3])
summary(lm(m.G26.090619.d$`15`$df[,2]~m.G26.090619.d$`15`$df[,3]))
plot(m.G26.090619.d$`20`$df[,2]~m.G26.090619.d$`20`$df[,3])
summary(lm(m.G26.090619.d$`20`$df[,2]~m.G26.090619.d$`20`$df[,3]))
plot(m.G26.090619.d$`25`$df[,2]~m.G26.090619.d$`25`$df[,3])
summary(lm(m.G26.090619.d$`25`$df[,2]~m.G26.090619.d$`25`$df[,3]))

#Ascending
plot(m.G26.090619.u$`25`$df[,2]~m.G26.090619.u$`25`$df[,3])
summary(lm(m.G26.090619.u$`25`$df[,2]~m.G26.090619.u$`25`$df[,3]))
plot(m.G26.090619.u$`30`$df[,2]~m.G26.090619.u$`30`$df[,3])
summary(lm(m.G26.090619.u$`30`$df[,2]~m.G26.090619.u$`30`$df[,3]))
plot(m.G26.090619.u$`35`$df[,2]~m.G26.090619.u$`35`$df[,3])
summary(lm(m.G26.090619.u$`35`$df[,2]~m.G26.090619.u$`35`$df[,3]))
plot(m.G26.090619.u$`40`$df[,2]~m.G26.090619.u$`40`$df[,3])
summary(lm(m.G26.090619.u$`40`$df[,2]~m.G26.090619.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26.100319$Curve))
for (i in 1:length(G26.100319$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26.100319$Curve[i])
}
G26.100319$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26.100319.d<-G26.100319[G26.100319$Direction == "down",]
G26.100319.u<-G26.100319[G26.100319$Direction == "up",]

#
#Fit curves
#

#Descending
m.G26.100319.d<-fitacis(G26.100319.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.100319.d)
coef(m.G26.100319.d)

#Ascending
m.G26.100319.u<-fitacis(G26.100319.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.100319.u)
coef(m.G26.100319.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G26.100319.d$`10`$df[,2]~m.G26.100319.d$`10`$df[,3])
summary(lm(m.G26.100319.d$`10`$df[,2]~m.G26.100319.d$`10`$df[,3]))
plot(m.G26.100319.d$`15`$df[,2]~m.G26.100319.d$`15`$df[,3])
summary(lm(m.G26.100319.d$`15`$df[,2]~m.G26.100319.d$`15`$df[,3]))
plot(m.G26.100319.d$`20`$df[,2]~m.G26.100319.d$`20`$df[,3])
summary(lm(m.G26.100319.d$`20`$df[,2]~m.G26.100319.d$`20`$df[,3]))
plot(m.G26.100319.d$`25`$df[,2]~m.G26.100319.d$`25`$df[,3])
summary(lm(m.G26.100319.d$`25`$df[,2]~m.G26.100319.d$`25`$df[,3]))

#Ascending
plot(m.G26.100319.u$`25`$df[,2]~m.G26.100319.u$`25`$df[,3])
summary(lm(m.G26.100319.u$`25`$df[,2]~m.G26.100319.u$`25`$df[,3]))
plot(m.G26.100319.u$`30`$df[,2]~m.G26.100319.u$`30`$df[,3])
summary(lm(m.G26.100319.u$`30`$df[,2]~m.G26.100319.u$`30`$df[,3]))
plot(m.G26.100319.u$`35`$df[,2]~m.G26.100319.u$`35`$df[,3])
summary(lm(m.G26.100319.u$`35`$df[,2]~m.G26.100319.u$`35`$df[,3]))
plot(m.G26.100319.u$`40`$df[,2]~m.G26.100319.u$`40`$df[,3])
summary(lm(m.G26.100319.u$`40`$df[,2]~m.G26.100319.u$`40`$df[,3]))

##
#Replicate 4
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26.110119$Curve))
for (i in 1:length(G26.110119$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26.110119$Curve[i])
}
G26.110119$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26.110119.u<-G26.110119[G26.110119$Direction == "up",]

#
#Fit curves
#

#Ascending
m.G26.110119.u<-fitacis(G26.110119.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26.110119.u)
coef(m.G26.110119.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Ascending
plot(m.G26.110119.u$`25`$df[,2]~m.G26.110119.u$`25`$df[,3])
summary(lm(m.G26.110119.u$`25`$df[,2]~m.G26.110119.u$`25`$df[,3]))
plot(m.G26.110119.u$`30`$df[,2]~m.G26.110119.u$`30`$df[,3])
summary(lm(m.G26.110119.u$`30`$df[,2]~m.G26.110119.u$`30`$df[,3]))
plot(m.G26.110119.u$`35`$df[,2]~m.G26.110119.u$`35`$df[,3])
summary(lm(m.G26.110119.u$`35`$df[,2]~m.G26.110119.u$`35`$df[,3]))
plot(m.G26.110119.u$`40`$df[,2]~m.G26.110119.u$`40`$df[,3])
summary(lm(m.G26.110119.u$`40`$df[,2]~m.G26.110119.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.G31.Rd<-rep(NA,length(G31.050418$Curve))
for (i in 1:length(G31.050418$Curve)){
  p.G31.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,G31.050418$Curve[i])
}
G31.050418$Rd<-p.G31.Rd

#Subset data depending on if they were collected in an ascending or descending order
G31.050418.d<-G31.050418[G31.050418$Direction == "down",]
G31.050418.u<-G31.050418[G31.050418$Direction == "up",]

#
#Fit curves
#

#Descending
m.G31.050418.d<-fitacis(G31.050418.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31.050418.d)
coef(m.G31.050418.d)

#Ascending
m.G31.050418.u<-fitacis(G31.050418.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31.050418.u)
coef(m.G31.050418.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G31.050418.d$`15`$df[,2]~m.G31.050418.d$`15`$df[,3])
summary(lm(m.G31.050418.d$`15`$df[,2]~m.G31.050418.d$`15`$df[,3]))
plot(m.G31.050418.d$`20`$df[,2]~m.G31.050418.d$`20`$df[,3])
summary(lm(m.G31.050418.d$`20`$df[,2]~m.G31.050418.d$`20`$df[,3]))
plot(m.G31.050418.d$`25`$df[,2]~m.G31.050418.d$`25`$df[,3])
summary(lm(m.G31.050418.d$`25`$df[,2]~m.G31.050418.d$`25`$df[,3]))
plot(m.G31.050418.d$`30`$df[,2]~m.G31.050418.d$`30`$df[,3])
summary(lm(m.G31.050418.d$`30`$df[,2]~m.G31.050418.d$`30`$df[,3]))

#Ascending
plot(m.G31.050418.u$`30`$df[,2]~m.G31.050418.u$`30`$df[,3])
summary(lm(m.G31.050418.u$`30`$df[,2]~m.G31.050418.u$`30`$df[,3]))
plot(m.G31.050418.u$`35`$df[,2]~m.G31.050418.u$`35`$df[,3])
summary(lm(m.G31.050418.u$`35`$df[,2]~m.G31.050418.u$`35`$df[,3]))
plot(m.G31.050418.u$`40`$df[,2]~m.G31.050418.u$`40`$df[,3])
summary(lm(m.G31.050418.u$`40`$df[,2]~m.G31.050418.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.G31.Rd<-rep(NA,length(G31.071118$Curve))
for (i in 1:length(G31.071118$Curve)){
  p.G31.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,G31.071118$Curve[i])
}
G31.071118$Rd<-p.G31.Rd

#Subset data depending on if they were collected in an ascending or descending order
G31.071118.d<-G31.071118[G31.071118$Direction == "down",]
G31.071118.u<-G31.071118[G31.071118$Direction == "up",]

#
#Fit curves
#

#Descending
m.G31.071118.d<-fitacis(G31.071118.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31.071118.d)
coef(m.G31.071118.d)

#Ascending
m.G31.071118.u<-fitacis(G31.071118.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31.071118.u)
coef(m.G31.071118.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G31.071118.d$`10`$df[,2]~m.G31.071118.d$`10`$df[,3])
summary(lm(m.G31.071118.d$`10`$df[,2]~m.G31.071118.d$`10`$df[,3]))
plot(m.G31.071118.d$`15`$df[,2]~m.G31.071118.d$`15`$df[,3])
summary(lm(m.G31.071118.d$`15`$df[,2]~m.G31.071118.d$`15`$df[,3]))
plot(m.G31.071118.d$`20`$df[,2]~m.G31.071118.d$`20`$df[,3])
summary(lm(m.G31.071118.d$`20`$df[,2]~m.G31.071118.d$`20`$df[,3]))
plot(m.G31.071118.d$`25`$df[,2]~m.G31.071118.d$`25`$df[,3])
summary(lm(m.G31.071118.d$`25`$df[,2]~m.G31.071118.d$`25`$df[,3]))
plot(m.G31.071118.d$`30`$df[,2]~m.G31.071118.d$`30`$df[,3])
summary(lm(m.G31.071118.d$`30`$df[,2]~m.G31.071118.d$`30`$df[,3]))

#Ascending
plot(m.G31.071118.u$`30`$df[,2]~m.G31.071118.u$`30`$df[,3])
summary(lm(m.G31.071118.u$`30`$df[,2]~m.G31.071118.u$`30`$df[,3]))
plot(m.G31.071118.u$`35`$df[,2]~m.G31.071118.u$`35`$df[,3])
summary(lm(m.G31.071118.u$`35`$df[,2]~m.G31.071118.u$`35`$df[,3]))
plot(m.G31.071118.u$`40`$df[,2]~m.G31.071118.u$`40`$df[,3])
summary(lm(m.G31.071118.u$`40`$df[,2]~m.G31.071118.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.G31.Rd<-rep(NA,length(G31.072218$Curve))
for (i in 1:length(G31.072218$Curve)){
  p.G31.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,G31.072218$Curve[i])
}
G31.072218$Rd<-p.G31.Rd

#Subset data depending on if they were collected in an ascending or descending order
G31.072218.d<-G31.072218[G31.072218$Direction == "down",]
G31.072218.u<-G31.072218[G31.072218$Direction == "up",]

#
#Fit curves
#

#Descending
m.G31.072218.d<-fitacis(G31.072218.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31.072218.d)
coef(m.G31.072218.d)

#Ascending
m.G31.072218.u<-fitacis(G31.072218.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31.072218.u)
coef(m.G31.072218.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G31.072218.d$`10`$df[,2]~m.G31.072218.d$`10`$df[,3])
summary(lm(m.G31.072218.d$`10`$df[,2]~m.G31.072218.d$`10`$df[,3]))
plot(m.G31.072218.d$`15`$df[,2]~m.G31.072218.d$`15`$df[,3])
summary(lm(m.G31.072218.d$`15`$df[,2]~m.G31.072218.d$`15`$df[,3]))
plot(m.G31.072218.d$`20`$df[,2]~m.G31.072218.d$`20`$df[,3])
summary(lm(m.G31.072218.d$`20`$df[,2]~m.G31.072218.d$`20`$df[,3]))
plot(m.G31.072218.d$`25`$df[,2]~m.G31.072218.d$`25`$df[,3])
summary(lm(m.G31.072218.d$`25`$df[,2]~m.G31.072218.d$`25`$df[,3]))
plot(m.G31.072218.d$`30`$df[,2]~m.G31.072218.d$`30`$df[,3])
summary(lm(m.G31.072218.d$`30`$df[,2]~m.G31.072218.d$`30`$df[,3]))
plot(m.G31.072218.d$`32`$df[,2]~m.G31.072218.d$`32`$df[,3])
summary(lm(m.G31.072218.d$`32`$df[,2]~m.G31.072218.d$`32`$df[,3]))

#Ascending
plot(m.G31.072218.u$`30`$df[,2]~m.G31.072218.u$`30`$df[,3])
summary(lm(m.G31.072218.u$`30`$df[,2]~m.G31.072218.u$`30`$df[,3]))
plot(m.G31.072218.u$`35`$df[,2]~m.G31.072218.u$`35`$df[,3])
summary(lm(m.G31.072218.u$`35`$df[,2]~m.G31.072218.u$`35`$df[,3]))
plot(m.G31.072218.u$`40`$df[,2]~m.G31.072218.u$`40`$df[,3])
summary(lm(m.G31.072218.u$`40`$df[,2]~m.G31.072218.u$`40`$df[,3]))

####
#Robinia
####

###
#21:15 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21.073118$Curve))
for (i in 1:length(R21.073118$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21.073118$Curve[i])
}
R21.073118$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21.073118.d<-R21.073118[R21.073118$Direction == "down",]
R21.073118.u<-R21.073118[R21.073118$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21.073118.d<-fitacis(R21.073118.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.073118.d)
coef(m.R21.073118.d)

#Ascending
m.R21.073118.u<-fitacis(R21.073118.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.073118.u)
coef(m.R21.073118.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21.073118.d$`10`$df[,2]~m.R21.073118.d$`10`$df[,3])
summary(lm(m.R21.073118.d$`10`$df[,2]~m.R21.073118.d$`10`$df[,3]))
plot(m.R21.073118.d$`15`$df[,2]~m.R21.073118.d$`15`$df[,3])
summary(lm(m.R21.073118.d$`15`$df[,2]~m.R21.073118.d$`15`$df[,3]))
plot(m.R21.073118.d$`20`$df[,2]~m.R21.073118.d$`20`$df[,3])
summary(lm(m.R21.073118.d$`20`$df[,2]~m.R21.073118.d$`20`$df[,3]))

#Ascending
plot(m.R21.073118.u$`20`$df[,2]~m.R21.073118.u$`20`$df[,3])
summary(lm(m.R21.073118.u$`20`$df[,2]~m.R21.073118.u$`20`$df[,3]))
plot(m.R21.073118.u$`25`$df[,2]~m.R21.073118.u$`25`$df[,3])
summary(lm(m.R21.073118.u$`25`$df[,2]~m.R21.073118.u$`25`$df[,3]))
plot(m.R21.073118.u$`30`$df[,2]~m.R21.073118.u$`30`$df[,3])
summary(lm(m.R21.073118.u$`30`$df[,2]~m.R21.073118.u$`30`$df[,3]))
plot(m.R21.073118.u$`35`$df[,2]~m.R21.073118.u$`35`$df[,3])
summary(lm(m.R21.073118.u$`35`$df[,2]~m.R21.073118.u$`35`$df[,3]))
plot(m.R21.073118.u$`40`$df[,2]~m.R21.073118.u$`40`$df[,3])
summary(lm(m.R21.073118.u$`40`$df[,2]~m.R21.073118.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21.100518$Curve))
for (i in 1:length(R21.100518$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21.100518$Curve[i])
}
R21.100518$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21.100518.d<-R21.100518[R21.100518$Direction == "down",]
R21.100518.u<-R21.100518[R21.100518$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21.100518.d<-fitacis(R21.100518.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.100518.d)
coef(m.R21.100518.d)

#Ascending
m.R21.100518.u<-fitacis(R21.100518.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.100518.u)
coef(m.R21.100518.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21.100518.d$`10`$df[,2]~m.R21.100518.d$`10`$df[,3])
summary(lm(m.R21.100518.d$`10`$df[,2]~m.R21.100518.d$`10`$df[,3]))
plot(m.R21.100518.d$`15`$df[,2]~m.R21.100518.d$`15`$df[,3])
summary(lm(m.R21.100518.d$`15`$df[,2]~m.R21.100518.d$`15`$df[,3]))
plot(m.R21.100518.d$`20`$df[,2]~m.R21.100518.d$`20`$df[,3])
summary(lm(m.R21.100518.d$`20`$df[,2]~m.R21.100518.d$`20`$df[,3]))

#Ascending
plot(m.R21.100518.u$`20`$df[,2]~m.R21.100518.u$`20`$df[,3])
summary(lm(m.R21.100518.u$`20`$df[,2]~m.R21.100518.u$`20`$df[,3]))
plot(m.R21.100518.u$`25`$df[,2]~m.R21.100518.u$`25`$df[,3])
summary(lm(m.R21.100518.u$`25`$df[,2]~m.R21.100518.u$`25`$df[,3]))
plot(m.R21.100518.u$`30`$df[,2]~m.R21.100518.u$`30`$df[,3])
summary(lm(m.R21.100518.u$`30`$df[,2]~m.R21.100518.u$`30`$df[,3]))
plot(m.R21.100518.u$`35`$df[,2]~m.R21.100518.u$`35`$df[,3])
summary(lm(m.R21.100518.u$`35`$df[,2]~m.R21.100518.u$`35`$df[,3]))
plot(m.R21.100518.u$`40`$df[,2]~m.R21.100518.u$`40`$df[,3])
summary(lm(m.R21.100518.u$`40`$df[,2]~m.R21.100518.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21.101018$Curve))
for (i in 1:length(R21.101018$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21.101018$Curve[i])
}
R21.101018$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21.101018.d<-R21.101018[R21.101018$Direction == "down",]
R21.101018.u<-R21.101018[R21.101018$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21.101018.d<-fitacis(R21.101018.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.101018.d)
coef(m.R21.101018.d)

#Ascending
m.R21.101018.u<-fitacis(R21.101018.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.101018.u)
coef(m.R21.101018.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21.101018.d$`10`$df[,2]~m.R21.101018.d$`10`$df[,3])
summary(lm(m.R21.101018.d$`10`$df[,2]~m.R21.101018.d$`10`$df[,3]))
plot(m.R21.101018.d$`15`$df[,2]~m.R21.101018.d$`15`$df[,3])
summary(lm(m.R21.101018.d$`15`$df[,2]~m.R21.101018.d$`15`$df[,3]))
plot(m.R21.101018.d$`20`$df[,2]~m.R21.101018.d$`20`$df[,3])
summary(lm(m.R21.101018.d$`20`$df[,2]~m.R21.101018.d$`20`$df[,3]))

#Ascending
plot(m.R21.101018.u$`20`$df[,2]~m.R21.101018.u$`20`$df[,3])
summary(lm(m.R21.101018.u$`20`$df[,2]~m.R21.101018.u$`20`$df[,3]))
plot(m.R21.101018.u$`25`$df[,2]~m.R21.101018.u$`25`$df[,3])
summary(lm(m.R21.101018.u$`25`$df[,2]~m.R21.101018.u$`25`$df[,3]))
plot(m.R21.101018.u$`30`$df[,2]~m.R21.101018.u$`30`$df[,3])
summary(lm(m.R21.101018.u$`30`$df[,2]~m.R21.101018.u$`30`$df[,3]))
plot(m.R21.101018.u$`35`$df[,2]~m.R21.101018.u$`35`$df[,3])
summary(lm(m.R21.101018.u$`35`$df[,2]~m.R21.101018.u$`35`$df[,3]))
plot(m.R21.101018.u$`40`$df[,2]~m.R21.101018.u$`40`$df[,3])
summary(lm(m.R21.101018.u$`40`$df[,2]~m.R21.101018.u$`40`$df[,3]))

##
#Replicate 4
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21.091719$Curve))
for (i in 1:length(R21.091719$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21.091719$Curve[i])
}
R21.091719$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21.091719.d<-R21.091719[R21.091719$Direction == "down",]
R21.091719.u<-R21.091719[R21.091719$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21.091719.d<-fitacis(R21.091719.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.091719.d)
coef(m.R21.091719.d)

#Ascending
m.R21.091719.u<-fitacis(R21.091719.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21.091719.u)
coef(m.R21.091719.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21.091719.d$`10`$df[,2]~m.R21.091719.d$`10`$df[,3])
summary(lm(m.R21.091719.d$`10`$df[,2]~m.R21.091719.d$`10`$df[,3]))
plot(m.R21.091719.d$`15`$df[,2]~m.R21.091719.d$`15`$df[,3])
summary(lm(m.R21.091719.d$`15`$df[,2]~m.R21.091719.d$`15`$df[,3]))
plot(m.R21.091719.d$`20`$df[,2]~m.R21.091719.d$`20`$df[,3])
summary(lm(m.R21.091719.d$`20`$df[,2]~m.R21.091719.d$`20`$df[,3]))

#Ascending
plot(m.R21.091719.u$`20`$df[,2]~m.R21.091719.u$`20`$df[,3])
summary(lm(m.R21.091719.u$`20`$df[,2]~m.R21.091719.u$`20`$df[,3]))
plot(m.R21.091719.u$`25`$df[,2]~m.R21.091719.u$`25`$df[,3])
summary(lm(m.R21.091719.u$`25`$df[,2]~m.R21.091719.u$`25`$df[,3]))
plot(m.R21.091719.u$`30`$df[,2]~m.R21.091719.u$`30`$df[,3])
summary(lm(m.R21.091719.u$`30`$df[,2]~m.R21.091719.u$`30`$df[,3]))
plot(m.R21.091719.u$`35`$df[,2]~m.R21.091719.u$`35`$df[,3])
summary(lm(m.R21.091719.u$`35`$df[,2]~m.R21.091719.u$`35`$df[,3]))
plot(m.R21.091719.u$`40`$df[,2]~m.R21.091719.u$`40`$df[,3])
summary(lm(m.R21.091719.u$`40`$df[,2]~m.R21.091719.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.R26.Rd<-rep(NA,length(R26.102219$Curve))
for (i in 1:length(R26.102219$Curve)){
  p.R26.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,R26.102219$Curve[i])
}
R26.102219$Rd<-p.R26.Rd

#Subset data depending on if they were collected in an ascending or descending order
R26.102219.d<-R26.102219[R26.102219$Direction == "down",]
R26.102219.u<-R26.102219[R26.102219$Direction == "up",]

#
#Fit curves
#

#Descending
m.R26.102219.d<-fitacis(R26.102219.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26.102219.d)
coef(m.R26.102219.d)

#Ascending
m.R26.102219.u<-fitacis(R26.102219.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26.102219.u)
coef(m.R26.102219.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R26.102219.d$`10`$df[,2]~m.R26.102219.d$`10`$df[,3])
summary(lm(m.R26.102219.d$`10`$df[,2]~m.R26.102219.d$`10`$df[,3]))
plot(m.R26.102219.d$`15`$df[,2]~m.R26.102219.d$`15`$df[,3])
summary(lm(m.R26.102219.d$`15`$df[,2]~m.R26.102219.d$`15`$df[,3]))
plot(m.R26.102219.d$`20`$df[,2]~m.R26.102219.d$`20`$df[,3])
summary(lm(m.R26.102219.d$`20`$df[,2]~m.R26.102219.d$`20`$df[,3]))
plot(m.R26.102219.d$`25`$df[,2]~m.R26.102219.d$`25`$df[,3])
summary(lm(m.R26.102219.d$`25`$df[,2]~m.R26.102219.d$`25`$df[,3]))

#Ascending
plot(m.R26.102219.u$`25`$df[,2]~m.R26.102219.u$`25`$df[,3])
summary(lm(m.R26.102219.u$`25`$df[,2]~m.R26.102219.u$`25`$df[,3]))
plot(m.R26.102219.u$`30`$df[,2]~m.R26.102219.u$`30`$df[,3])
summary(lm(m.R26.102219.u$`30`$df[,2]~m.R26.102219.u$`30`$df[,3]))
plot(m.R26.102219.u$`35`$df[,2]~m.R26.102219.u$`35`$df[,3])
summary(lm(m.R26.102219.u$`35`$df[,2]~m.R26.102219.u$`35`$df[,3]))
plot(m.R26.102219.u$`40`$df[,2]~m.R26.102219.u$`40`$df[,3])
summary(lm(m.R26.102219.u$`40`$df[,2]~m.R26.102219.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.R26.Rd<-rep(NA,length(R26.112219$Curve))
for (i in 1:length(R26.112219$Curve)){
  p.R26.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,R26.112219$Curve[i])
}
R26.112219$Rd<-p.R26.Rd

#Subset data depending on if they were collected in an ascending or descending order
R26.112219.d<-R26.112219[R26.112219$Direction == "down",]
R26.112219.u<-R26.112219[R26.112219$Direction == "up",]

#
#Fit curves
#

#Descending
m.R26.112219.d<-fitacis(R26.112219.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26.112219.d)
coef(m.R26.112219.d)

#Ascending
m.R26.112219.u<-fitacis(R26.112219.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26.112219.u)
coef(m.R26.112219.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R26.112219.d$`10`$df[,2]~m.R26.112219.d$`10`$df[,3])
summary(lm(m.R26.112219.d$`10`$df[,2]~m.R26.112219.d$`10`$df[,3]))
plot(m.R26.112219.d$`15`$df[,2]~m.R26.112219.d$`15`$df[,3])
summary(lm(m.R26.112219.d$`15`$df[,2]~m.R26.112219.d$`15`$df[,3]))
plot(m.R26.112219.d$`20`$df[,2]~m.R26.112219.d$`20`$df[,3])
summary(lm(m.R26.112219.d$`20`$df[,2]~m.R26.112219.d$`20`$df[,3]))
plot(m.R26.112219.d$`25`$df[,2]~m.R26.112219.d$`25`$df[,3])
summary(lm(m.R26.112219.d$`25`$df[,2]~m.R26.112219.d$`25`$df[,3]))

#Ascending
plot(m.R26.112219.u$`25`$df[,2]~m.R26.112219.u$`25`$df[,3])
summary(lm(m.R26.112219.u$`25`$df[,2]~m.R26.112219.u$`25`$df[,3]))
plot(m.R26.112219.u$`30`$df[,2]~m.R26.112219.u$`30`$df[,3])
summary(lm(m.R26.112219.u$`30`$df[,2]~m.R26.112219.u$`30`$df[,3]))
plot(m.R26.112219.u$`35`$df[,2]~m.R26.112219.u$`35`$df[,3])
summary(lm(m.R26.112219.u$`35`$df[,2]~m.R26.112219.u$`35`$df[,3]))
plot(m.R26.112219.u$`40`$df[,2]~m.R26.112219.u$`40`$df[,3])
summary(lm(m.R26.112219.u$`40`$df[,2]~m.R26.112219.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.R26.Rd<-rep(NA,length(R26.120319$Curve))
for (i in 1:length(R26.120319$Curve)){
  p.R26.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,R26.120319$Curve[i])
}
R26.120319$Rd<-p.R26.Rd

#Subset data depending on if they were collected in an ascending or descending order
R26.120319.d<-R26.120319[R26.120319$Direction == "down",] 
R26.120319.u<-R26.120319[R26.120319$Direction == "up",]

#
#Fit curves
#

#Descending
m.R26.120319.d<-fitacis(R26.120319.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26.120319.d)
coef(m.R26.120319.d) 

#Ascending
m.R26.120319.u<-fitacis(R26.120319.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26.120319.u)
coef(m.R26.120319.u) 

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R26.120319.d$`10`$df[,2]~m.R26.120319.d$`10`$df[,3])
summary(lm(m.R26.120319.d$`10`$df[,2]~m.R26.120319.d$`10`$df[,3]))
plot(m.R26.120319.d$`15`$df[,2]~m.R26.120319.d$`15`$df[,3])
summary(lm(m.R26.120319.d$`15`$df[,2]~m.R26.120319.d$`15`$df[,3]))
plot(m.R26.120319.d$`20`$df[,2]~m.R26.120319.d$`20`$df[,3])
summary(lm(m.R26.120319.d$`20`$df[,2]~m.R26.120319.d$`20`$df[,3]))
plot(m.R26.120319.d$`25`$df[,2]~m.R26.120319.d$`25`$df[,3])
summary(lm(m.R26.120319.d$`25`$df[,2]~m.R26.120319.d$`25`$df[,3]))

#Ascending
plot(m.R26.120319.u$`25`$df[,2]~m.R26.120319.u$`25`$df[,3])
summary(lm(m.R26.120319.u$`25`$df[,2]~m.R26.120319.u$`25`$df[,3]))
plot(m.R26.120319.u$`30`$df[,2]~m.R26.120319.u$`30`$df[,3])
summary(lm(m.R26.120319.u$`30`$df[,2]~m.R26.120319.u$`30`$df[,3]))
plot(m.R26.120319.u$`35`$df[,2]~m.R26.120319.u$`35`$df[,3])
summary(lm(m.R26.120319.u$`35`$df[,2]~m.R26.120319.u$`35`$df[,3]))
plot(m.R26.120319.u$`40`$df[,2]~m.R26.120319.u$`40`$df[,3])
summary(lm(m.R26.120319.u$`40`$df[,2]~m.R26.120319.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.R31.Rd<-rep(NA,length(R31.071318$Curve))
for (i in 1:length(R31.071318$Curve)){
  p.R31.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,R31.071318$Curve[i])
}
R31.071318$Rd<-p.R31.Rd

#Subset data depending on if they were collected in an ascending or descending order
R31.071318.d<-R31.071318[R31.071318$Direction == "down",]
R31.071318.u<-R31.071318[R31.071318$Direction == "up",]

#
#Fit curves
#

#Descending
m.R31.071318.d<-fitacis(R31.071318.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31.071318.d)
coef(m.R31.071318.d)

#Ascending
m.R31.071318.u<-fitacis(R31.071318.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31.071318.u)
coef(m.R31.071318.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R31.071318.d$`10`$df[,2]~m.R31.071318.d$`10`$df[,3])
summary(lm(m.R31.071318.d$`10`$df[,2]~m.R31.071318.d$`10`$df[,3]))
plot(m.R31.071318.d$`15`$df[,2]~m.R31.071318.d$`15`$df[,3])
summary(lm(m.R31.071318.d$`15`$df[,2]~m.R31.071318.d$`15`$df[,3]))
plot(m.R31.071318.d$`20`$df[,2]~m.R31.071318.d$`20`$df[,3])
summary(lm(m.R31.071318.d$`20`$df[,2]~m.R31.071318.d$`20`$df[,3]))
plot(m.R31.071318.d$`25`$df[,2]~m.R31.071318.d$`25`$df[,3])
summary(lm(m.R31.071318.d$`25`$df[,2]~m.R31.071318.d$`25`$df[,3]))
plot(m.R31.071318.d$`30`$df[,2]~m.R31.071318.d$`30`$df[,3])
summary(lm(m.R31.071318.d$`30`$df[,2]~m.R31.071318.d$`30`$df[,3]))

#Ascending
plot(m.R31.071318.u$`30`$df[,2]~m.R31.071318.u$`30`$df[,3])
summary(lm(m.R31.071318.u$`30`$df[,2]~m.R31.071318.u$`30`$df[,3]))
plot(m.R31.071318.u$`35`$df[,2]~m.R31.071318.u$`35`$df[,3])
summary(lm(m.R31.071318.u$`35`$df[,2]~m.R31.071318.u$`35`$df[,3]))
plot(m.R31.071318.u$`40`$df[,2]~m.R31.071318.u$`40`$df[,3])
summary(lm(m.R31.071318.u$`40`$df[,2]~m.R31.071318.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.R31.Rd<-rep(NA,length(R31.072318$Curve))
for (i in 1:length(R31.072318$Curve)){
  p.R31.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,R31.072318$Curve[i])
}
R31.072318$Rd<-p.R31.Rd

#Subset data depending on if they were collected in an ascending or descending order
R31.072318.d<-R31.072318[R31.072318$Direction == "down",]
R31.072318.u<-R31.072318[R31.072318$Direction == "up",]

#
#Fit curves
#

#Descending
m.R31.072318.d<-fitacis(R31.072318.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31.072318.d)
coef(m.R31.072318.d)

#Ascending
m.R31.072318.u<-fitacis(R31.072318.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31.072318.u)
coef(m.R31.072318.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R31.072318.d$`10`$df[,2]~m.R31.072318.d$`10`$df[,3])
summary(lm(m.R31.072318.d$`10`$df[,2]~m.R31.072318.d$`10`$df[,3]))
plot(m.R31.072318.d$`15`$df[,2]~m.R31.072318.d$`15`$df[,3])
summary(lm(m.R31.072318.d$`15`$df[,2]~m.R31.072318.d$`15`$df[,3]))
plot(m.R31.072318.d$`20`$df[,2]~m.R31.072318.d$`20`$df[,3])
summary(lm(m.R31.072318.d$`20`$df[,2]~m.R31.072318.d$`20`$df[,3]))
plot(m.R31.072318.d$`25`$df[,2]~m.R31.072318.d$`25`$df[,3])
summary(lm(m.R31.072318.d$`25`$df[,2]~m.R31.072318.d$`25`$df[,3]))
plot(m.R31.072318.d$`30`$df[,2]~m.R31.072318.d$`30`$df[,3])
summary(lm(m.R31.072318.d$`30`$df[,2]~m.R31.072318.d$`30`$df[,3]))

#Ascending
plot(m.R31.072318.u$`30`$df[,2]~m.R31.072318.u$`30`$df[,3])
summary(lm(m.R31.072318.u$`30`$df[,2]~m.R31.072318.u$`30`$df[,3]))
plot(m.R31.072318.u$`35`$df[,2]~m.R31.072318.u$`35`$df[,3])
summary(lm(m.R31.072318.u$`35`$df[,2]~m.R31.072318.u$`35`$df[,3]))
plot(m.R31.072318.u$`40`$df[,2]~m.R31.072318.u$`40`$df[,3])
summary(lm(m.R31.072318.u$`40`$df[,2]~m.R31.072318.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.R31.Rd<-rep(NA,length(R31.101418$Curve))
for (i in 1:length(R31.101418$Curve)){
  p.R31.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,R31.101418$Curve[i])
}
R31.101418$Rd<-p.R31.Rd

#Subset data depending on if they were collected in an ascending or descending order
R31.101418.d<-R31.101418[R31.101418$Direction == "down",]
R31.101418.u<-R31.101418[R31.101418$Direction == "up",]

#
#Fit curves
#

#Descending
m.R31.101418.d<-fitacis(R31.101418.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31.101418.d)
coef(m.R31.101418.d)

#Ascending
m.R31.101418.u<-fitacis(R31.101418.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31.101418.u)
coef(m.R31.101418.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R31.101418.d$`10`$df[,2]~m.R31.101418.d$`10`$df[,3])
summary(lm(m.R31.101418.d$`10`$df[,2]~m.R31.101418.d$`10`$df[,3]))
plot(m.R31.101418.d$`15`$df[,2]~m.R31.101418.d$`15`$df[,3])
summary(lm(m.R31.101418.d$`15`$df[,2]~m.R31.101418.d$`15`$df[,3]))
plot(m.R31.101418.d$`20`$df[,2]~m.R31.101418.d$`20`$df[,3])
summary(lm(m.R31.101418.d$`20`$df[,2]~m.R31.101418.d$`20`$df[,3]))
plot(m.R31.101418.d$`25`$df[,2]~m.R31.101418.d$`25`$df[,3])
summary(lm(m.R31.101418.d$`25`$df[,2]~m.R31.101418.d$`25`$df[,3]))
plot(m.R31.101418.d$`30`$df[,2]~m.R31.101418.d$`30`$df[,3])
summary(lm(m.R31.101418.d$`30`$df[,2]~m.R31.101418.d$`30`$df[,3]))

#Ascending
plot(m.R31.101418.u$`30`$df[,2]~m.R31.101418.u$`30`$df[,3])
summary(lm(m.R31.101418.u$`30`$df[,2]~m.R31.101418.u$`30`$df[,3]))
plot(m.R31.101418.u$`35`$df[,2]~m.R31.101418.u$`35`$df[,3])
summary(lm(m.R31.101418.u$`35`$df[,2]~m.R31.101418.u$`35`$df[,3]))
plot(m.R31.101418.u$`40`$df[,2]~m.R31.101418.u$`40`$df[,3])
summary(lm(m.R31.101418.u$`40`$df[,2]~m.R31.101418.u$`40`$df[,3]))

###
#Data which were removed
###

#Removed because R2 of fit <0.9
#m.A26.020820.u at 40 deg. C
#m.G21.082019.d all temperatures
#m.G26.110119.u 40 deg. C
#m.R21.073118.d at 10 deg. C
#m.R21.073118.u at 35 deg. C

#Removed based on Vcmax SE >40%
#m.M26.091119.u at 40 deg. C
#m.A26.020420.u at 40 deg. C
#m.A26.020820.u at 40 deg. C
#m.A31.071618.u at 40 deg. C
#m.G21.082019.u at 30 deg. C
#m.G26.110119.u at 40 deg. C (Vcmax SE is NaN)
#m.G31.071118.u at 40 deg. C

#Removed because they were not fit
#m.M26.091119.d Jmax at 10 deg. C
#m.M31.080218.u Jmax at 40 deg. C
#m.M31.120418.d Jmax at 10 deg. C
#m.A21.071218.d Jmax at 15 deg. C
#m.A21.071818.d Jmax at 15 deg. C
#m.G21.082019.u Jmax at 25 deg. C
#m.G26.090619.d Jmax at 10 deg. C
#m.G26.100319.d Jmax at 10 deg. C
#m.G31.072218.d Jmax at 10 and 15 deg. C
#m.R21.101018.d Jmax at 10 deg. C
#m.R26.102219.d Jmax at 10 deg. C
#m.R26.112219.d Jmax at 10 deg. C

###
#Extract A275 estimates from each A-Ci fit
###

##
#Morella
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
M21.071418.A275<-c(m.M21.071418.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.u$`20`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M21.071418.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
M21.120918.A275<-c(m.M21.120918.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.u$`20`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M21.120918.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
M21.012119.A275<-c(m.M21.012119.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.u$`20`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M21.012119.u$`40`$Photosyn(Ci=275)[[2]])

M21.A275<-c(M21.071418.A275,M21.120918.A275,M21.012119.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
M26.083019.A275<-c(m.M26.083019.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M26.083019.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
M26.091119.A275<-c(m.M26.091119.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M26.091119.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M26.091119.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M26.091119.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M26.091119.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M26.091119.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M26.091119.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 3
M26.092119.A275<-c(m.M26.092119.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M26.092119.u$`40`$Photosyn(Ci=275)[[2]])

M26.A275<-c(M26.083019.A275,M26.091119.A275,M26.092119.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
M31.080218.A275<-c(m.M31.080218.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.d$`30`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M31.080218.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
M31.120418.A275<-c(m.M31.120418.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M31.120418.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M31.120418.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M31.120418.d$`30`$Photosyn(Ci=275)[[2]],
                   m.M31.120418.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M31.120418.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M31.120418.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
M31.040919.A275<-c(m.M31.040919.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.d$`30`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M31.040919.u$`40`$Photosyn(Ci=275)[[2]])

M31.A275<-c(M31.080218.A275,M31.120418.A275,M31.040919.A275)

##
#Alnus
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
A21.071218.A275<-c(m.A21.071218.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.u$`20`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A21.071218.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
A21.071818.A275<-c(m.A21.071818.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.u$`20`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A21.071818.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
A21.101718.A275<-c(m.A21.101718.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.u$`20`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A21.101718.u$`40`$Photosyn(Ci=275)[[2]])

A21.A275<-c(A21.071218.A275,A21.071818.A275,A21.101718.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
A26.112419.A275<-c(m.A26.112419.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A26.112419.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
A26.020420.A275<-c(m.A26.020420.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A26.020420.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A26.020420.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A26.020420.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A26.020420.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A26.020420.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A26.020420.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 3
A26.020820.A275<-c(m.A26.020820.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A26.020820.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A26.020820.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A26.020820.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A26.020820.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A26.020820.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A26.020820.u$`35`$Photosyn(Ci=275)[[2]])

A26.A275<-c(A26.112419.A275,A26.020420.A275,A26.020820.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
A31.071618.A275<-c(m.A31.071618.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A31.071618.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A31.071618.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A31.071618.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A31.071618.d$`30`$Photosyn(Ci=275)[[2]],
                   m.A31.071618.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A31.071618.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 2
A31.092918.A275<-c(m.A31.092918.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.d$`30`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A31.092918.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
A31.112518.A275<-c(m.A31.112518.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.d$`30`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A31.112518.u$`40`$Photosyn(Ci=275)[[2]])

A31.A275<-c(A31.071618.A275,A31.092918.A275,A31.112518.A275)

##
#Gliricidia
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
G21.080318.A275<-c(m.G21.080318.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21.080318.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
G21.082019.A275<-c(m.G21.082019.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.082019.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21.082019.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
G21.082819.A275<-c(m.G21.082819.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21.082819.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 4
G21.112718.A275<-c(m.G21.112718.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21.112718.u$`40`$Photosyn(Ci=275)[[2]])

G21.A275<-c(G21.080318.A275,G21.112718.A275,G21.082019.A275,G21.082819.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
G26.082319.A275<-c(m.G26.082319.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G26.082319.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
G26.090619.A275<-c(m.G26.090619.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G26.090619.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
G26.100319.A275<-c(m.G26.100319.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G26.100319.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 4
G26.110119.A275<-c(m.G26.110119.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26.110119.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26.110119.u$`35`$Photosyn(Ci=275)[[2]])

G26.A275<-c(G26.082319.A275,G26.090619.A275,G26.100319.A275,G26.110119.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
G31.050418.A275<-c(m.G31.050418.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G31.050418.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G31.050418.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G31.050418.d$`30`$Photosyn(Ci=275)[[2]],
                   m.G31.050418.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G31.050418.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G31.050418.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
G31.071118.A275<-c(m.G31.071118.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G31.071118.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G31.071118.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G31.071118.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G31.071118.d$`30`$Photosyn(Ci=275)[[2]],
                   m.G31.071118.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G31.071118.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 3
G31.072218.A275<-c(m.G31.072218.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.d$`30`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.d$`32`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G31.072218.u$`40`$Photosyn(Ci=275)[[2]])

G31.A275<-c(G31.050418.A275,G31.071118.A275,G31.072218.A275)

##
#Robinia
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
R21.073118.A275<-c(m.R21.073118.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21.073118.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.073118.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.073118.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21.073118.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21.073118.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
R21.100518.A275<-c(m.R21.100518.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R21.100518.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
R21.101018.A275<-c(m.R21.101018.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R21.101018.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 4
R21.091719.A275<-c(m.R21.091719.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R21.091719.u$`40`$Photosyn(Ci=275)[[2]])

R21.A275<-c(R21.073118.A275,R21.100518.A275,R21.101018.A275,R21.091719.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
R26.102219.A275<-c(m.R26.102219.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R26.102219.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
R26.112219.A275<-c(m.R26.112219.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.d$`40`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R26.112219.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
R26.120319.A275<-c(m.R26.120319.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R26.120319.u$`40`$Photosyn(Ci=275)[[2]])

R26.A275<-c(R26.102219.A275,R26.112219.A275,R26.120319.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
R31.071318.A275<-c(m.R31.071318.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.d$`30`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R31.071318.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
R31.072318.A275<-c(m.R31.072318.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.d$`30`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R31.072318.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
R31.101418.A275<-c(m.R31.101418.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.d$`30`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R31.101418.u$`40`$Photosyn(Ci=275)[[2]])

R31.A275<-c(R31.071318.A275,R31.072318.A275,R31.101418.A275)

####
#The following assembles the data frames that are used to model A275, Vcmax, and Jmax as functions of temperature
####

###
#Morella
###

##
#21:15 deg. C growing temperature
##

dat.V.M21<-cbind(c(coef(m.M21.071418.d)[,2],coef(m.M21.071418.u)[,2],coef(m.M21.120918.d)[,2],
                   coef(m.M21.120918.u)[,2],coef(m.M21.012119.d)[,2],coef(m.M21.012119.u)[,2]))
dat.J.M21<-cbind(c(coef(m.M21.071418.d)[,3],coef(m.M21.071418.u)[,3],coef(m.M21.120918.d)[,3],
                   coef(m.M21.120918.u)[,3],coef(m.M21.012119.d)[,3],coef(m.M21.012119.u)[,3]))
dat.Tc.M21<-cbind(c(10,15,20,20,25,30,35,40,10,15,20,20,25,30,35,40,10,15,20,20,25,30,35,40))
dat.Dir.M21<-cbind(c("down","down","down","up","up","up","up","up","down","down","down","up","up","up","up","up",
                     "down","down","down","up","up","up","up","up"))
dat.Tr.M21<-rep("M21",24)
dat.Rep.M21<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
dat.L.M21<-c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6)

dat.M21<-as.data.frame(cbind(dat.Tr.M21,dat.Rep.M21,dat.Dir.M21,dat.L.M21,dat.Tc.M21,dat.V.M21,dat.J.M21,M21.A275))
colnames(dat.M21)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.M21,"ACi.dat.M21.csv")

##
#26:25 deg. C growing temperature
##

dat.V.M26<-cbind(c(coef(m.M26.083019.d)[,2],coef(m.M26.083019.u)[,2],
                   coef(m.M26.091119.d)[,2],coef(m.M26.091119.u)[1:3,2],
                   coef(m.M26.092119.d)[,2],coef(m.M26.092119.u)[,2]))
dat.J.M26<-cbind(c(coef(m.M26.083019.d)[,3],coef(m.M26.083019.u)[,3],
                   NA,coef(m.M26.091119.d)[2:4,3],coef(m.M26.091119.u)[1:3,3],
                   coef(m.M26.092119.d)[,3],coef(m.M26.092119.u)[,3]))
dat.Tc.M26<-cbind(c(10,15,20,25,25,30,35,40,
                    10,15,20,25,25,30,35,
                    10,15,20,25,25,30,35,40))
dat.Dir.M26<-cbind(c("down","down","down","down","up","up","up","up",
                     "down","down","down","down","up","up","up",
                     "down","down","down","down","up","up","up","up"))
dat.Tr.M26<-rep("M26",23)
dat.Rep.M26<-c(1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3)
dat.L.M26<-c(1,1,1,1,2,2,2,2,
             3,3,3,3,4,4,4,
             5,5,5,5,6,6,6,6)

dat.M26<-as.data.frame(cbind(dat.Tr.M26,dat.Rep.M26,dat.Dir.M26,dat.L.M26,dat.Tc.M26,dat.V.M26,dat.J.M26,M26.A275))
colnames(dat.M26)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.M26,"ACi.dat.M26.csv")

##
#31:25 deg. C growing temperature
##

dat.V.M31<-cbind(c(coef(m.M31.080218.d)[,2],coef(m.M31.080218.u)[,2],coef(m.M31.120418.d)[,2],
                   coef(m.M31.120418.u)[,2],coef(m.M31.040919.d)[,2],coef(m.M31.040919.u)[,2]))
dat.J.M31<-cbind(c(coef(m.M31.080218.d)[,3],coef(m.M31.080218.u)[1:2,3],NA,NA,coef(m.M31.120418.d)[2:4,3],
                   coef(m.M31.120418.u)[,3],coef(m.M31.040919.d)[,3],coef(m.M31.040919.u)[,3]))
dat.Tc.M31<-cbind(c(10,15,20,25,30,30,35,40,15,20,25,30,30,35,40,10,15,20,25,30,30,35,40))
dat.Dir.M31<-cbind(c("down","down","down","down","down","up","up","up","down","down","down","down","up","up","up",
                     "down","down","down","down","down","up","up","up"))
dat.Tr.M31<-rep("M31",23)
dat.Rep.M31<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
dat.L.M31<-c(1,1,1,1,1,2,2,2,3,3,3,3,4,4,4,5,5,5,5,5,6,6,6)

dat.M31<-as.data.frame(cbind(dat.Tr.M31,dat.Rep.M31,dat.Dir.M31,dat.L.M31,dat.Tc.M31,dat.V.M31,dat.J.M31,M31.A275))
colnames(dat.M31)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.M31,"ACi.dat.M31.csv")

###
#Alnus
###

##
#21:15 deg. C growing temperature
##

dat.V.A21<-cbind(c(coef(m.A21.071218.d)[,2],coef(m.A21.071218.u)[,2],coef(m.A21.071818.d)[,2],
                   coef(m.A21.071818.u)[,2],coef(m.A21.101718.d)[,2],coef(m.A21.101718.u)[,2]))
dat.J.A21<-cbind(c(coef(m.A21.071218.d)[1,3],NA,coef(m.A21.071218.d)[3,3],coef(m.A21.071218.u)[,3],NA,coef(m.A21.071818.d)[2:3,3],
                   coef(m.A21.071818.u)[,3],coef(m.A21.101718.d)[,3],coef(m.A21.101718.u)[,3]))
dat.Tc.A21<-cbind(c(10,15,20,20,25,30,35,40,10,15,20,20,25,30,35,40,10,15,20,20,25,30,35,40))
dat.Dir.A21<-cbind(c("down","down","down","up","up","up","up","up","down","down","down","up","up","up","up","up",
                     "down","down","down","up","up","up","up","up"))
dat.Tr.A21<-rep("A21",24)
dat.Rep.A21<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
dat.L.A21<-c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5,6,6,6,6,6)

dat.A21<-as.data.frame(cbind(dat.Tr.A21,dat.Rep.A21,dat.Dir.A21,dat.L.A21,dat.Tc.A21,dat.V.A21,dat.J.A21,A21.A275))
colnames(dat.A21)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.A21,"ACi.dat.A21.csv")

##
#26:20 deg. C growing temperature
##

dat.V.A26<-cbind(c(coef(m.A26.112419.d)[,2],coef(m.A26.112419.u)[,2],
                   coef(m.A26.020420.d)[,2],coef(m.A26.020420.u)[1:3,2],
                   coef(m.A26.020820.d)[,2],coef(m.A26.020820.u)[1:3,2]))
dat.J.A26<-cbind(c(coef(m.A26.112419.d)[,3],coef(m.A26.112419.u)[,3],
                   coef(m.A26.020420.d)[,3],coef(m.A26.020420.u)[1:3,3],
                   coef(m.A26.020820.d)[,3],coef(m.A26.020820.u)[1:3,3]))
dat.Tc.A26<-cbind(c(10,15,20,25,25,30,35,40,
                    10,15,20,25,25,30,35,
                    10,15,20,25,25,30,35))
dat.Dir.A26<-cbind(c("down","down","down","down","up","up","up","up",
                     "down","down","down","down","up","up","up",
                     "down","down","down","down","up","up","up"))
dat.Tr.A26<-rep("A26",22)
dat.Rep.A26<-c(1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,
               3,3,3,3,3,3,3)
dat.L.A26<-c(1,1,1,1,2,2,2,2,
             3,3,3,3,4,4,4,
             5,5,5,5,6,6,6)

dat.A26<-as.data.frame(cbind(dat.Tr.A26,dat.Rep.A26,dat.Dir.A26,dat.L.A26,dat.Tc.A26,dat.V.A26,dat.J.A26,A26.A275))
colnames(dat.A26)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.A26,"ACi.dat.A26.csv")

##
#31:25 deg. C growing temperature
##

dat.V.A31<-cbind(c(coef(m.A31.071618.d)[,2],coef(m.A31.071618.u)[1:2,2],coef(m.A31.092918.d)[,2],
                   coef(m.A31.092918.u)[,2],coef(m.A31.112518.d)[,2],coef(m.A31.112518.u)[,2]))
dat.J.A31<-cbind(c(coef(m.A31.071618.d)[,3],coef(m.A31.071618.u)[1:2,3],coef(m.A31.092918.d)[,3],
                   coef(m.A31.092918.u)[,3],coef(m.A31.112518.d)[,3],coef(m.A31.112518.u)[,3]))
dat.Tc.A31<-cbind(c(10,15,20,25,30,30,35,10,15,20,25,30,30,35,40,10,15,20,25,30,30,35,40))
dat.Dir.A31<-cbind(c("down","down","down","down","down","up","up","down","down","down","down","down","up","up","up",
                     "down","down","down","down","down","up","up","up"))
dat.Tr.A31<-rep("A31",23)
dat.Rep.A31<-c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3)
dat.L.A31<-c(1,1,1,1,1,2,2,3,3,3,3,3,4,4,4,5,5,5,5,5,6,6,6)

dat.A31<-as.data.frame(cbind(dat.Tr.A31,dat.Rep.A31,dat.Dir.A31,dat.L.A31,dat.Tc.A31,dat.V.A31,dat.J.A31,A31.A275))
colnames(dat.A31)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.A31,"ACi.dat.A31.csv")

###
#Gliricidia
###

##
#21:15 deg. C growing temperature
##

dat.V.G21<-cbind(c(coef(m.G21.080318.d)[,2],coef(m.G21.080318.u)[,2],coef(m.G21.112718.d)[,2],
                   coef(m.G21.112718.u)[,2],coef(m.G21.082019.u)[1,2],coef(m.G21.082019.u)[4:5,2],coef(m.G21.082819.d)[,2],coef(m.G21.082819.u)[,2]))
dat.J.G21<-cbind(c(coef(m.G21.080318.d)[,3],coef(m.G21.080318.u)[,3],coef(m.G21.112718.d)[,3],
                   coef(m.G21.112718.u)[,3],coef(m.G21.082019.u)[1,3],coef(m.G21.082019.u)[4:5,3],coef(m.G21.082819.d)[,3],coef(m.G21.082819.u)[,3]))
dat.Tc.G21<-cbind(c(10,15,20,20,25,30,35,40,10,15,20,20,25,30,35,40,20,35,40,10,15,20,20,25,30,35,40))
dat.Dir.G21<-cbind(c("down","down","down","up","up","up","up","up","down","down","down","up","up","up","up","up",
                     "up","up","up","down","down","down","up","up","up","up","up"))
dat.Tr.G21<-rep("G21",27)
dat.Rep.G21<-c(1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,3,3,3,4,4,4,4,4,4,4,4)
dat.L.G21<-c(1,1,1,2,2,2,2,2,3,3,3,4,4,4,4,4,5,5,5,6,6,6,7,7,7,7,7)

dat.G21<-as.data.frame(cbind(dat.Tr.G21,dat.Rep.G21,dat.Dir.G21,dat.L.G21,dat.Tc.G21,dat.V.G21,dat.J.G21,G21.A275))
colnames(dat.G21)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.G21,"ACi.dat.G21.csv")

##
#26:20 deg. C growing temperature
##

dat.V.G26<-cbind(c(coef(m.G26.082319.d)[,2],coef(m.G26.082319.u)[,2],
                   coef(m.G26.090619.d)[,2],coef(m.G26.090619.u)[,2],
                   coef(m.G26.100319.d)[,2],coef(m.G26.100319.u)[,2],
                   coef(m.G26.110119.u)[1:3,2]))
dat.J.G26<-cbind(c(coef(m.G26.082319.d)[,3],coef(m.G26.082319.u)[,3],
                   NA,coef(m.G26.090619.d)[2:4,3],coef(m.G26.090619.u)[,3],
                   NA,coef(m.G26.100319.d)[2:4,3],coef(m.G26.100319.u)[,3],
                   coef(m.G26.110119.u)[1:3,3]))
dat.Tc.G26<-cbind(c(10,15,20,25,25,30,35,40,
                    10,15,20,25,25,30,35,40,
                    10,15,20,25,25,30,35,40,
                    25,30,35))
dat.Dir.G26<-cbind(c("down","down","down","down","up","up","up","up",
                     "down","down","down","down","up","up","up","up",
                     "down","down","down","down","up","up","up","up",
                     "up","up","up"))
dat.Tr.G26<-rep("G26",27)
dat.Rep.G26<-c(1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3,
               4,4,4)
dat.L.G26<-c(1,1,1,1,2,2,2,2,
             3,3,3,3,4,4,4,4,
             5,5,5,5,6,6,6,6,
             7,7,7)

dat.G26<-as.data.frame(cbind(dat.Tr.G26,dat.Rep.G26,dat.Dir.G26,dat.L.G26,dat.Tc.G26,dat.V.G26,dat.J.G26,G26.A275))
colnames(dat.G26)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.G26,"ACi.dat.G26.csv")

##
#31:25 deg. C growing temperature
##

dat.V.G31<-cbind(c(coef(m.G31.050418.d)[,2],coef(m.G31.050418.u)[,2],
                   coef(m.G31.071118.d)[,2],coef(m.G31.071118.u)[1:2,2],
                   coef(m.G31.072218.d)[,2],coef(m.G31.072218.u)[,2]))
dat.J.G31<-cbind(c(coef(m.G31.050418.d)[,3],coef(m.G31.050418.u)[,3],
                   coef(m.G31.071118.d)[,3],coef(m.G31.071118.u)[1:2,3],
                   NA,NA,coef(m.G31.072218.d)[3:6,3],coef(m.G31.072218.u)[,3]))
dat.Tc.G31<-cbind(c(15,20,25,30,30,35,40,10,15,20,25,30,30,35,10,15,20,25,30,32.4,30,35,40))
dat.Dir.G31<-cbind(c("down","down","down","down","up","up","up","down","down","down","down","down","up","up",
                     "down","down","down","down","down","down","up","up","up"))
dat.Tr.G31<-rep("G31",23)
dat.Rep.G31<-c(1,1,1,1,1,1,1,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3)
dat.L.G31<-c(1,1,1,1,2,2,2,3,3,3,3,3,4,4,5,5,5,5,5,5,6,6,6)

dat.G31<-as.data.frame(cbind(dat.Tr.G31,dat.Rep.G31,dat.Dir.G31,dat.L.G31,dat.Tc.G31,dat.V.G31,dat.J.G31,G31.A275))
colnames(dat.G31)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.G31,"ACi.dat.G31.csv")

###
#Robinia
###

##
#21:15 deg. C growing temperature
##

dat.V.R21<-cbind(c(coef(m.R21.073118.d)[2:3,2],coef(m.R21.073118.u)[1:3,2],coef(m.R21.073118.u)[5,2],
                   coef(m.R21.100518.d)[,2],coef(m.R21.100518.u)[,2],
                   coef(m.R21.101018.d)[,2],coef(m.R21.101018.u)[,2],
                   coef(m.R21.091719.d)[,2],coef(m.R21.091719.u)[,2]))
dat.J.R21<-cbind(c(coef(m.R21.073118.d)[2:3,3],coef(m.R21.073118.u)[1:3,3],coef(m.R21.073118.u)[5,3],
                   coef(m.R21.100518.d)[,3],coef(m.R21.100518.u)[,3],
                   NA,coef(m.R21.101018.d)[2:3,2],coef(m.R21.101018.u)[,3],
                   coef(m.R21.091719.d)[,3],coef(m.R21.091719.u)[,3]))
dat.Tc.R21<-cbind(c(15,20,20,25,30,40,
                    10,15,20,20,25,30,35,40,
                    10,15,20,20,25,30,35,40,
                    10,15,20,20,25,30,35,40))
dat.Dir.R21<-cbind(c("down","down","up","up","up","up",
                     "down","down","down","up","up","up","up","up",
                     "down","down","down","up","up","up","up","up",
                     "down","down","down","up","up","up","up","up"))
dat.Tr.R21<-rep("R21",30)
dat.Rep.R21<-c(rep(1,6),rep(2,8),rep(3,8),rep(4,8))
dat.L.R21<-c(1,1,2,2,2,2,
             3,3,3,4,4,4,4,4,
             5,5,5,6,6,6,6,6,
             7,7,7,8,8,8,8,8)

dat.R21<-as.data.frame(cbind(dat.Tr.R21,dat.Rep.R21,dat.Dir.R21,dat.L.R21,dat.Tc.R21,dat.V.R21,dat.J.R21,R21.A275))
colnames(dat.R21)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.R21,"ACi.dat.R21.csv")

##
#26:20 deg. C growing temperature
##

dat.V.R26<-cbind(c(coef(m.R26.102219.d)[,2],coef(m.R26.102219.u)[,2],
                   coef(m.R26.112219.d)[,2],coef(m.R26.112219.u)[,2],
                   coef(m.R26.120319.d)[,2],coef(m.R26.120319.u)[,2]))
dat.J.R26<-cbind(c(NA,coef(m.R26.102219.d)[2:4,2],coef(m.R26.102219.u)[,3],
                   NA,coef(m.R26.112219.d)[2:5,3],coef(m.R26.112219.u)[,3],
                   coef(m.R26.120319.d)[,3],coef(m.R26.120319.u)[,3]))
dat.Tc.R26<-cbind(c(10,15,20,25,25,30,35,40,
                    10,15,20,25,40,25,30,35,40,
                    10,15,20,25,25,30,35,40))
dat.Dir.R26<-cbind(c("down","down","down","down","up","up","up","up",
                     "down","down","down","down","down","up","up","up","up",
                     "down","down","down","down","up","up","up","up"))
dat.Tr.R26<-rep("R26",25)
dat.Rep.R26<-c(1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3)
dat.L.R26<-c(1,1,1,1,2,2,2,2,
             3,3,3,3,3,4,4,4,4,
             5,5,5,5,6,6,6,6)

dat.R26<-as.data.frame(cbind(dat.Tr.R26,dat.Rep.R26,dat.Dir.R26,dat.L.R26,dat.Tc.R26,dat.V.R26,dat.J.R26,R26.A275))
colnames(dat.R26)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.R26,"ACi.dat.R26.csv")

##
#31:25 deg. C growing temperature
##

dat.V.R31<-cbind(c(coef(m.R31.071318.d)[,2],coef(m.R31.071318.u)[,2],coef(m.R31.072318.d)[,2],
                   coef(m.R31.072318.u)[,2],coef(m.R31.101418.d)[,2],coef(m.R31.101418.u)[,2]))
dat.J.R31<-cbind(c(coef(m.R31.071318.d)[,3],coef(m.R31.071318.u)[,3],coef(m.R31.072318.d)[,3],
                   coef(m.R31.072318.u)[,3],coef(m.R31.101418.d)[,3],coef(m.R31.101418.u)[,3]))
dat.Tc.R31<-cbind(c(10,15,20,25,30,30,35,40,
                    10,15,20,25,30,30,35,40,
                    10,15,20,25,30,30,35,40))
dat.Dir.R31<-cbind(c("down","down","down","down","down","up","up","up",
                     "down","down","down","down","down","up","up","up",
                     "down","down","down","down","down","up","up","up"))
dat.Tr.R31<-rep("R31",24)
dat.Rep.R31<-c(1,1,1,1,1,1,1,1,
               2,2,2,2,2,2,2,2,
               3,3,3,3,3,3,3,3)
dat.L.R31<-c(1,1,1,1,1,2,2,2,
             3,3,3,3,3,4,4,4,
             5,5,5,5,5,6,6,6)

dat.R31<-as.data.frame(cbind(dat.Tr.R31,dat.Rep.R31,dat.Dir.R31,dat.L.R31,dat.Tc.R31,dat.V.R31,dat.J.R31,R31.A275))
colnames(dat.R31)<-c("Tr","Plant","Dir","Leaf","Temp","Vcmax","Jmax","A275")
write.csv(dat.R31,"ACi.dat.R31.csv")