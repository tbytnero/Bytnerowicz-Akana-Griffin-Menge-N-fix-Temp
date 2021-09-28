###############################################################################################################
###############################################################################################################
#This script fits the A-Ci curves in order to calculate A275, Vcmax, and Jmax
###############################################################################################################
###############################################################################################################

#Run Resp_Calculations.R first

#Load necessary package
library(plantecophys)

#####
#Read in data
#####

M21a<-read.csv("MOCE21a_ACi.csv")
M21b<-read.csv("MOCE21b_ACi.csv")
M21c<-read.csv("MOCE21c_ACi.csv")
M26a<-read.csv("MOCE26a_ACi.csv")
M26b<-read.csv("MOCE26b_ACi.csv")
M26c<-read.csv("MOCE26c_ACi.csv")
M31a<-read.csv("MOCE31a_ACi.csv")
M31b<-read.csv("MOCE31b_ACi.csv")
M31c<-read.csv("MOCE31c_ACi.csv")
A21a<-read.csv("ALRU21a_ACi.csv")
A21b<-read.csv("ALRU21b_ACi.csv")
A21c<-read.csv("ALRU21c_ACi.csv")
A26a<-read.csv("ALRU26a_ACi.csv")
A26b<-read.csv("ALRU26b_ACi.csv")
A26c<-read.csv("ALRU26c_ACi.csv")
A31a<-read.csv("ALRU31a_ACi.csv")
A31b<-read.csv("ALRU31b_ACi.csv")
A31c<-read.csv("ALRU31c_ACi.csv")
G21a<-read.csv("GLSE21a_ACi.csv")
G21b<-read.csv("GLSE21b_ACi.csv")
G21c<-read.csv("GLSE21c_ACi.csv")
G21d<-read.csv("GLSE21d_ACi.csv")
G26a<-read.csv("GLSE26a_ACi.csv")
G26b<-read.csv("GLSE26b_ACi.csv")
G26c<-read.csv("GLSE26c_ACi.csv")
G26d<-read.csv("GLSE26d_ACi.csv")
G31a<-read.csv("GLSE31a_ACi.csv")
G31b<-read.csv("GLSE31b_ACi.csv")
G31c<-read.csv("GLSE31c_ACi.csv")
R21a<-read.csv("ROPS21a_ACi.csv")
R21b<-read.csv("ROPS21b_ACi.csv")
R21c<-read.csv("ROPS21c_ACi.csv")
R21d<-read.csv("ROPS21d_ACi.csv")
R26a<-read.csv("ROPS26a_ACi.csv")
R26b<-read.csv("ROPS26b_ACi.csv")
R26c<-read.csv("ROPS26c_ACi.csv")
R31a<-read.csv("ROPS31a_ACi.csv")
R31b<-read.csv("ROPS31b_ACi.csv")
R31c<-read.csv("ROPS31c_ACi.csv")

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
p.M21.Rd<-rep(NA,length(M21a$Curve))
for (i in 1:length(M21a$Curve)){
  p.M21.Rd[i]<-norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,M21a$Curve[i])
}
M21a$Rd<-p.M21.Rd

#Subset data depending on if they were collected in an ascending or descending order
M21a.d<-M21a[M21a$Direction == "down",]
M21a.u<-M21a[M21a$Direction == "up",]

#
#Fit curves
#

#Descending
m.M21a.d<-fitacis(M21a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21a.d)
coef(m.M21a.d)

#Ascending
m.M21a.u<-fitacis(M21a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21a.u)
coef(m.M21a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M21a.d$`10`$df[,2]~m.M21a.d$`10`$df[,3])
summary(lm(m.M21a.d$`10`$df[,2]~m.M21a.d$`10`$df[,3]))
plot(m.M21a.d$`15`$df[,2]~m.M21a.d$`15`$df[,3])
summary(lm(m.M21a.d$`15`$df[,2]~m.M21a.d$`15`$df[,3]))
plot(m.M21a.d$`20`$df[,2]~m.M21a.d$`20`$df[,3])
summary(lm(m.M21a.d$`20`$df[,2]~m.M21a.d$`20`$df[,3]))

#Ascending
plot(m.M21a.u$`20`$df[,2]~m.M21a.u$`20`$df[,3])
summary(lm(m.M21a.u$`20`$df[,2]~m.M21a.u$`20`$df[,3]))
plot(m.M21a.u$`25`$df[,2]~m.M21a.u$`25`$df[,3])
summary(lm(m.M21a.u$`25`$df[,2]~m.M21a.u$`25`$df[,3]))
plot(m.M21a.u$`30`$df[,2]~m.M21a.u$`30`$df[,3])
summary(lm(m.M21a.u$`30`$df[,2]~m.M21a.u$`30`$df[,3]))
plot(m.M21a.u$`35`$df[,2]~m.M21a.u$`35`$df[,3])
summary(lm(m.M21a.u$`35`$df[,2]~m.M21a.u$`35`$df[,3]))
plot(m.M21a.u$`40`$df[,2]~m.M21a.u$`40`$df[,3])
summary(lm(m.M21a.u$`40`$df[,2]~m.M21a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.M21.Rd<-rep(NA,length(M21b$Curve))
for (i in 1:length(M21b$Curve)){
  p.M21.Rd[i]<-norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,M21b$Curve[i])
}
M21b$Rd<-p.M21.Rd

#Subset data depending on if they were collected in an ascending or descending order
M21b.d<-M21b[M21b$Direction == "down",]
M21b.u<-M21b[M21b$Direction == "up",]

#
#Fit curves
#

#Descending
m.M21b.d<-fitacis(M21b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21b.d)
coef(m.M21b.d)

#Ascending
m.M21b.u<-fitacis(M21b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21b.u)
coef(m.M21b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M21b.d$`10`$df[,2]~m.M21b.d$`10`$df[,3])
summary(lm(m.M21b.d$`10`$df[,2]~m.M21b.d$`10`$df[,3]))
plot(m.M21b.d$`15`$df[,2]~m.M21b.d$`15`$df[,3])
summary(lm(m.M21b.d$`15`$df[,2]~m.M21b.d$`15`$df[,3]))
plot(m.M21b.d$`20`$df[,2]~m.M21b.d$`20`$df[,3])
summary(lm(m.M21b.d$`20`$df[,2]~m.M21b.d$`20`$df[,3]))

#Ascending
plot(m.M21b.u$`20`$df[,2]~m.M21b.u$`20`$df[,3])
summary(lm(m.M21b.u$`20`$df[,2]~m.M21b.u$`20`$df[,3]))
plot(m.M21b.u$`25`$df[,2]~m.M21b.u$`25`$df[,3])
summary(lm(m.M21b.u$`25`$df[,2]~m.M21b.u$`25`$df[,3]))
plot(m.M21b.u$`30`$df[,2]~m.M21b.u$`30`$df[,3])
summary(lm(m.M21b.u$`30`$df[,2]~m.M21b.u$`30`$df[,3]))
plot(m.M21b.u$`35`$df[,2]~m.M21b.u$`35`$df[,3])
summary(lm(m.M21b.u$`35`$df[,2]~m.M21b.u$`35`$df[,3]))
plot(m.M21b.u$`40`$df[,2]~m.M21b.u$`40`$df[,3])
summary(lm(m.M21b.u$`40`$df[,2]~m.M21b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.M21.Rd<-rep(NA,length(M21c$Curve))
for (i in 1:length(M21c$Curve)){
  p.M21.Rd[i]<-norm.Topt.s.lin(0.753/Rd.25.M21,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],18.5,M21c$Curve[i])
}
M21c$Rd<-p.M21.Rd

#Subset data depending on if they were collected in an ascending or descending order
M21c.d<-M21c[M21c$Direction == "down",]
M21c.u<-M21c[M21c$Direction == "up",]

#
#Fit curves
#

#Descending
m.M21c.d<-fitacis(M21c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21c.d)
coef(m.M21c.d)

#Ascending
m.M21c.u<-fitacis(M21c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M21c.u)
coef(m.M21c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M21c.d$`10`$df[,2]~m.M21c.d$`10`$df[,3])
summary(lm(m.M21c.d$`10`$df[,2]~m.M21c.d$`10`$df[,3]))
plot(m.M21c.d$`15`$df[,2]~m.M21c.d$`15`$df[,3])
summary(lm(m.M21c.d$`15`$df[,2]~m.M21c.d$`15`$df[,3]))
plot(m.M21c.d$`20`$df[,2]~m.M21c.d$`20`$df[,3])
summary(lm(m.M21c.d$`20`$df[,2]~m.M21c.d$`20`$df[,3]))

#Ascending
plot(m.M21c.u$`20`$df[,2]~m.M21c.u$`20`$df[,3])
summary(lm(m.M21c.u$`20`$df[,2]~m.M21c.u$`20`$df[,3]))
plot(m.M21c.u$`25`$df[,2]~m.M21c.u$`25`$df[,3])
summary(lm(m.M21c.u$`25`$df[,2]~m.M21c.u$`25`$df[,3]))
plot(m.M21c.u$`30`$df[,2]~m.M21c.u$`30`$df[,3])
summary(lm(m.M21c.u$`30`$df[,2]~m.M21c.u$`30`$df[,3]))
plot(m.M21c.u$`35`$df[,2]~m.M21c.u$`35`$df[,3])
summary(lm(m.M21c.u$`35`$df[,2]~m.M21c.u$`35`$df[,3]))
plot(m.M21c.u$`40`$df[,2]~m.M21c.u$`40`$df[,3])
summary(lm(m.M21c.u$`40`$df[,2]~m.M21c.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.M26.Rd<-rep(NA,length(M26a$Curve))
for (i in 1:length(M26a$Curve)){
  p.M26.Rd[i]<-norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,M26a$Curve[i])
}
M26a$Rd<-p.M26.Rd

#Subset data depending on if they were collected in an ascending or descending order
M26a.d<-M26a[M26a$Direction == "down",]
M26a.u<-M26a[M26a$Direction == "up",]

#
#Fit curves
#

#Descending
m.M26a.d<-fitacis(M26a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26a.d)
coef(m.M26a.d)

#Ascending
m.M26a.u<-fitacis(M26a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26a.u)
coef(m.M26a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M26a.d$`10`$df[,2]~m.M26a.d$`10`$df[,3])
summary(lm(m.M26a.d$`10`$df[,2]~m.M26a.d$`10`$df[,3]))
plot(m.M26a.d$`15`$df[,2]~m.M26a.d$`15`$df[,3])
summary(lm(m.M26a.d$`15`$df[,2]~m.M26a.d$`15`$df[,3]))
plot(m.M26a.d$`20`$df[,2]~m.M26a.d$`20`$df[,3])
summary(lm(m.M26a.d$`20`$df[,2]~m.M26a.d$`20`$df[,3]))
plot(m.M26a.d$`25`$df[,2]~m.M26a.d$`25`$df[,3])
summary(lm(m.M26a.d$`25`$df[,2]~m.M26a.d$`25`$df[,3]))

#Ascending
plot(m.M26a.u$`25`$df[,2]~m.M26a.u$`25`$df[,3])
summary(lm(m.M26a.u$`25`$df[,2]~m.M26a.u$`25`$df[,3]))
plot(m.M26a.u$`30`$df[,2]~m.M26a.u$`30`$df[,3])
summary(lm(m.M26a.u$`30`$df[,2]~m.M26a.u$`30`$df[,3]))
plot(m.M26a.u$`35`$df[,2]~m.M26a.u$`35`$df[,3])
summary(lm(m.M26a.u$`35`$df[,2]~m.M26a.u$`35`$df[,3]))
plot(m.M26a.u$`40`$df[,2]~m.M26a.u$`40`$df[,3])
summary(lm(m.M26a.u$`40`$df[,2]~m.M26a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.M26.Rd<-rep(NA,length(M26b$Curve))
for (i in 1:length(M26b$Curve)){
  p.M26.Rd[i]<-norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,M26b$Curve[i])
}
M26b$Rd<-p.M26.Rd

#Subset data depending on if they were collected in an ascending or descending order
M26b.d<-M26b[M26b$Direction == "down",]
M26b.u<-M26b[M26b$Direction == "up",]

#
#Fit curves
#

#Descending
m.M26b.d<-fitacis(M26b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26b.d)
coef(m.M26b.d)

#Ascending
m.M26b.u<-fitacis(M26b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26b.u)
coef(m.M26b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M26b.d$`10`$df[,2]~m.M26b.d$`10`$df[,3])
summary(lm(m.M26b.d$`10`$df[,2]~m.M26b.d$`10`$df[,3]))
plot(m.M26b.d$`15`$df[,2]~m.M26b.d$`15`$df[,3])
summary(lm(m.M26b.d$`15`$df[,2]~m.M26b.d$`15`$df[,3]))
plot(m.M26b.d$`20`$df[,2]~m.M26b.d$`20`$df[,3])
summary(lm(m.M26b.d$`20`$df[,2]~m.M26b.d$`20`$df[,3]))
plot(m.M26b.d$`25`$df[,2]~m.M26b.d$`25`$df[,3])
summary(lm(m.M26b.d$`25`$df[,2]~m.M26b.d$`25`$df[,3]))

#Ascending
plot(m.M26b.u$`25`$df[,2]~m.M26b.u$`25`$df[,3])
summary(lm(m.M26b.u$`25`$df[,2]~m.M26b.u$`25`$df[,3]))
plot(m.M26b.u$`30`$df[,2]~m.M26b.u$`30`$df[,3])
summary(lm(m.M26b.u$`30`$df[,2]~m.M26b.u$`30`$df[,3]))
plot(m.M26b.u$`35`$df[,2]~m.M26b.u$`35`$df[,3])
summary(lm(m.M26b.u$`35`$df[,2]~m.M26b.u$`35`$df[,3]))
plot(m.M26b.u$`40`$df[,2]~m.M26b.u$`40`$df[,3])
summary(lm(m.M26b.u$`40`$df[,2]~m.M26b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.M26.Rd<-rep(NA,length(M26c$Curve))
for (i in 1:length(M26c$Curve)){
  p.M26.Rd[i]<-norm.Topt.s.lin(0.491/Rd.25.M26,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],23.5,M26c$Curve[i])
}
M26c$Rd<-p.M26.Rd

#Subset data depending on if they were collected in an ascending or descending order
M26c.d<-M26c[M26c$Direction == "down",]
M26c.u<-M26c[M26c$Direction == "up",]

#
#Fit curves
#

#Descending
m.M26c.d<-fitacis(M26c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26c.d)
coef(m.M26c.d)

#Ascending
m.M26c.u<-fitacis(M26c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M26c.u)
coef(m.M26c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M26c.d$`10`$df[,2]~m.M26c.d$`10`$df[,3])
summary(lm(m.M26c.d$`10`$df[,2]~m.M26c.d$`10`$df[,3]))
plot(m.M26c.d$`15`$df[,2]~m.M26c.d$`15`$df[,3])
summary(lm(m.M26c.d$`15`$df[,2]~m.M26c.d$`15`$df[,3]))
plot(m.M26c.d$`20`$df[,2]~m.M26c.d$`20`$df[,3])
summary(lm(m.M26c.d$`20`$df[,2]~m.M26c.d$`20`$df[,3]))
plot(m.M26c.d$`25`$df[,2]~m.M26c.d$`25`$df[,3])
summary(lm(m.M26c.d$`25`$df[,2]~m.M26c.d$`25`$df[,3]))

#Ascending
plot(m.M26c.u$`25`$df[,2]~m.M26c.u$`25`$df[,3])
summary(lm(m.M26c.u$`25`$df[,2]~m.M26c.u$`25`$df[,3]))
plot(m.M26c.u$`30`$df[,2]~m.M26c.u$`30`$df[,3])
summary(lm(m.M26c.u$`30`$df[,2]~m.M26c.u$`30`$df[,3]))
plot(m.M26c.u$`35`$df[,2]~m.M26c.u$`35`$df[,3])
summary(lm(m.M26c.u$`35`$df[,2]~m.M26c.u$`35`$df[,3]))
plot(m.M26c.u$`40`$df[,2]~m.M26c.u$`40`$df[,3])
summary(lm(m.M26c.u$`40`$df[,2]~m.M26c.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.M31.Rd<-rep(NA,length(M31a$Curve))
for (i in 1:length(M31a$Curve)){
  p.M31.Rd[i]<-norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,M31a$Curve[i])
}
M31a$Rd<-p.M31.Rd

#Subset data depending on if they were collected in an ascending or descending order
M31a.d<-M31a[M31a$Direction == "down",]
M31a.u<-M31a[M31a$Direction == "up",]

#
#Fit curves
#

#Descending
m.M31a.d<-fitacis(M31a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31a.d)
coef(m.M31a.d)

#Ascending
m.M31a.u<-fitacis(M31a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31a.u)
coef(m.M31a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M31a.d$`10`$df[,2]~m.M31a.d$`10`$df[,3])
summary(lm(m.M31a.d$`10`$df[,2]~m.M31a.d$`10`$df[,3]))
plot(m.M31a.d$`15`$df[,2]~m.M31a.d$`15`$df[,3])
summary(lm(m.M31a.d$`15`$df[,2]~m.M31a.d$`15`$df[,3]))
plot(m.M31a.d$`20`$df[,2]~m.M31a.d$`20`$df[,3])
summary(lm(m.M31a.d$`20`$df[,2]~m.M31a.d$`20`$df[,3]))
plot(m.M31a.d$`25`$df[,2]~m.M31a.d$`25`$df[,3])
summary(lm(m.M31a.d$`25`$df[,2]~m.M31a.d$`25`$df[,3]))
plot(m.M31a.d$`30`$df[,2]~m.M31a.d$`30`$df[,3])
summary(lm(m.M31a.d$`30`$df[,2]~m.M31a.d$`30`$df[,3]))

#Ascending
plot(m.M31a.u$`30`$df[,2]~m.M31a.u$`30`$df[,3])
summary(lm(m.M31a.u$`30`$df[,2]~m.M31a.u$`30`$df[,3]))
plot(m.M31a.u$`35`$df[,2]~m.M31a.u$`35`$df[,3])
summary(lm(m.M31a.u$`35`$df[,2]~m.M31a.u$`35`$df[,3]))
plot(m.M31a.u$`40`$df[,2]~m.M31a.u$`40`$df[,3])
summary(lm(m.M31a.u$`40`$df[,2]~m.M31a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.M31.Rd<-rep(NA,length(M31b$Curve))
for (i in 1:length(M31b$Curve)){
  p.M31.Rd[i]<-norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,M31b$Curve[i])
}
M31b$Rd<-p.M31.Rd

#Subset data depending on if they were collected in an ascending or descending order
M31b.d<-M31b[M31b$Direction == "down",]
M31b.u<-M31b[M31b$Direction == "up",]

#
#Fit curves
#

#Descending
m.M31b.d<-fitacis(M31b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31b.d)
coef(m.M31b.d)

#Ascending
m.M31b.u<-fitacis(M31b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31b.u)
coef(m.M31b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M31b.d$`15`$df[,2]~m.M31b.d$`15`$df[,3])
summary(lm(m.M31b.d$`15`$df[,2]~m.M31b.d$`15`$df[,3]))
plot(m.M31b.d$`20`$df[,2]~m.M31b.d$`20`$df[,3])
summary(lm(m.M31b.d$`20`$df[,2]~m.M31b.d$`20`$df[,3]))
plot(m.M31b.d$`25`$df[,2]~m.M31b.d$`25`$df[,3])
summary(lm(m.M31b.d$`25`$df[,2]~m.M31b.d$`25`$df[,3]))
plot(m.M31b.d$`30`$df[,2]~m.M31b.d$`30`$df[,3])
summary(lm(m.M31b.d$`30`$df[,2]~m.M31b.d$`30`$df[,3]))

#Ascending
plot(m.M31b.u$`30`$df[,2]~m.M31b.u$`30`$df[,3])
summary(lm(m.M31b.u$`30`$df[,2]~m.M31b.u$`30`$df[,3]))
plot(m.M31b.u$`35`$df[,2]~m.M31b.u$`35`$df[,3])
summary(lm(m.M31b.u$`35`$df[,2]~m.M31b.u$`35`$df[,3]))
plot(m.M31b.u$`40`$df[,2]~m.M31b.u$`40`$df[,3])
summary(lm(m.M31b.u$`40`$df[,2]~m.M31b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.M31.Rd<-rep(NA,length(M31c$Curve))
for (i in 1:length(M31c$Curve)){
  p.M31.Rd[i]<-norm.Topt.s.lin(0.229/Rd.25.M31,ft.Topt.s.lin.M[11],ft.Topt.s.lin.M[12],ft.Topt.s.lin.M[13],ft.Topt.s.lin.M[14],28.5,M31c$Curve[i])
}
M31c$Rd<-p.M31.Rd

#Subset data depending on if they were collected in an ascending or descending order
M31c.d<-M31c[M31c$Direction == "down",]
M31c.u<-M31c[M31c$Direction == "up",]

#
#Fit curves
#

#Descending
m.M31c.d<-fitacis(M31c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31c.d)
coef(m.M31c.d)

#Ascending
m.M31c.u<-fitacis(M31c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.M31c.u)
coef(m.M31c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.M31c.d$`10`$df[,2]~m.M31c.d$`10`$df[,3])
summary(lm(m.M31c.d$`10`$df[,2]~m.M31c.d$`10`$df[,3]))
plot(m.M31c.d$`15`$df[,2]~m.M31c.d$`15`$df[,3])
summary(lm(m.M31c.d$`15`$df[,2]~m.M31c.d$`15`$df[,3]))
plot(m.M31c.d$`20`$df[,2]~m.M31c.d$`20`$df[,3])
summary(lm(m.M31c.d$`20`$df[,2]~m.M31c.d$`20`$df[,3]))
plot(m.M31c.d$`25`$df[,2]~m.M31c.d$`25`$df[,3])
summary(lm(m.M31c.d$`25`$df[,2]~m.M31c.d$`25`$df[,3]))
plot(m.M31c.d$`30`$df[,2]~m.M31c.d$`30`$df[,3])
summary(lm(m.M31c.d$`30`$df[,2]~m.M31c.d$`30`$df[,3]))

#Ascending
plot(m.M31c.u$`30`$df[,2]~m.M31c.u$`30`$df[,3])
summary(lm(m.M31c.u$`30`$df[,2]~m.M31c.u$`30`$df[,3]))
plot(m.M31c.u$`35`$df[,2]~m.M31c.u$`35`$df[,3])
summary(lm(m.M31c.u$`35`$df[,2]~m.M31c.u$`35`$df[,3]))
plot(m.M31c.u$`40`$df[,2]~m.M31c.u$`40`$df[,3])
summary(lm(m.M31c.u$`40`$df[,2]~m.M31c.u$`40`$df[,3]))

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
p.A21.Rd<-rep(NA,length(A21a$Curve))
for (i in 1:length(A21a$Curve)){
  p.A21.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,A21a$Curve[i])
}
A21a$Rd<-p.A21.Rd

#Subset data depending on if they were collected in an ascending or descending order
A21a.d<-A21a[A21a$Direction == "down",]
A21a.u<-A21a[A21a$Direction == "up",]

#
#Fit curves
#

#Descending
m.A21a.d<-fitacis(A21a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21a.d)
coef(m.A21a.d)

#Ascending
m.A21a.u<-fitacis(A21a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21a.u)
coef(m.A21a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A21a.d$`10`$df[,2]~m.A21a.d$`10`$df[,3])
summary(lm(m.A21a.d$`10`$df[,2]~m.A21a.d$`10`$df[,3]))
plot(m.A21a.d$`15`$df[,2]~m.A21a.d$`15`$df[,3])
summary(lm(m.A21a.d$`15`$df[,2]~m.A21a.d$`15`$df[,3]))
plot(m.A21a.d$`20`$df[,2]~m.A21a.d$`20`$df[,3])
summary(lm(m.A21a.d$`20`$df[,2]~m.A21a.d$`20`$df[,3]))

#Ascending
plot(m.A21a.u$`20`$df[,2]~m.A21a.u$`20`$df[,3])
summary(lm(m.A21a.u$`20`$df[,2]~m.A21a.u$`20`$df[,3]))
plot(m.A21a.u$`25`$df[,2]~m.A21a.u$`25`$df[,3])
summary(lm(m.A21a.u$`25`$df[,2]~m.A21a.u$`25`$df[,3]))
plot(m.A21a.u$`30`$df[,2]~m.A21a.u$`30`$df[,3])
summary(lm(m.A21a.u$`30`$df[,2]~m.A21a.u$`30`$df[,3]))
plot(m.A21a.u$`35`$df[,2]~m.A21a.u$`35`$df[,3])
summary(lm(m.A21a.u$`35`$df[,2]~m.A21a.u$`35`$df[,3]))
plot(m.A21a.u$`40`$df[,2]~m.A21a.u$`40`$df[,3])
summary(lm(m.A21a.u$`40`$df[,2]~m.A21a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.A21.Rd<-rep(NA,length(A21b$Curve))
for (i in 1:length(A21b$Curve)){
  p.A21.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,A21b$Curve[i])
}
A21b$Rd<-p.A21.Rd

#Subset data depending on if they were collected in an ascending or descending order
A21b.d<-A21b[A21b$Direction == "down",]
A21b.u<-A21b[A21b$Direction == "up",]

#
#Fit curves
#

#Descending
m.A21b.d<-fitacis(A21b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21b.d)
coef(m.A21b.d)

#Ascending
m.A21b.u<-fitacis(A21b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21b.u)
coef(m.A21b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A21b.d$`10`$df[,2]~m.A21b.d$`10`$df[,3])
summary(lm(m.A21b.d$`10`$df[,2]~m.A21b.d$`10`$df[,3]))
plot(m.A21b.d$`15`$df[,2]~m.A21b.d$`15`$df[,3])
summary(lm(m.A21b.d$`15`$df[,2]~m.A21b.d$`15`$df[,3]))
plot(m.A21b.d$`20`$df[,2]~m.A21b.d$`20`$df[,3])
summary(lm(m.A21b.d$`20`$df[,2]~m.A21b.d$`20`$df[,3]))

#Ascending
plot(m.A21b.u$`20`$df[,2]~m.A21b.u$`20`$df[,3])
summary(lm(m.A21b.u$`20`$df[,2]~m.A21b.u$`20`$df[,3]))
plot(m.A21b.u$`25`$df[,2]~m.A21b.u$`25`$df[,3])
summary(lm(m.A21b.u$`25`$df[,2]~m.A21b.u$`25`$df[,3]))
plot(m.A21b.u$`30`$df[,2]~m.A21b.u$`30`$df[,3])
summary(lm(m.A21b.u$`30`$df[,2]~m.A21b.u$`30`$df[,3]))
plot(m.A21b.u$`35`$df[,2]~m.A21b.u$`35`$df[,3])
summary(lm(m.A21b.u$`35`$df[,2]~m.A21b.u$`35`$df[,3]))
plot(m.A21b.u$`40`$df[,2]~m.A21b.u$`40`$df[,3])
summary(lm(m.A21b.u$`40`$df[,2]~m.A21b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.A21.Rd<-rep(NA,length(A21c$Curve))
for (i in 1:length(A21c$Curve)){
  p.A21.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A21,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],18.5,A21c$Curve[i])
}
A21c$Rd<-p.A21.Rd

#Subset data depending on if they were collected in an ascending or descending order
A21c.d<-A21c[A21c$Direction == "down",]
A21c.u<-A21c[A21c$Direction == "up",]

#
#Fit curves
#

#Descending
m.A21c.d<-fitacis(A21c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21c.d)
coef(m.A21c.d)

#Ascending
m.A21c.u<-fitacis(A21c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A21c.u)
coef(m.A21c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A21c.d$`10`$df[,2]~m.A21c.d$`10`$df[,3])
summary(lm(m.A21c.d$`10`$df[,2]~m.A21c.d$`10`$df[,3]))
plot(m.A21c.d$`15`$df[,2]~m.A21c.d$`15`$df[,3])
summary(lm(m.A21c.d$`15`$df[,2]~m.A21c.d$`15`$df[,3]))
plot(m.A21c.d$`20`$df[,2]~m.A21c.d$`20`$df[,3])
summary(lm(m.A21c.d$`20`$df[,2]~m.A21c.d$`20`$df[,3]))

#Ascending
plot(m.A21c.u$`20`$df[,2]~m.A21c.u$`20`$df[,3])
summary(lm(m.A21c.u$`20`$df[,2]~m.A21c.u$`20`$df[,3]))
plot(m.A21c.u$`25`$df[,2]~m.A21c.u$`25`$df[,3])
summary(lm(m.A21c.u$`25`$df[,2]~m.A21c.u$`25`$df[,3]))
plot(m.A21c.u$`30`$df[,2]~m.A21c.u$`30`$df[,3])
summary(lm(m.A21c.u$`30`$df[,2]~m.A21c.u$`30`$df[,3]))
plot(m.A21c.u$`35`$df[,2]~m.A21c.u$`35`$df[,3])
summary(lm(m.A21c.u$`35`$df[,2]~m.A21c.u$`35`$df[,3]))
plot(m.A21c.u$`40`$df[,2]~m.A21c.u$`40`$df[,3])
summary(lm(m.A21c.u$`40`$df[,2]~m.A21c.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.A26.Rd<-rep(NA,length(A26a$Curve))
for (i in 1:length(A26a$Curve)){
  p.A26.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,A26a$Curve[i])
}
A26a$Rd<-p.A26.Rd

#Subset data depending on if they were collected in an ascending or descending order
A26a.d<-A26a[A26a$Direction == "down",]
A26a.u<-A26a[A26a$Direction == "up",]

#
#Fit curves
#

#Descending
m.A26a.d<-fitacis(A26a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26a.d)
coef(m.A26a.d)

#Ascending
m.A26a.u<-fitacis(A26a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26a.u)
coef(m.A26a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A26a.d$`10`$df[,2]~m.A26a.d$`10`$df[,3])
summary(lm(m.A26a.d$`10`$df[,2]~m.A26a.d$`10`$df[,3]))
plot(m.A26a.d$`15`$df[,2]~m.A26a.d$`15`$df[,3])
summary(lm(m.A26a.d$`15`$df[,2]~m.A26a.d$`15`$df[,3]))
plot(m.A26a.d$`20`$df[,2]~m.A26a.d$`20`$df[,3])
summary(lm(m.A26a.d$`20`$df[,2]~m.A26a.d$`20`$df[,3]))
plot(m.A26a.d$`25`$df[,2]~m.A26a.d$`25`$df[,3])
summary(lm(m.A26a.d$`25`$df[,2]~m.A26a.d$`25`$df[,3]))

#Ascending
plot(m.A26a.u$`25`$df[,2]~m.A26a.u$`25`$df[,3])
summary(lm(m.A26a.u$`25`$df[,2]~m.A26a.u$`25`$df[,3]))
plot(m.A26a.u$`30`$df[,2]~m.A26a.u$`30`$df[,3])
summary(lm(m.A26a.u$`30`$df[,2]~m.A26a.u$`30`$df[,3]))
plot(m.A26a.u$`35`$df[,2]~m.A26a.u$`35`$df[,3])
summary(lm(m.A26a.u$`35`$df[,2]~m.A26a.u$`35`$df[,3]))
plot(m.A26a.u$`40`$df[,2]~m.A26a.u$`40`$df[,3])
summary(lm(m.A26a.u$`40`$df[,2]~m.A26a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.A26.Rd<-rep(NA,length(A26b$Curve))
for (i in 1:length(A26b$Curve)){
  p.A26.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,A26b$Curve[i])
}
A26b$Rd<-p.A26.Rd

#Subset data depending on if they were collected in an ascending or descending order
A26b.d<-A26b[A26b$Direction == "down",]
A26b.u<-A26b[A26b$Direction == "up",]

#
#Fit curves
#

#Descending
m.A26b.d<-fitacis(A26b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26b.d)
coef(m.A26b.d) 

#Ascending
m.A26b.u<-fitacis(A26b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26b.u)
coef(m.A26b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A26b.d$`10`$df[,2]~m.A26b.d$`10`$df[,3])
summary(lm(m.A26b.d$`10`$df[,2]~m.A26b.d$`10`$df[,3]))
plot(m.A26b.d$`15`$df[,2]~m.A26b.d$`15`$df[,3])
summary(lm(m.A26b.d$`15`$df[,2]~m.A26b.d$`15`$df[,3]))
plot(m.A26b.d$`20`$df[,2]~m.A26b.d$`20`$df[,3])
summary(lm(m.A26b.d$`20`$df[,2]~m.A26b.d$`20`$df[,3]))
plot(m.A26b.d$`25`$df[,2]~m.A26b.d$`25`$df[,3])
summary(lm(m.A26b.d$`25`$df[,2]~m.A26b.d$`25`$df[,3]))

#Ascending
plot(m.A26b.u$`25`$df[,2]~m.A26b.u$`25`$df[,3])
summary(lm(m.A26b.u$`25`$df[,2]~m.A26b.u$`25`$df[,3]))
plot(m.A26b.u$`30`$df[,2]~m.A26b.u$`30`$df[,3])
summary(lm(m.A26b.u$`30`$df[,2]~m.A26b.u$`30`$df[,3]))
plot(m.A26b.u$`35`$df[,2]~m.A26b.u$`35`$df[,3])
summary(lm(m.A26b.u$`35`$df[,2]~m.A26b.u$`35`$df[,3]))
plot(m.A26b.u$`40`$df[,2]~m.A26b.u$`40`$df[,3])
summary(lm(m.A26b.u$`40`$df[,2]~m.A26b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.A26.Rd<-rep(NA,length(A26c$Curve))
for (i in 1:length(A26c$Curve)){
  p.A26.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A26,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],23.5,A26c$Curve[i])
}
A26c$Rd<-p.A26.Rd

#Subset data depending on if they were collected in an ascending or descending order
A26c.d<-A26c[A26c$Direction == "down",]
A26c.u<-A26c[A26c$Direction == "up",]

#
#Fit curves
#

#Descending
m.A26c.d<-fitacis(A26c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26c.d)
coef(m.A26c.d) 

#Ascending
m.A26c.u<-fitacis(A26c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A26c.u)
coef(m.A26c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A26c.d$`10`$df[,2]~m.A26c.d$`10`$df[,3])
summary(lm(m.A26c.d$`10`$df[,2]~m.A26c.d$`10`$df[,3]))
plot(m.A26c.d$`15`$df[,2]~m.A26c.d$`15`$df[,3])
summary(lm(m.A26c.d$`15`$df[,2]~m.A26c.d$`15`$df[,3]))
plot(m.A26c.d$`20`$df[,2]~m.A26c.d$`20`$df[,3])
summary(lm(m.A26c.d$`20`$df[,2]~m.A26c.d$`20`$df[,3]))
plot(m.A26c.d$`25`$df[,2]~m.A26c.d$`25`$df[,3])
summary(lm(m.A26c.d$`25`$df[,2]~m.A26c.d$`25`$df[,3]))

#Ascending
plot(m.A26c.u$`25`$df[,2]~m.A26c.u$`25`$df[,3])
summary(lm(m.A26c.u$`25`$df[,2]~m.A26c.u$`25`$df[,3]))
plot(m.A26c.u$`30`$df[,2]~m.A26c.u$`30`$df[,3])
summary(lm(m.A26c.u$`30`$df[,2]~m.A26c.u$`30`$df[,3]))
plot(m.A26c.u$`35`$df[,2]~m.A26c.u$`35`$df[,3])
summary(lm(m.A26c.u$`35`$df[,2]~m.A26c.u$`35`$df[,3]))
plot(m.A26c.u$`40`$df[,2]~m.A26c.u$`40`$df[,3])
summary(lm(m.A26c.u$`40`$df[,2]~m.A26c.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.A31.Rd<-rep(NA,length(A31a$Curve))
for (i in 1:length(A31a$Curve)){
  p.A31.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,A31a$Curve[i])
}
A31a$Rd<-p.A31.Rd

#Subset data depending on if they were collected in an ascending or descending order
A31a.d<-A31a[A31a$Direction == "down",]
A31a.u<-A31a[A31a$Direction == "up",]

#
#Fit curves
#

#Descending
m.A31a.d<-fitacis(A31a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31a.d)
coef(m.A31a.d)

#Ascending
m.A31a.u<-fitacis(A31a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31a.u)
coef(m.A31a.u) 

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A31a.d$`10`$df[,2]~m.A31a.d$`10`$df[,3])
summary(lm(m.A31a.d$`10`$df[,2]~m.A31a.d$`10`$df[,3]))
plot(m.A31a.d$`15`$df[,2]~m.A31a.d$`15`$df[,3])
summary(lm(m.A31a.d$`15`$df[,2]~m.A31a.d$`15`$df[,3]))
plot(m.A31a.d$`20`$df[,2]~m.A31a.d$`20`$df[,3])
summary(lm(m.A31a.d$`20`$df[,2]~m.A31a.d$`20`$df[,3]))
plot(m.A31a.d$`25`$df[,2]~m.A31a.d$`25`$df[,3])
summary(lm(m.A31a.d$`25`$df[,2]~m.A31a.d$`25`$df[,3]))
plot(m.A31a.d$`30`$df[,2]~m.A31a.d$`30`$df[,3])
summary(lm(m.A31a.d$`30`$df[,2]~m.A31a.d$`30`$df[,3]))

#Ascending
plot(m.A31a.u$`30`$df[,2]~m.A31a.u$`30`$df[,3])
summary(lm(m.A31a.u$`30`$df[,2]~m.A31a.u$`30`$df[,3]))
plot(m.A31a.u$`35`$df[,2]~m.A31a.u$`35`$df[,3])
summary(lm(m.A31a.u$`35`$df[,2]~m.A31a.u$`35`$df[,3]))
plot(m.A31a.u$`40`$df[,2]~m.A31a.u$`40`$df[,3])
summary(lm(m.A31a.u$`40`$df[,2]~m.A31a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.A31.Rd<-rep(NA,length(A31b$Curve))
for (i in 1:length(A31b$Curve)){
  p.A31.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,A31b$Curve[i])
}
A31b$Rd<-p.A31.Rd

#Subset data depending on if they were collected in an ascending or descending order
A31b.d<-A31b[A31b$Direction == "down",]
A31b.u<-A31b[A31b$Direction == "up",]

#
#Fit curves
#

#Descending
m.A31b.d<-fitacis(A31b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31b.d)
coef(m.A31b.d)

#Ascending
m.A31b.u<-fitacis(A31b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31b.u)
coef(m.A31b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A31b.d$`10`$df[,2]~m.A31b.d$`10`$df[,3])
summary(lm(m.A31b.d$`10`$df[,2]~m.A31b.d$`10`$df[,3]))
plot(m.A31b.d$`15`$df[,2]~m.A31b.d$`15`$df[,3])
summary(lm(m.A31b.d$`15`$df[,2]~m.A31b.d$`15`$df[,3]))
plot(m.A31b.d$`20`$df[,2]~m.A31b.d$`20`$df[,3])
summary(lm(m.A31b.d$`20`$df[,2]~m.A31b.d$`20`$df[,3]))
plot(m.A31b.d$`25`$df[,2]~m.A31b.d$`25`$df[,3])
summary(lm(m.A31b.d$`25`$df[,2]~m.A31b.d$`25`$df[,3]))
plot(m.A31b.d$`30`$df[,2]~m.A31b.d$`30`$df[,3])
summary(lm(m.A31b.d$`30`$df[,2]~m.A31b.d$`30`$df[,3]))

#Ascending
plot(m.A31b.u$`30`$df[,2]~m.A31b.u$`30`$df[,3])
summary(lm(m.A31b.u$`30`$df[,2]~m.A31b.u$`30`$df[,3]))
plot(m.A31b.u$`35`$df[,2]~m.A31b.u$`35`$df[,3])
summary(lm(m.A31b.u$`35`$df[,2]~m.A31b.u$`35`$df[,3]))
plot(m.A31b.u$`40`$df[,2]~m.A31b.u$`40`$df[,3])
summary(lm(m.A31b.u$`40`$df[,2]~m.A31b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.A31.Rd<-rep(NA,length(A31c$Curve))
for (i in 1:length(A31c$Curve)){
  p.A31.Rd[i]<-norm.Topt.s.lin(0.938/Rd.25.A31,ft.Topt.s.lin.A[11],ft.Topt.s.lin.A[12],ft.Topt.s.lin.A[13],ft.Topt.s.lin.A[14],28.5,A31c$Curve[i])
}
A31c$Rd<-p.A31.Rd

#Subset data depending on if they were collected in an ascending or descending order
A31c.d<-A31c[A31c$Direction == "down",]
A31c.u<-A31c[A31c$Direction == "up",]

#
#Fit curves
#

#Descending
m.A31c.d<-fitacis(A31c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31c.d)
coef(m.A31c.d)

#Ascending
m.A31c.u<-fitacis(A31c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.A31c.u)
coef(m.A31c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.A31c.d$`10`$df[,2]~m.A31c.d$`10`$df[,3])
summary(lm(m.A31c.d$`10`$df[,2]~m.A31c.d$`10`$df[,3]))
plot(m.A31c.d$`15`$df[,2]~m.A31c.d$`15`$df[,3])
summary(lm(m.A31c.d$`15`$df[,2]~m.A31c.d$`15`$df[,3]))
plot(m.A31c.d$`20`$df[,2]~m.A31c.d$`20`$df[,3])
summary(lm(m.A31c.d$`20`$df[,2]~m.A31c.d$`20`$df[,3]))
plot(m.A31c.d$`25`$df[,2]~m.A31c.d$`25`$df[,3])
summary(lm(m.A31c.d$`25`$df[,2]~m.A31c.d$`25`$df[,3]))
plot(m.A31c.d$`30`$df[,2]~m.A31c.d$`30`$df[,3])
summary(lm(m.A31c.d$`30`$df[,2]~m.A31c.d$`30`$df[,3]))

#Ascending
plot(m.A31c.u$`30`$df[,2]~m.A31c.u$`30`$df[,3])
summary(lm(m.A31c.u$`30`$df[,2]~m.A31c.u$`30`$df[,3]))
plot(m.A31c.u$`35`$df[,2]~m.A31c.u$`35`$df[,3])
summary(lm(m.A31c.u$`35`$df[,2]~m.A31c.u$`35`$df[,3]))
plot(m.A31c.u$`40`$df[,2]~m.A31c.u$`40`$df[,3])
summary(lm(m.A31c.u$`40`$df[,2]~m.A31c.u$`40`$df[,3]))

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
p.G21.Rd<-rep(NA,length(G21a$Curve))
for (i in 1:length(G21a$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21a$Curve[i])
}
G21a$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21a.d<-G21a[G21a$Direction == "down",]
G21a.u<-G21a[G21a$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21a.d<-fitacis(G21a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21a.d)
coef(m.G21a.d)

#Ascending
m.G21a.u<-fitacis(G21a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21a.u)
coef(m.G21a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21a.d$`10`$df[,2]~m.G21a.d$`10`$df[,3])
summary(lm(m.G21a.d$`10`$df[,2]~m.G21a.d$`10`$df[,3]))
plot(m.G21a.d$`15`$df[,2]~m.G21a.d$`15`$df[,3])
summary(lm(m.G21a.d$`15`$df[,2]~m.G21a.d$`15`$df[,3]))
plot(m.G21a.d$`20`$df[,2]~m.G21a.d$`20`$df[,3])
summary(lm(m.G21a.d$`20`$df[,2]~m.G21a.d$`20`$df[,3]))

#Ascending
plot(m.G21a.u$`20`$df[,2]~m.G21a.u$`20`$df[,3])
summary(lm(m.G21a.u$`20`$df[,2]~m.G21a.u$`20`$df[,3]))
plot(m.G21a.u$`25`$df[,2]~m.G21a.u$`25`$df[,3])
summary(lm(m.G21a.u$`25`$df[,2]~m.G21a.u$`25`$df[,3]))
plot(m.G21a.u$`30`$df[,2]~m.G21a.u$`30`$df[,3])
summary(lm(m.G21a.u$`30`$df[,2]~m.G21a.u$`30`$df[,3]))
plot(m.G21a.u$`35`$df[,2]~m.G21a.u$`35`$df[,3])
summary(lm(m.G21a.u$`35`$df[,2]~m.G21a.u$`35`$df[,3]))
plot(m.G21a.u$`40`$df[,2]~m.G21a.u$`40`$df[,3])
summary(lm(m.G21a.u$`40`$df[,2]~m.G21a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21b$Curve))
for (i in 1:length(G21b$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21b$Curve[i])
}
G21b$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21b.d<-G21b[G21b$Direction == "down",]
G21b.u<-G21b[G21b$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21b.d<-fitacis(G21b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21b.d)
coef(m.G21b.d)

#Ascending
m.G21b.u<-fitacis(G21b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21b.u)
coef(m.G21b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21b.d$`10`$df[,2]~m.G21b.d$`10`$df[,3])
summary(lm(m.G21b.d$`10`$df[,2]~m.G21b.d$`10`$df[,3]))
plot(m.G21b.d$`15`$df[,2]~m.G21b.d$`15`$df[,3])
summary(lm(m.G21b.d$`15`$df[,2]~m.G21b.d$`15`$df[,3]))
plot(m.G21b.d$`20`$df[,2]~m.G21b.d$`20`$df[,3])
summary(lm(m.G21b.d$`20`$df[,2]~m.G21b.d$`20`$df[,3]))

#Ascending
plot(m.G21b.u$`20`$df[,2]~m.G21b.u$`20`$df[,3])
summary(lm(m.G21b.u$`20`$df[,2]~m.G21b.u$`20`$df[,3]))
plot(m.G21b.u$`25`$df[,2]~m.G21b.u$`25`$df[,3])
summary(lm(m.G21b.u$`25`$df[,2]~m.G21b.u$`25`$df[,3]))
plot(m.G21b.u$`30`$df[,2]~m.G21b.u$`30`$df[,3])
summary(lm(m.G21b.u$`30`$df[,2]~m.G21b.u$`30`$df[,3]))
plot(m.G21b.u$`35`$df[,2]~m.G21b.u$`35`$df[,3])
summary(lm(m.G21b.u$`35`$df[,2]~m.G21b.u$`35`$df[,3]))
plot(m.G21b.u$`40`$df[,2]~m.G21b.u$`40`$df[,3])
summary(lm(m.G21b.u$`40`$df[,2]~m.G21b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21c$Curve))
for (i in 1:length(G21c$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21c$Curve[i])
}
G21c$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21c.d<-G21c[G21c$Direction == "down",]
G21c.u<-G21c[G21c$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21c.d<-fitacis(G21c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21c.d)
coef(m.G21c.d)

#Ascending
m.G21c.u<-fitacis(G21c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21c.u)
coef(m.G21c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21c.d$`10`$df[,2]~m.G21c.d$`10`$df[,3])
summary(lm(m.G21c.d$`10`$df[,2]~m.G21c.d$`10`$df[,3]))
plot(m.G21c.d$`15`$df[,2]~m.G21c.d$`15`$df[,3])
summary(lm(m.G21c.d$`15`$df[,2]~m.G21c.d$`15`$df[,3]))
plot(m.G21c.d$`20`$df[,2]~m.G21c.d$`20`$df[,3])
summary(lm(m.G21c.d$`20`$df[,2]~m.G21c.d$`20`$df[,3]))

#Ascending
plot(m.G21c.u$`20`$df[,2]~m.G21c.u$`20`$df[,3])
summary(lm(m.G21c.u$`20`$df[,2]~m.G21c.u$`20`$df[,3]))
plot(m.G21c.u$`25`$df[,2]~m.G21c.u$`25`$df[,3])
summary(lm(m.G21c.u$`25`$df[,2]~m.G21c.u$`25`$df[,3]))
plot(m.G21c.u$`30`$df[,2]~m.G21c.u$`30`$df[,3])
summary(lm(m.G21c.u$`30`$df[,2]~m.G21c.u$`30`$df[,3]))
plot(m.G21c.u$`35`$df[,2]~m.G21c.u$`35`$df[,3])
summary(lm(m.G21c.u$`35`$df[,2]~m.G21c.u$`35`$df[,3]))
plot(m.G21c.u$`40`$df[,2]~m.G21c.u$`40`$df[,3])
summary(lm(m.G21c.u$`40`$df[,2]~m.G21c.u$`40`$df[,3]))

##
#Replicate 4
##

#Calculate RL for each measurement temperature
p.G21.Rd<-rep(NA,length(G21d$Curve))
for (i in 1:length(G21d$Curve)){
  p.G21.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G21,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],18.5,G21d$Curve[i])
}
G21d$Rd<-p.G21.Rd

#Subset data depending on if they were collected in an ascending or descending order
G21d.d<-G21d[G21d$Direction == "down",]
G21d.u<-G21d[G21d$Direction == "up",]

#
#Fit curves
#

#Descending
m.G21d.d<-fitacis(G21d.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21d.d)
coef(m.G21d.d)

#Ascending
m.G21d.u<-fitacis(G21d.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G21d.u)
coef(m.G21d.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G21d.d$`10`$df[,2]~m.G21d.d$`10`$df[,3])
summary(lm(m.G21d.d$`10`$df[,2]~m.G21d.d$`10`$df[,3]))
plot(m.G21d.d$`15`$df[,2]~m.G21d.d$`15`$df[,3])
summary(lm(m.G21d.d$`15`$df[,2]~m.G21d.d$`15`$df[,3]))
plot(m.G21d.d$`20`$df[,2]~m.G21d.d$`20`$df[,3])
summary(lm(m.G21d.d$`20`$df[,2]~m.G21d.d$`20`$df[,3]))

#Ascending
plot(m.G21d.u$`20`$df[,2]~m.G21d.u$`20`$df[,3])
summary(lm(m.G21d.u$`20`$df[,2]~m.G21d.u$`20`$df[,3]))
plot(m.G21d.u$`25`$df[,2]~m.G21d.u$`25`$df[,3])
summary(lm(m.G21d.u$`25`$df[,2]~m.G21d.u$`25`$df[,3]))
plot(m.G21d.u$`30`$df[,2]~m.G21d.u$`30`$df[,3])
summary(lm(m.G21d.u$`30`$df[,2]~m.G21d.u$`30`$df[,3]))
plot(m.G21d.u$`35`$df[,2]~m.G21d.u$`35`$df[,3])
summary(lm(m.G21d.u$`35`$df[,2]~m.G21d.u$`35`$df[,3]))
plot(m.G21d.u$`40`$df[,2]~m.G21d.u$`40`$df[,3])
summary(lm(m.G21d.u$`40`$df[,2]~m.G21d.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26a$Curve))
for (i in 1:length(G26a$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26a$Curve[i])
}
G26a$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26a.d<-G26a[G26a$Direction == "down",]
G26a.u<-G26a[G26a$Direction == "up",]

#
#Fit curves
#

#Descending
m.G26a.d<-fitacis(G26a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26a.d)
coef(m.G26a.d)

#Ascending
m.G26a.u<-fitacis(G26a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26a.u)
coef(m.G26a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G26a.d$`10`$df[,2]~m.G26a.d$`10`$df[,3])
summary(lm(m.G26a.d$`10`$df[,2]~m.G26a.d$`10`$df[,3]))
plot(m.G26a.d$`15`$df[,2]~m.G26a.d$`15`$df[,3])
summary(lm(m.G26a.d$`15`$df[,2]~m.G26a.d$`15`$df[,3]))
plot(m.G26a.d$`20`$df[,2]~m.G26a.d$`20`$df[,3])
summary(lm(m.G26a.d$`20`$df[,2]~m.G26a.d$`20`$df[,3]))
plot(m.G26a.d$`25`$df[,2]~m.G26a.d$`25`$df[,3])
summary(lm(m.G26a.d$`25`$df[,2]~m.G26a.d$`25`$df[,3]))

#Ascending
plot(m.G26a.u$`25`$df[,2]~m.G26a.u$`25`$df[,3])
summary(lm(m.G26a.u$`25`$df[,2]~m.G26a.u$`25`$df[,3]))
plot(m.G26a.u$`30`$df[,2]~m.G26a.u$`30`$df[,3])
summary(lm(m.G26a.u$`30`$df[,2]~m.G26a.u$`30`$df[,3]))
plot(m.G26a.u$`35`$df[,2]~m.G26a.u$`35`$df[,3])
summary(lm(m.G26a.u$`35`$df[,2]~m.G26a.u$`35`$df[,3]))
plot(m.G26a.u$`40`$df[,2]~m.G26a.u$`40`$df[,3])
summary(lm(m.G26a.u$`40`$df[,2]~m.G26a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26b$Curve))
for (i in 1:length(G26b$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26b$Curve[i])
}
G26b$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26b.d<-G26b[G26b$Direction == "down",]
G26b.u<-G26b[G26b$Direction == "up",]

#
#Fit curves
#

#Descending
m.G26b.d<-fitacis(G26b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26b.d)
coef(m.G26b.d)

#Ascending
m.G26b.u<-fitacis(G26b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26b.u)
coef(m.G26b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G26b.d$`10`$df[,2]~m.G26b.d$`10`$df[,3])
summary(lm(m.G26b.d$`10`$df[,2]~m.G26b.d$`10`$df[,3]))
plot(m.G26b.d$`15`$df[,2]~m.G26b.d$`15`$df[,3])
summary(lm(m.G26b.d$`15`$df[,2]~m.G26b.d$`15`$df[,3]))
plot(m.G26b.d$`20`$df[,2]~m.G26b.d$`20`$df[,3])
summary(lm(m.G26b.d$`20`$df[,2]~m.G26b.d$`20`$df[,3]))
plot(m.G26b.d$`25`$df[,2]~m.G26b.d$`25`$df[,3])
summary(lm(m.G26b.d$`25`$df[,2]~m.G26b.d$`25`$df[,3]))

#Ascending
plot(m.G26b.u$`25`$df[,2]~m.G26b.u$`25`$df[,3])
summary(lm(m.G26b.u$`25`$df[,2]~m.G26b.u$`25`$df[,3]))
plot(m.G26b.u$`30`$df[,2]~m.G26b.u$`30`$df[,3])
summary(lm(m.G26b.u$`30`$df[,2]~m.G26b.u$`30`$df[,3]))
plot(m.G26b.u$`35`$df[,2]~m.G26b.u$`35`$df[,3])
summary(lm(m.G26b.u$`35`$df[,2]~m.G26b.u$`35`$df[,3]))
plot(m.G26b.u$`40`$df[,2]~m.G26b.u$`40`$df[,3])
summary(lm(m.G26b.u$`40`$df[,2]~m.G26b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26c$Curve))
for (i in 1:length(G26c$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26c$Curve[i])
}
G26c$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26c.d<-G26c[G26c$Direction == "down",]
G26c.u<-G26c[G26c$Direction == "up",]

#
#Fit curves
#

#Descending
m.G26c.d<-fitacis(G26c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26c.d)
coef(m.G26c.d)

#Ascending
m.G26c.u<-fitacis(G26c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26c.u)
coef(m.G26c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G26c.d$`10`$df[,2]~m.G26c.d$`10`$df[,3])
summary(lm(m.G26c.d$`10`$df[,2]~m.G26c.d$`10`$df[,3]))
plot(m.G26c.d$`15`$df[,2]~m.G26c.d$`15`$df[,3])
summary(lm(m.G26c.d$`15`$df[,2]~m.G26c.d$`15`$df[,3]))
plot(m.G26c.d$`20`$df[,2]~m.G26c.d$`20`$df[,3])
summary(lm(m.G26c.d$`20`$df[,2]~m.G26c.d$`20`$df[,3]))
plot(m.G26c.d$`25`$df[,2]~m.G26c.d$`25`$df[,3])
summary(lm(m.G26c.d$`25`$df[,2]~m.G26c.d$`25`$df[,3]))

#Ascending
plot(m.G26c.u$`25`$df[,2]~m.G26c.u$`25`$df[,3])
summary(lm(m.G26c.u$`25`$df[,2]~m.G26c.u$`25`$df[,3]))
plot(m.G26c.u$`30`$df[,2]~m.G26c.u$`30`$df[,3])
summary(lm(m.G26c.u$`30`$df[,2]~m.G26c.u$`30`$df[,3]))
plot(m.G26c.u$`35`$df[,2]~m.G26c.u$`35`$df[,3])
summary(lm(m.G26c.u$`35`$df[,2]~m.G26c.u$`35`$df[,3]))
plot(m.G26c.u$`40`$df[,2]~m.G26c.u$`40`$df[,3])
summary(lm(m.G26c.u$`40`$df[,2]~m.G26c.u$`40`$df[,3]))

##
#Replicate 4
##

#Calculate RL for each measurement temperature
p.G26.Rd<-rep(NA,length(G26d$Curve))
for (i in 1:length(G26d$Curve)){
  p.G26.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G26,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],23.5,G26d$Curve[i])
}
G26d$Rd<-p.G26.Rd

#Subset data depending on if they were collected in an ascending or descending order
G26d.u<-G26d[G26d$Direction == "up",]

#
#Fit curves
#

#Ascending
m.G26d.u<-fitacis(G26d.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G26d.u)
coef(m.G26d.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Ascending
plot(m.G26d.u$`25`$df[,2]~m.G26d.u$`25`$df[,3])
summary(lm(m.G26d.u$`25`$df[,2]~m.G26d.u$`25`$df[,3]))
plot(m.G26d.u$`30`$df[,2]~m.G26d.u$`30`$df[,3])
summary(lm(m.G26d.u$`30`$df[,2]~m.G26d.u$`30`$df[,3]))
plot(m.G26d.u$`35`$df[,2]~m.G26d.u$`35`$df[,3])
summary(lm(m.G26d.u$`35`$df[,2]~m.G26d.u$`35`$df[,3]))
plot(m.G26d.u$`40`$df[,2]~m.G26d.u$`40`$df[,3])
summary(lm(m.G26d.u$`40`$df[,2]~m.G26d.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.G31.Rd<-rep(NA,length(G31a$Curve))
for (i in 1:length(G31a$Curve)){
  p.G31.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,G31a$Curve[i])
}
G31a$Rd<-p.G31.Rd

#Subset data depending on if they were collected in an ascending or descending order
G31a.d<-G31a[G31a$Direction == "down",]
G31a.u<-G31a[G31a$Direction == "up",]

#
#Fit curves
#

#Descending
m.G31a.d<-fitacis(G31a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31a.d)
coef(m.G31a.d)

#Ascending
m.G31a.u<-fitacis(G31a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31a.u)
coef(m.G31a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G31a.d$`15`$df[,2]~m.G31a.d$`15`$df[,3])
summary(lm(m.G31a.d$`15`$df[,2]~m.G31a.d$`15`$df[,3]))
plot(m.G31a.d$`20`$df[,2]~m.G31a.d$`20`$df[,3])
summary(lm(m.G31a.d$`20`$df[,2]~m.G31a.d$`20`$df[,3]))
plot(m.G31a.d$`25`$df[,2]~m.G31a.d$`25`$df[,3])
summary(lm(m.G31a.d$`25`$df[,2]~m.G31a.d$`25`$df[,3]))
plot(m.G31a.d$`30`$df[,2]~m.G31a.d$`30`$df[,3])
summary(lm(m.G31a.d$`30`$df[,2]~m.G31a.d$`30`$df[,3]))

#Ascending
plot(m.G31a.u$`30`$df[,2]~m.G31a.u$`30`$df[,3])
summary(lm(m.G31a.u$`30`$df[,2]~m.G31a.u$`30`$df[,3]))
plot(m.G31a.u$`35`$df[,2]~m.G31a.u$`35`$df[,3])
summary(lm(m.G31a.u$`35`$df[,2]~m.G31a.u$`35`$df[,3]))
plot(m.G31a.u$`40`$df[,2]~m.G31a.u$`40`$df[,3])
summary(lm(m.G31a.u$`40`$df[,2]~m.G31a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.G31.Rd<-rep(NA,length(G31b$Curve))
for (i in 1:length(G31b$Curve)){
  p.G31.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,G31b$Curve[i])
}
G31b$Rd<-p.G31.Rd

#Subset data depending on if they were collected in an ascending or descending order
G31b.d<-G31b[G31b$Direction == "down",]
G31b.u<-G31b[G31b$Direction == "up",]

#
#Fit curves
#

#Descending
m.G31b.d<-fitacis(G31b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31b.d)
coef(m.G31b.d)

#Ascending
m.G31b.u<-fitacis(G31b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31b.u)
coef(m.G31b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G31b.d$`10`$df[,2]~m.G31b.d$`10`$df[,3])
summary(lm(m.G31b.d$`10`$df[,2]~m.G31b.d$`10`$df[,3]))
plot(m.G31b.d$`15`$df[,2]~m.G31b.d$`15`$df[,3])
summary(lm(m.G31b.d$`15`$df[,2]~m.G31b.d$`15`$df[,3]))
plot(m.G31b.d$`20`$df[,2]~m.G31b.d$`20`$df[,3])
summary(lm(m.G31b.d$`20`$df[,2]~m.G31b.d$`20`$df[,3]))
plot(m.G31b.d$`25`$df[,2]~m.G31b.d$`25`$df[,3])
summary(lm(m.G31b.d$`25`$df[,2]~m.G31b.d$`25`$df[,3]))
plot(m.G31b.d$`30`$df[,2]~m.G31b.d$`30`$df[,3])
summary(lm(m.G31b.d$`30`$df[,2]~m.G31b.d$`30`$df[,3]))

#Ascending
plot(m.G31b.u$`30`$df[,2]~m.G31b.u$`30`$df[,3])
summary(lm(m.G31b.u$`30`$df[,2]~m.G31b.u$`30`$df[,3]))
plot(m.G31b.u$`35`$df[,2]~m.G31b.u$`35`$df[,3])
summary(lm(m.G31b.u$`35`$df[,2]~m.G31b.u$`35`$df[,3]))
plot(m.G31b.u$`40`$df[,2]~m.G31b.u$`40`$df[,3])
summary(lm(m.G31b.u$`40`$df[,2]~m.G31b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.G31.Rd<-rep(NA,length(G31c$Curve))
for (i in 1:length(G31c$Curve)){
  p.G31.Rd[i]<-norm.Topt.s.lin(0.839/Rd.25.G31,ft.Topt.s.lin.G[10],ft.Topt.s.lin.G[11],ft.Topt.s.lin.G[12],ft.Topt.s.lin.G[13],28.5,G31c$Curve[i])
}
G31c$Rd<-p.G31.Rd

#Subset data depending on if they were collected in an ascending or descending order
G31c.d<-G31c[G31c$Direction == "down",]
G31c.u<-G31c[G31c$Direction == "up",]

#
#Fit curves
#

#Descending
m.G31c.d<-fitacis(G31c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31c.d)
coef(m.G31c.d)

#Ascending
m.G31c.u<-fitacis(G31c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.G31c.u)
coef(m.G31c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.G31c.d$`10`$df[,2]~m.G31c.d$`10`$df[,3])
summary(lm(m.G31c.d$`10`$df[,2]~m.G31c.d$`10`$df[,3]))
plot(m.G31c.d$`15`$df[,2]~m.G31c.d$`15`$df[,3])
summary(lm(m.G31c.d$`15`$df[,2]~m.G31c.d$`15`$df[,3]))
plot(m.G31c.d$`20`$df[,2]~m.G31c.d$`20`$df[,3])
summary(lm(m.G31c.d$`20`$df[,2]~m.G31c.d$`20`$df[,3]))
plot(m.G31c.d$`25`$df[,2]~m.G31c.d$`25`$df[,3])
summary(lm(m.G31c.d$`25`$df[,2]~m.G31c.d$`25`$df[,3]))
plot(m.G31c.d$`30`$df[,2]~m.G31c.d$`30`$df[,3])
summary(lm(m.G31c.d$`30`$df[,2]~m.G31c.d$`30`$df[,3]))
plot(m.G31c.d$`32`$df[,2]~m.G31c.d$`32`$df[,3])
summary(lm(m.G31c.d$`32`$df[,2]~m.G31c.d$`32`$df[,3]))

#Ascending
plot(m.G31c.u$`30`$df[,2]~m.G31c.u$`30`$df[,3])
summary(lm(m.G31c.u$`30`$df[,2]~m.G31c.u$`30`$df[,3]))
plot(m.G31c.u$`35`$df[,2]~m.G31c.u$`35`$df[,3])
summary(lm(m.G31c.u$`35`$df[,2]~m.G31c.u$`35`$df[,3]))
plot(m.G31c.u$`40`$df[,2]~m.G31c.u$`40`$df[,3])
summary(lm(m.G31c.u$`40`$df[,2]~m.G31c.u$`40`$df[,3]))

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
p.R21.Rd<-rep(NA,length(R21a$Curve))
for (i in 1:length(R21a$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21a$Curve[i])
}
R21a$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21a.d<-R21a[R21a$Direction == "down",]
R21a.u<-R21a[R21a$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21a.d<-fitacis(R21a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21a.d)
coef(m.R21a.d)

#Ascending
m.R21a.u<-fitacis(R21a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21a.u)
coef(m.R21a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21a.d$`10`$df[,2]~m.R21a.d$`10`$df[,3])
summary(lm(m.R21a.d$`10`$df[,2]~m.R21a.d$`10`$df[,3]))
plot(m.R21a.d$`15`$df[,2]~m.R21a.d$`15`$df[,3])
summary(lm(m.R21a.d$`15`$df[,2]~m.R21a.d$`15`$df[,3]))
plot(m.R21a.d$`20`$df[,2]~m.R21a.d$`20`$df[,3])
summary(lm(m.R21a.d$`20`$df[,2]~m.R21a.d$`20`$df[,3]))

#Ascending
plot(m.R21a.u$`20`$df[,2]~m.R21a.u$`20`$df[,3])
summary(lm(m.R21a.u$`20`$df[,2]~m.R21a.u$`20`$df[,3]))
plot(m.R21a.u$`25`$df[,2]~m.R21a.u$`25`$df[,3])
summary(lm(m.R21a.u$`25`$df[,2]~m.R21a.u$`25`$df[,3]))
plot(m.R21a.u$`30`$df[,2]~m.R21a.u$`30`$df[,3])
summary(lm(m.R21a.u$`30`$df[,2]~m.R21a.u$`30`$df[,3]))
plot(m.R21a.u$`35`$df[,2]~m.R21a.u$`35`$df[,3])
summary(lm(m.R21a.u$`35`$df[,2]~m.R21a.u$`35`$df[,3]))
plot(m.R21a.u$`40`$df[,2]~m.R21a.u$`40`$df[,3])
summary(lm(m.R21a.u$`40`$df[,2]~m.R21a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21b$Curve))
for (i in 1:length(R21b$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21b$Curve[i])
}
R21b$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21b.d<-R21b[R21b$Direction == "down",]
R21b.u<-R21b[R21b$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21b.d<-fitacis(R21b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21b.d)
coef(m.R21b.d)

#Ascending
m.R21b.u<-fitacis(R21b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21b.u)
coef(m.R21b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21b.d$`10`$df[,2]~m.R21b.d$`10`$df[,3])
summary(lm(m.R21b.d$`10`$df[,2]~m.R21b.d$`10`$df[,3]))
plot(m.R21b.d$`15`$df[,2]~m.R21b.d$`15`$df[,3])
summary(lm(m.R21b.d$`15`$df[,2]~m.R21b.d$`15`$df[,3]))
plot(m.R21b.d$`20`$df[,2]~m.R21b.d$`20`$df[,3])
summary(lm(m.R21b.d$`20`$df[,2]~m.R21b.d$`20`$df[,3]))

#Ascending
plot(m.R21b.u$`20`$df[,2]~m.R21b.u$`20`$df[,3])
summary(lm(m.R21b.u$`20`$df[,2]~m.R21b.u$`20`$df[,3]))
plot(m.R21b.u$`25`$df[,2]~m.R21b.u$`25`$df[,3])
summary(lm(m.R21b.u$`25`$df[,2]~m.R21b.u$`25`$df[,3]))
plot(m.R21b.u$`30`$df[,2]~m.R21b.u$`30`$df[,3])
summary(lm(m.R21b.u$`30`$df[,2]~m.R21b.u$`30`$df[,3]))
plot(m.R21b.u$`35`$df[,2]~m.R21b.u$`35`$df[,3])
summary(lm(m.R21b.u$`35`$df[,2]~m.R21b.u$`35`$df[,3]))
plot(m.R21b.u$`40`$df[,2]~m.R21b.u$`40`$df[,3])
summary(lm(m.R21b.u$`40`$df[,2]~m.R21b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21c$Curve))
for (i in 1:length(R21c$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21c$Curve[i])
}
R21c$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21c.d<-R21c[R21c$Direction == "down",]
R21c.u<-R21c[R21c$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21c.d<-fitacis(R21c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21c.d)
coef(m.R21c.d)

#Ascending
m.R21c.u<-fitacis(R21c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21c.u)
coef(m.R21c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21c.d$`10`$df[,2]~m.R21c.d$`10`$df[,3])
summary(lm(m.R21c.d$`10`$df[,2]~m.R21c.d$`10`$df[,3]))
plot(m.R21c.d$`15`$df[,2]~m.R21c.d$`15`$df[,3])
summary(lm(m.R21c.d$`15`$df[,2]~m.R21c.d$`15`$df[,3]))
plot(m.R21c.d$`20`$df[,2]~m.R21c.d$`20`$df[,3])
summary(lm(m.R21c.d$`20`$df[,2]~m.R21c.d$`20`$df[,3]))

#Ascending
plot(m.R21c.u$`20`$df[,2]~m.R21c.u$`20`$df[,3])
summary(lm(m.R21c.u$`20`$df[,2]~m.R21c.u$`20`$df[,3]))
plot(m.R21c.u$`25`$df[,2]~m.R21c.u$`25`$df[,3])
summary(lm(m.R21c.u$`25`$df[,2]~m.R21c.u$`25`$df[,3]))
plot(m.R21c.u$`30`$df[,2]~m.R21c.u$`30`$df[,3])
summary(lm(m.R21c.u$`30`$df[,2]~m.R21c.u$`30`$df[,3]))
plot(m.R21c.u$`35`$df[,2]~m.R21c.u$`35`$df[,3])
summary(lm(m.R21c.u$`35`$df[,2]~m.R21c.u$`35`$df[,3]))
plot(m.R21c.u$`40`$df[,2]~m.R21c.u$`40`$df[,3])
summary(lm(m.R21c.u$`40`$df[,2]~m.R21c.u$`40`$df[,3]))

##
#Replicate 4
##

#Calculate RL for each measurement temperature
p.R21.Rd<-rep(NA,length(R21d$Curve))
for (i in 1:length(R21d$Curve)){
  p.R21.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R21,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],18.5,R21d$Curve[i])
}
R21d$Rd<-p.R21.Rd

#Subset data depending on if they were collected in an ascending or descending order
R21d.d<-R21d[R21d$Direction == "down",]
R21d.u<-R21d[R21d$Direction == "up",]

#
#Fit curves
#

#Descending
m.R21d.d<-fitacis(R21d.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21d.d)
coef(m.R21d.d)

#Ascending
m.R21d.u<-fitacis(R21d.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R21d.u)
coef(m.R21d.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R21d.d$`10`$df[,2]~m.R21d.d$`10`$df[,3])
summary(lm(m.R21d.d$`10`$df[,2]~m.R21d.d$`10`$df[,3]))
plot(m.R21d.d$`15`$df[,2]~m.R21d.d$`15`$df[,3])
summary(lm(m.R21d.d$`15`$df[,2]~m.R21d.d$`15`$df[,3]))
plot(m.R21d.d$`20`$df[,2]~m.R21d.d$`20`$df[,3])
summary(lm(m.R21d.d$`20`$df[,2]~m.R21d.d$`20`$df[,3]))

#Ascending
plot(m.R21d.u$`20`$df[,2]~m.R21d.u$`20`$df[,3])
summary(lm(m.R21d.u$`20`$df[,2]~m.R21d.u$`20`$df[,3]))
plot(m.R21d.u$`25`$df[,2]~m.R21d.u$`25`$df[,3])
summary(lm(m.R21d.u$`25`$df[,2]~m.R21d.u$`25`$df[,3]))
plot(m.R21d.u$`30`$df[,2]~m.R21d.u$`30`$df[,3])
summary(lm(m.R21d.u$`30`$df[,2]~m.R21d.u$`30`$df[,3]))
plot(m.R21d.u$`35`$df[,2]~m.R21d.u$`35`$df[,3])
summary(lm(m.R21d.u$`35`$df[,2]~m.R21d.u$`35`$df[,3]))
plot(m.R21d.u$`40`$df[,2]~m.R21d.u$`40`$df[,3])
summary(lm(m.R21d.u$`40`$df[,2]~m.R21d.u$`40`$df[,3]))

###
#26:20 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.R26.Rd<-rep(NA,length(R26a$Curve))
for (i in 1:length(R26a$Curve)){
  p.R26.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,R26a$Curve[i])
}
R26a$Rd<-p.R26.Rd

#Subset data depending on if they were collected in an ascending or descending order
R26a.d<-R26a[R26a$Direction == "down",]
R26a.u<-R26a[R26a$Direction == "up",]

#
#Fit curves
#

#Descending
m.R26a.d<-fitacis(R26a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26a.d)
coef(m.R26a.d)

#Ascending
m.R26a.u<-fitacis(R26a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26a.u)
coef(m.R26a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R26a.d$`10`$df[,2]~m.R26a.d$`10`$df[,3])
summary(lm(m.R26a.d$`10`$df[,2]~m.R26a.d$`10`$df[,3]))
plot(m.R26a.d$`15`$df[,2]~m.R26a.d$`15`$df[,3])
summary(lm(m.R26a.d$`15`$df[,2]~m.R26a.d$`15`$df[,3]))
plot(m.R26a.d$`20`$df[,2]~m.R26a.d$`20`$df[,3])
summary(lm(m.R26a.d$`20`$df[,2]~m.R26a.d$`20`$df[,3]))
plot(m.R26a.d$`25`$df[,2]~m.R26a.d$`25`$df[,3])
summary(lm(m.R26a.d$`25`$df[,2]~m.R26a.d$`25`$df[,3]))

#Ascending
plot(m.R26a.u$`25`$df[,2]~m.R26a.u$`25`$df[,3])
summary(lm(m.R26a.u$`25`$df[,2]~m.R26a.u$`25`$df[,3]))
plot(m.R26a.u$`30`$df[,2]~m.R26a.u$`30`$df[,3])
summary(lm(m.R26a.u$`30`$df[,2]~m.R26a.u$`30`$df[,3]))
plot(m.R26a.u$`35`$df[,2]~m.R26a.u$`35`$df[,3])
summary(lm(m.R26a.u$`35`$df[,2]~m.R26a.u$`35`$df[,3]))
plot(m.R26a.u$`40`$df[,2]~m.R26a.u$`40`$df[,3])
summary(lm(m.R26a.u$`40`$df[,2]~m.R26a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.R26.Rd<-rep(NA,length(R26b$Curve))
for (i in 1:length(R26b$Curve)){
  p.R26.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,R26b$Curve[i])
}
R26b$Rd<-p.R26.Rd

#Subset data depending on if they were collected in an ascending or descending order
R26b.d<-R26b[R26b$Direction == "down",]
R26b.u<-R26b[R26b$Direction == "up",]

#
#Fit curves
#

#Descending
m.R26b.d<-fitacis(R26b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26b.d)
coef(m.R26b.d)

#Ascending
m.R26b.u<-fitacis(R26b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26b.u)
coef(m.R26b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R26b.d$`10`$df[,2]~m.R26b.d$`10`$df[,3])
summary(lm(m.R26b.d$`10`$df[,2]~m.R26b.d$`10`$df[,3]))
plot(m.R26b.d$`15`$df[,2]~m.R26b.d$`15`$df[,3])
summary(lm(m.R26b.d$`15`$df[,2]~m.R26b.d$`15`$df[,3]))
plot(m.R26b.d$`20`$df[,2]~m.R26b.d$`20`$df[,3])
summary(lm(m.R26b.d$`20`$df[,2]~m.R26b.d$`20`$df[,3]))
plot(m.R26b.d$`25`$df[,2]~m.R26b.d$`25`$df[,3])
summary(lm(m.R26b.d$`25`$df[,2]~m.R26b.d$`25`$df[,3]))

#Ascending
plot(m.R26b.u$`25`$df[,2]~m.R26b.u$`25`$df[,3])
summary(lm(m.R26b.u$`25`$df[,2]~m.R26b.u$`25`$df[,3]))
plot(m.R26b.u$`30`$df[,2]~m.R26b.u$`30`$df[,3])
summary(lm(m.R26b.u$`30`$df[,2]~m.R26b.u$`30`$df[,3]))
plot(m.R26b.u$`35`$df[,2]~m.R26b.u$`35`$df[,3])
summary(lm(m.R26b.u$`35`$df[,2]~m.R26b.u$`35`$df[,3]))
plot(m.R26b.u$`40`$df[,2]~m.R26b.u$`40`$df[,3])
summary(lm(m.R26b.u$`40`$df[,2]~m.R26b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.R26.Rd<-rep(NA,length(R26c$Curve))
for (i in 1:length(R26c$Curve)){
  p.R26.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R26,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],23.5,R26c$Curve[i])
}
R26c$Rd<-p.R26.Rd

#Subset data depending on if they were collected in an ascending or descending order
R26c.d<-R26c[R26c$Direction == "down",] 
R26c.u<-R26c[R26c$Direction == "up",]

#
#Fit curves
#

#Descending
m.R26c.d<-fitacis(R26c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26c.d)
coef(m.R26c.d) 

#Ascending
m.R26c.u<-fitacis(R26c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R26c.u)
coef(m.R26c.u) 

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R26c.d$`10`$df[,2]~m.R26c.d$`10`$df[,3])
summary(lm(m.R26c.d$`10`$df[,2]~m.R26c.d$`10`$df[,3]))
plot(m.R26c.d$`15`$df[,2]~m.R26c.d$`15`$df[,3])
summary(lm(m.R26c.d$`15`$df[,2]~m.R26c.d$`15`$df[,3]))
plot(m.R26c.d$`20`$df[,2]~m.R26c.d$`20`$df[,3])
summary(lm(m.R26c.d$`20`$df[,2]~m.R26c.d$`20`$df[,3]))
plot(m.R26c.d$`25`$df[,2]~m.R26c.d$`25`$df[,3])
summary(lm(m.R26c.d$`25`$df[,2]~m.R26c.d$`25`$df[,3]))

#Ascending
plot(m.R26c.u$`25`$df[,2]~m.R26c.u$`25`$df[,3])
summary(lm(m.R26c.u$`25`$df[,2]~m.R26c.u$`25`$df[,3]))
plot(m.R26c.u$`30`$df[,2]~m.R26c.u$`30`$df[,3])
summary(lm(m.R26c.u$`30`$df[,2]~m.R26c.u$`30`$df[,3]))
plot(m.R26c.u$`35`$df[,2]~m.R26c.u$`35`$df[,3])
summary(lm(m.R26c.u$`35`$df[,2]~m.R26c.u$`35`$df[,3]))
plot(m.R26c.u$`40`$df[,2]~m.R26c.u$`40`$df[,3])
summary(lm(m.R26c.u$`40`$df[,2]~m.R26c.u$`40`$df[,3]))

###
#31:25 deg. C growing temperature
###

##
#Replicate 1
##

#Calculate RL for each measurement temperature
p.R31.Rd<-rep(NA,length(R31a$Curve))
for (i in 1:length(R31a$Curve)){
  p.R31.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,R31a$Curve[i])
}
R31a$Rd<-p.R31.Rd

#Subset data depending on if they were collected in an ascending or descending order
R31a.d<-R31a[R31a$Direction == "down",]
R31a.u<-R31a[R31a$Direction == "up",]

#
#Fit curves
#

#Descending
m.R31a.d<-fitacis(R31a.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31a.d)
coef(m.R31a.d)

#Ascending
m.R31a.u<-fitacis(R31a.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31a.u)
coef(m.R31a.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R31a.d$`10`$df[,2]~m.R31a.d$`10`$df[,3])
summary(lm(m.R31a.d$`10`$df[,2]~m.R31a.d$`10`$df[,3]))
plot(m.R31a.d$`15`$df[,2]~m.R31a.d$`15`$df[,3])
summary(lm(m.R31a.d$`15`$df[,2]~m.R31a.d$`15`$df[,3]))
plot(m.R31a.d$`20`$df[,2]~m.R31a.d$`20`$df[,3])
summary(lm(m.R31a.d$`20`$df[,2]~m.R31a.d$`20`$df[,3]))
plot(m.R31a.d$`25`$df[,2]~m.R31a.d$`25`$df[,3])
summary(lm(m.R31a.d$`25`$df[,2]~m.R31a.d$`25`$df[,3]))
plot(m.R31a.d$`30`$df[,2]~m.R31a.d$`30`$df[,3])
summary(lm(m.R31a.d$`30`$df[,2]~m.R31a.d$`30`$df[,3]))

#Ascending
plot(m.R31a.u$`30`$df[,2]~m.R31a.u$`30`$df[,3])
summary(lm(m.R31a.u$`30`$df[,2]~m.R31a.u$`30`$df[,3]))
plot(m.R31a.u$`35`$df[,2]~m.R31a.u$`35`$df[,3])
summary(lm(m.R31a.u$`35`$df[,2]~m.R31a.u$`35`$df[,3]))
plot(m.R31a.u$`40`$df[,2]~m.R31a.u$`40`$df[,3])
summary(lm(m.R31a.u$`40`$df[,2]~m.R31a.u$`40`$df[,3]))

##
#Replicate 2
##

#Calculate RL for each measurement temperature
p.R31.Rd<-rep(NA,length(R31b$Curve))
for (i in 1:length(R31b$Curve)){
  p.R31.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,R31b$Curve[i])
}
R31b$Rd<-p.R31.Rd

#Subset data depending on if they were collected in an ascending or descending order
R31b.d<-R31b[R31b$Direction == "down",]
R31b.u<-R31b[R31b$Direction == "up",]

#
#Fit curves
#

#Descending
m.R31b.d<-fitacis(R31b.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31b.d)
coef(m.R31b.d)

#Ascending
m.R31b.u<-fitacis(R31b.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31b.u)
coef(m.R31b.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R31b.d$`10`$df[,2]~m.R31b.d$`10`$df[,3])
summary(lm(m.R31b.d$`10`$df[,2]~m.R31b.d$`10`$df[,3]))
plot(m.R31b.d$`15`$df[,2]~m.R31b.d$`15`$df[,3])
summary(lm(m.R31b.d$`15`$df[,2]~m.R31b.d$`15`$df[,3]))
plot(m.R31b.d$`20`$df[,2]~m.R31b.d$`20`$df[,3])
summary(lm(m.R31b.d$`20`$df[,2]~m.R31b.d$`20`$df[,3]))
plot(m.R31b.d$`25`$df[,2]~m.R31b.d$`25`$df[,3])
summary(lm(m.R31b.d$`25`$df[,2]~m.R31b.d$`25`$df[,3]))
plot(m.R31b.d$`30`$df[,2]~m.R31b.d$`30`$df[,3])
summary(lm(m.R31b.d$`30`$df[,2]~m.R31b.d$`30`$df[,3]))

#Ascending
plot(m.R31b.u$`30`$df[,2]~m.R31b.u$`30`$df[,3])
summary(lm(m.R31b.u$`30`$df[,2]~m.R31b.u$`30`$df[,3]))
plot(m.R31b.u$`35`$df[,2]~m.R31b.u$`35`$df[,3])
summary(lm(m.R31b.u$`35`$df[,2]~m.R31b.u$`35`$df[,3]))
plot(m.R31b.u$`40`$df[,2]~m.R31b.u$`40`$df[,3])
summary(lm(m.R31b.u$`40`$df[,2]~m.R31b.u$`40`$df[,3]))

##
#Replicate 3
##

#Calculate RL for each measurement temperature
p.R31.Rd<-rep(NA,length(R31c$Curve))
for (i in 1:length(R31c$Curve)){
  p.R31.Rd[i]<-norm.Topt.s.lin(0.983/Rd.25.R31,ft.Topt.s.lin.R[10],ft.Topt.s.lin.R[11],ft.Topt.s.lin.R[12],ft.Topt.s.lin.R[13],28.5,R31c$Curve[i])
}
R31c$Rd<-p.R31.Rd

#Subset data depending on if they were collected in an ascending or descending order
R31c.d<-R31c[R31c$Direction == "down",]
R31c.u<-R31c[R31c$Direction == "up",]

#
#Fit curves
#

#Descending
m.R31c.d<-fitacis(R31c.d,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31c.d)
coef(m.R31c.d)

#Ascending
m.R31c.u<-fitacis(R31c.u,"Curve",useRd=TRUE,fitTPU=T,Tcorrect=F)
plot(m.R31c.u)
coef(m.R31c.u)

#
#Calculate R2 of fitted vs. measured assimilation values
#

#Descending
plot(m.R31c.d$`10`$df[,2]~m.R31c.d$`10`$df[,3])
summary(lm(m.R31c.d$`10`$df[,2]~m.R31c.d$`10`$df[,3]))
plot(m.R31c.d$`15`$df[,2]~m.R31c.d$`15`$df[,3])
summary(lm(m.R31c.d$`15`$df[,2]~m.R31c.d$`15`$df[,3]))
plot(m.R31c.d$`20`$df[,2]~m.R31c.d$`20`$df[,3])
summary(lm(m.R31c.d$`20`$df[,2]~m.R31c.d$`20`$df[,3]))
plot(m.R31c.d$`25`$df[,2]~m.R31c.d$`25`$df[,3])
summary(lm(m.R31c.d$`25`$df[,2]~m.R31c.d$`25`$df[,3]))
plot(m.R31c.d$`30`$df[,2]~m.R31c.d$`30`$df[,3])
summary(lm(m.R31c.d$`30`$df[,2]~m.R31c.d$`30`$df[,3]))

#Ascending
plot(m.R31c.u$`30`$df[,2]~m.R31c.u$`30`$df[,3])
summary(lm(m.R31c.u$`30`$df[,2]~m.R31c.u$`30`$df[,3]))
plot(m.R31c.u$`35`$df[,2]~m.R31c.u$`35`$df[,3])
summary(lm(m.R31c.u$`35`$df[,2]~m.R31c.u$`35`$df[,3]))
plot(m.R31c.u$`40`$df[,2]~m.R31c.u$`40`$df[,3])
summary(lm(m.R31c.u$`40`$df[,2]~m.R31c.u$`40`$df[,3]))

###
#Data which were removed
###

#Removed because R2 of fit <0.9
#m.A26c.u at 40 deg. C
#m.G21c.d all temperatures
#m.G26d.u 40 deg. C
#m.R21a.d at 10 deg. C
#m.R21a.u at 35 deg. C

#Removed based on Vcmax SE >40%
#m.M26b.u at 40 deg. C
#m.A26b.u at 40 deg. C
#m.A26c.u at 40 deg. C
#m.A31a.u at 40 deg. C
#m.G21c.u at 30 deg. C
#m.G26d.u at 40 deg. C (Vcmax SE is NaN)
#m.G31b.u at 40 deg. C

#Removed because they were not fit
#m.M26b.d Jmax at 10 deg. C
#m.M31a.u Jmax at 40 deg. C
#m.M31b.d Jmax at 10 deg. C
#m.A21a.d Jmax at 15 deg. C
#m.A21b.d Jmax at 15 deg. C
#m.G21c.u Jmax at 25 deg. C
#m.G26b.d Jmax at 10 deg. C
#m.G26c.d Jmax at 10 deg. C
#m.G31c.d Jmax at 10 and 15 deg. C
#m.R21c.d Jmax at 10 deg. C
#m.R26a.d Jmax at 10 deg. C
#m.R26b.d Jmax at 10 deg. C

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
M21a.A275<-c(m.M21a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M21a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M21a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M21a.u$`20`$Photosyn(Ci=275)[[2]],
                   m.M21a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M21a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M21a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M21a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
M21b.A275<-c(m.M21b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M21b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M21b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M21b.u$`20`$Photosyn(Ci=275)[[2]],
                   m.M21b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M21b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M21b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M21b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
M21c.A275<-c(m.M21c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M21c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M21c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M21c.u$`20`$Photosyn(Ci=275)[[2]],
                   m.M21c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M21c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M21c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M21c.u$`40`$Photosyn(Ci=275)[[2]])

M21.A275<-c(M21a.A275,M21b.A275,M21c.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
M26a.A275<-c(m.M26a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M26a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M26a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M26a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M26a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M26a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M26a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M26a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
M26b.A275<-c(m.M26b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M26b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M26b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M26b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M26b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M26b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M26b.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 3
M26c.A275<-c(m.M26c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M26c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M26c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M26c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M26c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.M26c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M26c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M26c.u$`40`$Photosyn(Ci=275)[[2]])

M26.A275<-c(M26a.A275,M26b.A275,M26c.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
M31a.A275<-c(m.M31a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M31a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M31a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M31a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M31a.d$`30`$Photosyn(Ci=275)[[2]],
                   m.M31a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M31a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M31a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
M31b.A275<-c(m.M31b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M31b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M31b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M31b.d$`30`$Photosyn(Ci=275)[[2]],
                   m.M31b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M31b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M31b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
M31c.A275<-c(m.M31c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.M31c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.M31c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.M31c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.M31c.d$`30`$Photosyn(Ci=275)[[2]],
                   m.M31c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.M31c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.M31c.u$`40`$Photosyn(Ci=275)[[2]])

M31.A275<-c(M31a.A275,M31b.A275,M31c.A275)

##
#Alnus
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
A21a.A275<-c(m.A21a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A21a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A21a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A21a.u$`20`$Photosyn(Ci=275)[[2]],
                   m.A21a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A21a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A21a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A21a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
A21b.A275<-c(m.A21b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A21b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A21b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A21b.u$`20`$Photosyn(Ci=275)[[2]],
                   m.A21b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A21b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A21b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A21b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
A21c.A275<-c(m.A21c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A21c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A21c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A21c.u$`20`$Photosyn(Ci=275)[[2]],
                   m.A21c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A21c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A21c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A21c.u$`40`$Photosyn(Ci=275)[[2]])

A21.A275<-c(A21a.A275,A21b.A275,A21c.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
A26a.A275<-c(m.A26a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A26a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A26a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A26a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A26a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A26a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A26a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A26a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
A26b.A275<-c(m.A26b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A26b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A26b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A26b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A26b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A26b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A26b.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 3
A26c.A275<-c(m.A26c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A26c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A26c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A26c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A26c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.A26c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A26c.u$`35`$Photosyn(Ci=275)[[2]])

A26.A275<-c(A26a.A275,A26b.A275,A26c.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
A31a.A275<-c(m.A31a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A31a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A31a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A31a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A31a.d$`30`$Photosyn(Ci=275)[[2]],
                   m.A31a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A31a.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 2
A31b.A275<-c(m.A31b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A31b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A31b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A31b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A31b.d$`30`$Photosyn(Ci=275)[[2]],
                   m.A31b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A31b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A31b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
A31c.A275<-c(m.A31c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.A31c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.A31c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.A31c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.A31c.d$`30`$Photosyn(Ci=275)[[2]],
                   m.A31c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.A31c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.A31c.u$`40`$Photosyn(Ci=275)[[2]])

A31.A275<-c(A31a.A275,A31b.A275,A31c.A275)

##
#Gliricidia
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
G21a.A275<-c(m.G21a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G21a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G21a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G21a.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G21a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G21a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
G21c.A275<-c(m.G21c.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21c.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
G21d.A275<-c(m.G21d.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G21d.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G21d.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G21d.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21d.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G21d.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G21d.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21d.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 4
G21b.A275<-c(m.G21b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G21b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G21b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G21b.u$`20`$Photosyn(Ci=275)[[2]],
                   m.G21b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G21b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G21b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G21b.u$`40`$Photosyn(Ci=275)[[2]])

G21.A275<-c(G21a.A275,G21b.A275,G21c.A275,G21d.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
G26a.A275<-c(m.G26a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G26a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G26a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G26a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G26a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G26a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
G26b.A275<-c(m.G26b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G26b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G26b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G26b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G26b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G26b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
G26c.A275<-c(m.G26c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G26c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G26c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G26c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G26c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G26c.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 4
G26d.A275<-c(m.G26d.u$`25`$Photosyn(Ci=275)[[2]],
                   m.G26d.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G26d.u$`35`$Photosyn(Ci=275)[[2]])

G26.A275<-c(G26a.A275,G26b.A275,G26c.A275,G26d.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
G31a.A275<-c(m.G31a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G31a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G31a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G31a.d$`30`$Photosyn(Ci=275)[[2]],
                   m.G31a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G31a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G31a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
G31b.A275<-c(m.G31b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G31b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G31b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G31b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G31b.d$`30`$Photosyn(Ci=275)[[2]],
                   m.G31b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G31b.u$`35`$Photosyn(Ci=275)[[2]])

#Replicate 3
G31c.A275<-c(m.G31c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.G31c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.G31c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.G31c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.G31c.d$`30`$Photosyn(Ci=275)[[2]],
                   m.G31c.d$`32`$Photosyn(Ci=275)[[2]],
                   m.G31c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.G31c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.G31c.u$`40`$Photosyn(Ci=275)[[2]])

G31.A275<-c(G31a.A275,G31b.A275,G31c.A275)

##
#Robinia
##

#
#21:15 deg. C growing temperature
#

#Replicate 1
R21a.A275<-c(m.R21a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21a.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
R21b.A275<-c(m.R21b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R21b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21b.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R21b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
R21c.A275<-c(m.R21c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R21c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21c.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R21c.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 4
R21d.A275<-c(m.R21d.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R21d.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R21d.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R21d.u$`20`$Photosyn(Ci=275)[[2]],
                   m.R21d.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R21d.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R21d.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R21d.u$`40`$Photosyn(Ci=275)[[2]])

R21.A275<-c(R21a.A275,R21b.A275,R21c.A275,R21d.A275)

#
#26:20 deg. C growing temperature
#

#Replicate 1
R26a.A275<-c(m.R26a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R26a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R26a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R26a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R26a.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R26a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R26a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R26a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
R26b.A275<-c(m.R26b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R26b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R26b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R26b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R26b.d$`40`$Photosyn(Ci=275)[[2]],
                   m.R26b.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R26b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R26b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R26b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
R26c.A275<-c(m.R26c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R26c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R26c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R26c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R26c.u$`25`$Photosyn(Ci=275)[[2]],
                   m.R26c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R26c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R26c.u$`40`$Photosyn(Ci=275)[[2]])

R26.A275<-c(R26a.A275,R26b.A275,R26c.A275)

#
#31:25 deg. C growing temperature
#

#Replicate 1
R31a.A275<-c(m.R31a.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R31a.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R31a.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R31a.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R31a.d$`30`$Photosyn(Ci=275)[[2]],
                   m.R31a.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R31a.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R31a.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 2
R31b.A275<-c(m.R31b.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R31b.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R31b.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R31b.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R31b.d$`30`$Photosyn(Ci=275)[[2]],
                   m.R31b.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R31b.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R31b.u$`40`$Photosyn(Ci=275)[[2]])

#Replicate 3
R31c.A275<-c(m.R31c.d$`10`$Photosyn(Ci=275)[[2]],
                   m.R31c.d$`15`$Photosyn(Ci=275)[[2]],
                   m.R31c.d$`20`$Photosyn(Ci=275)[[2]],
                   m.R31c.d$`25`$Photosyn(Ci=275)[[2]],
                   m.R31c.d$`30`$Photosyn(Ci=275)[[2]],
                   m.R31c.u$`30`$Photosyn(Ci=275)[[2]],
                   m.R31c.u$`35`$Photosyn(Ci=275)[[2]],
                   m.R31c.u$`40`$Photosyn(Ci=275)[[2]])

R31.A275<-c(R31a.A275,R31b.A275,R31c.A275)

####
#The following assembles the data frames that are used to model A275, Vcmax, and Jmax as functions of temperature
####

###
#Morella
###

##
#21:15 deg. C growing temperature
##

dat.V.M21<-cbind(c(coef(m.M21a.d)[,2],coef(m.M21a.u)[,2],coef(m.M21b.d)[,2],
                   coef(m.M21b.u)[,2],coef(m.M21c.d)[,2],coef(m.M21c.u)[,2]))
dat.J.M21<-cbind(c(coef(m.M21a.d)[,3],coef(m.M21a.u)[,3],coef(m.M21b.d)[,3],
                   coef(m.M21b.u)[,3],coef(m.M21c.d)[,3],coef(m.M21c.u)[,3]))
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

dat.V.M26<-cbind(c(coef(m.M26a.d)[,2],coef(m.M26a.u)[,2],
                   coef(m.M26b.d)[,2],coef(m.M26b.u)[1:3,2],
                   coef(m.M26c.d)[,2],coef(m.M26c.u)[,2]))
dat.J.M26<-cbind(c(coef(m.M26a.d)[,3],coef(m.M26a.u)[,3],
                   NA,coef(m.M26b.d)[2:4,3],coef(m.M26b.u)[1:3,3],
                   coef(m.M26c.d)[,3],coef(m.M26c.u)[,3]))
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

dat.V.M31<-cbind(c(coef(m.M31a.d)[,2],coef(m.M31a.u)[,2],coef(m.M31b.d)[,2],
                   coef(m.M31b.u)[,2],coef(m.M31c.d)[,2],coef(m.M31c.u)[,2]))
dat.J.M31<-cbind(c(coef(m.M31a.d)[,3],coef(m.M31a.u)[1:2,3],NA,NA,coef(m.M31b.d)[2:4,3],
                   coef(m.M31b.u)[,3],coef(m.M31c.d)[,3],coef(m.M31c.u)[,3]))
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

dat.V.A21<-cbind(c(coef(m.A21a.d)[,2],coef(m.A21a.u)[,2],coef(m.A21b.d)[,2],
                   coef(m.A21b.u)[,2],coef(m.A21c.d)[,2],coef(m.A21c.u)[,2]))
dat.J.A21<-cbind(c(coef(m.A21a.d)[1,3],NA,coef(m.A21a.d)[3,3],coef(m.A21a.u)[,3],NA,coef(m.A21b.d)[2:3,3],
                   coef(m.A21b.u)[,3],coef(m.A21c.d)[,3],coef(m.A21c.u)[,3]))
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

dat.V.A26<-cbind(c(coef(m.A26a.d)[,2],coef(m.A26a.u)[,2],
                   coef(m.A26b.d)[,2],coef(m.A26b.u)[1:3,2],
                   coef(m.A26c.d)[,2],coef(m.A26c.u)[1:3,2]))
dat.J.A26<-cbind(c(coef(m.A26a.d)[,3],coef(m.A26a.u)[,3],
                   coef(m.A26b.d)[,3],coef(m.A26b.u)[1:3,3],
                   coef(m.A26c.d)[,3],coef(m.A26c.u)[1:3,3]))
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

dat.V.A31<-cbind(c(coef(m.A31a.d)[,2],coef(m.A31a.u)[1:2,2],coef(m.A31b.d)[,2],
                   coef(m.A31b.u)[,2],coef(m.A31c.d)[,2],coef(m.A31c.u)[,2]))
dat.J.A31<-cbind(c(coef(m.A31a.d)[,3],coef(m.A31a.u)[1:2,3],coef(m.A31b.d)[,3],
                   coef(m.A31b.u)[,3],coef(m.A31c.d)[,3],coef(m.A31c.u)[,3]))
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

dat.V.G21<-cbind(c(coef(m.G21a.d)[,2],coef(m.G21a.u)[,2],coef(m.G21b.d)[,2],
                   coef(m.G21b.u)[,2],coef(m.G21c.u)[1,2],coef(m.G21c.u)[4:5,2],coef(m.G21d.d)[,2],coef(m.G21d.u)[,2]))
dat.J.G21<-cbind(c(coef(m.G21a.d)[,3],coef(m.G21a.u)[,3],coef(m.G21b.d)[,3],
                   coef(m.G21b.u)[,3],coef(m.G21c.u)[1,3],coef(m.G21c.u)[4:5,3],coef(m.G21d.d)[,3],coef(m.G21d.u)[,3]))
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

dat.V.G26<-cbind(c(coef(m.G26a.d)[,2],coef(m.G26a.u)[,2],
                   coef(m.G26b.d)[,2],coef(m.G26b.u)[,2],
                   coef(m.G26c.d)[,2],coef(m.G26c.u)[,2],
                   coef(m.G26d.u)[1:3,2]))
dat.J.G26<-cbind(c(coef(m.G26a.d)[,3],coef(m.G26a.u)[,3],
                   NA,coef(m.G26b.d)[2:4,3],coef(m.G26b.u)[,3],
                   NA,coef(m.G26c.d)[2:4,3],coef(m.G26c.u)[,3],
                   coef(m.G26d.u)[1:3,3]))
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

dat.V.G31<-cbind(c(coef(m.G31a.d)[,2],coef(m.G31a.u)[,2],
                   coef(m.G31b.d)[,2],coef(m.G31b.u)[1:2,2],
                   coef(m.G31c.d)[,2],coef(m.G31c.u)[,2]))
dat.J.G31<-cbind(c(coef(m.G31a.d)[,3],coef(m.G31a.u)[,3],
                   coef(m.G31b.d)[,3],coef(m.G31b.u)[1:2,3],
                   NA,NA,coef(m.G31c.d)[3:6,3],coef(m.G31c.u)[,3]))
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

dat.V.R21<-cbind(c(coef(m.R21a.d)[2:3,2],coef(m.R21a.u)[1:3,2],coef(m.R21a.u)[5,2],
                   coef(m.R21b.d)[,2],coef(m.R21b.u)[,2],
                   coef(m.R21c.d)[,2],coef(m.R21c.u)[,2],
                   coef(m.R21d.d)[,2],coef(m.R21d.u)[,2]))
dat.J.R21<-cbind(c(coef(m.R21a.d)[2:3,3],coef(m.R21a.u)[1:3,3],coef(m.R21a.u)[5,3],
                   coef(m.R21b.d)[,3],coef(m.R21b.u)[,3],
                   NA,coef(m.R21c.d)[2:3,2],coef(m.R21c.u)[,3],
                   coef(m.R21d.d)[,3],coef(m.R21d.u)[,3]))
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

dat.V.R26<-cbind(c(coef(m.R26a.d)[,2],coef(m.R26a.u)[,2],
                   coef(m.R26b.d)[,2],coef(m.R26b.u)[,2],
                   coef(m.R26c.d)[,2],coef(m.R26c.u)[,2]))
dat.J.R26<-cbind(c(NA,coef(m.R26a.d)[2:4,2],coef(m.R26a.u)[,3],
                   NA,coef(m.R26b.d)[2:5,3],coef(m.R26b.u)[,3],
                   coef(m.R26c.d)[,3],coef(m.R26c.u)[,3]))
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

dat.V.R31<-cbind(c(coef(m.R31a.d)[,2],coef(m.R31a.u)[,2],coef(m.R31b.d)[,2],
                   coef(m.R31b.u)[,2],coef(m.R31c.d)[,2],coef(m.R31c.u)[,2]))
dat.J.R31<-cbind(c(coef(m.R31a.d)[,3],coef(m.R31a.u)[,3],coef(m.R31b.d)[,3],
                   coef(m.R31b.u)[,3],coef(m.R31c.d)[,3],coef(m.R31c.u)[,3]))
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