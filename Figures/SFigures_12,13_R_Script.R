###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figures 12 and 13
###############################################################################################################
###############################################################################################################

####
#Run Temporal_Decline_EqS7_Parameter_Calculations.R prior to this script
####

####
#All figures were exported as PDFs from R Studio
#Dimensions for each figure are given with script specific to each figure
####

####
#Load Additional Package
####

library(MASS)

####
#Separate parameter estimates by treatment (species x growing temperature)
####

M21<-temp.parms.rm[temp.parms.rm$Tr=="M21",]
A21<-temp.parms.rm[temp.parms.rm$Tr=="A21",]
G21<-temp.parms.rm[temp.parms.rm$Tr=="G21",]
R21<-temp.parms.rm[temp.parms.rm$Tr=="R21",]
M26<-temp.parms.rm[temp.parms.rm$Tr=="M26",]
A26<-temp.parms.rm[temp.parms.rm$Tr=="A26",]
G26<-temp.parms.rm[temp.parms.rm$Tr=="G26",]
R26<-temp.parms.rm[temp.parms.rm$Tr=="R26",]
M31<-temp.parms.rm[temp.parms.rm$Tr=="M31",]
A31<-temp.parms.rm[temp.parms.rm$Tr=="A31",]
G31<-temp.parms.rm[temp.parms.rm$Tr=="G31",]
R31<-temp.parms.rm[temp.parms.rm$Tr=="R31",]

####
#Additional function
####

#Equation S7 for how N fixation drops as a function of exposure time and temperature
#(b and d are both a function of temperature; see functions below)
fall.norm<-function(b,c,d,t){
  y<-((1+b)/(1+b*d))*((1-d)/(1+b*exp(c*t))+d)
  y
}

####
#Vector of tau values
####
xv<-seq(-10,10,0.1)

####
#Calculate 95% CI of b and d as functions of tau
####

#b
vmat.b=mvrnorm(10000,mu=coef(b.mle2),Sigma=vcov(b.mle2))
low.b<-numeric(length(xv))
high.b<-numeric(length(xv))
dist.b=numeric(10000)
for(j in 1:201){
  for(i in 1:10000){
    dist.b[i]<-exp.mod(vmat.b[i,3],vmat.b[i,2],xv[j])
  }
  low.b[j]=quantile(dist.b,0.025)
  high.b[j]=quantile(dist.b,0.975)
}

#d
vmat.d=mvrnorm(10000,mu=coef(d.mle),Sigma=vcov(d.mle))
low.d<-numeric(length(xv))
high.d<-numeric(length(xv))
dist.d=numeric(10000)
for(j in 1:201){
  for(i in 1:10000){
    dist.d[i]<-nsig.mod(vmat.d[i,2],vmat.d[i,3],xv[j])
  }
  low.d[j]=quantile(dist.d,0.025)
  high.d[j]=quantile(dist.d,0.975)
}

####
#95% CI calculations for Supplementary Figure 13
#Skip this section if only plotting Supplementary Figure 12
####

#Time vector
time.sim=seq(0,6,0.01)

###
#21:15 deg. C growing temperature
###

##
#Morella
##

#Tau = -3.43 deg. C
low.n3.43<-numeric(length(time.sim))
high.n3.43<-numeric(length(time.sim))
dist.n3.43=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n3.43[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-3.43))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-3.43),time.sim[j])
  }
  low.n3.43[j]=quantile(na.omit(dist.n3.43),0.025)
  high.n3.43[j]=quantile(na.omit(dist.n3.43),0.975)
}

#Tau = 0.97 deg. C
low.p0.97<-numeric(length(time.sim))
high.p0.97<-numeric(length(time.sim))
dist.p0.97=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p0.97[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],0.97))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],0.97),time.sim[j])
  }
  low.p0.97[j]=quantile(na.omit(dist.p0.97),0.025)
  high.p0.97[j]=quantile(na.omit(dist.p0.97),0.975)
}

#Tau = 5.87 deg. C
low.p5.87<-numeric(length(time.sim))
high.p5.87<-numeric(length(time.sim))
dist.p5.87=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p5.87[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],5.87))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],5.87),time.sim[j])
  }
  low.p5.87[j]=quantile(na.omit(dist.p5.87),0.025)
  high.p5.87[j]=quantile(na.omit(dist.p5.87),0.975)
}

#Tau = 10.07 deg. C
low.p10.07<-numeric(length(time.sim))
high.p10.07<-numeric(length(time.sim))
dist.p10.07=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p10.07[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],10.07))-1,c.est,
                              nsig.mod(vmat.d[i,2],vmat.d[i,3],10.07),time.sim[j])
  }
  low.p10.07[j]=quantile(na.omit(dist.p10.07),0.025)
  high.p10.07[j]=quantile(na.omit(dist.p10.07),0.975)
}

##
#Alnus
##

#Tau = -7.72 deg. C
low.n7.72<-numeric(length(time.sim))
high.n7.72<-numeric(length(time.sim))
dist.n7.72=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n7.72[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-7.72))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-7.72),time.sim[j])
  }
  low.n7.72[j]=quantile(na.omit(dist.n7.72),0.025)
  high.n7.72[j]=quantile(na.omit(dist.n7.72),0.975)
}

#Tau = -1.92 deg. C
low.n1.92<-numeric(length(time.sim))
high.n1.92<-numeric(length(time.sim))
dist.n1.92=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n1.92[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-1.92))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-1.92),time.sim[j])
  }
  low.n1.92[j]=quantile(na.omit(dist.n1.92),0.025)
  high.n1.92[j]=quantile(na.omit(dist.n1.92),0.975)
}

#Tau = 2.48 deg. C
low.p2.48<-numeric(length(time.sim))
high.p2.48<-numeric(length(time.sim))
dist.p2.48=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.48[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.48))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.48),time.sim[j])
  }
  low.p2.48[j]=quantile(na.omit(dist.p2.48),0.025)
  high.p2.48[j]=quantile(na.omit(dist.p2.48),0.975)
}

#Tau = 7.58 deg. C
low.p7.58<-numeric(length(time.sim))
high.p7.58<-numeric(length(time.sim))
dist.p7.58=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p7.58[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],7.58))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],7.58),time.sim[j])
  }
  low.p7.58[j]=quantile(na.omit(dist.p7.58),0.025)
  high.p7.58[j]=quantile(na.omit(dist.p7.58),0.975)
}

##
#Gliricidia
##

#Tau = -6.08 deg. C
low.n6.08<-numeric(length(time.sim))
high.n6.08<-numeric(length(time.sim))
dist.n6.08=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n6.08[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-6.08))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-6.08),time.sim[j])
  }
  low.n6.08[j]=quantile(na.omit(dist.n6.08),0.025)
  high.n6.08[j]=quantile(na.omit(dist.n6.08),0.975)
}

#Tau = -1.08 deg. C
low.n1.08<-numeric(length(time.sim))
high.n1.08<-numeric(length(time.sim))
dist.n1.08=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n1.08[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-1.08))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-1.08),time.sim[j])
  }
  low.n1.08[j]=quantile(na.omit(dist.n1.08),0.025)
  high.n1.08[j]=quantile(na.omit(dist.n1.08),0.975)
}

#Tau = 2.92 deg. C
low.p2.92<-numeric(length(time.sim))
high.p2.92<-numeric(length(time.sim))
dist.p2.92=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.92[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.92))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.92),time.sim[j])
  }
  low.p2.92[j]=quantile(na.omit(dist.p2.92),0.025)
  high.p2.92[j]=quantile(na.omit(dist.p2.92),0.975)
}

#Tau = 7.72 deg. C
low.p7.72<-numeric(length(time.sim))
high.p7.72<-numeric(length(time.sim))
dist.p7.72=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p7.72[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],7.72))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],7.72),time.sim[j])
  }
  low.p7.72[j]=quantile(na.omit(dist.p7.72),0.025)
  high.p7.72[j]=quantile(na.omit(dist.p7.72),0.975)
}

##
#Robinia
##

#Tau = -6.69 deg. C
low.n6.69<-numeric(length(time.sim))
high.n6.69<-numeric(length(time.sim))
dist.n6.69=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n6.69[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-6.69))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-6.69),time.sim[j])
  }
  low.n6.69[j]=quantile(na.omit(dist.n6.69),0.025)
  high.n6.69[j]=quantile(na.omit(dist.n6.69),0.975)
}

#Tau = -1.29 deg. C
low.n1.29<-numeric(length(time.sim))
high.n1.29<-numeric(length(time.sim))
dist.n1.29=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n1.29[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-1.29))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-1.29),time.sim[j])
  }
  low.n1.29[j]=quantile(na.omit(dist.n1.29),0.025)
  high.n1.29[j]=quantile(na.omit(dist.n1.29),0.975)
}

#Tau = 3.11 deg. C
low.p3.11<-numeric(length(time.sim))
high.p3.11<-numeric(length(time.sim))
dist.p3.11=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p3.11[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],3.11))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],3.11),time.sim[j])
  }
  low.p3.11[j]=quantile(na.omit(dist.p3.11),0.025)
  high.p3.11[j]=quantile(na.omit(dist.p3.11),0.975)
}

#Tau = 7.91 deg. C
low.p7.91<-numeric(length(time.sim))
high.p7.91<-numeric(length(time.sim))
dist.p7.91=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p7.91[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],7.91))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],7.91),time.sim[j])
  }
  low.p7.91[j]=quantile(na.omit(dist.p7.91),0.025)
  high.p7.91[j]=quantile(na.omit(dist.p7.91),0.975)
}

###
#26:20 deg. C growing temperature
###

##
#Morella
##

#Tau = -6.84 deg. C
low.n6.84<-numeric(length(time.sim))
high.n6.84<-numeric(length(time.sim))
dist.n6.84=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n6.84[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-6.84))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-6.84),time.sim[j])
  }
  low.n6.84[j]=quantile(na.omit(dist.n6.84),0.025)
  high.n6.84[j]=quantile(na.omit(dist.n6.84),0.975)
}

#Tau = -2.34 deg. C
low.n2.34<-numeric(length(time.sim))
high.n2.34<-numeric(length(time.sim))
dist.n2.34=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n2.34[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-2.34))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-2.34),time.sim[j])
  }
  low.n2.34[j]=quantile(na.omit(dist.n2.34),0.025)
  high.n2.34[j]=quantile(na.omit(dist.n2.34),0.975)
}

#Tau = 2.46 deg. C
low.p2.46<-numeric(length(time.sim))
high.p2.46<-numeric(length(time.sim))
dist.p2.46=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.46[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.46))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.46),time.sim[j])
  }
  low.p2.46[j]=quantile(na.omit(dist.p2.46),0.025)
  high.p2.46[j]=quantile(na.omit(dist.p2.46),0.975)
}

#Tau = 6.06 deg. C
low.p6.06<-numeric(length(time.sim))
high.p6.06<-numeric(length(time.sim))
dist.p6.06=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p6.06[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],6.06))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],6.06),time.sim[j])
  }
  low.p6.06[j]=quantile(na.omit(dist.p6.06),0.025)
  high.p6.06[j]=quantile(na.omit(dist.p6.06),0.975)
}

##
#Alnus
##

#Tau = -6.84 deg. C
low.n6.84<-numeric(length(time.sim))
high.n6.84<-numeric(length(time.sim))
dist.n6.84=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n6.84[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-6.84))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-6.84),time.sim[j])
  }
  low.n6.84[j]=quantile(na.omit(dist.n6.84),0.025)
  high.n6.84[j]=quantile(na.omit(dist.n6.84),0.975)
}

#Tau = -2.54 deg. C
low.n2.54<-numeric(length(time.sim))
high.n2.54<-numeric(length(time.sim))
dist.n2.54=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n2.54[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-2.54))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-2.54),time.sim[j])
  }
  low.n2.54[j]=quantile(na.omit(dist.n2.54),0.025)
  high.n2.54[j]=quantile(na.omit(dist.n2.54),0.975)
}

#Tau = 2.66 deg. C
low.p2.66<-numeric(length(time.sim))
high.p2.66<-numeric(length(time.sim))
dist.p2.66=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.66[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.66))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.66),time.sim[j])
  }
  low.p2.66[j]=quantile(na.omit(dist.p2.66),0.025)
  high.p2.66[j]=quantile(na.omit(dist.p2.66),0.975)
}

#Tau = 6.86 deg. C
low.p6.86<-numeric(length(time.sim))
high.p6.86<-numeric(length(time.sim))
dist.p6.86=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p6.86[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],6.86))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],6.86),time.sim[j])
  }
  low.p6.86[j]=quantile(na.omit(dist.p6.86),0.025)
  high.p6.86[j]=quantile(na.omit(dist.p6.86),0.975)
}

##
#Gliricidia
##

#Tau = -8.26 deg. C
low.n8.26<-numeric(length(time.sim))
high.n8.26<-numeric(length(time.sim))
dist.n8.26=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n8.26[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-8.26))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-8.26),time.sim[j])
  }
  low.n8.26[j]=quantile(na.omit(dist.n8.26),0.025)
  high.n8.26[j]=quantile(na.omit(dist.n8.26),0.975)
}

#Tau = -3.26 deg. C
low.n3.26<-numeric(length(time.sim))
high.n3.26<-numeric(length(time.sim))
dist.n3.26=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n3.26[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-3.26))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-3.26),time.sim[j])
  }
  low.n3.26[j]=quantile(na.omit(dist.n3.26),0.025)
  high.n3.26[j]=quantile(na.omit(dist.n3.26),0.975)
}

#Tau = 1.54 deg. C
low.p1.54<-numeric(length(time.sim))
high.p1.54<-numeric(length(time.sim))
dist.p1.54=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p1.54[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],1.54))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],1.54),time.sim[j])
  }
  low.p1.54[j]=quantile(na.omit(dist.p1.54),0.025)
  high.p1.54[j]=quantile(na.omit(dist.p1.54),0.975)
}

#Tau = 6.14 deg. C
low.p6.14<-numeric(length(time.sim))
high.p6.14<-numeric(length(time.sim))
dist.p6.14=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p6.14[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],6.14))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],6.14),time.sim[j])
  }
  low.p6.14[j]=quantile(na.omit(dist.p6.14),0.025)
  high.p6.14[j]=quantile(na.omit(dist.p6.14),0.975)
}

##
#Robinia
##

#Tau = -6.32 deg. C
low.n6.32<-numeric(length(time.sim))
high.n6.32<-numeric(length(time.sim))
dist.n6.32=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n6.32[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-6.32))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-6.32),time.sim[j])
  }
  low.n6.32[j]=quantile(na.omit(dist.n6.32),0.025)
  high.n6.32[j]=quantile(na.omit(dist.n6.32),0.975)
}

#Tau = -2.62 deg. C
low.n2.62<-numeric(length(time.sim))
high.n2.62<-numeric(length(time.sim))
dist.n2.62=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n2.62[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-2.62))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-2.62),time.sim[j])
  }
  low.n2.62[j]=quantile(na.omit(dist.n2.62),0.025)
  high.n2.62[j]=quantile(na.omit(dist.n2.62),0.975)
}

#Tau = 3.18 deg. C
low.p3.18<-numeric(length(time.sim))
high.p3.18<-numeric(length(time.sim))
dist.p3.18=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p3.18[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],3.18))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],3.18),time.sim[j])
  }
  low.p3.18[j]=quantile(na.omit(dist.p3.18),0.025)
  high.p3.18[j]=quantile(na.omit(dist.p3.18),0.975)
}

#Tau = 7.38 deg. C
low.p7.38<-numeric(length(time.sim))
high.p7.38<-numeric(length(time.sim))
dist.p7.38=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p7.38[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],7.38))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],7.38),time.sim[j])
  }
  low.p7.38[j]=quantile(na.omit(dist.p7.38),0.025)
  high.p7.38[j]=quantile(na.omit(dist.p7.38),0.975)
}

###
#31:25 deg. C growing temperature
###

##
#Morella
##

#Tau = -8.26 deg. C
low.n8.26<-numeric(length(time.sim))
high.n8.26<-numeric(length(time.sim))
dist.n8.26=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n8.26[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-8.26))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-8.26),time.sim[j])
  }
  low.n8.26[j]=quantile(na.omit(dist.n8.26),0.025)
  high.n8.26[j]=quantile(na.omit(dist.n8.26),0.975)
}

#Tau = -2.26 deg. C
low.n2.26<-numeric(length(time.sim))
high.n2.26<-numeric(length(time.sim))
dist.n2.26=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n2.26[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-2.26))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-2.26),time.sim[j])
  }
  low.n2.26[j]=quantile(na.omit(dist.n2.26),0.025)
  high.n2.26[j]=quantile(na.omit(dist.n2.26),0.975)
}

#Tau = 2.34 deg. C
low.p2.34<-numeric(length(time.sim))
high.p2.34<-numeric(length(time.sim))
dist.p2.34=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.34[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.34))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.34),time.sim[j])
  }
  low.p2.34[j]=quantile(na.omit(dist.p2.34),0.025)
  high.p2.34[j]=quantile(na.omit(dist.p2.34),0.975)
}

##
#Alnus
##

#Tau = -4.17 deg. C
low.n4.17<-numeric(length(time.sim))
high.n4.17<-numeric(length(time.sim))
dist.n4.17=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n4.17[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-4.17))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-4.17),time.sim[j])
  }
  low.n4.17[j]=quantile(na.omit(dist.n4.17),0.025)
  high.n4.17[j]=quantile(na.omit(dist.n4.17),0.975)
}

#Tau = 2.13 deg. C
low.p2.13<-numeric(length(time.sim))
high.p2.13<-numeric(length(time.sim))
dist.p2.13=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.13[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.13))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.13),time.sim[j])
  }
  low.p2.13[j]=quantile(na.omit(dist.p2.13),0.025)
  high.p2.13[j]=quantile(na.omit(dist.p2.13),0.975)
}

#Tau = 6.83 deg. C
low.p6.83<-numeric(length(time.sim))
high.p6.83<-numeric(length(time.sim))
dist.p6.83=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p6.83[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],6.83))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],6.83),time.sim[j])
  }
  low.p6.83[j]=quantile(na.omit(dist.p6.83),0.025)
  high.p6.83[j]=quantile(na.omit(dist.p6.83),0.975)
}

##
#Gliricidia
##

#Tau = -6.74 deg. C
low.n6.74<-numeric(length(time.sim))
high.n6.74<-numeric(length(time.sim))
dist.n6.74=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n6.74[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-6.74))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-6.74),time.sim[j])
  }
  low.n6.74[j]=quantile(na.omit(dist.n6.74),0.025)
  high.n6.74[j]=quantile(na.omit(dist.n6.74),0.975)
}

#Tau = -0.14 deg. C
low.n0.14<-numeric(length(time.sim))
high.n0.14<-numeric(length(time.sim))
dist.n0.14=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n0.14[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-0.14))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-0.14),time.sim[j])
  }
  low.n0.14[j]=quantile(na.omit(dist.n0.14),0.025)
  high.n0.14[j]=quantile(na.omit(dist.n0.14),0.975)
}

#Tau = 3.56 deg. C
low.p3.56<-numeric(length(time.sim))
high.p3.56<-numeric(length(time.sim))
dist.p3.56=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p3.56[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],3.56))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],3.56),time.sim[j])
  }
  low.p3.56[j]=quantile(na.omit(dist.p3.56),0.025)
  high.p3.56[j]=quantile(na.omit(dist.p3.56),0.975)
}

##
#Robinia
##

#Tau = -7.94 deg. C
low.n7.94<-numeric(length(time.sim))
high.n7.94<-numeric(length(time.sim))
dist.n7.94=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n7.94[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-7.94))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-7.94),time.sim[j])
  }
  low.n7.94[j]=quantile(na.omit(dist.n7.94),0.025)
  high.n7.94[j]=quantile(na.omit(dist.n7.94),0.975)
}

#Tau = -1.14 deg. C
low.n1.14<-numeric(length(time.sim))
high.n1.14<-numeric(length(time.sim))
dist.n1.14=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.n1.14[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],-1.14))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],-1.14),time.sim[j])
  }
  low.n1.14[j]=quantile(na.omit(dist.n1.14),0.025)
  high.n1.14[j]=quantile(na.omit(dist.n1.14),0.975)
}

#Tau = 2.56 deg. C
low.p2.56<-numeric(length(time.sim))
high.p2.56<-numeric(length(time.sim))
dist.p2.56=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p2.56[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],2.56))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],2.56),time.sim[j])
  }
  low.p2.56[j]=quantile(na.omit(dist.p2.56),0.025)
  high.p2.56[j]=quantile(na.omit(dist.p2.56),0.975)
}

#Tau = 7.36 deg. C
low.p7.36<-numeric(length(time.sim))
high.p7.36<-numeric(length(time.sim))
dist.p7.36=numeric(10000)
for(j in 1:601){
  for(i in 1:10000){
    dist.p7.36[i]<-fall.norm(exp(exp.mod(vmat.b[i,3],vmat.b[i,2],7.36))-1,c.est,
                             nsig.mod(vmat.d[i,2],vmat.d[i,3],7.36),time.sim[j])
  }
  low.p7.36[j]=quantile(na.omit(dist.p7.36),0.025)
  high.p7.36[j]=quantile(na.omit(dist.p7.36),0.975)
}

###############################################################################################################
#Supplementary Figure 12
###############################################################################################################

#PDF dimension is 9x5 inches  

#Plotting Settings
nf<-layout(matrix(c(1,2,3),1,3,byrow=T),c(3,3,1.3),c(3,3,3),T)
layout.show(nf)
par(oma=c(2,1,1,1))
par(mar=c(4,6,4,1))
par(pty="s")

#S12a
plot(b~Tau,data=temp.parms.rm,col="white",xlab=NA,ylab=NA,xlim=c(-10,10),las=1,cex.axis=1.2)
mtext(text="a",side=3,cex=1.2,adj=0)
mtext(expression(italic(b)*'('*tau*')'),side=2,cex=1.2,line=3)
points(b~Tau,data=M21,col="darkorange1",pch=16,cex=1.5)
points(b~Tau,data=M21,col="dodgerblue1",pch=1,cex=1.5)
points(b~Tau,data=A21,col="darkturquoise",pch=16,cex=1.5)
points(b~Tau,data=A21,col="dodgerblue1",pch=1,cex=1.5)
points(b~Tau,data=G21,col="orangered2",pch=17,cex=1.5)
points(b~Tau,data=G21,col="dodgerblue1",pch=2,cex=1.5)
points(b~Tau,data=R21,col="dodgerblue3",pch=17,cex=1.5)
points(b~Tau,data=R21,col="dodgerblue1",pch=2,cex=1.5)
points(b~Tau,data=M26,col="darkorange1",pch=16,cex=1.5)
points(b~Tau,data=M26,col="gold1",pch=1,cex=1.5)
points(b~Tau,data=A26,col="darkturquoise",pch=16,cex=1.5)
points(b~Tau,data=A26,col="gold1",pch=1,cex=1.5)
points(b~Tau,data=G26,col="orangered2",pch=17,cex=1.5)
points(b~Tau,data=G26,col="gold1",pch=2,cex=1.5)
points(b~Tau,data=R26,col="dodgerblue3",pch=17,cex=1.5)
points(b~Tau,data=R26,col="gold1",pch=2,cex=1.5)
points(b~Tau,data=M31,col="darkorange1",pch=16,cex=1.5)
points(b~Tau,data=M31,col="orangered3",pch=1,cex=1.5)
points(b~Tau,data=A31,col="darkturquoise",pch=16,cex=1.5)
points(b~Tau,data=A31,col="orangered3",pch=1,cex=1.5)
points(b~Tau,data=G31,col="orangered2",pch=17,cex=1.5)
points(b~Tau,data=G31,col="orangered3",pch=2,cex=1.5)
points(b~Tau,data=R31,col="dodgerblue3",pch=17,cex=1.5)
points(b~Tau,data=R31,col="orangered3",pch=2,cex=1.5)
curve(exp(exp.mod(log(alpha.mle2),log(beta.mle2),x))-1,from=-10,to=10,lty=1,lwd=2,col="black",add=T)
points(exp(low.b)-1~xv,lty=2,type="l")
points(exp(high.b)-1~xv,lty=2,type="l")

#S12b
plot(d~Tau,data=temp.parms.rm,xlab=NA,ylab=NA,cex.axis=1.2,las=1,xlim=c(-10,10),col="white")
mtext(text="b",side=3,cex=1.2,adj=0)
mtext(expression(italic(d)*'('*tau*')'),side=2,cex=1.2,line=3)
points(d~Tau,data=M21,col="darkorange1",pch=16,cex=1.5)
points(d~Tau,data=M21,col="dodgerblue1",pch=1,cex=1.5)
points(d~Tau,data=A21,col="darkturquoise",pch=16,cex=1.5)
points(d~Tau,data=A21,col="dodgerblue1",pch=1,cex=1.5)
points(d~Tau,data=G21,col="orangered2",pch=17,cex=1.5)
points(d~Tau,data=G21,col="dodgerblue1",pch=2,cex=1.5)
points(d~Tau,data=R21,col="dodgerblue3",pch=17,cex=1.5)
points(d~Tau,data=R21,col="dodgerblue1",pch=2,cex=1.5)
points(d~Tau,data=M26,col="darkorange1",pch=16,cex=1.5)
points(d~Tau,data=M26,col="gold1",pch=1,cex=1.5)
points(d~Tau,data=A26,col="darkturquoise",pch=16,cex=1.5)
points(d~Tau,data=A26,col="gold1",pch=1,cex=1.5)
points(d~Tau,data=G26,col="orangered2",pch=17,cex=1.5)
points(d~Tau,data=G26,col="gold1",pch=2,cex=1.5)
points(d~Tau,data=R26,col="dodgerblue3",pch=17,cex=1.5)
points(d~Tau,data=R26,col="gold1",pch=2,cex=1.5)
points(d~Tau,data=M31,col="darkorange1",pch=16,cex=1.5)
points(d~Tau,data=M31,col="orangered3",pch=1,cex=1.5)
points(d~Tau,data=A31,col="darkturquoise",pch=16,cex=1.5)
points(d~Tau,data=A31,col="orangered3",pch=1,cex=1.5)
points(d~Tau,data=G31,col="orangered2",pch=17,cex=1.5)
points(d~Tau,data=G31,col="orangered3",pch=2,cex=1.5)
points(d~Tau,data=R31,col="dodgerblue3",pch=17,cex=1.5)
points(d~Tau,data=R31,col="orangered3",pch=2,cex=1.5)
curve(nsig.mod(log(p.mle),log(q.mle),x),from=-10,to=10,lty=1,lwd=2,col="black",add=T)
points(low.d~xv,lty=2,type="l")
points(high.d~xv,lty=2,type="l")

mtext(expression(tau*' = '*italic(T)[opt]*' '-' '*italic(T)[meas]*' ('*degree*'C)            '),side=1,line=-4,cex=1.2,outer=T)

#Legend
par(mar=c(0,0,0,0))
plot(0:10, 0:10, type='n', bty='n', xaxt='n', yaxt='n',xlab=NA,ylab=NA)
legend("left",c(expression(underline(bolditalic("Morella"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Alnus"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Gliricidia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Robinia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C')),
       col=c(NA,"darkorange1","darkorange1","darkorange1",NA,"darkturquoise","darkturquoise","darkturquoise",
             NA,"orangered2","orangered2","orangered2",NA,"dodgerblue3","dodgerblue3","dodgerblue3"),
       pch=c(NA,16,16,16,NA,16,16,16,NA,17,17,17,NA,17,17,17),bty="n",y.intersp = 0.8,cex=1,x.intersp = 0.7,pt.cex=1.5,xpd=T,
       text.col=c("white","white","white","white","white","white","white","white","white","white","white","white"))
legend("left",c(expression(underline(bolditalic("Morella"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Alnus"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Gliricidia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C'),
                expression(underline(bolditalic("Robinia"))),expression('21:15 '*degree*'C'),expression('26:20 '*degree*'C'),expression('31:25 '*degree*'C')),
       col=c(NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",NA,
             "dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3"),
       pch=c(NA,1,1,1,NA,1,1,1,NA,2,2,2,NA,2,2,2),bty="n",y.intersp = 0.8,cex=1,x.intersp = 0.7,pt.cex=1.5,xpd=T)

###############################################################################################################
#Supplementary Figure 13
###############################################################################################################

#Note: Data is multiplied by (1+b)/(1+b*d) to normalize to 1

#PDF dimension is 9x7 inches

#Plotting Settings
par(pty="s")
nf<-layout(matrix(seq(1,12,1),3,4,byrow=T),rep(3,12),rep(3,12),T)
layout.show(nf)
par(oma=c(5,5,4,4))
par(mar=c(0,0,0,0))

#S13a
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=c("0","20","40","60","80","100","120","140"),las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n3.43,rev(high.n3.43)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p0.97,rev(high.p0.97)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p5.87,rev(high.p5.87)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p10.07,rev(high.p10.07)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.92*max(M21_25C.cut$Vmax))~Time,data=M21_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.01*Vmax/(0.94*max(M21_30C.cut$Vmax))~Time,data=M21_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.04*Vmax/(1.28*max(M21_35C.cut$Vmax))~Time,data=M21_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.86*Vmax/(2.08*max(M21_40C.cut$Vmax))~Time,data=M21_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-3.43))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-3.43),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),0.97))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),0.97),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),5.87))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),5.87),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),10.07))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),10.07),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Mb21a),Mc21a,Md21a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Mb21b),Mc21b,Md21b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Mb21c),Mc21c,Md21c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Mb21d),Mc21d,Md21d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
mtext(expression(italic(Morella)),side=3,line=1,cex=1.5)
title(main=expression('  a'),cex.main=1.5,adj=0,line=-1)

#S13b
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n7.72,rev(high.n7.72)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n1.92,rev(high.n1.92)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.48,rev(high.p2.48)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p7.58,rev(high.p7.58)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.92*max(A21_25C.cut$Vmax))~Time,data=A21_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.12*Vmax/(1.15*max(A21_30C.cut$Vmax))~Time,data=A21_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.05*Vmax/(0.96*max(A21_35C.cut$Vmax))~Time,data=A21_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.28*Vmax/(1.36*max(A21_40C.cut$Vmax))~Time,data=A21_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-7.72))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-7.72),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-1.92))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-1.92),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.48))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.48),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),7.58))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),7.58),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Ab21a),Ac21a,Ad21a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Ab21b),Ac21b,Ad21b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Ab21c),Ac21c,Ad21c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Ab21d),Ac21d,Ad21d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
mtext(expression(italic(Alnus)),side=3,line=1,cex=1.5)
title(main=expression('  b'),cex.main=1.5,adj=0,line=-1)

#S13c
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n6.08,rev(high.n6.08)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n1.08,rev(high.n1.08)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.92,rev(high.p2.92)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p7.72,rev(high.p7.72)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.35*0.97*max(G21_25C.cut$Vmax))~Time,data=G21_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.04*Vmax/(0.94*max(G21_30C.cut$Vmax))~Time,data=G21_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(Vmax/(0.7*0.90*max(G21_35C.cut$Vmax))~Time,data=G21_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.19*Vmax/(1.15*max(G21_40C.cut$Vmax))~Time,data=G21_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-6.08))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-6.08),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-1.08))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-1.08),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.92))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.92),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),7.72))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),7.72),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Gb21a),Gc21a,Gd21a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Gb21b),Gc21b,Gd21b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Gb21c),Gc21c,Gd21c,x),from=0,to=6,lty=3,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Gb21d),Gc21d,Gd21d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
mtext(expression(italic(Gliricidia)),side=3,line=1,cex=1.5)
title(main=expression('  c'),cex.main=1.5,adj=0,line=-1)

#S13d
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n6.69,rev(high.n6.69)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n1.29,rev(high.n1.29)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p3.11,rev(high.p3.11)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p7.91,rev(high.p7.91)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.91*max(R21_25C.cut$Vmax))~Time,data=R21_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.03*Vmax/(1.03*max(R21_30C.cut$Vmax))~Time,data=R21_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.03*Vmax/(1.01*max(R21_35C.cut$Vmax))~Time,data=R21_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.04*Vmax/(1.25*max(R21_40C.cut$Vmax))~Time,data=R21_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-6.69))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-6.69),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-1.29))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-1.29),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),3.11))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),3.11),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),7.91))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),7.91),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Rb21a),Rc21a,Rd21a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Rb21b),Rc21b,Rd21b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Rb21c),Rc21c,Rd21c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Rb21d),Rc21d,Rd21d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
mtext(expression(italic(Robinia)),side=3,line=1,cex=1.5)
mtext(expression('21:15 '*degree*'C'),side=4,line=1,cex=1.5)
title(main=expression('  d'),cex.main=1.5,adj=0,line=-1)

#S13e
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=c("0","20","40","60","80","100","120","140"),las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n6.84,rev(high.n6.84)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n2.34,rev(high.n2.34)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.46,rev(high.p2.46)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p6.06,rev(high.p6.06)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(2.64*Vmax/(2.67*max(M26_25C.cut$Vmax))~Time,data=M26_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.01*Vmax/(0.95*max(M26_30C.cut$Vmax))~Time,data=M26_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(Vmax/(0.89*max(M26_35C.cut$Vmax))~Time,data=M26_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.55*Vmax/(1.67*max(M26_40C.cut$Vmax))~Time,data=M26_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-6.84))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-6.84),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-2.34))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-2.34),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.46))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.46),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),6.06))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),6.06),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Mb26a),Mc26a,Md26a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Mb26b),Mc26b,Md26b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Mb26c),Mc26c,Md26c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Mb26d),Mc26d,Md26d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
title(main=expression('  e'),cex.main=1.5,adj=0,line=-1)

#S13f
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n6.84,rev(high.n6.84)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n2.54,rev(high.n2.54)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.66,rev(high.p2.66)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p6.86,rev(high.p6.86)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.85*max(A26_25C.cut$Vmax))~Time,data=A26_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.05*Vmax/(1.02*max(A26_30C.cut$Vmax))~Time,data=A26_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(Vmax/(0.25*1.07*max(A26_35C.cut$Vmax))~Time,data=A26_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.18*Vmax/(1.25*max(A26_40C.cut$Vmax))~Time,data=A26_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-6.84))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-6.84),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-2.54))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-2.54),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.66))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.66),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),6.86))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),6.86),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Ab26a),Ac26a,Ad26a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Ab26b),Ac26b,Ad26b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Ab26c),Ac26c,Ad26c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Ab26d),Ac26d,Ad26d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
title(main=expression('  f'),cex.main=1.5,adj=0,line=-1)

#S13g
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n8.26,rev(high.n8.26)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n3.26,rev(high.n3.26)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p1.54,rev(high.p1.54)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p6.14,rev(high.p6.14)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.74*max(G26_25C.cut$Vmax))~Time,data=G26_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.03*Vmax/(0.95*max(G26_30C.cut$Vmax))~Time,data=G26_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.06*Vmax/(1.02*max(G26_35C.cut$Vmax))~Time,data=G26_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.77*Vmax/(1.87*max(G26_40C.cut$Vmax))~Time,data=G26_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-8.26))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-8.26),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-3.26))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-3.26),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),1.54))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),1.54),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),6.14))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),6.14),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Gb26a),Gc26a,Gd26a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Gb26b),Gc26b,Gd26b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Gb26c),Gc26c,Gd26c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Gb26d),Gc26d,Gd26d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
title(main=expression('  g'),cex.main=1.5,adj=0,line=-1)

#S13h
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=F,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n6.32,rev(high.n6.32)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n2.62,rev(high.n2.62)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p3.18,rev(high.p3.18)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p7.38,rev(high.p7.38)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.89*max(R26_25C.cut$Vmax))~Time,data=R26_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(1.02*Vmax/(0.97*max(R26_30C.cut$Vmax))~Time,data=R26_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(6.11*Vmax/(6.19*max(R26_35C.cut$Vmax))~Time,data=R26_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(8.73*Vmax/(8.67*max(R26_40C.cut$Vmax))~Time,data=R26_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-6.32))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-6.32),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-2.62))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-2.62),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),3.18))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),3.18),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),7.38))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),7.38),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Rb26a),Rc26a,Rd26a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Rb26b),Rc26b,Rd26b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Rb26c),Rc26c,Rd26c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Rb26d),Rc26d,Rd26d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
mtext(expression('26:20 '*degree*'C'),side=4,line=1,cex=1.5)
title(main=expression('  h'),cex.main=1.5,adj=0,line=-1)

#S13i
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=T,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=c("0","20","40","60","80","100","120","140"),las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n8.26,rev(high.n8.26)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n2.26,rev(high.n2.26)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.34,rev(high.p2.34)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.9*0.95*max(M31_30C.cut$Vmax))~Time,data=M31_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.01*Vmax/(0.94*max(M31_35C.cut$Vmax))~Time,data=M31_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(Vmax/(0.91*max(M31_40C.cut$Vmax))~Time,data=M31_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-8.26))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-8.26),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-2.26))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-2.26),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.34))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.34),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Mb31a),Mc31a,Md31a,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Mb31b),Mc31b,Md31b,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Mb31c),Mc31c,Md31c,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
title(main=expression('  i'),cex.main=1.5,adj=0,line=-1)

#S13j
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=T,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n4.17,rev(high.n4.17)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.13,rev(high.p2.13)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p6.83,rev(high.p6.83)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.8*0.97*max(A31_30C.cut$Vmax))~Time,data=A31_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.02*Vmax/(0.99*max(A31_35C.cut$Vmax))~Time,data=A31_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.12*Vmax/(1.34*max(A31_40C.cut$Vmax))~Time,data=A31_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-4.17))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-4.17),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.13))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.13),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),6.83))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),6.83),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Ab31a),Ac31a,Ad31a,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Ab31b),Ac31b,Ad31b,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Ab31c),Ac31c,Ad31c,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
title(main=expression('  j'),cex.main=1.5,adj=0,line=-1)

#S13k
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=T,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n6.74,rev(high.n6.74)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n0.14,rev(high.n0.14)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p3.56,rev(high.p3.56)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.8*0.98*max(G31_30C.cut$Vmax))~Time,data=G31_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(Vmax/(0.89*max(G31_35C.cut$Vmax))~Time,data=G31_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.02*Vmax/(0.91*max(G31_40C.cut$Vmax))~Time,data=G31_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-6.74))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-6.74),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-0.14))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-0.14),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),3.56))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),3.56),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Gb31a),Gc31a,Gd31a,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Gb31b),Gc31b,Gd31b,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(Gb31c),Gc31c,Gd31c,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
title(main=expression('  k'),cex.main=1.5,adj=0,line=-1)

#S13l
plot(0:10,0:10,col="white",xlim=c(0,6),ylim=c(0,1.4),xlab=NA,ylab=NA,xaxt="n",yaxt="n")
axis(1,at=c(0,1,2,3,4,5,6),labels=T,cex.axis=1.5)
axis(2,at=c(0,.2,.4,.6,.8,1,1.2,1.4),labels=F,las=1,cex.axis=1.5)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n7.94,rev(high.n7.94)),
        col=adjustcolor("dodgerblue3", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.n1.14,rev(high.n1.14)),
        col=adjustcolor("orange", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p2.56,rev(high.p2.56)),
        col=adjustcolor("red", alpha.f = 0.10), border = NA)
polygon(x=c(time.sim,rev(time.sim)),y=c(low.p7.36,rev(high.p7.36)),
        col=adjustcolor("darkred", alpha.f = 0.10), border = NA)
points(Vmax/(0.93*max(R31_25C.cut$Vmax))~Time,data=R31_25C,pch=20,
       col=adjustcolor("dodgerblue3",alpha.f = 0.05),cex=0.01)
points(Vmax/(0.91*max(R31_30C.cut$Vmax))~Time,data=R31_30C,pch=20,
       col=adjustcolor("orange",alpha.f = 0.05),cex=0.01)
points(1.08*Vmax/(1.07*max(R31_35C.cut$Vmax))~Time,data=R31_35C,pch=20,
       col=adjustcolor("red",alpha.f = 0.05),cex=0.01)
points(1.11*Vmax/(1.08*max(R31_40C.cut$Vmax))~Time,data=R31_40C,pch=20,
       col=adjustcolor("darkred",alpha.f = 0.05),cex=0.01)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-7.94))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-7.94),x),
      from=0,to=6,lty=1,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),-1.14))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),-1.14),x),
      from=0,to=6,lty=1,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),2.56))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),2.56),x),
      from=0,to=6,lty=1,lwd=1.5,col="red",add=T)
curve(fall.norm(exp(exp.mod(log(alpha.mle2),log(beta.mle2),7.36))-1,c.est,
                nsig.mod(log(p.mle),log(q.mle),7.36),x),
      from=0,to=6,lty=1,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Rb31a),Rc31a,Rd31a,x),from=0,to=6,lty=2,lwd=1.5,col="dodgerblue3",add=T)
curve(fall.norm(exp(Rb31b),Rc31b,Rd31b,x),from=0,to=6,lty=2,lwd=1.5,col="orange",add=T)
curve(fall.norm(exp(Rb31d),Rc31d,Rd31d,x),from=0,to=6,lty=2,lwd=1.5,col="darkred",add=T)
curve(fall.norm(exp(Rb31c),Rc31c,Rd31c,x),from=0,to=6,lty=2,lwd=1.5,col="red",add=T)
mtext(expression('31:25 '*degree*'C'),side=4,line=1,cex=1.5)
title(main=expression('  l'),cex.main=1.5,adj=0,line=-1)

mtext(expression('Nitrogen fixation (% of max)'),side=2,line=2.9,cex=1.5,outer=T)
mtext(expression('Time (hr)'),side=1,line=3,cex=1.5,outer=T)