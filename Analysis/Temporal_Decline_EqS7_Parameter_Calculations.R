###############################################################################################################
###############################################################################################################
#This script fits estimates the parameter values in Equation S7
###############################################################################################################
###############################################################################################################

###This script focuses on the falling portion of Equation S4 that occurs above the breakpoint calculated in
#Piecewise_Reg_Constant_Temp_Fits.R

###First, the modified negative sigmoid function is Equation S4 is fit to estimate parameters "b", "c", and "d"
#for each continuous N-fixation measurement at constant temperatures. Next, the median value of "c" is
#calculated and "b" and "d" are regressed across tau (Equations S5 and S6) to parameterize Equation S7

####
#Load Necessary Package
####

library(bbmle)

####
#Define functions
####

#How N-fixation changes as a function of time for each individual (falling portion of Equation S4)
fall<-function(a,b,c,d,t){
  y<-a*((1-d)/(1+b*exp(c*t))+d)
  y
}

#How b changes as a function of tau 
exp.mod <- function(alpha,beta,tau){
  y <- exp(alpha)*exp(exp(beta)*tau)
  y
}

#How d changes as a function of tau
nsig.mod <- function(p,q,tau){
  y <-1/(1+exp(p)*exp(exp(q)*tau))
  y
}

###############################################################################################################
#Subset data from breakpoint to 6 hr mark
#Vmax data are approximately normalized to 1 in order to make it easier to fit falling function
###############################################################################################################

###
#Morella
###

##
#21:15 deg. C growing temperature
##

#25 deg. C
M21_25C<-read.csv("MOCE21_long25.csv")
#Breakpoint was visual estimated due to unequilibrated rates with fit breakpoint
print(st<-which(M21_25C$Time==M21_25C$Time[M21_25C$Time>1.4][1]))
print(end<-which(M21_25C$Time==M21_25C$Time[M21_25C$Time>6][1]))
M21_25C.cut<-M21_25C[st:end,]
M21_25C.cut$Vmax.norm<-M21_25C.cut$Vmax/max(M21_25C.cut$Vmax)

#30 deg. C
M21_30C<-read.csv("MOCE21_long30.csv")
#Breakpoint was visual estimated due to unequilibrated rates with fit breakpoint
print(st<-which(M21_30C$Time==M21_30C$Time[M21_30C$Time>1.4][1]))
print(end<-which(M21_30C$Time==M21_30C$Time[M21_30C$Time>6][1]))
M21_30C.cut<-M21_30C[st:end,]
M21_30C.cut$Vmax.norm<-M21_30C.cut$Vmax/max(M21_30C.cut$Vmax)

#35 deg. C
M21_35C<-read.csv("MOCE21_long35.csv")
print(st<-which(M21_35C$Time==M21_35C$Time[M21_35C$Time>2.44][1]))
print(end<-which(M21_35C$Time==M21_35C$Time[M21_35C$Time>6][1]))
M21_35C.cut<-M21_35C[st:end,]
M21_35C.cut$Vmax.norm<-M21_35C.cut$Vmax/max(M21_35C.cut$Vmax)

#40 deg. C
M21_40C<-read.csv("MOCE21_long40.csv")
print(st<-which(M21_40C$Time==M21_40C$Time[M21_40C$Time>0.15][1]))
print(end<-which(M21_40C$Time==M21_40C$Time[M21_40C$Time>6][1]))
M21_40C.cut<-M21_40C[st:end,]
M21_40C.cut$Vmax.norm<-M21_40C.cut$Vmax/max(M21_40C.cut$Vmax)

##
#26:20 deg. C growing temperature
##

#25 deg. C
M26_25C<-read.csv("MOCE26_long25.csv")
#Breakpoint was visual estimated due to unequilibrated rates with fit breakpoint
print(st<-which(M26_25C$Time==M26_25C$Time[M26_25C$Time>0.6][1]))
print(end<-length(M26_25C$Time))
M26_25C.cut<-M26_25C[st:end,]
M26_25C.cut$Vmax.norm<-M26_25C.cut$Vmax/max(M26_25C.cut$Vmax)

#30 deg. C
M26_30C<-read.csv("MOCE26_long30.csv")
print(st<-which(M26_30C$Time==M26_30C$Time[M26_30C$Time>0.78][1]))
print(end<-length(M26_30C$Time))
M26_30C.cut<-M26_30C[st:end,]
M26_30C.cut$Vmax.norm<-M26_30C.cut$Vmax/max(M26_30C.cut$Vmax)

#35 deg. C
M26_35C<-read.csv("MOCE26_long35.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(M26_35C$Time==M26_35C$Time[M26_35C$Time>1][1]))
print(end<-length(M26_35C$Time))
M26_35C.cut<-M26_35C[st:end,]
M26_35C.cut$Vmax.norm<-M26_35C.cut$Vmax/max(M26_35C.cut$Vmax)

#40 deg. C
M26_40C<-read.csv("MOCE26_long40.csv")
print(st<-which(M26_40C$Time==M26_40C$Time[M26_40C$Time>0.32][1]))
print(end<-length(M26_40C$Time))
M26_40C.cut<-M26_40C[st:end,]
M26_40C.cut$Vmax.norm<-M26_40C.cut$Vmax/max(M26_40C.cut$Vmax)

##
#31:25 deg. C growing temperature
##

#30 deg. C
M31_30C<-read.csv("MOCE31_long30.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(M31_30C$Time==M31_30C$Time[M31_30C$Time>1][1]))
print(end<-length(M31_30C$Time))
M31_30C.cut<-M31_30C[st:end,]
M31_30C.cut$Vmax.norm<-M31_30C.cut$Vmax/max(M31_30C.cut$Vmax)
M31_30C.cut$Vmax.norm<-M31_30C.cut$Vmax.norm/0.9

#35 deg. C
M31_35C<-read.csv("MOCE31_long35.csv")
print(st<-which(M31_35C$Time==M31_35C$Time[M31_35C$Time>3.9][1]))
print(end<-which(M31_35C$Time==M31_35C$Time[M31_35C$Time>6][1]))
M31_35C.cut<-M31_35C[st:end,]
M31_35C.cut$Vmax.norm<-M31_35C.cut$Vmax/max(M31_35C.cut$Vmax)

#40 deg. C
M31_40C<-read.csv("MOCE31_long40.csv")
#Breakpoint was visual estimated due to unequilibrated rates with fit breakpoint
print(st<-which(M31_40C$Time==M31_40C$Time[M31_40C$Time>3][1]))
print(end<-length(M31_40C$Time))
M31_40C.cut<-M31_40C[st:end,]
M31_40C.cut$Vmax.norm<-M31_40C.cut$Vmax/max(M31_40C.cut$Vmax)

###
#Alnus
###

##
#21:15 deg. C growing temperature
##

#25 deg. C
A21_25C<-read.csv("ALRU21_long25.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(A21_25C$Time==A21_25C$Time[A21_25C$Time>4][1]))
print(end<-which(A21_25C$Time==A21_25C$Time[A21_25C$Time>6][1]))
A21_25C.cut<-A21_25C[st:end,]
A21_25C.cut$Vmax.norm<-A21_25C.cut$Vmax/max(A21_25C.cut$Vmax)

#30 deg. C
A21_30C<-read.csv("ALRU21_long30.csv")
print(st<-which(A21_30C$Time==A21_30C$Time[A21_30C$Time>2.98][1]))
print(end<-which(A21_30C$Time==A21_30C$Time[A21_30C$Time>6][1]))
A21_30C.cut<-A21_30C[st:end,]
A21_30C.cut$Vmax.norm<-A21_30C.cut$Vmax/max(A21_30C.cut$Vmax)

#35 deg. C
A21_35C<-read.csv("ALRU21_long35.csv")
print(st<-which(A21_35C$Time==A21_35C$Time[A21_35C$Time>0.51][1]))
print(end<-which(A21_35C$Time==A21_35C$Time[A21_35C$Time>6][1]))
A21_35C.cut<-A21_35C[st:end,]
A21_35C.cut$Vmax.norm<-A21_35C.cut$Vmax/max(A21_35C.cut$Vmax)

#40 deg. C
A21_40C<-read.csv("ALRU21_long40.csv")
print(st<-which(A21_40C$Time==A21_40C$Time[A21_40C$Time>0.21][1]))
print(end<-which(A21_40C$Time==A21_40C$Time[A21_40C$Time>6][1]))
A21_40C.cut<-A21_40C[st:end,]
A21_40C.cut$Vmax.norm<-A21_40C.cut$Vmax/max(A21_40C.cut$Vmax)

##
#26:20 deg. C growing temperature
##

#25 deg. C
A26_25C<-read.csv("ALRU26_long25.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(A26_25C$Time==A26_25C$Time[A26_25C$Time>1][1]))
print(end<-length(A26_25C$Time))
A26_25C.cut<-A26_25C[st:end,]
A26_25C.cut$Vmax.norm<-A26_25C.cut$Vmax/max(A26_25C.cut$Vmax)

#30 deg. C
A26_30C<-read.csv("ALRU26_long30.csv")
print(st<-which(A26_30C$Time==A26_30C$Time[A26_30C$Time>2.88][1]))
print(end<-length(A26_30C$Time))
A26_30C.cut<-A26_30C[st:end,]
A26_30C.cut$Vmax.norm<-A26_30C.cut$Vmax/max(A26_30C.cut$Vmax)

#35 deg. C
A26_35C<-read.csv("ALRU26_long35.csv")
print(st<-which(A26_35C$Time==A26_35C$Time[A26_35C$Time>1.83][1]))
print(end<-length(A26_35C$Time))
A26_35C.cut<-A26_35C[st:end,]
A26_35C.cut$Vmax.norm<-4*A26_35C.cut$Vmax/max(A26_35C.cut$Vmax)

#40 deg. C
A26_40C<-read.csv("ALRU26_long40.csv")
print(st<-which(A26_40C$Time==A26_40C$Time[A26_40C$Time>0.07][1]))
print(end<-length(A26_40C$Time))
A26_40C.cut<-A26_40C[st:end,]
A26_40C.cut$Vmax.norm<-A26_40C.cut$Vmax/max(A26_40C.cut$Vmax)

##
#31:25 deg. C growing temperature
##

#30 deg. C
A31_30C<-read.csv("ALRU31_long30.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(A31_30C$Time==A31_30C$Time[A31_30C$Time>2][1]))
print(end<-which(A31_30C$Time==A31_30C$Time[A31_30C$Time>6][1]))
A31_30C.cut<-A31_30C[st:end,]
A31_30C.cut$Vmax.norm<-A31_30C.cut$Vmax/max(A31_30C.cut$Vmax)
A31_30C.cut$Vmax.norm<-A31_30C.cut$Vmax.norm/0.8

#35 deg. C
A31_35C<-read.csv("ALRU31_long35.csv")
print(st<-which(A31_35C$Time==A31_35C$Time[A31_35C$Time>0.12][1]))
print(end<-which(A31_35C$Time==A31_35C$Time[A31_35C$Time>6][1]))
A31_35C.cut<-A31_35C[st:end,]
A31_35C.cut$Vmax.norm<-A31_35C.cut$Vmax/max(A31_35C.cut$Vmax)

#40 deg. C
A31_40C<-read.csv("ALRU31_long40.csv")
print(st<-which(A31_40C$Time==A31_40C$Time[A31_40C$Time>0.35][1]))
print(end<-length(A31_40C$Time))
A31_40C.cut<-A31_40C[st:end,]
A31_40C.cut$Vmax.norm<-A31_40C.cut$Vmax/max(A31_40C.cut$Vmax)

###
#Gliricidia
###

##
#21:15 deg. C growing temperature
##

#25 deg. C
G21_25C<-read.csv("GLSE21_long25.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(G21_25C$Time==G21_25C$Time[G21_25C$Time>2][1]))
print(end<-which(G21_25C$Time==G21_25C$Time[G21_25C$Time>6][1]))
G21_25C.cut<-G21_25C[st:end,]
G21_25C.cut$Vmax.norm<-G21_25C.cut$Vmax/max(G21_25C.cut$Vmax)
G21_25C.cut$Vmax.norm<-G21_25C.cut$Vmax.norm/0.35

#30 deg. C
G21_30C<-read.csv("GLSE21_long30.csv")
print(st<-which(G21_30C$Time==G21_30C$Time[G21_30C$Time>1.5][1]))
print(end<-which(G21_30C$Time==G21_30C$Time[G21_30C$Time>6][1]))
G21_30C.cut<-G21_30C[st:end,]
G21_30C.cut$Vmax.norm<-G21_30C.cut$Vmax/max(G21_30C.cut$Vmax)

#35 deg. C
G21_35C<-read.csv("GLSE21_long35.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(G21_35C$Time==G21_35C$Time[G21_35C$Time>2][1]))
print(end<-which(G21_35C$Time==G21_35C$Time[G21_35C$Time>6][1]))
G21_35C.cut<-G21_35C[st:end,]
G21_35C.cut$Vmax.norm<-G21_35C.cut$Vmax/max(G21_35C.cut$Vmax)
G21_35C.cut$Vmax.norm<-G21_35C.cut$Vmax.norm/0.7

#40 deg. C
G21_40C<-read.csv("GLSE21_long40.csv")
print(st<-which(G21_40C$Time==G21_40C$Time[G21_40C$Time>0.03][1]))
print(end<-which(G21_40C$Time==G21_40C$Time[G21_40C$Time>6][1]))
G21_40C.cut<-G21_40C[st:end,]
G21_40C.cut$Vmax.norm<-G21_40C.cut$Vmax/max(G21_40C.cut$Vmax)

##
#26:20 deg. C growing temperature
##

#25 deg. C
G26_25C<-read.csv("GLSE26_long25.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(G26_25C$Time==G26_25C$Time[G26_25C$Time>1][1]))
print(end<-length(G26_25C$Time))
G26_25C.cut<-G26_25C[st:end,]
G26_25C.cut$Vmax.norm<-G26_25C.cut$Vmax/max(G26_25C.cut$Vmax)

#30 deg. C
G26_30C<-read.csv("GLSE26_long30.csv")
print(st<-which(G26_30C$Time==G26_30C$Time[G26_30C$Time>0.24][1]))
print(end<-length(G26_30C$Time))
G26_30C.cut<-G26_30C[st:end,]
G26_30C.cut$Vmax.norm<-G26_30C.cut$Vmax/max(G26_30C.cut$Vmax)

#35 deg. C
G26_35C<-read.csv("GLSE26_long35.csv")
print(st<-which(G26_35C$Time==G26_35C$Time[G26_35C$Time>0.34][1]))
print(end<-length(G26_35C$Time))
G26_35C.cut<-G26_35C[st:end,]
G26_35C.cut$Vmax.norm<-G26_35C.cut$Vmax/max(G26_35C.cut$Vmax)

#40 deg. C
G26_40C<-read.csv("GLSE26_long40.csv")
print(st<-which(G26_40C$Time==G26_40C$Time[G26_40C$Time>0.07][1]))
print(end<-length(G26_40C$Time))
G26_40C.cut<-G26_40C[st:end,]
G26_40C.cut$Vmax.norm<-G26_40C.cut$Vmax/max(G26_40C.cut$Vmax)

##
#31:25 deg. C growing temperature
##

#30 deg. C
G31_30C<-read.csv("GLSE31_long30.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(G31_30C$Time==G31_30C$Time[G31_30C$Time>1.5][1]))
print(end<-which(G31_30C$Time==G31_30C$Time[G31_30C$Time>6][1]))
G31_30C.cut<-G31_30C[st:end,]
G31_30C.cut$Vmax.norm<-G31_30C.cut$Vmax/max(G31_30C.cut$Vmax)
G31_30C.cut$Vmax.norm<-G31_30C.cut$Vmax.norm/0.8

#35 deg. C
G31_35C<-read.csv("GLSE31_long35.csv")
print(st<-which(G31_35C$Time==G31_35C$Time[G31_35C$Time>3.9][1]))
print(end<-which(G31_35C$Time==G31_35C$Time[G31_35C$Time>6][1]))
G31_35C.cut<-G31_35C[st:end,]
G31_35C.cut$Vmax.norm<-G31_35C.cut$Vmax/max(G31_35C.cut$Vmax)

#40 deg. C
G31_40C<-read.csv("GLSE31_long40.csv")
print(st<-which(G31_40C$Time==G31_40C$Time[G31_40C$Time>0.18][1]))
print(end<-which(G31_40C$Time==G31_40C$Time[G31_40C$Time>6][1]))
G31_40C.cut<-G31_40C[st:end,]
G31_40C.cut$Vmax.norm<-G31_40C.cut$Vmax/max(G31_40C.cut$Vmax)

###
#Robinia
###

##
#21:15 deg. C growing temperature
##

#25 deg. C
R21_25C<-read.csv("ROPS21_long25.csv")
print(st<-which(R21_25C$Time==R21_25C$Time[R21_25C$Time>3.93][1]))
print(end<-which(R21_25C$Time==R21_25C$Time[R21_25C$Time>6][1]))
R21_25C.cut<-R21_25C[st:end,]
R21_25C.cut$Vmax.norm<-R21_25C.cut$Vmax/max(R21_25C.cut$Vmax)

#30 deg. C
R21_30C<-read.csv("ROPS21_long30.csv")
print(st<-which(R21_30C$Time==R21_30C$Time[R21_30C$Time>1.13][1]))
print(end<-which(R21_30C$Time==R21_30C$Time[R21_30C$Time>6][1]))
R21_30C.cut<-R21_30C[st:end,]
R21_30C.cut$Vmax.norm<-R21_30C.cut$Vmax/max(R21_30C.cut$Vmax)

#35 deg. C
R21_35C<-read.csv("ROPS21_long35.csv")
print(st<-which(R21_35C$Time==R21_35C$Time[R21_35C$Time>0.15][1]))
print(end<-which(R21_35C$Time==R21_35C$Time[R21_35C$Time>6][1]))
R21_35C.cut<-R21_35C[st:end,]
R21_35C.cut$Vmax.norm<-R21_35C.cut$Vmax/max(R21_35C.cut$Vmax)

#40 deg. C
R21_40C<-read.csv("ROPS21_long40.csv")
#Data were cut after 0.5 hr because N-fixation rates were essentially zero after that
#point and convergence was poor if model fit included those data
print(st<-which(R21_40C$Time==R21_40C$Time[R21_40C$Time>0.08][1]))
print(end<-which(R21_40C$Time==R21_40C$Time[R21_40C$Time>0.5][1]))
R21_40C.cut<-R21_40C[st:end,]
R21_40C.cut$Vmax.norm<-R21_40C.cut$Vmax/max(R21_40C.cut$Vmax)

##
#26:20 deg. C growing temperature
##

#25 deg. C
R26_25C<-read.csv("ROPS26_long25.csv")
#There was no decline so breakpoint was visually estimated
print(st<-which(R26_25C$Time==R26_25C$Time[R26_25C$Time>0.5][1]))
print(end<-length(R26_25C$Time))
R26_25C.cut<-R26_25C[st:end,]
R26_25C.cut$Vmax.norm<-R26_25C.cut$Vmax/max(R26_25C.cut$Vmax)

#30 deg. C
R26_30C<-read.csv("ROPS26_long30.csv")
print(st<-which(R26_30C$Time==R26_30C$Time[R26_30C$Time>0.98][1]))
print(end<-length(R26_30C$Time))
R26_30C.cut<-R26_30C[st:end,]
R26_30C.cut$Vmax.norm<-R26_30C.cut$Vmax/max(R26_30C.cut$Vmax)

#35 deg. C
R26_35C<-read.csv("ROPS26_long35.csv")
#Breakpoint was visual estimated due to unequilibrated rates with fit breakpoint
print(st<-which(R26_35C$Time==R26_35C$Time[R26_35C$Time>0.28][1]))
print(end<-length(R26_35C$Time))
R26_35C.cut<-R26_35C[st:end,]
R26_35C.cut$Vmax.norm<-R26_35C.cut$Vmax/max(R26_35C.cut$Vmax)

#40 deg. C
R26_40C<-read.csv("ROPS26_long40.csv")
print(st<-which(R26_40C$Time==R26_40C$Time[R26_40C$Time>0.34][1]))
print(end<-length(R26_40C$Time))
R26_40C.cut<-R26_40C[st:end,]
R26_40C.cut$Vmax.norm<-R26_40C.cut$Vmax/max(R26_40C.cut$Vmax)

##
#31:25 deg. C growing temperature
##

#25 deg. C
R31_25C<-read.csv("ROPS31_long25.csv")
#Breakpoint was visual estimated due to unequilibrated rates with fit breakpoint
print(st<-which(R31_25C$Time==R31_25C$Time[R31_25C$Time>3][1]))
print(end<-which(R31_25C$Time==R31_25C$Time[R31_25C$Time>6][1]))
R31_25C.cut<-R31_25C[st:end,]
R31_25C.cut$Vmax.norm<-R31_25C.cut$Vmax/max(R31_25C.cut$Vmax)

#30 deg. C
R31_30C<-read.csv("ROPS31_long30.csv")
print(st<-which(R31_30C$Time==R31_30C$Time[R31_30C$Time>4.49][1]))
print(end<-which(R31_30C$Time==R31_30C$Time[R31_30C$Time>6][1]))
R31_30C.cut<-R31_30C[st:end,]
R31_30C.cut$Vmax.norm<-R31_30C.cut$Vmax/max(R31_30C.cut$Vmax)

#35 deg. C
R31_35C<-read.csv("ROPS31_long35.csv")
print(st<-which(R31_35C$Time==R31_35C$Time[R31_35C$Time>0.15][1]))
print(end<-which(R31_35C$Time==R31_35C$Time[R31_35C$Time>6][1]))
R31_35C.cut<-R31_35C[st:end,]
R31_35C.cut$Vmax.norm<-R31_35C.cut$Vmax/max(R31_35C.cut$Vmax)

#40 deg. C
R31_40C<-read.csv("ROPS31_long40.csv")
print(st<-which(R31_40C$Time==R31_40C$Time[R31_40C$Time>0.08][1]))
print(end<-which(R31_40C$Time==R31_40C$Time[R31_40C$Time>6][1]))
R31_40C.cut<-R31_40C[st:end,]
R31_40C.cut$Vmax.norm<-R31_40C.cut$Vmax/max(R31_40C.cut$Vmax)

###############################################################################################################
#Now fit falling function from Equation S4 to each individual run
###############################################################################################################

####
#Negative log-likelihood function
####

fall_one_normNLL <- function(sdNase,a1,b1,c1,d1,t1,Nasedat){
  Nasemean<-fall(a1,exp(b1),c1,d1,t1)
  -(sum(dnorm(Nasedat,mean=Nasemean,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

####
#Maximum likelihood fits
####

###
#Morella
###

##
#21:15 deg. C growing temperature
##

##25 deg. C
M21a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.05),c1=1.3,d1=0.9),
                            data=list(t1=M21_25C.cut$Time,Nasedat=M21_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M21a_fall_one_normNLL)
#Save estimated coefficients
ft<-coef(M21a_fall_one_normNLL)
Mb21a<-ft[3]
Md21a<-ft[5]
Mc21a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M21_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
M21b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.001),c1=1.3,d1=0.8),
                            data=list(t1=M21_30C.cut$Time,Nasedat=M21_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M21b_fall_one_normNLL)
#Save coefficients
ft<-coef(M21b_fall_one_normNLL)
Mb21b<-ft[3]
Md21b<-ft[5]
Mc21b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M21_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
M21c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.005),c1=1.3,d1=0.6),
                            data=list(t1=M21_35C.cut$Time,Nasedat=M21_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M21c_fall_one_normNLL)
#Save coefficients
ft<-coef(M21c_fall_one_normNLL)
Mb21c<-ft[3]
Md21c<-ft[5]
Mc21c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M21_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
M21d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.7),c1=1.3,d1=0.04),
                            data=list(t1=M21_40C.cut$Time,Nasedat=M21_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M21d_fall_one_normNLL)
#Save coefficients
ft<-coef(M21d_fall_one_normNLL)
Mb21d<-ft[3]
Md21d<-ft[5]
Mc21d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M21_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#26:20 deg. C growing temperature
##

##25 deg. C
M26a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.7),
                            data=list(t1=M26_25C.cut$Time,Nasedat=M26_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M26a_fall_one_normNLL)
#Save coefficients
ft<-coef(M26a_fall_one_normNLL)
Mb26a<-ft[3]
Md26a<-ft[5]
Mc26a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M26_25C.cut,xlim=c(0,6),ylim=c(0,1.2))
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
M26b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.6),
                            data=list(t1=M26_30C.cut$Time,Nasedat=M26_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M26b_fall_one_normNLL)
#Save coefficients
ft<-coef(M26b_fall_one_normNLL)
Mb26b<-ft[3]
Md26b<-ft[5]
Mc26b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M26_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
M26c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.95),
                            data=list(t1=M26_35C.cut$Time,Nasedat=M26_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M26c_fall_one_normNLL)
#Save coefficients
ft<-coef(M26c_fall_one_normNLL)
Mb26c<-ft[3]
Md26c<-ft[5]
Mc26c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M26_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
M26d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.4),
                            data=list(t1=M26_40C.cut$Time,Nasedat=M26_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M26d_fall_one_normNLL)
#Save coefficients
ft<-coef(M26d_fall_one_normNLL)
Mb26d<-ft[3]
Md26d<-ft[5]
Mc26d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M26_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#31:25 deg. C growing temperature
##

##30 deg. C
M31a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=M31_30C.cut$Time,Nasedat=M31_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M31a_fall_one_normNLL)
#Save coefficients
ft<-coef(M31a_fall_one_normNLL)
Mb31a<-ft[3]
Md31a<-ft[5]
Mc31a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M31_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
M31b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.1),c1=1.3,d1=0.7),
                            data=list(t1=M31_35C.cut$Time,Nasedat=M31_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M31b_fall_one_normNLL)
#Save coefficients
ft<-coef(M31b_fall_one_normNLL)
Mb31b<-ft[3]
Md31b<-ft[5]
Mc31b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M31_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
M31c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.01),c1=1.3,d1=0.8),
                            data=list(t1=M31_40C.cut$Time,Nasedat=M31_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(M31c_fall_one_normNLL)
#Save coefficients
ft<-coef(M31c_fall_one_normNLL)
Mb31c<-ft[3]
Md31c<-ft[5]
Mc31c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=M31_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

###
#Alnus
###

##
#21:15 deg. C growing temperature
##

##25 deg. C
A21a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=A21_25C.cut$Time,Nasedat=A21_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A21a_fall_one_normNLL)
#Save coefficients
ft<-coef(A21a_fall_one_normNLL)
Ab21a<-ft[3]
Ad21a<-ft[5]
Ac21a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A21_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
A21b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.006),c1=1.3,d1=0.8),
                            data=list(t1=A21_30C.cut$Time,Nasedat=A21_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A21b_fall_one_normNLL)
#Save coefficients
ft<-coef(A21b_fall_one_normNLL)
Ab21b<-ft[3]
Ad21b<-ft[5]
Ac21b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A21_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
A21c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.009),c1=1.3,d1=0.15),
                            data=list(t1=A21_35C.cut$Time,Nasedat=A21_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A21c_fall_one_normNLL)
#Save coefficients
ft<-coef(A21c_fall_one_normNLL)
Ab21c<-ft[3]
Ad21c<-ft[5]
Ac21c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A21_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
A21d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.2),c1=1.3,d1=0.05),
                            data=list(t1=A21_40C.cut$Time,Nasedat=A21_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A21d_fall_one_normNLL)
#Save coefficients
ft<-coef(A21d_fall_one_normNLL)
Ab21d<-ft[3]
Ad21d<-ft[5]
Ac21d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A21_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#26:20 deg. C growing temperature
##

##25 deg. C
A26a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=A26_25C.cut$Time,Nasedat=A26_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A26a_fall_one_normNLL)
#Save coefficients
ft<-coef(A26a_fall_one_normNLL)
Ab26a<-ft[3]
Ad26a<-ft[5]
Ac26a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A26_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
A26b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.7),
                            data=list(t1=A26_30C.cut$Time,Nasedat=A26_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A26b_fall_one_normNLL)
#Save coefficients
ft<-coef(A26b_fall_one_normNLL)
Ab26b<-ft[3]
Ad26b<-ft[5]
Ac26b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A26_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
A26c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.9),
                            data=list(t1=A26_35C.cut$Time,Nasedat=A26_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A26c_fall_one_normNLL)
#Save coefficients
ft<-coef(A26c_fall_one_normNLL)
Ab26c<-ft[3]
Ad26c<-ft[5]
Ac26c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A26_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
A26d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.3),c1=1.3,d1=0.1),
                            data=list(t1=A26_40C.cut$Time,Nasedat=A26_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A26d_fall_one_normNLL)
#Save coefficients
ft<-coef(A26d_fall_one_normNLL)
Ab26d<-ft[3]
Ad26d<-ft[5]
Ac26d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A26_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#31:25 deg. C growing temperature
##

##30 deg. C
A31a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=A31_30C.cut$Time,Nasedat=A31_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A31a_fall_one_normNLL)
#Save coefficients
ft<-coef(A31a_fall_one_normNLL)
Ab31a<-ft[3]
Ad31a<-ft[5]
Ac31a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A31_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
A31b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.2),
                            data=list(t1=A31_35C.cut$Time,Nasedat=A31_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A31b_fall_one_normNLL)
#Save coefficients
ft<-coef(A31b_fall_one_normNLL)
Ab31b<-ft[3]
Ad31b<-ft[5]
Ac31b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A31_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

#40 deg. C
A31c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.04),c1=1.3,d1=0.2),
                            data=list(t1=A31_40C.cut$Time,Nasedat=A31_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(A31c_fall_one_normNLL)
#Save coefficients
ft<-coef(A31c_fall_one_normNLL)
Ab31c<-ft[3]
Ad31c<-ft[5]
Ac31c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=A31_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

###
#Gliricidia
###

##
#21:15 deg. C growing temperature
##

##25 deg. C
G21a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=G21_25C.cut$Time,Nasedat=G21_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G21a_fall_one_normNLL)
#Save coefficients
ft<-coef(G21a_fall_one_normNLL)
Gb21a<-ft[3]
Gd21a<-ft[5]
Gc21a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G21_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
G21b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.02),c1=1.3,d1=0.6),
                            data=list(t1=G21_30C.cut$Time,Nasedat=G21_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G21b_fall_one_normNLL)
#Save coefficients
ft<-coef(G21b_fall_one_normNLL)
Gb21b<-ft[3]
Gd21b<-ft[5]
Gc21b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G21_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
G21c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=G21_35C.cut$Time,Nasedat=G21_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G21c_fall_one_normNLL)
#Save coefficients
ft<-coef(G21c_fall_one_normNLL)
Gb21c<-ft[3]
Gd21c<-ft[5]
Gc21c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G21_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
G21d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.3),c1=1.3,d1=0.4),
                            data=list(t1=G21_40C.cut$Time,Nasedat=G21_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G21d_fall_one_normNLL)
#Save coefficients
ft<-coef(G21d_fall_one_normNLL)
Gb21d<-ft[3]
Gd21d<-ft[5]
Gc21d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G21_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#26:20 deg. C growing temperature
##

##25 deg. C
G26a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=G26_25C.cut$Time,Nasedat=G26_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G26a_fall_one_normNLL)
#Save coefficients
ft<-coef(G26a_fall_one_normNLL)
Gb26a<-ft[3]
Gd26a<-ft[5]
Gc26a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G26_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
G26b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.5),
                            data=list(t1=G26_30C.cut$Time,Nasedat=G26_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G26b_fall_one_normNLL)
#Save coefficients
ft<-coef(G26b_fall_one_normNLL)
Gb26b<-ft[3]
Gd26b<-ft[5]
Gc26b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G26_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
G26c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.3),
                            data=list(t1=G26_35C.cut$Time,Nasedat=G26_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G26c_fall_one_normNLL)
#Save coefficients
ft<-coef(G26c_fall_one_normNLL)
Gb26c<-ft[3]
Gd26c<-ft[5]
Gc26c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G26_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
G26d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.05),
                            data=list(t1=G26_40C.cut$Time,Nasedat=G26_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G26d_fall_one_normNLL)
#Save coefficients
ft<-coef(G26d_fall_one_normNLL)
Gb26d<-ft[3]
Gd26d<-ft[5]
Gc26d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G26_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#31:25 deg. C growing temperature
##

##30 deg. C
G31a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=G31_30C.cut$Time,Nasedat=G31_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G31a_fall_one_normNLL)
#Save coefficients
ft<-coef(G31a_fall_one_normNLL)
Gb31a<-ft[3]
Gd31a<-ft[5]
Gc31a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G31_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
G31b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.05),c1=1.3,d1=0.75),
                            data=list(t1=G31_35C.cut$Time,Nasedat=G31_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G31b_fall_one_normNLL)
#Save coefficients
ft<-coef(G31b_fall_one_normNLL)
Gb31b<-ft[3]
Gd31b<-ft[5]
Gc31b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G31_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
G31c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.25),c1=1.3,d1=0.5),
                            data=list(t1=G31_40C.cut$Time,Nasedat=G31_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(G31c_fall_one_normNLL)
#Save coefficients
ft<-coef(G31c_fall_one_normNLL)
Gb31c<-ft[3]
Gd31c<-ft[5]
Gc31c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=G31_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

###
#Robinia
###

##
#21:15 deg. C growing temperature
##

##25 deg. C
R21a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.01),c1=1.3,d1=0.7),
                            data=list(t1=R21_25C.cut$Time,Nasedat=R21_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R21a_fall_one_normNLL)
#Save coefficients
ft<-coef(R21a_fall_one_normNLL)
Rb21a<-ft[3]
Rd21a<-ft[5]
Rc21a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R21_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
R21b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.02),c1=1.3,d1=0.3),
                            data=list(t1=R21_30C.cut$Time,Nasedat=R21_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R21b_fall_one_normNLL)
#Save coefficients
ft<-coef(R21b_fall_one_normNLL)
Rb21b<-ft[3]
Rd21b<-ft[5]
Rc21b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R21_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
R21c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.1),
                            data=list(t1=R21_35C.cut$Time,Nasedat=R21_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R21c_fall_one_normNLL)
#Save coefficients
ft<-coef(R21c_fall_one_normNLL)
Rb21c<-ft[3]
Rd21c<-ft[5]
Rc21c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R21_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
R21d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1.3,b1=log(1.01),c1=1.3,d1=0.001),
                            data=list(t1=R21_40C.cut$Time,Nasedat=R21_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R21d_fall_one_normNLL)
#Save coefficients
ft<-coef(R21d_fall_one_normNLL)
Rb21d<-ft[3]
Rd21d<-ft[5]
Rc21d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R21_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#26:20 deg. C growing temperature
##

##25 deg. C
R26a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=-1000,c1=1.3,d1=1),
                            data=list(t1=R26_25C.cut$Time,Nasedat=R26_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R26a_fall_one_normNLL)
#Save coefficients
ft<-coef(R26a_fall_one_normNLL)
Rb26a<-ft[3]
Rd26a<-ft[5]
Rc26a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R26_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
R26b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.6),
                            data=list(t1=R26_30C.cut$Time,Nasedat=R26_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R26b_fall_one_normNLL)
#Save coefficients
ft<-coef(R26b_fall_one_normNLL)
Rb26b<-ft[3]
Rd26b<-ft[5]
Rc26b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R26_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
R26c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.3),
                            data=list(t1=R26_35C.cut$Time,Nasedat=R26_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R26c_fall_one_normNLL)
#Save coefficients
ft<-coef(R26c_fall_one_normNLL)
Rb26c<-ft[3]
Rd26c<-ft[5]
Rc26c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R26_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
R26d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.03),c1=1.3,d1=0.1),
                            data=list(t1=R26_40C.cut$Time,Nasedat=R26_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R26d_fall_one_normNLL)
#Save coefficients
ft<-coef(R26d_fall_one_normNLL)
Rb26d<-ft[3]
Rd26d<-ft[5]
Rc26d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R26_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##
#31:25 deg. C growing temperature
##

##25 deg. C
R31a_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.1),c1=1.3,d1=0.97),
                            data=list(t1=R31_25C.cut$Time,Nasedat=R31_25C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R31a_fall_one_normNLL)
#Save coefficients
ft<-coef(R31a_fall_one_normNLL)
Rb31a<-ft[3]
Rd31a<-ft[5]
Rc31a<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R31_25C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##30 deg. C
R31b_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.05),c1=1.3,d1=0.8),
                            data=list(t1=R31_30C.cut$Time,Nasedat=R31_30C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R31b_fall_one_normNLL)
#Save coefficients
ft<-coef(R31b_fall_one_normNLL)
Rb31b<-ft[3]
Rd31b<-ft[5]
Rc31b<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R31_30C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##35 deg. C
R31c_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.1),c1=1.3,d1=0.2),
                            data=list(t1=R31_35C.cut$Time,Nasedat=R31_35C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R31c_fall_one_normNLL)
#Save coefficients
ft<-coef(R31c_fall_one_normNLL)
Rb31c<-ft[3]
Rd31c<-ft[5]
Rc31c<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R31_35C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

##40 deg. C
R31d_fall_one_normNLL<-mle2(fall_one_normNLL,start=list(sdNase=-1,a1=1,b1=log(0.01),c1=1.3,d1=0.08),
                            data=list(t1=R31_40C.cut$Time,Nasedat=R31_40C.cut$Vmax.norm),
                            control=list(maxit=20000))
summary(R31d_fall_one_normNLL)
#Save coefficients
ft<-coef(R31d_fall_one_normNLL)
Rb31d<-ft[3]
Rd31d<-ft[5]
Rc31d<-ft[4]
#Plot data and curve
plot(Vmax.norm~Time,cex=0.01,data=R31_40C.cut)
curve(fall(ft[2],exp(ft[3]),ft[4],ft[5],x),from=0,to=6,col="orange",lwd=2,add=T)

####
#Compile estimated parameters into a table
####

#Make species-level vectors of measured temperatures
T.a.M<-c(25.6,30.0,34.9,39.1,26.1,30.6,35.4,39.0,28.6,34.6,39.2) #Morella
T.a.A<-c(24.7,30.5,34.9,40.0,25.9,30.2,35.4,39.6,28.9,35.2,39.9) #Alnus
T.a.G<-c(25.5,30.5,34.5,39.3,25.4,30.4,35.2,39.8,29.0,35.6,39.3) #Gliricidia
T.a.R<-c(25.2,30.6,35.0,39.8,25.7,29.4,35.2,39.4,24.2,31.0,34.7,39.5) #Robinia

#Redefine temperature as tau
T.d.M<-c(c(T.a.M[1:4]-29.03),c(T.a.M[5:8]-32.94),c(T.a.M[9:11]-36.86)) #Morella
T.d.A<-c(c(T.a.A[1:4]-32.42),c(T.a.A[5:8]-32.74),c(T.a.A[9:11]-33.07)) #Alnus
T.d.G<-c(c(T.a.G[1:4]-31.58),c(T.a.G[5:8]-33.66),c(T.a.G[9:11]-35.74)) #Gliricidia
T.d.R<-c(c(T.a.R[1:4]-31.89),c(T.a.R[5:8]-32.02),c(T.a.R[9:12]-32.14)) #Robinia

b.all<-c(Ab21a,Ab21b,Ab21c,Ab21d,Ab26a,Ab26b,Ab26c,Ab26d,Ab31a,Ab31b,Ab31c,
         Mb21a,Mb21b,Mb21c,Mb21d,Mb26a,Mb26b,Mb26c,Mb26d,Mb31a,Mb31b,Mb31c,
         Gb21a,Gb21b,Gb21c,Gb21d,Gb26a,Gb26b,Gb26c,Gb26d,Gb31a,Gb31b,Gb31c,
         Rb21a,Rb21b,Rb21c,Rb21d,Rb26a,Rb26b,Rb26c,Rb26d,Rb31a,Rb31b,Rb31c,Rb31d) #b (must be exponentiated)
c.all<-c(Ac21a,Ac21b,Ac21c,Ac21d,Ac26a,Ac26b,Ac26c,Ac26d,Ac31a,Ac31b,Ac31c,
         Mc21a,Mc21b,Mc21c,Mc21d,Mc26a,Mc26b,Mc26c,Mc26d,Mc31a,Mc31b,Mc31c,
         Gc21a,Gc21b,Gc21c,Gc21d,Gc26a,Gc26b,Gc26c,Gc26d,Gc31a,Gc31b,Gc31c,
         Rc21a,Rc21b,Rc21c,Rc21d,Rc26a,Rc26b,Rc26c,Rc26d,Rc31a,Rc31b,Rc31c,Rc31d) #c
d.all<-c(Ad21a,Ad21b,Ad21c,Ad21d,Ad26a,Ad26b,Ad26c,Ad26d,Ad31a,Ad31b,Ad31c,
         Md21a,Md21b,Md21c,Md21d,Md26a,Md26b,Md26c,Md26d,Md31a,Md31b,Md31c,
         Gd21a,Gd21b,Gd21c,Gd21d,Gd26a,Gd26b,Gd26c,Gd26d,Gd31a,Gd31b,Gd31c,
         Rd21a,Rd21b,Rd21c,Rd21d,Rd26a,Rd26b,Rd26c,Rd26d,Rd31a,Rd31b,Rd31c,Rd31d) #d
Tr.all<-c("A21","A21","A21","A21","A26","A26","A26","A26","A31","A31","A31",
          "M21","M21","M21","M21","M26","M26","M26","M26","M31","M31","M31",
          "G21","G21","G21","G21","G26","G26","G26","G26","G31","G31","G31",
          "R21","R21","R21","R21","R26","R26","R26","R26","R31","R31","R31","R31") #Treatment
Tau.all<-c(T.d.A,T.d.M,T.d.G,T.d.R) #Tau
Temp.all<-c(T.a.A,T.a.M,T.a.G,T.a.R) #Measurement temperature

#Save estimates as a table
write.csv(cbind(Tr=Tr.all,Temp=Temp.all,Tau=Tau.all,b=exp(b.all),c=c.all,d=d.all),"Temp_Decline_Parms.csv")

##
#Outliers (Tr,Temp,Tau,b,c,d)
##

#M26,26.1,-6.84,0.62,1.15,0.67
#R26,35.2,3.18,6.42,0.31,0.03
#R26,39.4,7.38,5.33,0.46,0.01

###############################################################################################################
#Solve for median value of "c" and fit Equations S5 and S6 to "b" and "d", respectively
###############################################################################################################

###
#Set "c" estimates to NA if there was no decline and for outliers in order to calculate median "c" value
###

c.all.rm<-rep(NA,length(c.all))
for(i in 1:length(c.all)){
  if(c.all[i]==1.3){
    c.all.rm[i]<-NA
  } else {
    c.all.rm[i]<-c.all[i]
  }
}
c.all.rm[16]<-NA
c.all.rm[40:41]<-NA

print(c.est<-median(na.omit(c.all.rm))) #0.9572269

#95% CI
quantile(na.omit(c.all.rm),0.025)
quantile(na.omit(c.all.rm),0.975)

###
#Model "b" and "d" as functions of tau
###

#Read in table (created above)
temp.parms<-read.csv("Temp_Decline_Parms.csv")

#Remove outliers
temp.parms.rm<-rbind(temp.parms[1:15,],temp.parms[17:39,],temp.parms[42:45,])

###
#NLL functions
###

#Function for modeling "b"
b_normNLL <- function(sdNase,beta,alpha,x,bdat){
  b.nll <- exp.mod(alpha,beta,x)
  -(sum(dnorm(bdat,mean=b.nll,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

#Function for modeling "d"
d_normNLL <- function(sdNase,p,q,x,ddat){
  d.nll <- nsig.mod(p,q,x)
  -(sum(dnorm(ddat,mean=d.nll,sd=exp(sdNase),log=TRUE),na.rm=TRUE))
}

###
#Fit NLL functions with maximum likelihood
###

##
#b
##

#Starting values for modeling b as a function of tau
alpha.est<-0.02114
beta.est<-0.33945

#Fit regular exponential function
b.mle<-mle2(b_normNLL,start=list(sdNase=-1,beta=log(beta.est),alpha=log(alpha.est)),
            data=list(x=temp.parms.rm$Tau,bdat=temp.parms.rm$b),
            control=list(maxit=20000))

#Fit Equation S5 (notice transformed data)
b.mle2<-mle2(b_normNLL,start=list(sdNase=-1,beta=log(beta.est),alpha=log(alpha.est)),
             data=list(x=temp.parms.rm$Tau,bdat=log(temp.parms.rm$b+1)),
             control=list(maxit=20000))
b.mle2.ci<-confint(b.mle2)

#Compare models for b~tau
AICctab(b.mle,b.mle2,nobs=42) #b.mle2 (Equation S5) is best

#Estimated parameters
print(alpha.mle2<-exp(coef(b.mle2)[[3]])) #alpha (this is equivalent to alpha_b in the manuscript)
print(beta.mle2<-exp(coef(b.mle2)[[2]])) #beta (this is equivalent to beta_b in the manuscript)

#95% CI
exp(b.mle2.ci[3,1]);exp(b.mle2.ci[3,2]) #alpha
exp(b.mle2.ci[2,1]);exp(b.mle2.ci[2,2]) #beta

##
#d
##

#Starting values for modeling d as a function of tau
p.est<-0.7022
q.est<-0.2528

#Fit Equation S6
d.mle<-mle2(d_normNLL,start=list(sdNase=-1,p=log(p.est),q=log(q.est)),
            data=list(x=temp.parms.rm$Tau,ddat=temp.parms.rm$d),
            control=list(maxit=20000))
d.mle.ci<-confint(d.mle)

#Estimated parameters
print(p.mle<-exp(coef(d.mle)[[2]])) #p (this is equivalent to alpha_d in the manuscript
print(q.mle<-exp(coef(d.mle)[[3]])) #q (this is equivalent to beta_d in the manuscript)

#95% CI
exp(d.mle.ci[2,1]);exp(d.mle.ci[2,2]) #p
exp(d.mle.ci[3,1]);exp(d.mle.ci[3,2]) #q
