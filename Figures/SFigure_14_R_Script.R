###############################################################################################################
###############################################################################################################
#This script generates Supplementary Figure 14
###############################################################################################################
###############################################################################################################

####
#Figure were exported as PDFs from R Studio
#PDF dimension is 10.5 x 3.5 inches
####

#Read in data
x1.dat<-read.csv("x1.Topt.r15.r40.csv") #simulated data at measurement heating rate
og.dat<-read.csv("og.Topt.r15.r40.csv") #measured data (see Nfix_Model_Comparison.R for how these values were calculated)

#Convert rates from proportion to percent for plotting
x1.dat$r15<-x1.dat$r15*100
x1.dat$r40<-x1.dat$r40*100
og.dat$r15<-og.dat$r15*100
og.dat$r40<-og.dat$r40*100

#Plotting settings
nf<-layout(matrix(c(1,2,3,4),1,4,byrow=T),c(4,4,4,2),c(4,4,4,4),T)
layout.show(nf)
par(mar=c(4,6,4,1))
par(oma=c(0,0,0,0))
par(pty="s")

#S14a
plot(25:40,25:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(0,100),ylim=c(0,100),cex.lab=1.5,cex.axis=1.2)
mtext(expression('N fixation at 15 '*degree*'C'),side=2,cex=1.2,line=4.75)
mtext(expression('(simulated; % of max)'),side=2,cex=1.2,line=3)
mtext(expression('N fixation at 15 '*degree*'C'),side=1,cex=1.2,line=3)
mtext(expression('(measured; % of max)'),side=1,cex=1.2,line=4.75)
abline(a=0,b=1,lty=2)
points(c(x1.dat$r15)~c(og.dat$r15),pch=c(16,16,16,16,16,16,17,17,17,17,17,17),
       col=c("darkorange1","darkorange1","darkorange1","darkturquoise","darkturquoise","darkturquoise",
             "orangered2","orangered2","orangered2","dodgerblue3","dodgerblue3","dodgerblue3"),cex=1.5)
points(c(x1.dat$r15)~c(og.dat$r15),pch=c(1,1,1,1,1,1,2,2,2,2,2,2),
       col=c("dodgerblue1","gold1","orangered3","dodgerblue1","gold1","orangered3",
             "dodgerblue1","gold1","orangered3","dodgerblue1","gold1","orangered3"),cex=1.5)
points(c(x1.dat$r15[10])~c(og.dat$r15[10]),pch=17,col="dodgerblue3",cex=1.5)
points(c(x1.dat$r15[10])~c(og.dat$r15[10]),pch=2,col="dodgerblue1",cex=1.5)
points(c(x1.dat$r15[11])~c(og.dat$r15[11]),pch=17,col="dodgerblue3",cex=1.5)
points(c(x1.dat$r15[11])~c(og.dat$r15[11]),pch=2,col="gold1",cex=1.5)
points(c(x1.dat$r15[12])~c(og.dat$r15[12]),pch=17,col="dodgerblue3",cex=1.5)
points(c(x1.dat$r15[12])~c(og.dat$r15[12]),pch=2,col="orangered3",cex=1.5)
mtext(text="a",side=3,cex=1.2,adj=0)

#S14b
plot(25:40,25:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(25,40),ylim=c(25,40),cex.lab=1.5,cex.axis=1.2)
mtext(expression(italic('T')[opt]),side=2,cex=1.2,line=4.75)
mtext(expression('(simulated; '*degree*'C)'),side=2,cex=1.2,line=3)
mtext(expression(italic('T')[opt]),side=1,cex=1.2,line=3)
mtext(expression('(measured; '*degree*'C)'),side=1,cex=1.2,line=4.75)
abline(a=0,b=1,lty=2)
points(x1.dat$Topt~og.dat$Topt,pch=c(16,16,16,16,16,16,17,17,17,17,17,17),
       col=c("darkorange1","darkorange1","darkorange1","darkturquoise","darkturquoise","darkturquoise",
             "orangered2","orangered2","orangered2","dodgerblue3","dodgerblue3","dodgerblue3"),cex=1.5)
points(x1.dat$Topt~og.dat$Topt,pch=c(1,1,1,1,1,1,2,2,2,2,2,2),
       col=c("dodgerblue1","gold1","orangered3","dodgerblue1","gold1","orangered3",
             "dodgerblue1","gold1","orangered3","dodgerblue1","gold1","orangered3"),cex=1.5)
points(x1.dat$Topt[11]~og.dat$Topt[11],pch=17,col="dodgerblue3",cex=1.5)
points(x1.dat$Topt[11]~og.dat$Topt[11],pch=2,col="gold1",cex=1.5)
points(x1.dat$Topt[12]~og.dat$Topt[12],pch=17,col="dodgerblue3",cex=1.5)
points(x1.dat$Topt[12]~og.dat$Topt[12],pch=2,col="orangered3",cex=1.5)
points(x1.dat$Topt[5]~og.dat$Topt[5],pch=16,col="darkturquoise",cex=1.5)
points(x1.dat$Topt[5]~og.dat$Topt[5],pch=1,col="gold1",cex=1.5)
points(x1.dat$Topt[6]~og.dat$Topt[6],pch=16,col="darkturquoise",cex=1.5)
points(x1.dat$Topt[6]~og.dat$Topt[6],pch=1,col="orangered3",cex=1.5)
mtext(text="b",side=3,cex=1.2,adj=0)

#S14c
plot(25:40,25:40,col="white",xlab=NA,ylab=NA,las=1,xlim=c(0,100),ylim=c(0,100),cex.lab=1.5,cex.axis=1.2)
mtext(expression('N fixation at 40 '*degree*'C'),side=2,cex=1.2,line=4.75)
mtext(expression('(simulated; % of max)'),side=2,cex=1.2,line=3)
mtext(expression('N fixation at 40 '*degree*'C'),side=1,cex=1.2,line=3)
mtext(expression('(measured; % of max)'),side=1,cex=1.2,line=4.75)
abline(a=0,b=1,lty=2)
points(c(x1.dat$r40)~c(og.dat$r40),pch=c(16,16,16,16,16,16,17,17,17,17,17,17),
       col=c("darkorange1","darkorange1","darkorange1","darkturquoise","darkturquoise","darkturquoise",
             "orangered2","orangered2","orangered2","dodgerblue3","dodgerblue3","dodgerblue3"),cex=1.5)
points(c(x1.dat$r40)~c(og.dat$r40),pch=c(1,1,1,1,1,1,2,2,2,2,2,2),
       col=c("dodgerblue1","gold1","orangered3","dodgerblue1","gold1","orangered3",
             "dodgerblue1","gold1","orangered3","dodgerblue1","gold1","orangered3"),cex=1.5)
points(c(x1.dat$r40[10])~c(og.dat$r40[10]),pch=17,col="dodgerblue3",cex=1.5)
points(c(x1.dat$r40[10])~c(og.dat$r40[10]),pch=2,col="dodgerblue1",cex=1.5)
points(c(x1.dat$r40[11])~c(og.dat$r40[11]),pch=17,col="dodgerblue3",cex=1.5)
points(c(x1.dat$r40[11])~c(og.dat$r40[11]),pch=2,col="gold1",cex=1.5)
points(c(x1.dat$r40[12])~c(og.dat$r40[12]),pch=17,col="dodgerblue3",cex=1.5)
points(c(x1.dat$r40[12])~c(og.dat$r40[12]),pch=2,col="orangered3",cex=1.5)
points(c(x1.dat$r40[1])~c(og.dat$r40[1]),pch=16,col="darkorange1",cex=1.5)
points(c(x1.dat$r40[1])~c(og.dat$r40[1]),pch=1,col="dodgerblue1",cex=1.5)
mtext(text="c",side=3,cex=1.2,adj=0)

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
       col=c(NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3",
             NA,"dodgerblue1","gold1","orangered3",NA,"dodgerblue1","gold1","orangered3"),
       pch=c(NA,1,1,1,NA,1,1,1,NA,2,2,2,NA,2,2,2),bty="n",y.intersp = 0.8,cex=1,x.intersp = 0.7,pt.cex=1.5,xpd=T)