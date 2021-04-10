## note - need to first create/load data1 and data2 - use data_mgt.R


#############################################
load("mcmcSDE_set1.RData"); #provides myresGp

myresGpA<-window(myresGp,start=2000);
rm(myresGp);

load("mcmcSDE_set2.RData"); #provides myresGp

myresGpB<-window(myresGp,start=2000);
rm(myresGp);


myFcastA<-summary(myresGpA[,c(8:27),]);# 8:27, actual, 28:47 = latent, 48:57, fits
lowerA<-myFcastA$quantiles[,"2.5%"]
midA<-myFcastA$quantiles[,"50%"]
upperA<-myFcastA$quantiles[,"97.5%"]

myFitsA<-summary(myresGpA[,c(48:57),]);# 8:27, actual, 28:47 = latent, 48:57, fits
midFitA<-myFitsA$quantiles[,"50%"]

myFcastB<-summary(myresGpB[,c(8:27),]);# 8:27, actual, 28:47 = latent, 48:57, fits
lowerB<-myFcastB$quantiles[,"2.5%"]
midB<-myFcastB$quantiles[,"50%"]
upperB<-myFcastB$quantiles[,"97.5%"]

myFitsB<-summary(myresGpB[,c(48:57),]);# 8:27, actual, 28:47 = latent, 48:57, fits
midFitB<-myFitsB$quantiles[,"50%"]


####################################################
############# plot
png("plotNLSDE1.png",width=480*1.5,height=480*1.5,pointsize=14);
# Setup plotting environment
mycol<-"black";
par(mar=c(5,6,5,5));
par(cex.axis=2.0);par(cex.lab=2.5);par(bg="white");par(fg=mycol);par(col.axis=mycol);par(col=mycol);par(col.main=mycol);
par(cex.main=2.0);par(col.lab=mycol);par(las=1);par(xaxs="r");par(yaxs="r");
par(mfrow=c(1,1));
# end setup


plot(data1A$Time2[1:231],data1A$Observation2[1:231],pch=21,type="n",axes=FALSE,ylab="",xlab="",
     xlim=c(0,25),ylim=c(100,1200))
rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col = "grey90")
grid(col="grey60")

lines(data1A$Time2[1:231],data1A$Observation2[1:231],lwd=2,lty=1,col="skyblue");
points(data1A$Time2[1:231],data1A$Observation2[1:231],pch=23,col="skyblue",bg="darkblue");

lines(data1B$Time2[1:231],data1B$Observation2[1:231],lwd=2,lty=1,col="green");
points(data1B$Time2[1:231],data1B$Observation2[1:231],pch=23,col="green",bg="darkgreen");


axis(1,seq(0,25,5),padj=0.3,cex.axis=1.5);
mtext("Time (days)",1,line=3.5,cex=2);par(las=1);
axis(2,cex.axis=1.5);par(las=3);
mtext(expression('Volume mm' ^3), 2, line=4,cex=2);par(las=1);
axis(4,padj=0.3,cex.axis=1.5);par(las=3);
mtext("",4,line=3,cex=2);
title("Bayesian diffusion forecasting",line=2);par(las=1);

x<-seq(23.1,25.0,by=0.1)

#points(c(x[1],x[1]),c(lower[1],upper[1]),pch=22,col="red",bg="black")
for(i in 1:length(lowerA)){
lines(c(x[i],x[i]),c(lowerA[i],upperA[i]),pch=22,col="grey60",bg="grey60")
  points(c(x[i]),c(midA[i]),pch="+",col="black",bg="black")

  lines(c(x[i],x[i]),c(lowerB[i],upperB[i]),pch=22,col="grey60",bg="grey60")
  points(c(x[i]),c(midB[i]),pch="+",col="black",bg="black")
  
  }

x<-data1A$Time1[251:260]
#points(c(x[1],x[1]),c(lower[1],upper[1]),pch=22,col="red",bg="black")
for(i in 1:length(lowerA)){
  points(c(x[i]),c(midFitA[i]),pch=24,col="blue",bg="magenta",cex=1.0)
  points(c(x[i]),c(midFitB[i]),pch=24,col="blue",bg="magenta",cex=1.0)

  }
box();


dev.off()








