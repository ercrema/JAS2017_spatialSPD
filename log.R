######################
### Load Libraries ###
######################

library(R.utils)
library(sp)
library(doParallel)
library(devtools)
library(zoo)

################################
### Install and Load rcarbon ###
################################

# The following installs the specific github commit used in the paper
# An up-to-date version of rcarbon can be installed with the following command:
# install_github("ahb108/rcarbon")
# Note that the updated version may not work or generate slightly different result. 

install_github("ahb108/rcarbon@1530ae2")
library(rcarbon)


########################################################
####### Case Study N.1: Generate Artifical Data ########
########################################################



# The following script generate a sample of radiocarbon dates from a spatially varying demographic process.
# The key objects that will be generated are two data.frames, one with the full data set (samplePoints) and the latter with a simulated effect of unven sampling across space (samplePointsThinned). Both data.frames contain the following columns:
# x ... x coordinate of each sample
# y ... y coordinate of each sample
# t ... "true" calendar date (in BP) of the sample
# zone ... zone of each sample. One between A,B,C, and D. Each zone are assumed to have distinct population trajectory (except for A and D which are assumed to be identical
# c14age ... back-calibrated C14 Age
# c14error ... error of the back-calibrated C14 Age
# SiteID ... unique identifier of the site where the samples were recovered. Samples with identical spatial coordinates were assumed to be coming from the same site/location



## Generate intensity curve over time

XX=7000:3001 
AD=data.frame(x=c(7000,5000,4000,3000),y=c(0.1,0.6,0.7,0.7))
fitAD=lm(y~poly(x,3),data=AD)
predAD=as.numeric(predict(fitAD,list(x=XX)))

B=data.frame(x=c(7000,6000,4000,3500,3000),y=c(0.1,0.9,0.4,0.3,0.2))
fitB=lm(y~poly(x,4),data=B)
predB=as.numeric(predict(fitB,list(x=XX)))

C=data.frame(x=c(7000,6500,6000,4000,3500,3000),y=c(0.1,0.15,0.16,0.3,0.9,0.95))
fitC=lm(y~poly(x,3),data=C)
predC=as.numeric(predict(fitC,list(x=XX)))


## Generate probability vectors in radiocarbon age

predADdat=data.frame(calBP=XX, PrDens=predAD)
predADdat$PrDens <- predADdat$PrDens/sum(predADdat$PrDens)
class(predADdat)="CalGrid"
predADdat <- uncalibrate(predADdat,compact=FALSE)


predBdat=data.frame(calBP=XX, PrDens=predB)
predBdat$PrDens <- predBdat$PrDens/sum(predBdat$PrDens)
class(predBdat)="CalGrid"
predBdat <- uncalibrate(predBdat,compact=FALSE)


predCdat=data.frame(calBP=XX, PrDens=predC)
predCdat$PrDens <- predCdat$PrDens/sum(predCdat$PrDens)
class(predCdat)="CalGrid"
predCdat <- uncalibrate(predCdat,compact=FALSE)



## Copy intensities in the baseMatrix array: 

baseMatrix<-array(NA,dim=c(46331,40,40))
baseMatrix[,1:20,1:20]=predADdat$PrDens
baseMatrix[,21:40,21:40]=predADdat$PrDens
baseMatrix[,1:20,21:40]=predCdat$PrDens
baseMatrix[,21:40,1:20]=predBdat$PrDens


## Samples radiocarbon dates from from the baseMatrix array

set.seed(12345) 
n=5000 #number of 14C dates to be sampled

samplePoints<-sample(1:length(baseMatrix),replace=TRUE,size=n,prob=as.numeric(baseMatrix))

x=numeric(length=n)  
y=numeric(length=n)  
t=numeric(length=n)

for (i in 1:n)
    {
        tmp=arrayInd(samplePoints[i],.dim=dim(baseMatrix))
        t[i]=tmp[1]
        x[i]=tmp[2]
        y[i]=tmp[3]
    }

samplePoints=data.frame(x=x,y=y,t=t)
samplePoints$y=c(1:40)[match(samplePoints$y,c(40:1))]


## Assign Zones to sample points

samplePoints$zone=NA
samplePoints$zone[which(samplePoints$x<=20&samplePoints$y>=20)]="A"
samplePoints$zone[which(samplePoints$x<=20&samplePoints$y<=20)]="C"
samplePoints$zone[which(samplePoints$x>=20&samplePoints$y>=20)]="B"
samplePoints$zone[which(samplePoints$x>=20&samplePoints$y<=20)]="D"

samplePoints$t=predCdat$CRA[samplePoints$t]
samplePoints$c14age=samplePoints$t
samplePoints$c14error=sample(20:60,size=nrow(samplePoints),replace=TRUE)

## Simulate research intensity bias (reduce 50% of the samples in quadrants C & D)

index=which(samplePoints$x>20)
indexToRemove=sample(index,size=round(length(index)*0.50),replace=FALSE) 
samplePointsThinned=samplePoints[-indexToRemove,]

## Create "sites" for full data and thinned data

sites=unique(data.frame(x=samplePoints$x,y=samplePoints$y))
sites$SiteID=paste("S",1:nrow(sites),sep="")
samplePoints$SiteID=NA

for (s in 1:nrow(sites))
    {
        index=which(samplePoints$x==sites$x[s]&samplePoints$y==sites$y[s])
        samplePoints$SiteID[index]=sites$SiteID[s]
    }


sites=unique(data.frame(x=samplePointsThinned$x,y=samplePointsThinned$y))
sites$SiteID=paste("S",1:nrow(sites),sep="")
samplePointsThinned$SiteID=NA

for (s in 1:nrow(sites))
    {
        index=which(samplePointsThinned$x==sites$x[s]&samplePointsThinned$y==sites$y[s])
        samplePointsThinned$SiteID[index]=sites$SiteID[s]
    }


#########################################
####### Case Study N.1: Analysis ########
#########################################


## Compute summed probability distributions for each zone for both full and thinned datasets 

zones <- c("A","B","C","D")
timeRange=c(7000,3000) #time range of analysis

for (z in 1:length(zones))
{
tmp=subset(samplePoints,zone==zones[z])
bin=binPrep(sites=tmp$SiteID,ages=tmp$c14age,h=200)
cal=calibrate(x=tmp$c14age,errors=as.integer(tmp$c14error),timeRange=timeRange,normalised=FALSE)
assign(paste("spd.",zones[z],sep=""),spd(x=cal,timeRange=timeRange,bins=bin,spdnormalised=TRUE))


tmp=subset(samplePointsThinned,zone==zones[z])
bin=binPrep(sites=tmp$SiteID,ages=tmp$c14age,h=200)
cal=calibrate(x=tmp$c14age,errors=as.integer(tmp$c14error),timeRange=timeRange,normalised=FALSE)
assign(paste("spdT.",zones[z],sep=""),spd(x=cal,timeRange=timeRange,bins=bin,spdnormalised=TRUE))
}

## Compute "sampled" Geometric Growth Rates

breakSize=500
breaks=seq(7000,3000,-breakSize)
nBreaks=length(breaks)-1


ABCD = matrix(NA,nrow=4,ncol=nBreaks)
ABCDRoca = matrix(NA,nrow=4,ncol=nBreaks-1)

ABCD.t = ABCD
ABCDRoca.t = ABCDRoca

for (z in 1:length(zones))
	{
	tmp <- get(paste("spd.",zones[z],sep=""))
	for (i in 1:nBreaks)
		{
		ABCD[z,i]=sum(tmp$grid$PrDens[tmp$grid$calBP<=breaks[i]&tmp$grid$calBP>=breaks[i+1]])
		}		
	for (i in 1:c(nBreaks-1))
		{
		ABCDRoca[z,i]=(ABCD[z,i+1]/ABCD[z,i])^(1/breakSize)-1
		}
	}


for (z in 1:length(zones))
	{
	tmp <- get(paste("spdT.",zones[z],sep=""))
	for (i in 1:nBreaks)
		{
		ABCD.t [z,i]=sum(tmp$grid$PrDens[tmp$grid$calBP<=breaks[i]&tmp$grid$calBP>=breaks[i+1]])
		}		
	for (i in 1:c(nBreaks-1))
		{
		ABCDRoca.t [z,i]=(ABCD[z,i+1]/ABCD[z,i])^(1/breakSize)-1
		}
	}

# Compute "true" Geometric Growth Rates

breaksReal=c(0,500,1000,1500,2000,2500,3000,3500,4000)

ADsum=numeric()
Bsum=numeric()
Csum=numeric()


ADsum.r=numeric()
Bsum.r=numeric()
Csum.r=numeric()


for (i in 1:8)
		{
		ADsum[i]=sum(predAD[c(breaksReal[i]+1):breaksReal[i+1]])
		Bsum[i]=sum(predB[c(breaksReal[i]+1):breaksReal[i+1]])
		Csum[i]=sum(predC[c(breaksReal[i]+1):breaksReal[i+1]])
		}		
	for (i in 1:7)
		{
		ADsum.r[i]=(ADsum[i+1]/ADsum[i])^(1/breakSize)-1
		Bsum.r[i]=(Bsum[i+1]/Bsum[i])^(1/breakSize)-1
		Csum.r[i]=(Csum[i+1]/Csum[i])^(1/breakSize)-1

		}





### Spatial Permutation Test  ###

##Define Sample Sites and generate SpatialPoints class object
samplePointsSITES=unique(data.frame(SiteID=samplePoints$SiteID,x=samplePoints$x,y=samplePoints$y))
samplePointsThinnedSITES=unique(data.frame(SiteID=samplePointsThinned$SiteID,x=samplePointsThinned$x,y=samplePointsThinned$y))
samplePointsSITES.sp=samplePointsSITES
samplePointsThinnedSITES.sp=samplePointsThinnedSITES
row.names(samplePointsSITES.sp)=samplePointsSITES.sp$SiteID
row.names(samplePointsThinnedSITES.sp)=samplePointsThinnedSITES.sp$SiteID
samplePointsSITES.sp=samplePointsSITES.sp[,-1]
samplePointsThinnedSITES.sp=samplePointsThinnedSITES.sp[,-1]
coordinates(samplePointsSITES.sp)<-c("x","y")
coordinates(samplePointsThinnedSITES.sp)<-c("x","y")

##Compute Distance Matrix
distSamples=spDists(samplePointsSITES.sp,samplePointsSITES.sp,longlat=FALSE)
distSamplesT=spDists(samplePointsThinnedSITES.sp,samplePointsThinnedSITES.sp,longlat=FALSE)

##Compute Spatial Weights (using a bandwidth of 6 units)
weights6=spweights(distSamples,h=6) 
weightsT6=spweights(distSamplesT,h=6)


## calibrate all dates and generate bins 
#generate bins
bin=binPrep(sites=samplePoints$SiteID,ages=samplePoints$c14age,h=200)
binThinned=binPrep(sites=samplePointsThinned$SiteID,ages=samplePointsThinned$c14age,h=200)
#calibrate dates
cal=calibrate(x=samplePoints$c14age,errors=as.integer(samplePoints$c14error),timeRange=c(7000,3000),normalised=FALSE)  #calibrate dates
calThinned=calibrate(x=samplePointsThinned$c14age,errors=as.integer(samplePointsThinned$c14error),timeRange=c(7000,3000),normalised=FALSE)

## Spatial permutation tests
result.locations=SPpermTest(calDates=cal,bins=bin,timeRange=c(7000,3000),locations=samplePointsSITES.sp,spatialweights=weights6,nsim=10000,breaks=seq(7000,3000,-500),permute="locations",ncores=3)
resultThinned.locations=SPpermTest(calThinned,bins=binThinned,timeRange=c(7000,3000),locations=samplePointsThinnedSITES.sp,spatialweights=weightsT6,nsim=10000,breaks=seq(7000,3000,-500),permute="locations",ncores=3)


#########################################
####### Case Study N.1: Figures  ########
#########################################



### Code for Figure 1 ##

breakSize=500
breaks=seq(7000,3000,-breakSize)
nBreaks=length(breaks)-1

par(mfrow=c(2,2))
par(mar=c(3,3,2,1))
plot(0,xlim=c(7000,3000),ylim=c(0,1),xlab="cal BP",ylab="Relative Intensity",main="a")
lines(XX,predAD)
lines(XX,predB,col="blue")
lines(XX,predC,col="red")
legend(x=4550,y=1,legend=c("A&D","B","C"),lwd=2,col=c("black","blue","red"),bg="white",cex=0.7)



plot(spd.A$grid$calBP, spd.A$grid$SPD,col=rgb(0,0,0,0.5),lwd=0.5,type="n",xlim=c(7000,3000),main="b",ylim=c(0,0.00068),xlab="cal BP",ylab="SPD")

for (x in c(1,3,5,7))
{
rect(xleft=breaks[x],xright=breaks[x+1],ybottom=-100,ytop=100,border=NA,col=rgb(0,0,0,0.2))
}

for (x in c(2,4,6,8))
{
rect(xleft=breaks[x],xright=breaks[x+1],ybottom=-100,ytop=100,border=NA,col=rgb(0,0,0,0.05))
}

for (x in 1:c(length(breaks)-2))
{
if (x%%2==0)
{	
arrows(x0=breaks[x]-250,x1=breaks[x+1]-250,y0=0.00066,y1=0.00066,code=3,length=0.05)
text(x=breaks[x+1],y=0.00068,labels=as.roman(x))
} else {
arrows(x0=breaks[x]-250,x1=breaks[x+1]-250,y0=0.0006,y1=0.0006,code=3,length=0.05)
text(x=breaks[x+1],y=0.00062,labels=as.roman(x))
}
}

lines(spd.A$grid$calBP, spd.A$grid$PrDens,col=rgb(0,0,0,0.5),lwd=0.5)
lines(spd.D$grid$calBP, spd.D$grid$PrDens,col=rgb(0,0,0,0.5),lwd=0.5)
lines(spd.B$grid$calBP, spd.B$grid$PrDens,col=rgb(0,0,1,0.5),lwd=0.5)
lines(spd.C$grid$calBP, spd.C$grid$PrDens,col=rgb(1,0,0,0.5),lwd=0.5)

lines(spd.A$grid$calBP,rollmean(spd.A$grid$PrDens,k=200,fill=NA),lwd=2)
lines(spd.D$grid$calBP,rollmean(spd.D$grid$PrDens,k=200,fill=NA),lwd=2)
lines(spd.B$grid$calBP,rollmean(spd.B$grid$PrDens,k=200,fill=NA),lwd=2,col="blue")
lines(spd.C$grid$calBP,rollmean(spd.C$grid$PrDens,k=200,fill=NA),lwd=2,col="red")


plot(ABCDRoca[1,],ylim=range(ABCDRoca),xlab="Transitions",type="n",ylab="Rate of Growth",axes=F,main="c")
lines(ADsum.r)
points(ADsum.r,pch=20)

lines(Bsum.r,col="blue")
points(Bsum.r,pch=20,col="blue")


lines(Csum.r,col="red")
points(Csum.r,pch=20,col="red")

abline(h=0,lty=2,col="darkgreen")

axis(side=1,at=1:7,labels=as.roman(1:7))
axis(side=2)
box()


plot(ABCDRoca[1,],ylim=range(ABCDRoca),xlab="Transitions",type="l",ylab="Rate of Growth",axes=F,main="d")
points(ABCDRoca[1,],pch=20)
lines(ABCDRoca[4,],pch=20)
points(ABCDRoca[4,],,pch=20)
lines(ABCDRoca[2,],pch=20,col="blue")
points(ABCDRoca[2,],pch=20,col="blue")
lines(ABCDRoca[3,],pch=20,col="red")
points(ABCDRoca[3,],pch=20,col="red")
abline(h=0,lty=2,col="darkgreen")
axis(side=1,at=1:7,labels=as.roman(1:7))
axis(side=2)

box()


dev.print(device=pdf,"figure.1.pdf")




### Code for Figure 2 ##


layout(matrix(c(1:5,1,6:9,10:14,10,15:18),4,5,T),width=c(0.3,1,1,1,1))
par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,ylim=c(0,1),xlim=c(0,1))
mtext(side=2,"Full Data",line=-3)


par(mar=c(0,0,1,0))
for (i in 1:c(7))
{
plot(result.locations,index=i,option="test")
if (i==1){
text(10,30,"A",cex=2)
text(30,30,"B",cex=2)
text(30,10,"D",cex=2)
text(10,10,"C",cex=2)
}
mtext(side=3,as.roman(i),line=-0.2)
}
plot(1,type="n",axes=F)

par(mar=c(0,0,0,0))
plot(1,type="n",axes=F,ylim=c(0,1),xlim=c(0,1))
mtext(side=2,"Thinned Data",line=-3)

par(mar=c(0,0,1,0))

for (i in 1:7)
{
plot(resultThinned.locations,index=i,option="test")
mtext(side=3,as.roman(i))
}

plot(1,type="n",axes=F)
par(mar=c(0,0,0,0))
legend("center",title="",legend=c("p<0.05 (Positive Deviation)","q<0.05 (Positive Deviation)","p<0.05 (Negative Deviation)","q<0.05 (Negative Deviation)")
       ,pch=20,col=c("orange","red","cornflowerblue","darkblue"),bg="white",cex=1.2,bty="n")

dev.print(device=pdf,"figure2.pdf")






#########################################
####### Case Study N.2: Analysis ########
#########################################


# Load Euroevol Dataset
data(euroevol)


### Subset radiocarbon dates for the interval 8000 to 5000 Cal BP (c7200-4200 C14BP) ###
rangeStart=8000
rangeEnd=5000
edge=800
yearRange=c(rangeStart,rangeEnd)
euroevol=subset(euroevol,C14Age<=c(rangeStart-edge)&C14Age>=c(rangeEnd-edge))


#Create vector of cutpoints for the chronological blocks
breaksize=500
breaks=seq(rangeStart,rangeEnd,-breaksize) #cutpoints of the chronological blocks


#Create a SpatialPoints class object
SITES = unique(data.frame(SiteID=euroevol$SiteID,Longitude=euroevol$Longitude,Latitude=euroevol$Latitude)) #extrapolate sites
locations=data.frame(Longitude=SITES$Longitude,Latitude=SITES$Latitude)
rownames(locations)=SITES$SiteID
coordinates(locations)<-c("Longitude","Latitude")
proj4string(locations)<- CRS("+proj=longlat +datum=WGS84")


## Compute Distance and Spatial Weights ##
# Compute great-arc distance
distSamples=spDists(locations,locations,longlat = TRUE)
# Compute distance based weights  (using a 100km bandwidth)
spatialweights=spweights(distSamples,h=100)


### Binning and Calibration###
bins=binPrep(sites=euroevol$SiteID,ages=euroevol$C14Age,h=200)  
calDates=calibrate(x=euroevol$C14Age,errors=euroevol$C14SD,timeRange=yearRange,normalised=FALSE)


### Compute SPD and geometric growth ratee ###
spd=spd(x=calDates,timeRange=yearRange,bins=bins,datenormalised=FALSE,spdnormalised=TRUE) #SPD

nBreaks=length(breaks)-1
spd.blocksum=numeric()
spd.roc=numeric()

for (i in 1:nBreaks)
{
	spd.blocksum[i]=sum(spd$grid$PrDens[spd$grid$calBP<=breaks[i]&spd$grid$calBP>=breaks[i+1]])
}

for (i in 1:c(nBreaks-1))
{
spd.roc[i]=(spd.blocksum[i+1]/spd.blocksum[i])^(1/breaksize)-1
}


### Spatial Permutation Test ###

res.locations=SPpermTest(calDates,timeRange=yearRange,bins=bins,
			 locations=locations,spatialweights=spatialweights,
			 breaks=breaks,ncores=3,nsim=10000,permute="locations",datenormalised=FALSE)




#########################################
####### Case Study N.2: Figures  ########
#########################################




## Code for Figure 3 ##



par(mfrow=c(1,2))
plot(spd$grid$calBP, spd$grid$PrDens,col=rgb(0,0,0,0.5),lwd=0.5,type="n",xlim=c(8000,5000),main="a",ylim=c(0,max(spd$grid$PrDens)+0.0002),xlab="cal BP",ylab="SPD")



for (x in c(1,3,5,7))
{
rect(xleft=breaks[x],xright=breaks[x+1],ybottom=-100,ytop=100,border=NA,col=rgb(0,0,0,0.2))
}

for (x in c(2,4,6,8))
{
rect(xleft=breaks[x],xright=breaks[x+1],ybottom=-100,ytop=100,border=NA,col=rgb(0,0,0,0.05))
}

for (x in 1:c(length(breaks)-2))
{
if (x%%2==0)
{	
arrows(x0=breaks[x]-250,x1=breaks[x+1]-250,y0=0.00070,y1=0.00070,code=3,length=0.05)
text(x=breaks[x+1],y=0.00073,labels=as.roman(x))
} else {
arrows(x0=breaks[x]-250,x1=breaks[x+1]-250,y0=0.00067,y1=0.00067,code=3,length=0.05)
text(x=breaks[x+1],y=0.00070,labels=as.roman(x))
}
}

lines(spd$grid$calBP, spd$grid$PrDens,col=rgb(0,0,0,0.5),lwd=0.5)
lines(spd$grid$calBP,rollmean(spd$grid$PrDens,k=200,fill=NA),lwd=2)


plot(spd.roc,ylim=range(-0.0005,0.0018),xlim=c(0.5,5.5),xlab="Transitions",type="l",ylab="Geometric Growth Rate",axes=F,main="b")
points(spd.roc,pch=20)
abline(h=0,lty=2,col="red")
axis(side=1,at=1:5,labels=as.roman(1:5))
axis(side=2)
dev.print(device=pdf,"figure3.pdf")


## Code for Figure 4 ##


library(rworldmap)
base=getMap(resolution="low")
xrange=bbox(res.locations$locations)[1,]
yrange=bbox(res.locations$locations)[2,]


par(mfrow=c(2,3))

for (i in 1:5)
{
par(mar=c(0.1,0.1,0,0.5))
plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange)
plot(res.locations,index=i,add=TRUE,option="raw",breakRange=c(-0.005,0.005))
legend("topleft",legend=c(NA),border=NA,title=as.roman(i),cex=2,bty="n")
}

plot(res.locations,option="rawlegend",breakRange=c(-0.005,0.005),rd=3,legSize=1.6)

dev.print(device=pdf,"figure4.pdf")



## Code for Figure 5 ##



par(mfrow=c(2,3))

for (i in 1:5)
{
par(mar=c(0.1,0.1,0,0.5))	
plot(base,col="antiquewhite3",border="antiquewhite3",xlim=xrange,ylim=yrange)
plot(res.locations,index=i,add=TRUE,option="test")
legend("topleft",legend=c(NA),border=NA,title=as.roman(i),cex=2,bty="n")
}

plot(res.locations,option="testlegend",legSize=2)

dev.print(device=pdf,"figure5.pdf")


