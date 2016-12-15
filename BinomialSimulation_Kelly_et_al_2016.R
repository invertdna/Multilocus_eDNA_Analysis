#5-locus simulation; eDNA binomial occupancy probability assuming variable detection rate (per-locus, per-species) and a false positive detection rate of 0. See Kelly et al. for details.
library(vegan); library(grid); library(Hmisc) #load libraries
set.seed(95) #for reproducibility
setwd("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Figures") #folder for writing figure

Nsimulations=1000 #number of simulations to carry out
Nspecies = 100 #number of species present in hypothetical community. This is the theoretical maximum number of species to be detected.
Nreplicates = 5 #max number of biological replicates (e.g., simulated individual water samples) in simulated survey

#Assign detection probabilities
p11_A<-rbeta(Nspecies, 1,5)  #assign p11_A (i.e., the probability of true detection using locus A) for each of Nspecies species, such that most species have a low probability of detection, and about a third of them have a zero probability of dection with locus A
p11_B<-p11_B<-rbeta(Nspecies, 1.5,10)  #the distribution of probs is skewed left, but less so
p11_C<-rbeta(Nspecies, 1.5,5) #about 20% of them have a zero probability of dection with locus C, but most species are still below .2 probability
p11_D<-rbeta(Nspecies, .01,.04) #more than 70% of them have a zero probability of dection with locus D, several are nearly 100%
p11_E<-rbeta(Nspecies, .2,10) #about 90% of them have a zero probability of dection with locus E, and a very few have nearly 100% detection probability



####Multilocus Case (figure left panel)
#######################################
singleLocus.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol= Nreplicates))
	for (j in 1: Nreplicates){
		for (i in 1:Nsimulations){singleLocus.df[i,j]<-specnumber(data.frame(rbinom(Nspecies,j,p11_A),rbinom(Nspecies,j,p11_B),rbinom(Nspecies,j,p11_C),rbinom(Nspecies,j,p11_D),rbinom(Nspecies,j,p11_E))[,sample(1:5,1)])
			}	
		}
twoLocus.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
	for (j in 1:Nreplicates){
		for (i in 1:Nsimulations){twoLocus.df[i,j]<-specnumber(rowSums(data.frame(rbinom(Nspecies,j,p11_A),rbinom(Nspecies,j,p11_B),rbinom(Nspecies,j,p11_C),rbinom(Nspecies,j,p11_D),rbinom(Nspecies,j,p11_E))[,sample(1:5,2)]))
			}	
		}
threeLocus.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
for (j in 1:Nreplicates){
	for (i in 1:Nsimulations){threeLocus.df[i,j]<-specnumber(rowSums(data.frame(rbinom(Nspecies,j,p11_A),rbinom(Nspecies,j,p11_B),rbinom(Nspecies,j,p11_C),rbinom(Nspecies,j,p11_D),rbinom(Nspecies,j,p11_E))[,sample(1:5,3)]))
		}	
	}
fourLocus.df <-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
for (j in 1:Nreplicates){
	for (i in 1:Nsimulations){fourLocus.df[i,j]<-specnumber(rowSums(data.frame(rbinom(Nspecies,j,p11_A),rbinom(Nspecies,j,p11_B),rbinom(Nspecies,j,p11_C),rbinom(Nspecies,j,p11_D),rbinom(Nspecies,j,p11_E))[,sample(1:5,4)]))
		}	
	}
fiveLocus.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
for (j in 1:Nreplicates){
	for (i in 1:Nsimulations){fiveLocus.df[i,j]<-specnumber(rowSums(data.frame(rbinom(Nspecies,j,p11_A),rbinom(Nspecies,j,p11_B),rbinom(Nspecies,j,p11_C),rbinom(Nspecies,j,p11_D),rbinom(Nspecies,j,p11_E))))
		}	
	}


####Single Case (figure right panel)
#######################################
locusA.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
	for (j in 1:Nreplicates){
		for (i in 1:Nsimulations){locusA.df[i,j]<-specnumber(rbinom(Nspecies,j,p11_A))
			}
	}
locusB.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
	for (j in 1:Nreplicates){
		for (i in 1:Nsimulations){locusB.df[i,j]<-specnumber(rbinom(Nspecies,j,p11_B))
			}
	}
locusC.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
	for (j in 1:Nreplicates){
		for (i in 1:Nsimulations){locusC.df[i,j]<-specnumber(rbinom(Nspecies,j,p11_C))
			}
	}
locusD.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
	for (j in 1:Nreplicates){
		for (i in 1:Nsimulations){locusD.df[i,j]<-specnumber(rbinom(Nspecies,j,p11_D))
			}
	}
locusE.df<-as.data.frame(matrix(NA, nrow=Nsimulations, ncol=Nreplicates))
	for (j in 1:Nreplicates){
		for (i in 1:Nsimulations){locusE.df[i,j]<-specnumber(rbinom(Nspecies,j,p11_E))
			}
	}


##PLOTTING
####PLOTTING
pdf("BinomSimulations.pdf",width=12, height=5)
par(mfrow=c(1,2))
par(mar=c(3,4,1,0))

##multi-locus case
shift=.07
boxplot(data.frame(fiveLocus.df), ylim=c(0,100), xaxt='n', xlab="Replicate Samples", ylab="Species Detected", col="grey20",at=c(1:5)-2*shift, boxwex=.2, staplewex=.2, outline=F, range=0)
	axis(1,1:5,labels=1:5)
			ydata= apply(data.frame(fiveLocus.df),2, median)
			xdata= c(1:5)-2*shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey20", lwd=2)

boxplot(data.frame(fourLocus.df), add=T, col='grey40', xaxt='n', at=c(1:5)-shift, boxwex=.2, staplewex=.2, outline=F, range=0)
			ydata= apply(data.frame(fourLocus.df),2, median)
			xdata= c(1:5)-shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey40", lwd=2)

boxplot(data.frame(threeLocus.df), add=T, col='grey60', xaxt='n', at=c(1:5), boxwex=.2, staplewex=.2, outline=F, range=0)
			ydata= apply(data.frame(threeLocus.df),2, median)
			xdata= c(1:5)
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey60", lwd=2)

boxplot(data.frame(twoLocus.df), add=T, col='grey80', xaxt='n', at=c(1:5)+shift, boxwex=.2, staplewex=.2, outline=F, range=0)
			ydata= apply(data.frame(twoLocus.df),2, median)
			xdata= c(1:5)+shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey80", lwd=2)

boxplot(data.frame(singleLocus.df), add=T, col='grey100', xaxt='n', at=c(1:5)+2*shift, boxwex=.2, staplewex=.2, outline=F, range=0)
			ydata= apply(data.frame(singleLocus.df),2, median)
			xdata= c(1:5)+2*shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey90", lwd=2)

legend(0.15,105,c('5','4', '3', '2', '1'), fill=c("grey20","grey40","grey60","grey80","grey100"), title="N Loci", bty='n')
text(.35, 0, "A", cex=2)

par(mar=c(3,1,1,1))
##single-locus case
shift=0.1
boxplot(data.frame(locusA.df), ylim=c(0,100), xaxt='n', yaxt="n", ylab="", xlab="Replicate Samples", col="grey20", at=c(1:5)-2*shift ,boxwex=.2, staplewex=.2, outline=F, range=0)
	axis(1,1:5,labels=1:5)
			ydata= apply(data.frame(locusA.df),2, median)
			xdata= c(1:5)-2*shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey20", lwd=2)

boxplot(data.frame(locusB.df), add=T, boxwex=.2, staplewex=.2, at=c(1:5)-shift, xaxt='n', yaxt="n", col="grey40", outline=F, range=0)
			ydata= apply(data.frame(locusB.df),2, median)
			xdata= c(1:5)-shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey40", lwd=2)
			
boxplot(data.frame(locusC.df), add=T,boxwex=.2, staplewex=.2, at=c(1:5), xaxt='n', yaxt="n", col="grey60", outline=F, range=0)
			ydata= apply(data.frame(locusC.df),2, median)
			xdata= c(1:5)
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey60", lwd=2)
			
boxplot(data.frame(locusD.df), add=T,boxwex=.2, staplewex=.2, at=c(1:5)+shift, yaxt="n", xaxt='n', col="grey80", outline=F)
			ydata= apply(data.frame(locusD.df),2, median)
			xdata= c(1:5)+shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey80", lwd=2)
			
boxplot(data.frame(locusE.df), add=T,boxwex=.2, staplewex=.2, at=c(1:5)+2*shift, yaxt="n", xaxt='n', col="grey100", outline=F, range=0)
			ydata= apply(data.frame(locusE.df),2, median)
			xdata= c(1:5)+2*shift
			nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
				a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
			curve(a*log(x)+b, xlim=c(.5,8), add=T, col="grey90", lwd=2)
			
legend(0,105,c('Locus A','Locus B', 'Locus C', 'Locus D', 'Locus E'), fill=c("grey20","grey40","grey60","grey80","grey100"), bty='n', y.intersp=1.2)
text(.3, 0, "B", cex=2)
mtext("Replicate Samples", side=1, line=2, outer=F, at=c(-0.08))

#add distributions of p11 values, if desired
 subplot(hist(p11_A,xaxt='n', yaxt='n', main="", ylab="", xlab=""), 1.4,99, size=c(.5,.2))
 subplot(hist(p11_B,xaxt='n', yaxt='n', main="", ylab="", xlab=""), 1.4,93, size=c(.5,.2))
 subplot(hist(p11_C,xaxt='n', yaxt='n', main="", ylab="", xlab=""), 1.4,87, size=c(.5,.2))
 subplot(hist(p11_D,xaxt='n', yaxt='n', main="", ylab="", xlab=""), 1.4,81, size=c(.5,.2))
 subplot(hist(p11_E,xaxt='n', yaxt='n', main="", ylab="", xlab=""), 1.4,75, size=c(.5,.2))
dev.off()


##Effect of going from 1 to 2 loci, single replicate
summary(singleLocus.df[,1]) #one locus, one replicate
summary(twoLocus.df[,1]) #two loci, one replicate

##Effect of going from 1 to 2 replicates, single locus  #Depends strongly on locus.
(median(locusC.df[,2])-median(locusC.df[,1]))/median(locusC.df[,1]) #best-case 68% increase
(median(locusD.df[,2])-median(locusD.df[,1]))/median(locusD.df[,1]) #worst-case 0% increase
