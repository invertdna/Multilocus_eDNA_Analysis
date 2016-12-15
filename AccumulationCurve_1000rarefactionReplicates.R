######accumulation curve for Multilocus analysis, using 1000 rarefaction draws created/stored in /Data/Rarefaction_1000_18s_families.RData, etc.
library(wesanderson); library(vegan)

load("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Data/ThreeGenes_Rarefied.rData")
load("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Data/Rarefaction_1000_18s_families_4sites.RData") #incorporates 16s and COI as well, because of the way those files were saved in series.

#how variable are rarefaction replicates?

Nfams=NA; for(i in 1:1000){Nfams[i]<-specnumber(colSums(rareDataList16s[[i]]))}; sd(Nfams)/mean(Nfams)
Nfams=NA; for(i in 1:1000){Nfams[i]<-specnumber(colSums(rareDataList18s[[i]]))}; sd(Nfams)/mean(Nfams)
Nfams=NA; for(i in 1:1000){Nfams[i]<-specnumber(colSums(rareDataListCOI[[i]]))}; sd(Nfams)/mean(Nfams)


#create manual count dataset
famsManual_Sites<-aggregate(familiesManual, list(substr(row.names(familiesManual),1,3)), sum); row.names(famsManual_Sites)<-famsManual_Sites[,1]; famsManual_Sites<-famsManual_Sites[,-1]
	famsManual_Sites[famsManual_Sites>0]<-1 


resultsList=list() #create empty list; each element of which will be a dataframe summarizing the Families found for each data type in each rarefaction draw
NSites=1:nrow(Site_total_16s)  #here, we have a max of 4 geographic sites we can sample
NsamplesPerSite=1000 #how many different draws do you want to sample?  each will draw a different set of N sites, and each will contain a different draw from the pool of 1000 rarefaction replicates
index=1
for (n in NSites){
	for (i in 1: NsamplesPerSite){
		focalSites=sample(unique(row.names(Site_total_16s)),n) #choose your set of n focal sample sites at random
					#for rarefaction draw 'i' (rareDataList16s is 1000 different rarefaction draws from the overall 16s dataset) from set of sites 'n', calculate how many unique families occur within the set of focal sites.  Create a dataframe in which to store the results of Families detected at those sites for each of the four different data types
if(n==1)
{temp16s<-aggregate(rareDataList16s[[i]][row.names(Site_total_16s)%in%focalSites,], by=list(master16s[which(master16s $OTU_name_swarm%in%colnames(rareDataList16s[[i]])),"OTU_taxon_family"]), sum)
	names(temp16s)<-c("Family", "Count")
temp18s<-aggregate(rareDataList18s[[i]][row.names(Site_total_18s)%in%focalSites,], by=list(master18s[which(master18s $OTU_name_swarm%in%colnames(rareDataList18s[[i]])),"OTU_taxon_family"]), sum)
	names(temp18s)<-c("Family", "Count")
tempCOI<-aggregate(rareDataListCOI[[i]][row.names(Site_total_COI)%in%focalSites,], by=list(masterCOI[which(masterCOI $OTU_name_swarm%in%colnames(rareDataListCOI[[i]])),"OTU_taxon_family"]), sum)
	names(tempCOI)<-c("Family", "Count")
tempManual<-aggregate(t(famsManual_Sites[row.names(famsManual_Sites)%in%focalSites,]), list(colnames(famsManual_Sites)), sum)
	names(tempManual)<-c("Family", "Count")
	} else 
{temp16s<-aggregate(t(rareDataList16s[[i]][row.names(Site_total_16s)%in%focalSites,]), list(master16s[which(master16s $OTU_name_swarm%in%colnames(rareDataList16s[[i]])),"OTU_taxon_family"]), sum)
	temp16s<-data.frame(temp16s[,1], rowSums(temp16s[,2:ncol(temp16s)])); names(temp16s)<-c("Family", "Count")
temp18s<-aggregate(t(rareDataList18s[[i]][row.names(Site_total_18s)%in%focalSites,]), list(master18s[which(master18s $OTU_name_swarm%in%colnames(rareDataList18s[[i]])),"OTU_taxon_family"]), sum)
	temp18s<-data.frame(temp18s[,1], rowSums(temp18s[,2:ncol(temp18s)])); names(temp18s)<-c("Family", "Count")
tempCOI<-aggregate(t(rareDataListCOI[[i]][row.names(Site_total_COI)%in%focalSites,]), list(masterCOI[which(masterCOI $OTU_name_swarm%in%colnames(rareDataListCOI[[i]])),"OTU_taxon_family"]), sum)
	tempCOI<-data.frame(tempCOI[,1], rowSums(tempCOI[,2:ncol(tempCOI)])); names(tempCOI)<-c("Family", "Count")
tempManual<-aggregate(t(famsManual_Sites[row.names(famsManual_Sites)%in%focalSites,]), list(colnames(famsManual_Sites)), sum)
	tempManual <-data.frame(tempManual[,1], rowSums(tempManual[,2:ncol(tempManual)])); names(tempManual)<-c("Family", "Count")
	}
	
suppressWarnings(allData<-merge(merge(merge(temp16s, temp18s, by="Family", all=T), tempCOI, by="Family", all=T), tempManual, by="Family", all=T))
row.names(allData)<-allData[,1]; allData<-allData[,-1]; names(allData)<-c("temp16s", "temp18s", "tempCOI", "tempManual")
	allData[is.na(allData)]<-0; allData[allData>0]<-1 
resultsList[[index]]<-allData

		index=index+1
	}
}

FamilyRichness16s<-NA; for(i in 1:(max(NSites)*NsamplesPerSite)){FamilyRichness16s[i]<-sum(resultsList[[i]][,"temp16s"])}
FamilyRichness18s<-NA; for(i in 1:(max(NSites)*NsamplesPerSite)){FamilyRichness18s[i]<-sum(resultsList[[i]][,"temp18s"])}
FamilyRichnessCOI<-NA; for(i in 1:(max(NSites)*NsamplesPerSite)){FamilyRichnessCOI[i]<-sum(resultsList[[i]][,"tempCOI"])}
FamilyRichnessManual<-NA; for(i in 1:(max(NSites)*NsamplesPerSite)){FamilyRichnessManual[i]<-sum(resultsList[[i]][,"tempManual"])}
FamilyRichnessAllLoci<-NA; for(i in 1:(max(NSites)*NsamplesPerSite)){FamilyRichnessAllLoci[i]<-specnumber(rowSums(resultsList[[i]][,1:3]))}
FamilyRichnessAllData<-NA; for(i in 1:(max(NSites)*NsamplesPerSite)){FamilyRichnessAllData[i]<-specnumber(rowSums(resultsList[[i]][,1:4]))}



pdf(paste0("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Figures/AccumulationCurve_Families_",Sys.Date(),".pdf"), width=8, height=8)
shift=0.05
			ydata= FamilyRichnessAllData; xdata= c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites))))-2*shift
	boxplot(ydata ~ xdata, ylim=c(15,410), xlim=c(0.6,4.5),col ="dodgerblue", xlab="Sampled Sites", ylab="Unique Taxonomic Families", boxwex=.2, xaxt='n', at=(1:max(NSites))-2*shift)
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col="dodgerblue", lwd=2)

			ydata= FamilyRichnessAllLoci; xdata= c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites))))-shift
	boxplot(FamilyRichnessAllLoci ~ c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites)))), col =wes_palette("Darjeeling")[5], add=T, boxwex=.2, xaxt='n', at=(1:max(NSites))-shift)
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col=wes_palette("Darjeeling")[5], lwd=2)

			ydata= FamilyRichness18s; xdata= c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites))))
	boxplot(FamilyRichness18s ~ c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites)))), col =wes_palette("Darjeeling")[2], add=T, boxwex=.2, xaxt='n', at=(1:max(NSites)))
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=1, b=1))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col=wes_palette("Darjeeling")[2], lwd=2)

			ydata= FamilyRichnessCOI; xdata= c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites))))+shift
	boxplot(FamilyRichnessCOI ~ c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites)))), col =wes_palette("Darjeeling")[3], add=T, boxwex=.2, xaxt='n', at=(1:max(NSites))+shift)
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col=wes_palette("Darjeeling")[3], lwd=2)

			ydata= FamilyRichness16s; xdata= c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites))))+2*shift
	boxplot(FamilyRichness16s ~ c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites)))) , col =wes_palette("Darjeeling")[1], add=T, boxwex=.2, xaxt='n', at=(1:max(NSites))+2*shift)
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col=wes_palette("Darjeeling")[1], lwd=2)

			ydata= FamilyRichnessManual; xdata= c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites))))+3*shift		
	boxplot(FamilyRichnessManual ~ c(rep(1:max(NSites), rep(NsamplesPerSite,max(NSites)))), col =wes_palette("Darjeeling")[4], add=T, boxwex=.2, xaxt='n', at=(1:max(NSites))+3*shift)
		nls.mod=nls(ydata~a*log(xdata)+b, start=list(a=600, b=300))
			a=coefficients(nls.mod)[1] ; b=coefficients(nls.mod)[2]
		curve(a*log(x)+b, xlim=c(1,8), add=T, col=wes_palette("Darjeeling")[4], lwd=2)
	
	axis(1, 1:max(NSites))
	
	legend("topleft", fill=c("dodgerblue", wes_palette("Darjeeling")[c(5,1,2,3,4)]), legend=c("All data","3 loci","16s", "18s", "COI", "Manual Survey"), cex=0.8)
dev.off()	
	

