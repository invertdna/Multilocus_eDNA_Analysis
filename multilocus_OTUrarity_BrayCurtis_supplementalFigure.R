#plot mean among-PCR-replicate Bray-Curtis distances as a function of OTU rarity for each of 3 loci

####16s#####
		binNo<-cut(1:nrow(rareData16s), 10, labels=F)  #10 equal-sized groups
		MeanBC<-list()
		replicates=substr(names(rareData16s), 1,5) #create vector of replicate names
		for (j in 1:length(unique(binNo))){
			distmat=vegdist(t(rareData16s[binNo==j,]), upper=T, diag=T) #calculate bray-curtis distances
			distList=list(NA) ; distList.tri=list(NA); index=1  #create empty lists to fill in, and then create list of matrices, each of which is distances among replicates
				index=1
				for (i in unique(replicates)){
					rowMatch<-replicates%in%i
					distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
							distList.tri[[index]]<-	distList[[index]][upper.tri(	distList[[index]])]
					index=index+1
				}
			MeanBC[[j]]<-unlist(lapply(distList, mean))
		}
		xdata16s<-rep(1:length(unique(binNo)), each=length(out16s)/length(unique(binNo)))
		out16s<-unlist(MeanBC)
		#dim(out16s)<-c(length(out16s)/length(unique(binNo)),length(unique(binNo))) #reshape

#####COI
		binNo<-cut(1:nrow(rareDataCOI), 10, labels=F)  #10 equal-sized groups
		MeanBC<-list()
		replicates=substr(names(rareDataCOI), 1,5) #create vector of replicate names
		for (j in 1:length(unique(binNo))){
			distmat=vegdist(t(rareDataCOI[binNo==j,]), upper=T, diag=T) #calculate bray-curtis distances
			distList=list(NA) ; distList.tri=list(NA); index=1  #create empty lists to fill in, and then create list of matrices, each of which is distances among replicates
				index=1
				for (i in unique(replicates)){
					rowMatch<-replicates%in%i
					distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
							distList.tri[[index]]<-	distList[[index]][upper.tri(	distList[[index]])]
					index=index+1
				}
			MeanBC[[j]]<-unlist(lapply(distList, mean))
		}
		xdataCOI<-rep(1:length(unique(binNo)), each=length(outCOI)/length(unique(binNo)))
		outCOI<-unlist(MeanBC)
		#dim(outCOI)<-c(length(outCOI)/length(unique(binNo)),length(unique(binNo))) #reshape

####18s
		binNo<-cut(1:nrow(rareData18s), 10, labels=F)  #10 equal-sized groups
		MeanBC<-list()
		replicates=substr(names(rareData18s), 1,5) #create vector of replicate names
		for (j in 1:length(unique(binNo))){
			distmat=vegdist(t(rareData18s[binNo==j,]), upper=T, diag=T) #calculate bray-curtis distances
			distList=list(NA) ; distList.tri=list(NA); index=1  #create empty lists to fill in, and then create list of matrices, each of which is distances among replicates
				index=1
				for (i in unique(replicates)){
					rowMatch<-replicates%in%i
					distList[[index]]<-as.matrix(distmat)[rowMatch, rowMatch]
							distList.tri[[index]]<-	distList[[index]][upper.tri(	distList[[index]])]
					index=index+1
				}
			MeanBC[[j]]<-unlist(lapply(distList, mean))
		}
		xdata18s<-rep(1:length(unique(binNo)), each=length(out18s)/length(unique(binNo)))
		out18s<-unlist(MeanBC)
		#dim(out18s)<-c(length(out18s)/length(unique(binNo)),length(unique(binNo))) #reshape


#####PLOTTING####
pdf("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Manuscript/v2_Frontiers_revisions_Dec2016/SupplementalFigure_BrayCurtis_OTUrarity.pdf", width=9, height=7)
par(mar=c(5,4,2,2))
boxplot(out16s~xdata16s, ylim=c(0,0.8), xaxt="n", ylab="Bray-Curtis Distances Among Technical Replicates", col =wes_palette("Darjeeling")[1], boxwex=.2)
	#axis(1, 1:10)
	axis(1, 1:10, line=-.5, labels=c("Top 10%\nofOTUs", "10-20%", "20-30%", "30-40%", "40-50%", "50-60%", "60-70%", "70-80%", "80-90%", "Bottom 10%\nofOTUs"), lty=0, las=2, cex.axis=0.8)
	
boxplot(outCOI~xdataCOI, add=T, col =wes_palette("Darjeeling")[3], boxwex=.2, at= unique(xdataCOI)+.2, xaxt="n")
boxplot(out18s~xdata18s, add=T, col =wes_palette("Darjeeling")[2], boxwex=.2, at= unique(xdataCOI)+.4, xaxt="n")
legend(.1,0.835,legend=c("16s","COI","18s"), fill=wes_palette("Darjeeling")[c(1,3,2)])
dev.off()


#mtext("")




