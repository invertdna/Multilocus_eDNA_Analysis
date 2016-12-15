#Analysis script for Multilocus manuscript; Kelly et al.; Frontiers in Marine Science; largely written June 2016

#README: this is the main analysis script used for the multilocus analysis in Kelly et al. (Frontiers in Marine Science, submitted 2016). It loads the datafile ThreeGenes_Rarefied.rData, and calls supplementary scripts, which depend on objects created by running this main analysis script. So run this one first. Remember to enter a path to a working directory for this script, on line 23.



##########LOAD LIBRARIES and FUNCTIONS#############
libs=c("vegan","MASS","DESeq2","data.table", "gdata","lattice","plyr","dplyr", "gridExtra", "ggplot2", "taxize", "lmtest", "wesanderson", "xtable")
	lapply(libs, require, character.only = TRUE)  #"gsubfn"

bit=function(x){print(x[c(1:10),c(1:10)]); print(dim(x))}  #shows top-left corner of a data.frame/matrix; useful for large datasets

ProbOcc=function(x, psi, p11, p10, K){
	(psi*(p11^x)*(1-p11)^(K-x)) / ((psi*(p11^x)*(1-p11)^(K-x))+(((1-psi)*(p10^x))*((1-p10)^(K-x))))
	}   #from Lahoz-Monfort, Mol Ecol Res 2015

gg_colors <- function(n) {
	  hues = seq(15, 375, length=n+1)
	  hcl(h=hues, l=65, c=100)[1:n]
	}

###################################################
setwd("") #set working directory; insert relevant path here.

####Load data and libraries
load("ThreeGenes_Rarefied.rData")  #load rarefied data.

###Source taxon comparison; also used to create VENN DIAGRAM (diagram written to Figures folder) if part of script is uncommented; can also do subset for just animals, etc
source("Analysis/VennDiagram_ComparingFamiliesMultilocus.R")  #note that for all data types, greatest number of taxa are unique to that type




###Accumulation Curve, if desired.  Will create a pdf in the Figures folder.  Uses 1000 draws from 1000 rarefaction replicates to create accumulation curve of Families with each locus/data type
#source("AccumulationCurve_1000rarefactionReplicates.R")

###Binomial Simulations.
#source("Analysis/BinomialSimulation_Kelly_et_al_2016.R")
	#will run simulations and write PDF to file


#write table res.df for supplementary info
# res.df[is.na(res.df)]<-"-"
# colnames(res.df)[6:9]<-c("18S","16S","COI","Manual")
# tab.suppl<- xtable(res.df, auto=T, type="latex")
	# align(tab.suppl)<-c("|c","|c","|c","|c","|c","|c","|c","|c","|c","|c|")
	# print(tab.suppl, file="Figures/table_suppl.tex", include.rownames=F, scalebox='0.75')


###RESULTS section 1
#Our tow-net surveys recovered 49 distinct taxa, of which 36 could be annotated to Family level. 
length(names(envir)[15:59]) #total taxa
length(grep("idae",names(envir)[15:59])) #annotated to Family

#eDNA surveys with PCR primers targeting three gene regions—16S, COI, and 18S—detected a total of 1.8x104 unique operational taxonomic units (OTUs), of which a total of 1.38x104 could be classified by matching known sequences in Genbank (Table 1).
read.summary<-data.frame(
c(	dim(rareData16s)[1], #total OTUs
	sum(rowSums(!is.na(master16s[master16s$OTU_name_swarm%in%row.names(rareData16s),5:9]))>0), #N OTUs annotated to any taxonomic level
	sum(!is.na(master16s[master16s$OTU_name_swarm%in%row.names(rareData16s),"OTU_taxon_family"])),	#N OTUs annotated to Family
	length(unique(master16s[master16s$OTU_name_swarm%in%row.names(rareData16s),"OTU_taxon_family"]))-1 #N unique Families)
	),

c(	dim(rareDataCOI)[1], #total OTUs
	sum(rowSums(!is.na(masterCOI[masterCOI$OTU_name_swarm%in%row.names(rareDataCOI),5:9]))>0), #N OTUs annotated to any taxonomic level
	sum(!is.na(masterCOI[masterCOI$OTU_name_swarm%in%row.names(rareDataCOI),"OTU_taxon_family"])),	#N OTUs annotated to Family
	length(unique(masterCOI[masterCOI$OTU_name_swarm%in%row.names(rareDataCOI),"OTU_taxon_family"]))-1 #N unique Families)
	),

c(	dim(rareData18s)[1], #total OTUs
	sum(rowSums(!is.na(master18s[master18s$OTU_name_swarm%in%row.names(rareData18s),5:9]))>0), #N OTUs annotated to any taxonomic level
	sum(!is.na(master18s[master18s$OTU_name_swarm%in%row.names(rareData18s),"OTU_taxon_family"])),	#N OTUs annotated to Family
	length(unique(master18s[master18s$OTU_name_swarm%in%row.names(rareData18s),"OTU_taxon_family"]))-1 #N unique Families)
	)
) ; colnames(read.summary)<-c("16s", "COI", "18s") ; row.names(read.summary)<-c("Total OTUs", "Annotated OTUs", "OTUs Annotated to Family", "Unique Families")

#tab1<- xtable(read.summary, auto=T, type="latex")
#	align(tab1)<-c("|c","|c","|c","|c|")
#	print(tab1, file="Figures/table1.tex")
		#never got this to work; need to paste latex preamble
	# fileConn<-file("Figures/table1.tex")	
		# writeLines("\\documentclass[11pt,letterpaper]{article}\n\\begin{document}\n", fileConn)
		# print(tab1, file=fileConn, append=T)
		# write("\n\\end{document}", fileConn, append=T)
	# close(fileConn)


#Of those that could be annotated, the samples represented a total 366 unique taxonomic Families across the Eukarya, including representatives from 33 Phyla from 9 Kingdoms. 

dim(res.df)[1]  #N unique Families
length(unique(res.df$Phylum)) #N unique Phyla
length(unique(res.df$Kingdom)) #N unique Kingdoms

#The identities of taxa detected varied dramatically across genes, and between genetic vs. manual methods of sampling (Figure 1)
#[Venn diagram, generated above]

#No Family was detected by all four data types, while [10] Families were detected by 3 data types, and [80] were detected by two data types. Hence, many of the Families detected by a given data type were unique to that data type (i.e.,  [333] Families detected only by a single data type), underscoring the differences in suites of taxa sampled by the different methods. 

	table(rowSums(res.df[,6:9]))
	head(res.df)

#	Each sampling method showed a distinct phylogenetic bias, consistent with the idea that primer-site mismatches are likely to be driving the observed patterns (Figure 2).
#OTUsByKingdomByDataType.pdf, generated below

#Because each data type reveals a nearly non-overlapping set of taxa, the accumulation curves reflect the idea that adding data types dramatically increases taxonomic coverage (Fig 3).
#[accumulation curve, generated above]


#biases within primer sets; proportion of lineages detected belonging to a particular group (e.g., Kingdom)
data.frame(names(table(lookup18s$Kingdom)), table(lookup18s$Kingdom)/length(lookup18s$Kingdom)) 

	###18s
		#unique lineages recovered by 18s in rarefied dataset:
		subLookup18s<-unique(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),5:9])	
		#table of OTUs by Kingdom
		m18s<-table(master18s$OTU_taxon_kingdom[match(row.names(rareData18s), master18s$OTU_name_swarm)]); m18s <- m18s/sum(m18s)
		#table of unique lineages (to Family) by Kingdom
		m<-table(subLookup18s$OTU_taxon_kingdom); m/sum(m)

	###COI
		#unique lineages recovered by COI in rarefied dataset:
		subLookupCOI<-unique(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),5:9])	
		#table of OTUs by Kingdom
		mCOI<-table(masterCOI$OTU_taxon_kingdom[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm)]); mCOI <- mCOI/sum(mCOI)
		#table of unique lineages (to Family) by Kingdom
		m<-table(subLookupCOI$OTU_taxon_kingdom); m/sum(m)

	###16s
		#unique lineages recovered by 16s in rarefied dataset:
		subLookup16s<-unique(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),5:9])	
		#table of OTUs by Kingdom
		m16s<-table(master16s$OTU_taxon_kingdom[match(row.names(rareData16s), master16s$OTU_name_swarm)]); m16s <- m16s/sum(m16s)
		#table of unique lineages (to Family) by Kingdom
		m<-table(subLookup16s$OTU_taxon_kingdom); m/sum(m)

KingdomsComparison<-merge(merge(merge(data.frame(m18s),data.frame(mCOI), by="Var1", all=T),data.frame(m16s), by="Var1", all=T), data.frame(table(lookupManual$Kingdom)/length(lookupManual$Kingdom)),by="Var1", all=T); row.names(KingdomsComparison)<-KingdomsComparison[,1]; KingdomsComparison= KingdomsComparison[,-1]; names(KingdomsComparison)<-c("18s", "COI", "16s", "Manual")

#18s and COI get similar complements of taxa, looking at the level of Kingdom...
# pdf("Figures/OTUsByKingdomByDataType.pdf", width=7)
# par(mfrow=c(4,1)); par(mar=c(2,3,1,1))
# for (i in 1:4){barplot(KingdomsComparison[,i], col=wes_palette("Darjeeling")[i], ylim=c(0,1)); 
	# text(.5, .9, names(KingdomsComparison)[i], font=2, cex=1.3)
	# }
	# axis(1, at=seq(0.5,7.5,1)*1.22, labels=row.names(KingdomsComparison), las=2, tick=F, line=-5.5, cex.axis=0.9)
# dev.off()
	
plot(KingdomsComparison[,1:2], main="Kingdoms\n(proportion of OTUs detected)"); text(KingdomsComparison[,1], KingdomsComparison[,2]-.01, row.names(KingdomsComparison[,1:2]), cex=.5); curve(1*x, add=T, lty=2)


#the counts of sequence reads look similar to OTUs recovered
reads18s_kingdom<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_kingdom"]), sum)
readsCOI_kingdom<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_kingdom"]), sum)
readsKingdomsComparisons<-merge(reads18s_kingdom, readsCOI_kingdom, by="Group.1", all=T); names(readsKingdomsComparisons)<-c("Kingdom", "18s", "COI")
plot(readsKingdomsComparisons[,2:3]); text(readsKingdomsComparisons[,2], readsKingdomsComparisons[,3]-5000, readsKingdomsComparisons[,1], cex=.5)
plot(log(readsKingdomsComparisons[,2:3]))


#looking at the level of Order, 16s continues to be quite different from 18s/COI
reads16s_order<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_order"]), sum)
reads18s_order<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_order"]), sum)
readsCOI_order<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_order"]), sum)
readsOrderComparisons<-merge(merge(reads18s_order, readsCOI_order, by="Group.1", all=T), reads16s_order, by="Group.1", all=T); names(readsOrderComparisons)<-c("Order", "18s", "COI", "16s")


OTUs16s_order<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_order"]), specnumber)
OTUsCOI_order<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_order"]), specnumber)
OTUs18s_order<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_order"]), specnumber)
OTUsOrderComparisons<-merge(merge(OTUs18s_order, OTUsCOI_order, by="Group.1", all=T), OTUs16s_order, by="Group.1", all=T); names(OTUsOrderComparisons)<-c("Order", "18s", "COI", "16s"); OTUsOrderComparisons[is.na(OTUsOrderComparisons)]<-0

#looking at the level of family, 16s continues to be quite different from 18s/COI
reads16s_family<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_family"]), sum)
reads18s_family<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_family"]), sum)
readsCOI_family<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_family"]), sum)
readsfamilyComparisons<-merge(merge(reads18s_family, readsCOI_family, by="Group.1", all=T), reads16s_family, by="Group.1", all=T); names(readsfamilyComparisons)<-c("family", "18s", "COI", "16s")


OTUs16s_family<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_family"]), specnumber)
OTUsCOI_family<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_family"]), specnumber)
OTUs18s_family<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_family"]), specnumber)
OTUsfamilyComparisons<-merge(merge(OTUs18s_family, OTUsCOI_family, by="Group.1", all=T), OTUs16s_family, by="Group.1", all=T); names(OTUsfamilyComparisons)<-c("family", "18s", "COI", "16s") ; OTUsfamilyComparisons[is.na(OTUsfamilyComparisons)]<-0


#Where a taxon was detected by more than one eDNA marker, we asked whether the different eDNA markers yielded consistent estimates of diversity within the taxon. The only pair of markers that yielded consistent results among the taxa shared was COI-18S. Families detected by both 18S and COI had correlated numbers of OTUs.

df=OTUsfamilyComparisons[specnumber(OTUsfamilyComparisons[,2:3])==2,2:3] #18s/COI
	df=log(df) ; dim(df)
	summary(lm(df[,1]~df[,2])) #linear regression
	cor.test(df[,1],df[,2], method="spearman")
	#negative binomial regression
	summary(glm.nb(df[,1]~df[,2]))
	
df=OTUsfamilyComparisons[specnumber(OTUsfamilyComparisons[,c(2,4)])==2,c(2,4)] #18s/16s
	df=log(df) ; dim(df)
	summary(lm(df[,1]~df[,2]))
	cor.test(df[,1],df[,2], method="spearman")
	
df=OTUsfamilyComparisons[specnumber(OTUsfamilyComparisons[,c(3,4)])==2,c(3,4)] #COI/16s
	df=log(df) ; dim(df)
	summary(lm(df[,1]~df[,2]))
	cor.test(df[,1],df[,2], method="spearman")
	
	
#	Between-marker correlations in read-counts-per-taxon mirrored the OTUs-per-Family results, with the 18S and COI primer sets reflecting correlated numbers of reads for the Phyla, Classes, Orders, and Families detected at both loci.
 

	#Phylum
reads16s_phylum<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_phylum"]), sum)
reads18s_phylum<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_phylum"]), sum)
readsCOI_phylum<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_phylum"]), sum)
readsphylumComparisons<-merge(merge(reads18s_phylum, readsCOI_phylum, by="Group.1", all=T), reads16s_phylum, by="Group.1", all=T); names(readsphylumComparisons)<-c("phylum", "18s", "COI", "16s")
pairs(log(readsphylumComparisons[,2:4]))
summary(lm(log(readsphylumComparisons[,2])~ log(readsphylumComparisons[,3]))) #18s/COI
summary(lm(log(readsphylumComparisons[,2])~ log(readsphylumComparisons[,4]))) #18s/16s
summary(lm(log(readsphylumComparisons[,3])~ log(readsphylumComparisons[,4]))) #COI/16s

cor.test(log(readsphylumComparisons[,2]), log(readsphylumComparisons[,3]), method="spearman")
cor.test(log(readsphylumComparisons[,2]), log(readsphylumComparisons[,4]), method="spearman")
cor.test(log(readsphylumComparisons[,3]), log(readsphylumComparisons[,4]), method="spearman")


	#class
reads16s_class<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_class"]), sum)
reads18s_class<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_class"]), sum)
readsCOI_class<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_class"]), sum)
readsclassComparisons<-merge(merge(reads18s_class, readsCOI_class, by="Group.1", all=T), reads16s_class, by="Group.1", all=T); names(readsclassComparisons)<-c("class", "18s", "COI", "16s")
pairs(log(readsclassComparisons[,2:4]))
summary(lm(log(readsclassComparisons[,2])~ log(readsclassComparisons[,3]))) #18s/COI
summary(lm(log(readsclassComparisons[,2])~ log(readsclassComparisons[,4]))) #18s/16s
summary(lm(log(readsclassComparisons[,3])~ log(readsclassComparisons[,4]))) #COI/16s

cor.test(log(readsclassComparisons[,2]), log(readsclassComparisons[,3]), method="spearman")
cor.test(log(readsclassComparisons[,2]), log(readsclassComparisons[,4]), method="spearman")
cor.test(log(readsclassComparisons[,3]), log(readsclassComparisons[,4]), method="spearman")



	#order
reads16s_order<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_order"]), sum)
reads18s_order<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_order"]), sum)
readsCOI_order<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_order"]), sum)
readsorderComparisons<-merge(merge(reads18s_order, readsCOI_order, by="Group.1", all=T), reads16s_order, by="Group.1", all=T); names(readsorderComparisons)<-c("order", "18s", "COI", "16s")
pairs(log(readsorderComparisons[,2:4]))
summary(lm(log(readsorderComparisons[,2])~ log(readsorderComparisons[,3]))) #18s/COI
summary(lm(log(readsorderComparisons[,2])~ log(readsorderComparisons[,4]))) #18s/16s
summary(lm(log(readsorderComparisons[,3])~ log(readsorderComparisons[,4]))) #COI/16s

cor.test(log(readsorderComparisons[,2]), log(readsorderComparisons[,3]), method="spearman")
cor.test(log(readsorderComparisons[,2]), log(readsorderComparisons[,4]), method="spearman")
cor.test(log(readsorderComparisons[,3]), log(readsorderComparisons[,4]), method="spearman")



	#family
reads16s_family<-aggregate(rowSums(rareData16s), list(master16s[match(row.names(rareData16s), master16s$OTU_name_swarm),"OTU_taxon_family"]), sum)
reads18s_family<-aggregate(rowSums(rareData18s), list(master18s[match(row.names(rareData18s), master18s$OTU_name_swarm),"OTU_taxon_family"]), sum)
readsCOI_family<-aggregate(rowSums(rareDataCOI), list(masterCOI[match(row.names(rareDataCOI), masterCOI$OTU_name_swarm),"OTU_taxon_family"]), sum)
readsfamilyComparisons<-merge(merge(reads18s_family, readsCOI_family, by="Group.1", all=T), reads16s_family, by="Group.1", all=T); names(readsfamilyComparisons)<-c("family", "18s", "COI", "16s")
pairs(log(readsfamilyComparisons[,2:4]))
summary(lm(log(readsfamilyComparisons[,2])~ log(readsfamilyComparisons[,3]))) #18s/COI
summary(lm(log(readsfamilyComparisons[,2])~ log(readsfamilyComparisons[,4]))) #18s/16s
summary(lm(log(readsfamilyComparisons[,3])~ log(readsfamilyComparisons[,4]))) #COI/16s

cor.test(log(readsfamilyComparisons[,2]), log(readsfamilyComparisons[,3]), method="spearman")
cor.test(log(readsfamilyComparisons[,2]), log(readsfamilyComparisons[,4]), method="spearman")
cor.test(log(readsfamilyComparisons[,3]), log(readsfamilyComparisons[,4]), method="spearman")


####PERMANOVA (with raw, unrarefied data)
#	We used a permutation-based analysis of variance (PERMANOVA) to apportion the variance in Jaccard distances among our geographic sites (N = 4), among our biological samples within sites (3/ geographic sites) and among PCR replicates within each of those biological samples (3 or 4 PCRs / 

sites=substr(colnames(OTUsDecontam16s),1,3)
replicates=substr(colnames(OTUsDecontam16s),1,5)
presAbs16s<-t(OTUsDecontam16s) ; presAbs16s[presAbs16s>0]<-1
			distmat=vegdist(presAbs16s, method="jaccard")
				var16s<-adonis(distmat~sites+replicates); var16s



sites=substr(colnames(OTUsDecontam18s),1,3)
replicates=substr(colnames(OTUsDecontam18s),1,5)
presAbs18s<-t(OTUsDecontam18s) ; presAbs18s[presAbs18s>0]<-1
			distmat=vegdist(presAbs18s, method="jaccard")
			var18s<-adonis(distmat~sites+replicates); var18s



sites=substr(colnames(OTUsDecontamCOI),1,3)
replicates=substr(colnames(OTUsDecontamCOI),1,5)
presAbsCOI<-t(OTUsDecontamCOI) ; presAbsCOI[presAbsCOI>0]<-1
			distmat=vegdist(presAbsCOI, method="jaccard")
				varCOI<-adonis(distmat~sites+replicates); varCOI



sites= envir$Site.code
presAbsManual<-envir[,15:59] ; presAbsManual[presAbsManual>0]<-1
			distmat=vegdist(presAbsManual, method="jaccard")
				varManual<-adonis(distmat~sites); varManual
				
##create table of variance				
varTable<-data.frame(round(var16s$aov.tab[1:3,5],3),round(varCOI$aov.tab[1:3,5],3),round(var18s$aov.tab[1:3,5],3),c(round(varManual$aov.tab[1:2,5],3), "-"))
	row.names(varTable)<-c("Among Sites", "Among Water Samples\nWithin Sites", "Among PCR Replicates")
	colnames(varTable)<-c("16s","COI","18s","Manual")

#Table 2
# tab2<- xtable(varTable, auto=T, type="latex")
	# align(tab2)<-c("|c","|c","|c","|c","|c|")
	# outfile<-"table2.tex"
	# print(tab2, file=outfile)


#PCA plots for supplement

pdf("SupplFigure1_PCAplots.pdf")
par(mfrow=c(2,2))
 pca.16s<-princomp(vegdist(presAbs16s, method="jaccard"))
 	sites=substr(colnames(OTUsDecontam16s),1,3)
			ordiplot(pca.16s, display="sites", main="16s")
				ordihull(pca.16s,groups= sites,draw="polygon",col="grey90",label=F)
				orditorp(pca.16s,display="sites", labels=sites)
 pca.COI<-princomp(vegdist(presAbsCOI, method="jaccard"))
 	sites=substr(colnames(OTUsDecontamCOI),1,3)
			ordiplot(pca.COI, display="sites", main="COI")
				ordihull(pca.COI,groups= sites,draw="polygon",col="grey90",label=F)
				orditorp(pca.COI,display="sites", labels=sites)
 pca.18s<-princomp(vegdist(presAbs18s, method="jaccard"))
 	sites=substr(colnames(OTUsDecontam18s),1,3)
			ordiplot(pca.18s, display="sites", main="18s")
				ordihull(pca.18s,groups= sites,draw="polygon",col="grey90",label=F)
				orditorp(pca.18s,display="sites", labels=sites)
 pca.Manual<-princomp(vegdist(presAbsManual, method="jaccard"))
 	sites=substr(envir[,1],1,3)
			ordiplot(pca.Manual, display="sites", main="Manual")
				ordihull(pca.Manual,groups= sites,draw="polygon",col="grey90",label=F)
				orditorp(pca.Manual,display="sites", labels=sites)
dev.off()




