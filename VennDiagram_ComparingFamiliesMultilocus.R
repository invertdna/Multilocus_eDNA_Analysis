library(VennDiagram); library(wesanderson)

##create combined lookup table for lineages of all molecular data
combinedLookup=rbind(lookup18s, lookupCOI, lookup16s, lookupManual)
	combinedLookup <- unique(combinedLookup)

##look at intersections
totalFams=union(manualFams, AllFamilies_3genes) #note this is derived from rarefied data, so it's a fair comparison across loci
intersect(manualFams, AllFamilies_3genes)
#length(intersect(manualFams, AllFamilies_3genes))/length(manualFams)

#create table of pres/abs for each taxon from each data type
res.df=data.frame(totalFams, rep(0, length(totalFams)), rep(0, length(totalFams)), rep(0, length(totalFams)), rep(0, length(totalFams))) ; names(res.df)=c("Taxon", "locus18s", "locus16s", "locusCOI", "Manual")
	res.df[match(colnames(families18s), totalFams)[!is.na(match(colnames(families18s), totalFams))],2]<-1
	res.df[match(colnames(families16s), totalFams)[!is.na(match(colnames(families16s), totalFams))],3]<-1
	res.df[match(colnames(familiesCOI), totalFams)[!is.na(match(colnames(familiesCOI), totalFams))],4]<-1
	res.df[match(manualFams, totalFams)[!is.na(match(manualFams, totalFams))],5]<-1
row.names(res.df)=res.df[,1]; res.df=res.df[,-1]
res.df=res.df[order(rowSums(res.df), decreasing=T),]
head(res.df, 20)

#summary table; how many data types detect each family?
table(rowSums(res.df))

#add full taxonomy to results
res.df=data.frame(combinedLookup[match(row.names(res.df), combinedLookup$Family),], res.df)

#look at just animals, if desired
animals=res.df[!is.na(res.df$Kingdom)&res.df$Kingdom=="Metazoa",]
animals=animals[order(animals$Phylum,animals$Class,animals$Order),]

fams16s<-colnames(families16s)
fams18s<-colnames(families18s)
famsCOI<-colnames(familiesCOI)

#create venn diagram

# pdf(paste0("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Figures/VennDiagram_",Sys.Date(),".pdf"))
    # # Reference four-set diagram
    # venn.plot <- draw.quad.venn(
    # area1 = length(fams16s),
    # area2 = length(fams18s),
    # area3 = length(famsCOI),
    # area4 = length(manualFams),
    # n12 = length(intersect(fams16s, fams18s)),
    # n13 = length(intersect(fams16s, famsCOI)),
    # n14 = length(intersect(fams16s, manualFams)),
    # n23 = length(intersect(fams18s, famsCOI)),
    # n24 = length(intersect(fams18s, manualFams)),
    # n34 = length(intersect(famsCOI, manualFams)),
    # n123 = length(Reduce(intersect, list(fams16s, fams18s, famsCOI))),
    # n124 = length(Reduce(intersect, list(fams16s, fams18s, manualFams))),
    # n134 = length(Reduce(intersect, list(fams16s, famsCOI, manualFams))),
    # n234 = length(Reduce(intersect, list(fams18s, famsCOI, manualFams))),
    # n1234 = 0,
    # category = c("16S", "18S", "COI", "Manual"),
    # fill = wes_palette("Darjeeling")[1:4],
    # lty = "dashed",
    # cex = 2,
    # cat.cex = 2,
    # cat.col = wes_palette("Darjeeling")[1:4]
    # )
# dev.off()



fams16s<-res.df[res.df$locus16s&res.df$Kingdom=="Metazoa",5]
fams18s<-res.df[res.df$locus18s&res.df$Kingdom=="Metazoa",5]
famsCOI<-res.df[res.df$locusCOI&res.df$Kingdom=="Metazoa",5]
manFams<-res.df[res.df$Manual&res.df$Kingdom=="Metazoa",5]

#create venn diagram

# pdf(paste0("/Users/rpk/GoogleDrive/Kelly_Lab/Projects/MultiLocus_PugetSound/Figures/VennDiagram_ANIMALS_",Sys.Date(),".pdf"))
    # # Reference four-set diagram
    # venn.plot <- draw.quad.venn(
    # area1 = length(fams16s),
    # area2 = length(fams18s),
    # area3 = length(famsCOI),
    # area4 = length(manFams),
    # n12 = length(intersect(fams16s, fams18s)),
    # n13 = length(intersect(fams16s, famsCOI)),
    # n14 = length(intersect(fams16s, manFams)),
    # n23 = length(intersect(fams18s, famsCOI)),
    # n24 = length(intersect(fams18s, manFams)),
    # n34 = length(intersect(famsCOI, manFams)),
    # n123 = length(Reduce(intersect, list(fams16s, fams18s, famsCOI))),
    # n124 = length(Reduce(intersect, list(fams16s, fams18s, manFams))),
    # n134 = length(Reduce(intersect, list(fams16s, famsCOI, manFams))),
    # n234 = length(Reduce(intersect, list(fams18s, famsCOI, manFams))),
    # n1234 = 0,
    # category = c("16S", "18S", "COI", "Manual"),
    # fill = wes_palette("Darjeeling")[1:4],
    # lty = "dashed",
    # cex = 2,
    # cat.cex = 2,
    # cat.col = wes_palette("Darjeeling")[1:4]
    # )
# dev.off()

