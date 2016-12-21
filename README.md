# Multilocus_eDNA_Analysis
Companion scripts and data for Kelly et al. (Frontiers in Marine Science) multilocus eDNA paper. doi: 10.3389/fmars.2016.00283.

#Notes on data and metadata:
COI and 18s data were sequenced on Illumina MiSeq at Nickerson Lab (Northwest Genomics Center) at University of Washington on March 24-28, 2016.

16s data were a subset of those published in Kelly et al. https://peerj.com/articles/2444/, which had been sequenced in June 2015 on an Illumina NextSeq at Stanford University.

For the COI and 18S data, Jimmy Kralj had prepared 8 libraries in total, pooled into 2 pools (4 libraries/pool, 1 pool/MiSeq run).  Each library contained a mix of 18s and COI amplicons (see sequencing pool metadata sheet).  Each set of amplicons, in turn, had an oligo tag associated with it via second-round PCR. The 16s data had been prepared in the same way, as reflected in the PeerJ publication.

In addition to the Rdata file and scripts using those processed data, we also provide the fasta files for the OTUs analyzed at each locus, the OTU table (sample-by-OTU abundances) and metadata sheets for the relevant sequencing runs. We intend this repository to be more useful than those that Genbank provides, which don't allow this kind of associated data in any obvious way. You would use the metadata sheets attached here to work with the raw fastq files archived in genbank's Sequence Read Archive. 

 
