rm(list=ls())
getwd()
setwd('~/Desktop/Harris Lab/Tadpole_Sequences/OTUTable_Tadpoles_97uclust/')
getwd()


require("biom")
require("qiimer")
#require("phyloseq")
#source("http://bioconductor.org/biocLite.R")
#biocLite("phyloseq")
require("phyloseq")
require("otu2ot")


otufile <- "OTUTable.biom"
mapfile <- "Mapping_Tadpoles.txt"
trefile <- "rep_set.tre"
envir <- import_qiime_sample_data(mapfile)
myData <- import_biom(otufile, trefile)
myData <- merge_phyloseq(myData,envir)

colnames(tax_table(myData)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                 o = "Order", f = "Family", g = "Genus", s = "Species")
ntaxa(myData)

sample_variables(myData)
rank_names(myData)

# Look at abundances

S.obs <- function(x = ""){
  rowSums(x > 0) * 1
}
OTUtable <- otu_table(myData)
S.obs(OTUtable)

library("vegan")

veganotu <- function(physeq) {
  require("vegan")
  OTU = otu_table(physeq)
  if (taxa_are_rows(OTU)) {
    OTU = t(OTU)
  }
  return(as(OTU, "matrix"))
}

SbyS <- veganotu(myData)
site1 <- SbyS[1,]
dim(site1)
# Observed Richness
Tad.Obs <- S.obs(SbyS)

# Good's coverage
C <- function(x = ""){
  #1 - (sum(x == 1) / rowSums(x))
  #rowSums(x)
  
}
Tad.GC <- C(SbyS)
Tad.GC

# Subset Pseudomonas
GP.chl <- subset_taxa(myData, Genus == "g__Pseudomonas")

# using taxa_names(GP.chl), select pseudomonas reads from the fastq, 
# then oligotype

PseudoRow.names <- rownames(otu_table(GP.chl))
Pseudo.names <- taxa_names(GP.chl)

# Subset the FASTA file for reads with tags that correspond to the OTUs
# That have been assigned Pseudomonas as the taxonomy
library("seqinr")
fastafile<- read.fasta(file = "./extractedseqs250_rep_set_aligned_pfiltered_Reformated.fasta", 
                       seqtype = "DNA",as.string = TRUE, set.attributes = FALSE)

head(fastafile)
fasta.names <- names(fastafile)
split.names <- strsplit(fasta.names, "_")
split.OTUs <- lapply(split.names, `[[`, 3)
test <- fastafile[lapply(strsplit(names(fastafile), "_"), `[[`, 3) %in% Pseudo.names]


# Export file 
write.fasta(sequences = test, nbchar = 60, names = names(fastafile),
            file.out = "./extractedseqs250_rep_set_aligned_pfiltered_Pseudo_only.fasta", open = "w")
# Now test oligotyping
File <- "./extractedseqs250_rep_set_aligned_pfiltered_Pseudo_only.fasta"
Pseudo.OT <- MED(File, minseq = 21, entropymin = 0.6, Plot = T)
head(Pseudo.OT)

GetEnvironmentDatafromFileR(File, Start = 2, Stop = 9, test = T)

ENV <- GetEnvironmentDatafromFileR(File, Start = 2, Stop =9, test = F)

#Table0 <- SampleXOT_Table(
#  Pseudo.OT,
#  ENV = ENV,
#  mosaicPlot = T,
#  filterByMinAbund = 0
#)

#Table10 <- SampleXOT_Table(
#  Pseudo.OT,
#  ENV = ENV,
#  mosaicPlot = T,
#  filterByMinAbund = 10
#)
library("ggplot2")
library("stringr")
library("reshape2")
library("grid")

Pseudo.SbyS <- otu_table(GP.chl)
s_num<-dim(Pseudo.SbyS)[2]
totals<-colSums(Pseudo.SbyS)


s_num<-dim(Pseudo.SbyS)[2]
otu.perc<-matrix(rep("NA",times=(dim(Pseudo.SbyS)[1]*s_num)),nrow=dim(Pseudo.SbyS)[1],ncol=s_num)
rownames(otu.perc)<-rownames(Pseudo.SbyS)
colnames(otu.perc)=colnames(Pseudo.SbyS)
totals<-colSums(Pseudo.SbyS)

top.otus = as.data.frame(head(sort(rowSums(Pseudo.SbyS), decreasing = T),20))


for(s in c(1:s_num)){
  vec<-(Pseudo.SbyS[,s]/totals[s])*100
  otu.perc[,s]<-vec
  rm(vec)
}

write.table(otu.perc, "otu_perc.txt", col.names = NA, row.names = T, sep = "\t")
otu.perc = read.table("otu_perc.txt", sep = "\t", row.names=1, header=T, check.names=F)

otu.20.perc = otu.perc[rownames(top.otus),]
otu.20.melt.perc = melt(cbind(rownames(otu.20.perc),otu.20.perc))

colnames(otu.20.melt.perc)[1] = "OTU"
colnames(otu.20.melt.perc)[2] = "Sample"
colnames(otu.20.melt.perc)[3] = "Abundance"

#Plot it - lots of playing around with theme elements from ggplot possible
pdf("my_plot.pdf", height = 8, width = 8)
otu.plot.data.perc = ggplot(otu.20.melt.perc, aes(x=Sample, y = Abundance, fill = OTU))
otu.barplot.perc = otu.plot.data.perc + geom_bar(stat = "identity", position = "stack")
otu.barplot.perc + facet_grid(OTU ~.) + ggtitle("My Plot") + 
  theme(legend.position = "null", strip.text.y = element_text(size = 10, angle = 0), 
        axis.text.x = element_text(size = 8 , angle = 90), axis.text.y = element_text(size = 6), 
        panel.grid.major = element_line(size = 0.1)) + scale_y_continuous(breaks = c(0,50))
dev.off()


# So let's look at the two dominant OTUs, 400315 and 269930

# Dominant OTUs correspond to only one sequence, so oligotyping was unnecessary

