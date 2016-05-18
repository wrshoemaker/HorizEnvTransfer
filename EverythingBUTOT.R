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


biom_file <- system.file("OTUTable.biom", package = "biom")
x <- read_biom(biom_file)

otufile = "OTUTable.biom"
mapfile = "Mapping_Tadpoles.txt"
trefile = "rep_set.tre"
envir=import_qiime_sample_data(mapfile)
myData = import_biom(otufile, trefile)
myData = merge_phyloseq(myData,envir)

colnames(tax_table(myData)) <- c(k = "Kingdom", p = "Phylum", c = "Class", 
                                 o = "Order", f = "Family", g = "Genus", s = "Species")
ntaxa(myData)

sample_variables(myData)
rank_names(myData)
get_taxa_unique(myData, "Rank6")

# Subset Pseudomonas
GP.chl = subset_taxa(myData, Genus == "g__Pseudomonas")

# using taxa_names(GP.chl), select pseudomonas reads from the fastq, 
# then oligotype

PseudoRow.names <- rownames(otu_table(GP.chl))
Pseduo.names <- taxa_names(GP.chl)
