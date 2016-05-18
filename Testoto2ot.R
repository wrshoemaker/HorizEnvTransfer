rm(list=ls())
getwd()
setwd('~/github/HorizEnvTransfer/')
getwd()
#install.packages("seqinr")
#install.packages("otu2ot_1.4.tar.gz", repos = NULL, type = "source")

File = "./Tutorial_otu2ot/SupplementaryMaterial/HGB_0013_GXJPMPL01A3OQX.fasta"

OT.seq.concat <- MED(File, minseq = 21, entropymin = 0.6, Plot = T)

head(OTs.seq.concat)

# OR we could import the sequences and then apply MED
OT.seq.concat <- MEDMat(AlignedSequences = Sequences, minseq=21, entropymin = 0.6, Plot = T)

# Retrieve sample information from the FASTA header
# Test length using 
GetEnvironmentDatafromFileR(File, Start = 2, Stop = 9, test = T)

ENV <- GetEnvironmentDatafromFileR(File, Start = 2, Stop =9, test = F)

Table0 <- SampleXOT_Table(
  OT.seq.concat,
  ENV = ENV,
  mosaicPlot = T,
  filterByMinAbund = 0
  )

Table0 <- SampleXOT_Table(
  OT.seq.concat,
  ENV = ENV,
  mosaicPlot = T,
  filterByMinAbund = 10
)
