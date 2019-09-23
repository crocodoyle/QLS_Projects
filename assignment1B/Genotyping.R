# Assignment 1 B

# Install packages and require libraries ####
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Rsamtools")
library('Rsamtools')

BiocManager::install("GenomicTools")

library("seqinr")

# Define directories
getwd()
# datadir<-"C:/Users/User/Documents/GitHub/QLS_Projects/assignment1B"
datadir<-"C:/Users/Rodrigo Migueles/Documents/QLS_Projects/assignment1B"
setwd(datadir)

# Read BAM file and define IGV calls -> C and reads -> D ####
bamFile=scanBam('17.41000000-42000000.HG00096.wgs.ILLUMINA.bwa.GBR.high_cov_pcr_free.20140203.bam')
names(bamFile[[1]]) # Labels or names of the columns in the bam
lst <- lapply(names(bamFile[[1]]), function(elt) {
  do.call(c, unname(lapply(bamFile, "[[", elt)))
})
names(lst) <- names(bamFile[[1]])
bamFile_df=do.call("DataFrame", lst)

NR<-length(bamFile[[1]]$pos) # Number of reads
RP<-bamFile[[1]]$pos # Read positions
RL<-bamFile[[1]]$qwidth # Read length
ARP<-RP+RL #Absolute read position
ReadSequence<-NULL # Assign Read Sequences:
k<-1
for (k in (1:NR)) {
  ReadSequence[k]<-as.character(bamFile[[1]]$seq[k])  
}

# Read reference genome chr 17 -> G  ####
Ref<-read.fasta(file = "chr17-Copy.fasta")
RefLength<-length(Ref[[1]])

# Build huge alignment
for (h in seq(41000000,42000000,1)) {
  
}


# Define parameters
n<-20# Number of reads that cover a given position

# Create G table
G<-data.frame(matrix(0,n,10))
names(G) = c("AA","AC","AG","AT","CC","CG","CT","GG","GT","TT")
## Define nucleotide position
p<-NULL # Nucleotide position
# Start big for loop that reads one nucleotide position at the time here.
PofD<-data.frame(matrix(0,4,4))
names(PofD) = c("nuc","counts", "probability", "nuc found?")
PofD[,1]<-c("A","C","G","T")
# Fill in the PofD table:
for (j in 1:n) {
  PofD[,2]<-c(7,7,0,0) # input counts: 
  # count incidences of each nucleotide at position p
  # how to count how many ACGTs from the reads data?
  j<-j+1
}
EffReads<-sum(PofD[,"counts"])#Number of reads that effectively covered p with a non-zero nucleotide incidence
n==EffReads # Make sure all reads provide a nucleotide
# Calculate the probability of each nucleotide from the data 
# in the PofD table and identify undetected nucleotides:
for (i in 1:4) {
  PofD[i,3]<-PofD[i,2]/EffReads
  if (PofD[i,2]==0) {
    PofD[i,4]<-F
  }else {
    PofD[i,4]<-T
  }
  i<-i+1
}
PofD # See the table, deactivate when in a loop
# Evaluate the probability distribution and determine the genotype:
if (sum(PofD[,4])==1) {
  print("Homozygous")
} else if (sum(PofD[,4])==2) {
    if (diff(PofD[which(PofD[,2]!=0),3])==0) {
    print(paste("Heterozygous:",
                PofD[which(PofD[,2]!=0),1]))
      # Alele<-which(PofD[,2]!=0)
      # G[Alele[1],]<-(PofD[which(PofD[,2]!=0),3][1])/4
      # G[,Alele[2]]<-G+(PofD[which(PofD[,2]!=0),3][2])/4
    }else{
      print(paste("Two options:",
                  PofD[which(PofD[,2]!=0),1],
                  "with a probability of",
                  PofD[which(PofD[,2]!=0),3]))
      #Threshold strategy should go here
        }
  }else{
    print(paste("Inconclusive (multiallelic), best guess is",
                PofD[which(PofD[,3]==max(PofD[,3])),1],
                "with a probability of", max(PofD[,3])))
}

## Pr(G) ####
# G[2:4,1] = 0
# G[3:4,2] = 0
# G[4,3] = 0
sum(G)==1
G
