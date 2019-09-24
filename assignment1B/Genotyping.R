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

# Set parameters: ####
NR<-length(bamFile[[1]]$pos) # Number of reads
NR<-1000
RP<-bamFile[[1]]$pos[1:NR] # Read positions
RL<-bamFile[[1]]$qwidth[1:NR] # Read lengths
Begin<-RP[1]-1 # Relative position of the beginning of the analysis

# Assign Read Sequences: ####
ReadSequence<-NULL
ReadCoverage<-list()
ReadEnd<-NULL
k<-1
for (k in (1:NR)) {
  ReadSequence[k]<-as.character(bamFile[[1]]$seq[k])  # Extract the sequence
  ReadEnd[k]<-RP[k]+RL[k] # Read end position
  # ReadCoverage[[k]]<-seq((RP[k]-Begin),(ReadEnd[k]-Begin),1)
}

# Read reference genome chr 17 -> G  ####
Ref<-read.fasta(file = "chr17-Copy.fasta")
RefLength<-length(Ref[[1]])

# Split reads calls ####
ReadSplit<-list()
k<-1
for (k in (1:NR)) {
  ReadSplit[[k]]<-strsplit(ReadSequence[k],split='')
  k<-k+1
}


# Chose reads to be considered for each site ####
# Lower<-which(RP>=40999780)
# Higher<-which(ReadEnd<=41000040)
# S<-unique(Lower,Higher)
# ChosenReads<-list()
# k<-1
# h<-40999768
# for (h in (RP[1]:ReadEnd[NR])) {
#   for (k in (1:NR)) {
#     if (h %in% NP[[k]]) {
#       ChosenReads[[h]]<-c(ChosenReads[[h]],k)
#     }
#     k<-k+1
#   }
#   h<-h+1
# }


# Assign absolute nucleotide position ####
NP<-list()
k<-1
for (k in (1:NR)) {
  p<-RP[k]:(RP[k]+249)
  NP[[k]]<-p # Nucleotide position
  k<-k+1
}
# Ignore Error in RP[k]:(RP[k] + 249) : NA/NaN argument

# Assign nucleotide calls
Calls<-list()
k<-1
for (k in (1:NR)) {
  Calls[[k]]<-ReadSplit[[k]] # Nucleotide calls
  k<-k+1
}


# Count nucleotide incidences ####
incidences<-array(data = 0,dim = c(5,(ReadEnd[NR]-Begin),NR))
h<-40999768
k<-1
for (h in (RP[1]:ReadEnd[NR])) { # For all nucleotide positions
  RNP<-(h-Begin) # Assign the relative nucleotide position
  print(RNP)
  for (k in (1:NR)) { # Then for all reads 
    print(k)
    if (h %in% NP[[k]]) { # If the read covers the position, then:
      N<-Calls[[k]][[1]][match(x = h, table = NP[[1]])] # Set the 
      # nucleotide the read calls at that position
      # and count it depending what it is (A,C,G,T):
      if (is.na(N)) {# If there is no nucleotide, count as an exception to go on
        incidences[5,RNP,k]=incidences[1,RNP,k]+1
      }else{
        if (N=="A") {
          incidences[1,RNP,k]=incidences[1,RNP,k]+1
        }
        if (N=="C") {
          incidences[2,RNP,k]=incidences[2,RNP,k]+1
        }
        if (N=="G") {
          incidences[3,RNP,k]=incidences[3,RNP,k]+1
        }
        if (N=="T") {
          incidences[4,RNP,k]=incidences[4,RNP,k]+1
        }
      }
      k<-k+1
    }else{
      k<-k+1
    }
    # k<-k+1
  }
  h<-h+1
}
incidences[,180:200,2]
votes<-sum(incidences[4,45,])


# AL<-seq((40999768-Begin),(40999798-Begin),1)
# AL<-seq((40999768-Begin),(ReadEnd[NR]-Begin),1) # Relative Alignment Length
# h<-RP[1]-Begin


# Build huge alignment

# Create G table
G<-data.frame(matrix(0,n,10))
names(G) = c("AA","AC","AG","AT","CC","CG","CT","GG","GT","TT")
## Define nucleotide position
p<-NULL # Nucleotide position
# Start big for loop that reads one nucleotide position at the time here.

PofD<-data.frame(matrix(0,5,4))
names(PofD) = c("nuc","counts", "probability", "nuc found?")
PofD[,1]<-c("A","C","G","T","NA")
# Fill in the PofD table:
Genotype<-data.frame(matrix(0,(ReadEnd[NR]-Begin),5))
names(Genotype) = c("Homo/Hetero/Uneven/Multiallelic",
                "Best guess prob",
                "3",
                "4",
                "5")
h<-40999768
for (h in (RP[1]:ReadEnd[NR])) {
  RNP<-(h-Begin) # Assign the relative nucleotide position
  print(RNP)
  PofD[,2]<-c(sum(incidences[1,RNP,]),
              sum(incidences[2,RNP,]),
              sum(incidences[3,RNP,]),
              sum(incidences[4,RNP,]),
              sum(incidences[5,RNP,])) 
  # input counts: 
  # count incidences of each nucleotide at position p
  # how to count how many ACGTs from the reads data?
  EffReads<-sum(PofD[,2])#Number of reads that effectively covered p with a non-zero nucleotide incidence
  # Calculate the probability of each nucleotide from the data 
  # in the PofD table and identify undetected nucleotides:
    for (i in 1:4) { # Calculate probabilities
      PofD[i,3]<-PofD[i,2]/EffReads
      if (PofD[i,2]==0) { # Detected nucleotide?
        PofD[i,4]<-F
      }else {
        PofD[i,4]<-T
      }
      i<-i+1
    }
  # Evaluate the probability distribution and determine the genotype:
  if (sum(PofD[,4])==1) {
    # print("Homozygous")
    Genotype[RNP,1]<-1
    Genotype[RNP,2]<-PofD[which.max(PofD[,2]),3]
  } else if (sum(PofD[,4])==2) {
    if (diff(PofD[which(PofD[,2]!=0),3])==0) {
      # print(paste("Heterozygous:",
      #             PofD[which(PofD[,2]!=0),1]))
      Genotype[RNP,1]<-2
      Genotype[RNP,2]<-max(PofD[which(PofD[,2]!=0),3])
      # Alele<-which(PofD[,2]!=0)
      # G[Alele[1],]<-(PofD[which(PofD[,2]!=0),3][1])/4
      # G[,Alele[2]]<-G+(PofD[which(PofD[,2]!=0),3][2])/4
    }else{
      # print(paste("Two options:",
      #             PofD[which(PofD[,2]!=0),1],
      #             "with a probability of",
      #             PofD[which(PofD[,2]!=0),3]))
      Genotype[RNP,1]<-3
      Genotype[RNP,2]<-PofD[which.max(PofD[,2]),3]
      #Threshold strategy should go here
    }
  }else{
    # print(paste("Inconclusive (multiallelic), best guess is",
    #             PofD[which(PofD[,3]==max(PofD[,3])),1],
    #             "with a probability of", max(PofD[,3])))
    Genotype[RNP,1]<-4
    Genotype[RNP,2]<-PofD[which.max(PofD[,2]),3]
  }
  h<-h+1
}
# PofD # See the table, deactivate when in a loop
Genotype[245,]

## Pr(G) ####
# G[2:4,1] = 0
# G[3:4,2] = 0
# G[4,3] = 0
# sum(G)==1
# G
