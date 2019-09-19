## Genome Assembly in R 
## Rodrigo Migueles
## QLS FOundations
## Install packages
install.packages("seqinr")
library(seqinr)
## Read Fasta File ####
G<-read.fasta(file = "C:/Users/Rodrigo.TARANIS/OneDrive - McGill University/Courses/1 - 1F - Foudations of Quantitative Life Sciences I/randomGenome.fasta",as.string = F)# In TARANIS
G<-read.fasta(file = "C:/Users/User/OneDrive - McGill University/Courses/1 - 1F - Foudations of Quantitative Life Sciences I/randomGenome.fasta",as.string = F)
G<-read.fasta(file = "C:/Users/Rodrigo Migueles/OneDrive - McGill University/Courses/1 - 1F - Foudations of Quantitative Life Sciences I/randomGenome.fasta",as.string = F)
## Get Genome Length ####
GL<-length(G[[1]])
GL<-10000 # Genome length
#G<-paste0(G[[1]][1:GL],collapse = "") # Select a subpart of RG and turn this list into a character
NR<-GL/10;NR # Number of reads
RP<-sample(1:GL, NR, replace=T);RP # Read Positions
RLAve<-350 # Read length average
RLVar<-30 # Read length variance
RL<-round(rnorm(NR,RLAve,RLVar), 0);RL # Read lengths
C<-(sum(RL)/GL);C # Coverage
## Alternative: Number of reads depending on the desired coverage ##
#C<-5 # Desired Coverage
#NR<-GL*C/RL;NR # Number of reads
G[[1]][1:7]#[RP[i]:RP[i]+RL[i]]
Reads<-list()# Create list of reads
i<-1
for(i in 1:NR) {
  RE<-RP[i]+RL[i]
  Reads[[i]]<-paste0(G[[1]][RP[i]:RE],collapse = "")
}
# Pick a read 
#Reads[[5]]
ReadIndex<-list()
ReadIndex<-sample(1:length(Reads), 1, replace=F)
Seed<-Reads[ReadIndex]
Seed
# Set pattern length
W<-6 # Choose optimal length based on probability and GL
W<-floor(logb((GL/C),base = 4))
# Read pattern
Pattern<-Seed[[1]][1:W];Pattern
# Scan
gregexpr(Pattern,G)
b <- as.numeric(invisible(gregexpr('gacc',G)[[1]]))
# Here's where I'm stuck: I don't know how to look for 
# the pattern inside the reads. Note that I dont want to look for the 
# pattern in the beginning of the reads only but in the whole body of the reads.

#find(Pattern)
match(x = Pattern,table = GP)
Inc<-array()
R<-3 #Number of times the pattern was found among the reads
# If the number of repeats is less or equal to Coverage, 
# then assign the reads containing a repeat to an 
# Overlapping Reads Sequence
if(R <= C){
  ORS<-c(Read1,Read2) #Overlapping Read Sequences
}
else
{
  PL<-PL+1
}
# Remove overlapping duplicates from assembly to create a single contig
CL<-list() # Contig List
CL[[i]]<-unique(ORS)



