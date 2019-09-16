## Genome Assembly in R 
## Rodrigo Migueles
## QLS FOundations
## Install packages
install.packages("seqinr")
library(seqinr)
## Read Fasta File ####
#RG<-read.fasta(file = "C:/Users/Rodrigo.TARANIS/Documents/randomGenome.fasta",as.string = F)# In TARANIS
RG<-read.fasta(file = "C:/Users/User/OneDrive - McGill University/Courses/1 - 1F - Foudations of Quantitative Life Sciences I/randomGenome.fasta",as.string = F)
## Get Genome Length ####
GL<-length(RG[[1]])
GL<-100 # Genome length
G<-RG[[1]][1:GL] # Select a subpart of RG 
GP<-G[1]
i<-1

for (i in 1:(GL-1)){
  GP<-paste(GP,G[i+1],sep = "")
  i<-i+1
  #print(GP)
}
# and turn this list into a character
NR<-GL/10;NR # Number of reads
RP<-sample(1:GL, NR, replace=T);RP # Read Positions
RLAve<-15 # Read length average
RLVar<-2 # Read length variance
RL<-round(rnorm(NR,RLAve,RLVar), 0);RL # Read lengths
C<-(sum(RL)/GL);C # Coverage
## Alternative: Number of reads depending on the desired coverage ##
#C<-5 # Desired Coverage
#NR<-GL*C/RL;NR # Number of reads
Reads<-list()# Create list of reads
for(i in 1:NR) {
  Reads[[i]]<-G[RP[i]:(RP[i]+RL[i])]
}
# Pick a read 
Reads[[1]][1:2]
Seed<-Reads[sample(1:length(Reads), 1, replace=F)]
Seed
# Set pattern length
PL<-6 # Choose optimal length based on probability and GL
# Read pattern
Pattern<-Seed[[1]][1:PL];Pattern
# Scan
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
