#Assignment 4 - Furkan Guvenc, 999800979

#Part II -------------------------------------------
install.packages("stringr")
install.packages("dplyr")
library(stringr)
library(dplyr)

motifs <- read.table("yeast_tf_motifs.txt", header = F, 
                     col.names = c("Transcriptionfactor", "Motif"))

#reading the fasta file into R

fasta <- readLines("orf_genomic_1000_554.fasta", warn = F)

#Index of yeast gene names
s <- grep(">", fasta)
#creating dataframe listing the range of the length of the sequences of each gene
genes <- data.frame(index = s, from = s+1, to = c((s-1)[-1], length(fasta)))

#constructing the genes into the dataframe
sequences <- rep(NA, length(s))
for(i in 1:length(s)) {
  sequences[i]<-paste(c(str_sub(paste(fasta[genes$from[i]:genes$to[i]], collapse=""),0,1000), 
                        str_sub(paste(fasta[genes$from[i]:genes$to[i]], collapse=""),-1000,-1)), 
                      collapse = "")
}

#putting together the ORF name and 2k sequences
sequenceDatabase <- data.frame(sequenceNames=substr(fasta[s],2,8), 
                               sequence = sequences)

#determining the occurrence of each motif within the genes


#first 1000 bases of the sequences

sequenceDatabase$promoter <- substr(sequenceDatabase$sequence, 0,1000)

#count the incidence of each pattern within the promoter sequence of genes

for (i in 1:nrow(motifs)){
  motifs$totalcount[i]<-sum(str_count(sequenceDatabase$promoter, motifs$Motif[i]))
}

#or counting the number of genes that have at least one match to the sequence motif

for (i in 1:nrow(motifs)) { 
  motifs$count[i] <- sum((grepl(motifs$Motif[i], sequenceDatabase$promoter)), na.rm = T)
}

#Calculate fraction of occurrence of each pattern over the gene list

for (i in 1:nrow(motifs)) {
  motifs$totalcount_fraction[i]<-(motifs$totalcount[i]/nrow(sequenceDatabase)*100)
  motifs$count_fraction[i]<-(motifs$count[i]/nrow(sequenceDatabase)*100)
}



#output motif and their frequency of occurence

print(select(motifs, Motif, totalcount, count, totalcount_fraction, count_fraction))


