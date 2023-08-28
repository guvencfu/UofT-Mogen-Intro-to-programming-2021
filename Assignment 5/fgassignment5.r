#Assignment 5 - Furkan Guvenc, 999800979

install.packages("dplyr")
install.packages("stringr")
install.packages("ggplot2")

library(dplyr)
library(stringr)
library(ggplot2)

#reading FASTA file, using function written in previous assignment

fasta <- readLines("yeastChromosome16Genes-DNA.fasta", warn = F)

read_fasta = function(fasta){
  #Index of yeast gene names
  s <- grep(">", fasta)
  #creating dataframe listing the range of the length of the sequences of each gene
  genes <- data.frame(index = s, from = s+1, to = c((s-1)[-1], length(fasta)))
  
  #constructing the genes into the dataframe
  sequences <- rep(NA, length(s))
  for(i in 1:length(s)) {
    sequences[i]<-paste(fasta[genes$from[i]:genes$to[i]], collapse="")
  }
  
  #putting together the ORF name and 2k sequences
  sequenceDatabase <- na.omit(data.frame(SEQNAME=str_extract(fasta[s], "Y[A-P](R|L)[0-9]{3}(W|C)"), 
                                 SEQUENCE = sequences))
}

seqs <- read_fasta(fasta)

#the fasta file contains gene names that do not fall under the yeast ORF convention.
#I have omitted these to design primers specifically for the yeast ORFs


#creating the primer3 input
primer3_input <- list()
for (i in 1:nrow(seqs)){
  primer3_input[[i]] <- c(paste0("PRIMER_SEQUENCE_ID=",seqs$SEQNAME[i]),
                        paste0("SEQUENCE=",toupper(seqs$SEQUENCE[i])),
                        "PRIMER_FILE_FLAG=0",
                        "PRIMER_EXPLAIN_FLAG=1",
                        "=")
}

proper_primer3_input <- unlist(primer3_input)


#create the primer3 input file
cat(paste(proper_primer3_input, collapse = "\n"), file="primer_input.txt")

#generate primer3 output file
system2("primer3-1.1.4-WINXP/bin/primer3_core.exe",
        stdin="primer_input.txt",
        stdout="primer_output.txt",
        stderr=T)



#primer3 output file read in R
output_primer3<-scan("primer_output.txt", character(), sep = "\n")

#remove the lines from the primer3 output file that did not give results

output_primer3 <- output_primer3[-(49896:49909)]


#best primer pairs (the first one provided by primer3)
{
gene_name<- c(gsub("^.*=", "", (str_subset(output_primer3, 
                                               "Y[A-P](R|L)[0-9]{3}(W|C)"))))
primer_left<-c(gsub("^.*=", "", (str_subset(output_primer3, 
                                                "PRIMER_LEFT_SEQUENCE="))))
primer_right<-c(gsub("^.*=", "", (str_subset(output_primer3, 
                                                 "PRIMER_RIGHT_SEQUENCE="))))
df_primers <- data.frame(ORF = gene_name, left= primer_left, right = primer_right)

head(df_primers)

}


#Building df containing the primer characteristics
{
tm <- as.numeric(c(gsub("^.*=", "", (str_subset(output_primer3, "PRIMER_LEFT_TM=")))))
gc <- as.numeric(c(gsub("^.*=", "", (str_subset(output_primer3, "PRIMER_LEFT_GC_PERCENT=")))))
left_length <- c(nchar(gsub("^.*=", "", (str_subset(output_primer3, "PRIMER_LEFT_SEQUENCE=")))))
df_characteristics <- data.frame(ORF = gene_name, tm = tm, gc = gc, length = left_length)
write.table(df_characteristics, "all_primers_yChr16.txt")
}

#graphical outputs

{
  #mean primer length, tm, %gc
  mean_length <- mean(df_characteristics$length)
  mean_tm <- mean(df_characteristics$tm)
  mean_gc <- mean(df_characteristics$gc)
  sd_length <- sd(df_characteristics$length)
  sd_tm <- sd(df_characteristics$tm)
  sd_gc <- sd(df_characteristics$gc)
  stats <- matrix(c(mean_gc, sd_gc, mean_tm, sd_tm, mean_length, sd_length), ncol = 3)
  colnames(stats) <- c("GC", "Tm", "Primer_length")
  rownames(stats)<- c("Mean", "St.Dev")
  stats
}

{
  #plotting the tm vs gc data
  cor_Test <- cor.test(df_characteristics$gc, df_characteristics$tm, method = "pearson", conf.level = 0.95)
  cor.val <- round(cor_Test$estimate, 3)
  p.val <- round(cor_Test$p.value, 5)
  cor.label <- paste0("R = ", cor.val)
  p.label <- paste0("p = ", p.val)
  ggplot(df_characteristics, aes(gc, tm)) +
  geom_point() +
  annotate(x = 59, y = 59.5,  geom = "text", 
           label = cor.label, size = 4) +
    annotate(x = 59, y = 59.4,  geom = "text", 
             label = p.label, size = 4)
}


{
  #histogram for the TM values
  ggplot(df_characteristics, aes(x=tm))+
    geom_histogram()
}

#boxplot of TMs with GC <50% or >50%


{
greater_50<-filter(df_characteristics, gc>50)
greater_50$gc <- ">50"
less_50<-filter(df_characteristics, gc<50)
less_50$gc <- "<50"
df_boxplot<- union(less_50, greater_50)

#plotting the boxplot
significance <- t.test(greater_50$tm, less_50$tm)
p.value <- round(significance$p.value, 5)
p.label <- paste0("p = ", p.value)
ggplot(df_boxplot, aes(gc, tm)) + 
  geom_boxplot()+
  annotate(x = 1.5, y = 60.75, geom = "text",
           label = p.label, size = 4)
}
