#Assignment 2 - Furkan Guvenc, 999800979

#Part 1 -------

fileName = "seq.fasta"  #Change this to use your own file.

fileContents = read.table(fileName, sep="\t", header=T, stringsAsFactors=F, check.names=F)

sequence = paste(fileContents[[1]], sep="", collapse="")

cat ("Loaded sequence", names(fileContents),"\n")
cat ("Sequence:\n", sequence, "\n")

compSequence = chartr("ATCG", "TAGC", sequence)


cat ("The complement sequence is:\n", compSequence, "\n\n")

#split the string into an array of characters
sequenceVector = strsplit(compSequence, NULL)[[1]]

aCount = 0;
cCount = 0;
gCount = 0;
tCount = 0;
sequenceLength = 0;
for(base in sequenceVector)
{
  if(base == "A") {
    aCount = aCount+1
  }
  else if(base == "C") {
    cCount = cCount+1
  }
  else if(base == "G") {
    gCount = gCount+1
  }
  else if(base == "T") {
    tCount = tCount+1
  }
  sequenceLength = sequenceLength + 1
}
cat ("Sequence statistics:\n")
cat ("Fraction of A's : ",aCount/sequenceLength,"\n")
cat ("Fraction of T's : ",tCount/sequenceLength,"\n")
cat ("Fraction of G's : ",gCount/sequenceLength,"\n")
cat ("Fraction of C's : ",cCount/sequenceLength,"\n")
cat ("------------------------------\n")

#Output of the code
#Sequence statistics:
# Fraction of A's :  0.239978 
# Fraction of T's :  0.2246019 
# Fraction of G's :  0.2311917 
# Fraction of C's :  0.3042284 


#Explanation of the code: This script aims to intake a particular FASTA sequence, count number
#of nucleic acid bases and calculate sequence statisticcs. The output of the input sequence,
#reverse complement of the sequence as well as percentage of each nucleotide is provided.

#Part1 output


# Loaded sequence >gi|68248531|gb|DQ088379.1| Homo sapiens synuclein, alpha (non A4 component of amyloid precursor) 
#(SNCA) gene, complete cds 
# Sequence:
#   CATCCGAGATAGGGACGAGGAGCACGCTGCAGGGAAAGCAGCGAGCGCCGGGAGAGGGGCGGGCAGAAGCGCTGACAAATCAGCGGTGGGGGCGGAGAGCCGAGGAG
# AAGGAGAAGGAGGAGGACTAGGAGGAGGAGGACGGCGACGACCAGAAGGGGCCCAAGAGAGGGGGCGAGCGACCGAGCGCCGCGACGCGGAAGTGAGGTGCGTGCGGGCTGC
# AGCGCAGACCCCGGCCCGGCCCCTCCGAGAGCGTCCTGGGCGCTCCCTCACGCCTTGCCTTCAAGCCTTCTGCCTTTCCACCCTCGTGAGCGGAGAACTGGGAGTGGCCATT
# CGACGACAGGTTAGCGGGTTTGCCTCCCACTCCCCCAGCCTCGCGTCGCCGGCTCACAGCGGCCTCCTCTGGGGACAGTCCCCCCCGGGTGCCGCCTCCGCCCTTCCTGTGC
# GCTCCTTTTCCTTCTTCTTTCCTATTAAATATTATTTGGGAATTGTTTAAATTTTTTTTTTAAAAAAAGAGAGAGGCGGGGAGGAGTCGGAGTTGTGGAGAAGCAGAGGGAC
# TCAGGTAAGTACCTGTGGATCTAAACGGGCGTCTTTGGAAATCCTGGAGAACGCCGGGTGGGAGACGAATGGTCGTGGGCACCGGGAGGGGGTGGTGCTGCCATGAGGACCC
# GCTGGGCCAGGTCTCTGGGAGGTGAGTACTTGTCCCTTTGGGGAGCCTAAGGAAAGAGACTTGACCTGGCTTTCGTCCTGCTTCTGATATTCCCTTCTCCACAAGGGCTGAG
# AGATTAGGCTGCTTCTCCGGGATCCGCTTTTCCCCGGGAAACGCGAGGATGCTCCATGGAGCGTGAGCATCCAACTTTTCTCTCACATAAAATCTGTCTGCCCGCTCTCTTG
# GTTTTTCTCTGTAAAGTAAGCAAGCTGCGTTTGGCAAATAATGAAATGGAAGTGCAGGGAGGCCAAGTCAACAGGTGGTAACGGGTTAACAAGTGCTGGCGCGGGGTCCGCT
# AGGGTGGAGGCTGAGAACGCCCCCTCGGGTGGCTGGCGCGGGGTTGGAGACGGCCGGCGAGTGTGAGCGGCGCCTGCTCAGGGTAGATAGCTGAGGGCGGGGGTGGATGTTG
# GATGGATTAGAACCATCACACTTGGGCCCGCTGTTTGCCTGAGGTTGAACCACACCCCGAGTGAGTAGTTAGTTCTGTTGCCTACGCCTTTCCACCATCAACCTGTTAGCCT
# TCTTCTGGGATTCATGTTAAGGATACCCCTGACCCTAAGCCTCCAGCTTCCATGCTTCTAACTCATACTGTTACCCTTTAGACCCCGGGAATTTAAAAAAGGGGTTAATCTT
# TTCATGCAACTCCACTTCTGAAATGCAGTAATAACAACTCAGAGGATTCATCCTAATCCGTGGTTAGGTGGCTAGACTTTTACTAGCCAAGATGGATGGGAGATGCTAAATT
# TTTAATGCCAGAGCTAAAAATGTCTGCTTTGTCCAATGGTTAAATGAGTGTACACTTAAAAGAGTCTCACACTTTGGAGGGTTTCTCATGATTTTTCAGTGTTTTTTGTTTA
# TTTTTCCCCGAAAGTTCTCATTCAAAGTGTATTTTATGTTTTCCAGTGTGGTGTAAAGGAATTCATTAGCCATGGATGTATTCATGAAAGGACTTTCAAAGGCCAAGGAGGG
# AGTTGTGGCTGCTGCTGAGAAAACCAAACAGGGTGTGGCAGAAGCAGCAGGAAAGACAAAAGAGGGTGTTCTCTATGTAGGTAGGTAAACCCCAAATGTCAGTTTGGTGCTT
# GTTCATGAGTGATGGGTTAGGATAATCAATACTT 

# The complement sequence is:
#   GTAGGCTCTATCCCTGCTCCTCGTGCGACGTCCCTTTCGTCGCTCGCGGCCCTCTCCCCGCCCGTCTTCGCGACTGTTTAGTCGCCACCCCCGCCTCTCGGCTCCTCTTC
# CTCTTCCTCCTCCTGATCCTCCTCCTCCTGCCGCTGCTGGTCTTCCCCGGGTTCTCTCCCCCGCTCGCTGGCTCGCGGCGCTGCGCCTTCACTCCACGCACGCCCGACGTCGCG
# TCTGGGGCCGGGCCGGGGAGGCTCTCGCAGGACCCGCGAGGGAGTGCGGAACGGAAGTTCGGAAGACGGAAAGGTGGGAGCACTCGCCTCTTGACCCTCACCGGTAAGCTGCTG
# TCCAATCGCCCAAACGGAGGGTGAGGGGGTCGGAGCGCAGCGGCCGAGTGTCGCCGGAGGAGACCCCTGTCAGGGGGGGCCCACGGCGGAGGCGGGAAGGACACGCGAGGAAAA
# GGAAGAAGAAAGGATAATTTATAATAAACCCTTAACAAATTTAAAAAAAAAATTTTTTTCTCTCTCCGCCCCTCCTCAGCCTCAACACCTCTTCGTCTCCCTGAGTCCATTCAT
# GGACACCTAGATTTGCCCGCAGAAACCTTTAGGACCTCTTGCGGCCCACCCTCTGCTTACCAGCACCCGTGGCCCTCCCCCACCACGACGGTACTCCTGGGCGACCCGGTCCAG
# AGACCCTCCACTCATGAACAGGGAAACCCCTCGGATTCCTTTCTCTGAACTGGACCGAAAGCAGGACGAAGACTATAAGGGAAGAGGTGTTCCCGACTCTCTAATCCGACGAAG
# AGGCCCTAGGCGAAAAGGGGCCCTTTGCGCTCCTACGAGGTACCTCGCACTCGTAGGTTGAAAAGAGAGTGTATTTTAGACAGACGGGCGAGAGAACCAAAAAGAGACATTTCA
# TTCGTTCGACGCAAACCGTTTATTACTTTACCTTCACGTCCCTCCGGTTCAGTTGTCCACCATTGCCCAATTGTTCACGACCGCGCCCCAGGCGATCCCACCTCCGACTCTTGC
# GGGGGAGCCCACCGACCGCGCCCCAACCTCTGCCGGCCGCTCACACTCGCCGCGGACGAGTCCCATCTATCGACTCCCGCCCCCACCTACAACCTACCTAATCTTGGTAGTGTG
# AACCCGGGCGACAAACGGACTCCAACTTGGTGTGGGGCTCACTCATCAATCAAGACAACGGATGCGGAAAGGTGGTAGTTGGACAATCGGAAGAAGACCCTAAGTACAATTCCT
# ATGGGGACTGGGATTCGGAGGTCGAAGGTACGAAGATTGAGTATGACAATGGGAAATCTGGGGCCCTTAAATTTTTTCCCCAATTAGAAAAGTACGTTGAGGTGAAGACTTTAC
# GTCATTATTGTTGAGTCTCCTAAGTAGGATTAGGCACCAATCCACCGATCTGAAAATGATCGGTTCTACCTACCCTCTACGATTTAAAAATTACGGTCTCGATTTTTACAGACG
# AAACAGGTTACCAATTTACTCACATGTGAATTTTCTCAGAGTGTGAAACCTCCCAAAGAGTACTAAAAAGTCACAAAAAACAAATAAAAAGGGGCTTTCAAGAGTAAGTTTCAC
# ATAAAATACAAAAGGTCACACCACATTTCCTTAAGTAATCGGTACCTACATAAGTACTTTCCTGAAAGTTTCCGGTTCCTCCCTCAACACCGACGACGACTCTTTTGGTTTGTC
# CCACACCGTCTTCGTCGTCCTTTCTGTTTTCTCCCACAAGAGATACATCCATCCATTTGGGGTTTACAGTCAAACCACGAACAAGTACTCACTACCCAATCCTATTAGTTATGAA 
# 
# Sequence statistics:
# Fraction of A's :  0.239978 
# Fraction of T's :  0.2246019 
# Fraction of G's :  0.2311917 
# Fraction of C's :  0.3042284



#Part 2 --------------------------------------------------------------------------

cat ("this Seems to be ok \n")

sum = 0
sum = sum + 1
sum = sum + 1
sum = sum + 2
sum = sum + 2
cat (sum,"\n")

sequence1 = "AAAA"
sequence2 = "TTTT"
sequenceSum = paste0(sequence1, sequence2)
paste(sequenceSum)

cat ("\n-------------------\n")

email.theU = "myemail.utoronto.ca"
cat (email.theU, "\n")

sequence1 = "What's going on here?"
paste(sequence1)

#Part 3 --------------------------------------------------------------------------

numbers = c(1, 2.5, 3.4, 7.8)

#1. Compute and print the mean and standard deviation of these four numbers
mean(numbers)
sd(numbers)

cat ("\n-------------------\n")

cat (numbers - 5, "\n")
cat (numbers[1:4] - 5)
#2. Explain what how the above two lines are different

#In the first line, the value 5 is subtracted from all of the values in the vector
#in the second line, the specific values of the matrix (1 to 4) are recalled and subtracted 5.


cat ("\n-------------------\n")

sequence = "ATCAAATTCCGCGCATTCTTAGGCCATTCTT"
cat (sequence,"\n")
#3. Mask the subsequence "ATTC" with "xxxx"  and print the new sequence

masked_sequence = chartr("ATTC", "xxxx", sequence)
paste(masked_sequence)
cat ("\n-------------------\n")

cat (sequence,"\n")
#4. Remove the x's from the sequence
removed_x = gsub('x', '', masked_sequence)
removed_x
cat ("\n-------------------\n")

string = "This is a kind of a long string"
#5. Use the substr function to create a string with only the last word this sentence 
#and print it out
#you may try this: what happens when you put a negative number for the start index?
#This works, but there are easier, more conventional ways of solving this problem.

last_word = substr(string, 26, 31)
print(last_word)

cat ("\n-------------------\n")

#Part 4 ---------------------------------------------------------------------------

# Script: translate_sequence.r
# Problem description: read in sequence from file, compute stats, translate RNA into
# 	       amino acid sequence

nuc_sequence = read.table("assignment2_sequence.txt", header = F,
                          sep = "\t", check.names = F, stringsAsFactors = F)

seq1 = paste(nuc_sequence[[1]][3:7], sep = "", collapse = "")

cat("Length of the sequence is:", nchar(seq1), "\n")

sequence = substr(seq1, start = 61, stop = 150)
cat("Coding sequence: \n", sequence, "\n")

aasequence = ""

index=1

while(index < nchar(sequence)){
  
  codon = substr(sequence, index, index+2)
  
  codon = gsub("AUG", "Start", codon)
  codon = gsub("UUC", "Phe", codon)
  codon = gsub("CUU", "Leu", codon)
  codon = gsub("CCA", "Pro", codon)
  codon = gsub("ACC", "Thr", codon)
  codon = gsub("GAU", "Asp", codon)
  codon = gsub("CGA", "Arg", codon)
  codon = gsub("UGA", "Stop", codon)
  
  aasequence = append(aasequence, codon)
  
  index = index + 3
  
}

cat("The amino acid sequence is:", aasequence, "\n")
cat(paste(aasequence, collapse = ""), file = "aaoutput.txt")


