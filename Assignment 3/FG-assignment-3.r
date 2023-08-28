#Assignment 3 - Furkan Guvenc, 999800979

#Part1 -----------------------------------------------------------------------

#Yeast genome analysis script

db <- read.table("yeast_orflist.txt", header = T, sep = "\t", comment.char = "")


#using R built-in functions to return shortest and longest genes

#Longest
print(db[which.max(db$Length),])
#Shortest
print(db[which.min(db$Length),])


#using a for loop to determine longest and shortest gene


for (i in 1:nrow(db)) {
  if (db[i, "Length"] == max(db$Length)) {
  print(db[i,])
  }  
  else if (db[i, "Length"] == min(db$Length)) {
  print(db[i,])
  }  
}

#Shortest gene YJR151W-A is a putative protein of unknown function
#Longest gene YLR106C is an ATPase involved in ATP-dependent remodelling of the pre-60S particle, export from the nucleus.

#Output Watson and Crick strand genes from the database

x=1

watson <- c()
crick <- c()


while (x <= nrow(db)){
  if (grepl("W", db$ORF[x]) == T) {
  watson = append(watson, db$ORF[x])
  }  
  else if (grepl("C", db$ORF[x]) == T) {
  crick = append(crick, db$ORF[x])
  }
  x = x + 1
}

write.table(watson, "yeast_watsonstrand_genes.txt", quote = F, sep = "\t")
write.table(crick, "yeast_crickstrand_genes.txt", quote = F, sep = "\t")

#mean ORF length

print(mean(db$Length))

#Part2------------------------------------------------------------------------------

#A) Open the phenotype file for reading

phen_assignment <- read.delim("phenotype_data.tab", header = T, sep = "\t", as.is = T)
interactions <- read.table("yeast_interactions.txt", header = T, sep = "\t", comment.char = "")


total<- length(phen_assignment$PMID)
print(total)

unique_pmid <- length(unique(phen_assignment$PMID))
print(unique_pmid)

#140638 total PMID entries
#5880 unique PMID entries

#B) Creating a list for the inviable genes

inviable_genes = c()
viable_genes = c()
n = 1

while (n <= nrow(phen_assignment)){
  if (phen_assignment[n, "mitochondrial.genome.maintenance..abnormal"] == "inviable") {
  inviable_genes <- unique(append(inviable_genes, phen_assignment$IMI1[n]))
  }
  n = n + 1
}

length(inviable_genes)


#C) Number of interactors for a given protein


interactors <- as.data.frame((table(interactions$INTERACTOR_A)))
interactors

#D) Number of interactors for genes that are listed as inviable

inviable_interactions <- interactors[match(inviable_genes_list$inviable_genes, interactors$Var1, nomatch = 0), ]

#Average number of interactors of genes annotated a "inviable"

mean(inviable_interactions$Freq)

#Average number of interactors of all genes in the interaction list

mean(interactors$Freq)

#on average, the genes that have been annotated as inviable (essential) (~32.9) have a greater
#number of interactions compared to all the genes combined (~23.5). 

#Part3 ------------------------------------------------------------------------------------------

go_mapping <- as.data.frame(read.table("go_slim_mapping.tab", header = T, sep = "\t"))

bud <- c()
x <- 1

while (x <= nrow(go_mapping)){
  if (go_mapping[x, "cellular_component"] == "cellular bud"){
    bud <- append(bud, go_mapping$HRA1[x])
  }
  x <- x + 1
}
  
print(bud)   


