#Assignment 6 - Furkan Guvenc, 99980979

#Part 1: loading and manipulating data

{
expression_data <- read.table("assignment6_data.txt", header = T)
cancer_labels <- read.table("assignment6_sample_labels.txt")


cancer_labels$V1 <- gsub("1", "non-metastatic", cancer_labels$V1)
cancer_labels$V1 <- gsub("2", "metastatic", cancer_labels$V1)

colnames(expression_data) <- cancer_labels$V1
cancer_data_matrix <- data.matrix(expression_data)
}


#print information on BRCA1

{
brca1 <- cancer_data_matrix["BRCA1",]
brca1
}
#mean & stdev of BRCA1 expression profile across all samples

{
brca1 <- as.array(brca1)

brca1_stats <- c(mean(brca1), sd(brca1))
brca1_cols <- c("mean", "stdev")
brca1_data <- as.data.frame(brca1_stats, brca1_cols)
brca1_data
}
#mean and stdev of all the genes in the first metastatic group

{
first_metast_col<- cancer_data_matrix[, 45]
metast_stats <- c(mean(first_metast_col), sd(first_metast_col))
metast_cols<- c("mean", "stdev")
metast_data <- as.data.frame(metast_stats, metast_cols)
metast_data
}

#column indices of all metastatic samples

met_samples <- grep("^metastatic$", colnames(cancer_data_matrix))

#column indices of all non-metastatic samples

nonmet_samples <- grep("^non-metastatic$", colnames(cancer_data_matrix))

#using previous vectors, create vector containing expr differences in genes of the og matrix
#find the mean of all of the metastatic and the non-metastatic genes, then find the difference

for (i in nrow(cancer_data_matrix)){
 expr_difference <- rowMeans(cancer_data_matrix[, met_samples]) - rowMeans(cancer_data_matrix[, nonmet_samples])
}

#print the genes with the largest positive difference and the largest negative difference

print(expr_difference[which.min(expr_difference)]) #minimum
print(expr_difference[which.max(expr_difference)]) #maximum

#Part 2: Statistical analysis of gene expression data

#t test for the differences in means of metastatic and non metastatic samples

{
pValues <- apply(cancer_data_matrix, 1, function (x)
  t.test(x[met_samples], x[nonmet_samples])$p.value)

tStat <- apply(cancer_data_matrix, 1, function (x)
  t.test(x[met_samples], x[nonmet_samples])$statistic)

 
stat_df <- data.frame(pValues, tStat)
stat_df<-na.omit(stat_df)
}

#list of significantly differentially expressed genes, sorted by their order of significance.
{
library(dplyr)

diff_exp <- arrange(stat_df, pValues)
head(diff_exp)
}

#how many genes are significantly differentially expressed in the list

sum(stat_df$pValues < 0.05) #2513 genes are differentially-significantly expressed

#how many genes are significantly under-expressed in metastatic vs. non-metastatic tumors

with(stat_df, c(
  sum(pValues < 0.05 & tStat < 1))) #1331 genes

#How many genes are significantly over-expressed in metastatic vs. non metastatic tumors

with(stat_df, c(
  sum(pValues < 0.05 & tStat > 1))) #1182 genes
  