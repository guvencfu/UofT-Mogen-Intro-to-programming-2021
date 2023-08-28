#q1 -----

# There are 9 columns.
# 
# GeneID - NCBI gene ID
# AssociatedGenes - Symbol corresponding to GeneID 
# RelatedGenes - genes related to the current gene
# ConceptID - Identifier associated with disorder associated with this gene.
# DiseaseName - Full name of the condition associated with the gene.
# SourceName - Sources that use this name
# SourceID - identifier used by this source
# DiseaseMIM - MIM number used for this condition
   #MIM number - Number denoting the gene in the Mendelian Inheritance in Man
# LastUpdated - date the database was updated by NCBI.

#q2 -----

cut -f 1 gene_condition_source_id.txt| sort | uniq | wc -l #4504 unique terms
cut -f 2 gene_condition_source_id.txt | sort | uniq | wc -l #4505 unique terms
cut -f 3 gene_condition_source_id.txt | sort | uniq | wc -l #1408 unique terms
cut -f 4 gene_condition_source_id.txt | sort | uniq | wc -l #6404 unique terms
cut -f 5 gene_condition_source_id.txt | sort | uniq | wc -l #6405 unique terms
cut -f 6 gene_condition_source_id.txt | sort | uniq | wc -l #9 unique terms
cut -f 7 gene_condition_source_id.txt | sort | uniq | wc -l #3761 unique terms
cut -f 8 gene_condition_source_id.txt | sort | uniq | wc -l #5891 unique terms
cut -f 9 gene_condition_source_id.txt | sort | uniq | wc -l #538 unique terms

#q3 -----
grep cancer -w gene_condition_source_id.txt | wc -l
#161

#q4 -----
grep cancer gene_condition_source_id.txt | sort | uniq | wc -l
  #145

#q5 -----
grep TP53 -w gene_condition_source_id.txt | sort | uniq | wc -l
#20

#q6 -----
#There are 20 unique entries for TP53 in ClinVar database
#There are 111 unique entries for TP53 in Human Ontology Database (HOD)
#ClinVar exlusively cancer, HOD contains cancer and other afflictions.
#Better to use Human Ontology Database; more comprehensive and diverse.

