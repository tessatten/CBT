library("devtools")
library("ggplot2")
library("vegan")
library("ineq")
library("tidyverse")
library("praise")
library("phyloseq")

#data needs to be a  a matrix
otu_data_table <- read_csv(file = "alpha_otu_table_all_noNA.csv") #import otu table
#otu_data_table <- read_csv(file = "CVG_otutable_foranosim.csv") #import otu table
otu_data_matrix <- data.matrix(otu_data_table) #make the otu table into a matrix (it has to be a matrix)
rownames(otu_data_matrix) <- otu_data_table$`CDR3 AA sequence` #relabel the column names at the OTU numbers
otu_data_matrix <- otu_data_matrix[,-1] #remove null column
otu_data_matrix <- otu_data_matrix[,-1] #remove 2nd null column

UCB_alpha_otu_table <- otu_table(otu_data_matrix, taxa_are_rows = T) #turn OTU table into a phyloseq object

tax_matrix <- as.matrix(distinctCDR3s) #make the taxa table into a matrix (it has to be a matrix)
rownames(tax_matrix) <- distinctCDR3s$`CDR3 AA sequence` #relabel the column names to be the OTU numbers

UCB_alpha_tax_table <- tax_table(tax_matrix)

#sample_data_table <- read.csv(file = "CVG_metadata_JULY.csv") #import metadata for samples (all samples)
sample_data_table <- read.csv(file = "TCR_statistics_Jan.csv") #import metadata for samples (samples that worked)


sample_data_table <- subset(sample_data_table, subset=(sample_data_table$Chain=="alpha"))
#sample_data_table <- subset(sample_data_table, subset=(sample_data_table$Patient.or.control=="Patient"))

sample_matrix <- as.matrix(sample_data_table) #make it a matrix

rownames(sample_matrix) <- sample_data_table$Sample_ID_phyloseq #change column names, 
#sample_matrix <- sample_matrix[,-1]  #remove null column
sample_matrix <- as.data.frame(sample_matrix) #make it a data frame again
#sam_data is depreciated so switched to sample_data
#required a data frame so switched it back to a data frame


UCB_sample_table <- sample_data(sample_matrix) #import to phyloseq


#phyloseq_Tcell_alpha <- merge_phyloseq(UCB_otu_table, UCB_sample_table) #merge phyloseq objects
phyloseq_Tcell_alpha <- merge_phyloseq(UCB_alpha_otu_table, UCB_sample_table, UCB_alpha_tax_table) #merge phyloseq objects
#sample_names(UCB_alpha_otu_table)
#sample_names(UCB_sample_table)

GPr  = transform_sample_counts(GlobalPatterns, function(x) x / sum(x) )
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE)

##################################
#data needs to be a  a matrix
otu_data_table <- read_csv(file = "beta_otu_table_all_noNA.csv") #import otu table
#otu_data_table <- read_csv(file = "CVG_otutable_foranosim.csv") #import otu table
otu_data_matrix <- data.matrix(otu_data_table) #make the otu table into a matrix (it has to be a matrix)
rownames(otu_data_matrix) <- otu_data_table$`CDR3 AA sequence` #relabel the column names at the OTU numbers
otu_data_matrix <- otu_data_matrix[,-1] #remove null column
otu_data_matrix <- otu_data_matrix[,-1] #remove 2nd null column

UCB_beta_otu_table <- otu_table(otu_data_matrix, taxa_are_rows = T) #turn OTU table into a phyloseq object

#sample_data_table <- read.csv(file = "CVG_metadata_JULY.csv") #import metadata for samples (all samples)
sample_data_table <- read.csv(file = "TCR_statistics_Jan.csv") #import metadata for samples (samples that worked)
sample_data_table$Month.after.transplant2 <- factor(sample_data_table$Month.after.transplant, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))
sample_data_table$Month.after.transplant.numeric <- NULL

sample_data_table <- subset(sample_data_table, subset=(sample_data_table$Chain=="beta"))
#sample_data_table <- subset(sample_data_table, subset=(sample_data_table$Patient.or.control=="Patient"))

sample_matrix <- as.matrix(sample_data_table) #make it a matrix

rownames(sample_matrix) <- sample_data_table$Sample_ID_phyloseq #change column names, 
#sample_matrix <- sample_matrix[,-1]  #remove null column
sample_matrix <- as.data.frame(sample_matrix) #make it a data frame again
#sam_data is depreciated so switched to sample_data
#required a data frame so switched it back to a data frame


UCB_sample_table <- sample_data(sample_matrix) #import to phyloseq

tax_matrix <- as.matrix(distinctCDR3s) #make the taxa table into a matrix (it has to be a matrix)
rownames(tax_matrix) <- distinctCDR3s$`CDR3 AA sequence` #relabel the column names to be the OTU numbers

UCB_beta_tax_table <- tax_table(tax_matrix)


#phyloseq_Tcell_alpha <- merge_phyloseq(UCB_otu_table, UCB_sample_table) #merge phyloseq objects
phyloseq_Tcell_beta <- merge_phyloseq(UCB_beta_otu_table, UCB_sample_table, UCB_beta_tax_table) #merge phyloseq objects
#sample_names(UCB_otu_table)
#sample_names(UCB_sample_table)