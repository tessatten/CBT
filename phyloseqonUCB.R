library("ape")
library("ggplot2")
library("plyr")
library("scales")
library("vegan")
library("knitr")
library("dplyr")
library("praise")
library("tidyverse")
library("phyloseq")
library("Bios2cor")
library("grid")
library("ggpubr")
library(devtools)
library(nlcor)
library('stringi')
library("microbiome")

#here is am going to do phyloseq analysis on UCB

sample_data(phyloseq_Tcell_alpha)$Month.after.transplant2 = factor(sample_data(phyloseq_Tcell_alpha)$Month.after.transplant2, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

p = plot_richness(phyloseq_Tcell_alpha, x="Month.after.transplant2", color="Individual", measures=c("Shannon", "Observed", "Chao1"))
p + #geom_boxplot(data = p$data, aes(x = Sample_bin, y = value, color = NULL), 
                 #alpha = 0.1)+
  geom_point(size = 10)+ 
  labs(title = "", x = "Month after transplant", y = "Alpha Diversity Measure", color = " ") + theme(
                   plot.title = element_text(color="black", size=28, face="bold"),
                   axis.title.x = element_text(color="black", size=28, face="bold"),
                   axis.title.y = element_text(color="black", size=28, face="bold"),
                   axis.text=element_text(size=20),
                   legend.title = element_text(color = "black", size = 26),
                   legend.text = element_text(color = "black", size = 18))
ggsave("diversityInPatientsOverTimeall3alpha.pdf", width = 47, height = 40, units = c('cm'))

sample_data(phyloseq_Tcell_beta)$Month.after.transplant2 = factor(sample_data(phyloseq_Tcell_beta)$Month.after.transplant2, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

p = plot_richness(phyloseq_Tcell_beta, x="Month.after.transplant2", color="Individual", measures=c("Shannon", "Observed", "Chao1"))
p + #geom_boxplot(data = p$data, aes(x = Sample_bin, y = value, color = NULL), 
  #beta = 0.1)+
  geom_point(size = 10)+ 
  labs(title = "", x = "Month after transplant", y = "Beta Diversity Measure", color = " ") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 18))
ggsave("diversityInPatientsOverTimeall3beta.pdf", width = 47, height = 40, units = c('cm'))

p = plot_richness(phyloseq_Tcell_alpha, x="Sample.bin", color="Individual", measures=c("Shannon", "Observed", "Chao1"))
p + geom_boxplot(data = p$data, aes(x = Sample.bin, y = value, color = NULL), 
  alpha = 0.1)+
  geom_point(size = 10)+ 
  labs(title = "", x = "Month after transplant", y = "Alpha Diversity Measure", color = " ") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 18))
ggsave("diversityInPatientsOverTimeall3alphaboxplot.pdf", width = 47, height = 40, units = c('cm'))


p = plot_richness(phyloseq_Tcell_beta, x="Sample.bin", color="Individual", measures=c("Shannon", "Observed", "Chao1"))
p + geom_boxplot(data = p$data, aes(x = Sample.bin, y = value, color = NULL), 
  alpha = 0.1)+
  geom_point(size = 10)+ 
  labs(title = "", x = "Month after transplant", y = "Beta Diversity Measure", color = " ") + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=20),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 18))
ggsave("diversityInPatientsOverTimeall3betaboxplot.pdf", width = 47, height = 40, units = c('cm'))


#can i heatmap???
gpt <- prune_taxa(names(sort(taxa_sums(phyloseq_Tcell_alpha),TRUE)[1:100]), phyloseq_Tcell_alpha)
plot_heatmap(gpt, sample.label="Sample_ID_phyloseq")
ggsave("heatmap100alpha.pdf", width = 47, height = 40, units = c('cm'))

gpt <- prune_taxa(names(sort(taxa_sums(phyloseq_Tcell_beta),TRUE)[1:100]), phyloseq_Tcell_beta)
plot_heatmap(gpt, sample.label="Sample_ID_phyloseq")
ggsave("heatmap100beta.pdf", width = 47, height = 40, units = c('cm'))


topN <- 20
most_abundant_taxa <- sort(taxa_sums(phyloseq_Tcell_alpha), TRUE)[1:100]
print(most_abundant_taxa)
GP20 <- prune_taxa(names(most_abundant_taxa), phyloseq_Tcell_alpha)

plot_heatmap(GP20, sample.label="Sample_ID_phyloseq")


#can i make the heatmap smaller?
gpt <- prune_taxa(names(sort(taxa_sums(phyloseq_Tcell_alpha),TRUE)[1:100]), phyloseq_Tcell_alpha)
gpt

gpt = prune_taxa(taxa_sums(phyloseq_Tcell_alpha) > 100, phyloseq_Tcell_alpha)
gpt <- prune_taxa(names(sort(taxa_sums(gpt),TRUE)[1:100]), gpt)
plot_heatmap(gpt, "NMDS", "bray", "Individual")
ggsave("heatmap50alpha.pdf", width = 47, height = 40, units = c('cm'))

#phyloseq_Tcell_alpha.rel <- microbiome::transform(phyloseq_Tcell_alpha, "compositional") #transforms into rel abundance

#i'd like to transform my phyloseq data into relative abundance
#extract cdrs in all tthe samples
GPr  = transform_sample_counts(phyloseq_Tcell_alpha, function(x) x / sum(x) ) #transform to relative abundace
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE) #only CDR3s with a mean greater than 10^-5 are kept

gpt <- prune_taxa(names(sort(taxa_sums(GPfr),TRUE)[1:100]), GPfr) #pcks out the 100 most abundant CDR3s
set.seed(67)
heatmapAlpha = plot_heatmap(gpt, "NMDS", "bray", "Individual", low="#000033", high="#FF3300")
heatmapAlpha$scales$scales[[2]]$name <- "Abundance"
heatmapAlpha + theme (axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))+
  ylab('CDR3 sequence')+theme(axis.text=element_text(size=36),
                              axis.title=element_text(size=36,face="bold"))+ theme(legend.text = element_text(size=25), legend.title = element_text(size=25,face="bold"))
ggsave("heatmapalphaCBT.png", width = 60, height = 80, units = c('cm'), dpi =400 )


#this is the bit I used
gpt <- prune_taxa(names(sort(taxa_sums(phyloseq_Tcell_beta),TRUE)[1:10000]), phyloseq_Tcell_beta) #pcks out the 100 most abundant CDR3s

GPr  = transform_sample_counts(gpt, function(x) x / sum(x) ) #transform to relative abundace
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE) #only CDR3s with a mean greater than 10^-5 are kept
gpt <- prune_taxa(names(sort(taxa_sums(GPfr),TRUE)[1:20]), GPfr) #pcks out the 100 most abundant CDR3s
otusBETA <- otu_table(gpt)

write.csv(otusBETA, file='otusBETA.csv')
gpt <- prune_taxa(names(sort(taxa_sums(GPfr),TRUE)[1:100]), GPfr) #pcks out the 100 most abundant CDR3s
set.seed(68)
heatmapBeta = plot_heatmap(gpt, "NMDS", "bray", "Individual", low="#000033", high="#FF3300")
heatmapBeta$scales$scales[[2]]$name <- "Abundance"
heatmapBeta + theme (axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))+
  ylab('CDR3 sequence')+theme(axis.text=element_text(size=36),
                              axis.title=element_text(size=36,face="bold"))+ theme(legend.text = element_text(size=25), legend.title = element_text(size=25,face="bold"))
ggsave("heatmapbetaCBT.png", width = 60, height = 80, units = c('cm'), dpi =400 )

#ordinate

GP.ord <- ordinate(phyloseq_Tcell_alpha, "NMDS", "bray")
p1 = plot_ordination(phyloseq_Tcell_alpha, GP.ord, type="samples", color="Individual", title="taxa") + geom_point(size = 10)
print(p1)


#test

gpt <- prune_taxa(names(sort(taxa_sums(phyloseq_Tcell_alpha),TRUE)[1:10000]), phyloseq_Tcell_alpha) #pcks out the 100 most abundant CDR3s

GPr  = transform_sample_counts(gpt, function(x) x / sum(x) ) #transform to relative abundace
GPfr = filter_taxa(GPr, function(x) mean(x) > 1e-5, TRUE) #only CDR3s with a mean greater than 10^-5 are kept
gpt <- prune_taxa(names(sort(taxa_sums(GPfr),TRUE)[1:20]), GPfr) #pcks out the 100 most abundant CDR3s
otusALPHA <- otu_table(gpt)

write.csv(otusALPHA, file='otusALPHA.csv')
gpt <- prune_taxa(names(sort(taxa_sums(GPfr),TRUE)[1:100]), GPfr) #pcks out the 100 most abundant CDR3s
set.seed(68)
heatmapAlpha = plot_heatmap(gpt, "NMDS", "bray", "Individual", low="#000033", high="#FF3300")
heatmapAlpha$scales$scales[[2]]$name <- "Abundance"
heatmapAlpha + theme (axis.text.x = element_text(size=20), axis.text.y = element_text(size=20))+
  ylab('CDR3 sequence')+theme(axis.text=element_text(size=36),
                              axis.title=element_text(size=36,face="bold"))+ theme(legend.text = element_text(size=25), legend.title = element_text(size=25,face="bold"))
ggsave("heatmapalphaCBT.png", width = 60, height = 80, units = c('cm'), dpi =400 )

######
#subset patient and control
control.beta.phy = subset_samples(phyloseq_Tcell_beta, Patient.or.control=="Control")
patient.beta.phy = subset_samples(phyloseq_Tcell_beta, Patient.or.control=="Patient")
