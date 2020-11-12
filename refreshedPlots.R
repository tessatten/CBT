

library("devtools")
library("ggplot2")
library("vegan")
library("ineq")
library("tidyverse")
library("praise")
library("patchwork")

#this is a file for refreshed plots

tcrstats_new$Month.after.transplant2 <- factor(tcrstats_new$Month.after.transplant, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))
alpha$Month.after.transplant2 <- factor(alpha$Month.after.transplant, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))
beta$Month.after.transplant2 <- factor(beta$Month.after.transplant, levels = c("0", "0.75",  "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

noUCBalpha <- subset(noUCBs, subset=(noUCBs$Chain=="alpha"))
noUCBbeta <- subset(noUCBs, subset=(noUCBs$Chain=="beta"))

#########################################
#number of unique clonotypes in all samples
ggplot(tcrstats_new, aes(x=tcrstats_new$Sample.ID, y=tcrstats_new$Number.of.unique.sequences, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Sample ID") + 
  ylab("Number of Unique Clonotypes") +
  guides(fill=FALSE)+ theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=12),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=8))

ggsave("justonceclonos.pdf", width = 57, height = 35, units = c('cm'))

#########################################
#number of vj rearrangement reads in all samples
noUCBs <- subset(tcrstats_new, subset=(tcrstats_new$Month.after.transplant !="UCB"))

ggplot(noUCBs, aes(x=noUCBs$Sample.ID, y=noUCBs$Number.of.VJ.rearrangements, color=noUCBs$Individual, fill=noUCBs$Individual)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Sample ID") + 
  ylab("Number of Reads") +
  guides(fill=FALSE) + theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=12),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=8))

ggsave("numberofReads.pdf", width = 57, height = 35, units = c('cm'))

#########################################

alpha <- subset(tcrstats_new, subset=(tcrstats_new$Chain=="alpha"))
beta <- subset(tcrstats_new, subset=(tcrstats_new$Chain=="beta"))

noUCBalpha <- subset(noUCBs, subset=(noUCBs$Chain=="alpha"))
noUCBbeta <- subset(noUCBs, subset=(noUCBs$Chain=="beta"))

#########################################
#number of vj rearrangement reads in alpha samples
ggplot(noUCBalpha, aes(x=noUCBalpha$Sample.ID, y=noUCBalpha$Percentage.of.rearrangements.found, color=noUCBalpha$Individual, fill=noUCBalpha$Individual)) +
  geom_bar(colour="black", stat="identity") + 
  xlab("Sample ID") + 
  ylab("Percentage of reads for which the clonotype was identified (%)") +
  guides(fill=FALSE)+ theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=12),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=8))+ coord_cartesian(ylim = c(0, 100))

ggsave("percofreadsAlpha.pdf", width = 40, height = 40, units = c('cm'))

#########################################

#number of vj rearrangement reads in beta samples
ggplot(noUCBbeta, aes(x=noUCBbeta$Sample.ID, y=noUCBbeta$Percentage.of.rearrangements.found, color=noUCBbeta$Individual, fill=noUCBbeta$Individual)) +
  geom_bar(colour="black", stat="identity") +
  xlab("Sample ID") + 
  ylab("Percentage of reads for which the clonotype was identified (%)") +
  guides(fill=FALSE)+ theme_bw() + theme(
    plot.title = element_text(color="black", size=20, face="bold"),
    axis.title.x = element_text(color="black", size=20, face="bold"),
    axis.title.y = element_text(color="black", size=20, face="bold"),
    axis.text=element_text(size=12),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=8))+ coord_cartesian(ylim = c(0, 100))

ggsave("percofreadsBeta.pdf", width = 40, height = 40, units = c('cm'))

#########################################
#subset data with lab and clinical info
clinical = subset(tcrstats_new, subset=(tcrstats_new$Patient.or.control =="Patient" & tcrstats_new$CD4.count > 0.01))
#clinical = subset(tcrstats_new, subset=(tcrstats_new$Patient.or.control =="Patient" & tcrstats_new$CD4.count > 0.01 & tcrstats_new$Sample.ID != "D6CD4-a" & tcrstats_new$Sample.ID != "D6CD4-b"))
clinical$Month.after.transplant3 <- factor(clinical$Month.after.transplant, levels = c("1", "2", "3", "6", "12"))

#########################################

#TRECs vs shannon
ggplot(clinical, aes(x=clinical$TRECs, y=clinical$Shannon.entropy, color=clinical$Individual, fill=clinical$Individual)) +
  geom_point(size=10) + labs(title = " ", x = "TRECs", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("TRECSvsShannonless.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#TRECs vs gini
ggplot(clinical, aes(x=clinical$TRECs, y=clinical$Gini.coefficient, color=clinical$Individual, fill=clinical$Individual)) +
  geom_point(size=10) + labs(title = " ", x = "TRECs", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("TRECSvsGini.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#TCD4 count over time
ggplot(clinical, aes(x=clinical$Month.after.transplant3, y=clinical$CD4.count, color=clinical$Individual, fill=clinical$Individual)) +
  geom_point(size=16) + labs(title = " ", x = "Month after transplant", y = "CD4 count", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("CD4andMonths.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#shannon vs cd4 count by month
ggplot(clinical, aes(x=clinical$CD4.count, y=clinical$Shannon.entropy, color=clinical$Month.after.transplant3, fill=clinical$Month.after.transplant3)) +
  geom_point(size=16) + labs(title = " ", x = "CD4 count", y = "Shannon index", color = "Month after transplant") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("CD4vsShannonMonths.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#shannon vs cd4 count by month
ggplot(clinical, aes(x=clinical$CD4.count, y=clinical$Shannon.entropy, color=clinical$Individual, fill=clinical$Individual)) +
  geom_point(size=16) + labs(title = " ", x = "CD4 count", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("CD4vsShannonPatient.pdf", width = 40, height = 30, units = c('cm'))
############################################
alphaclinical = subset(clinical, subset=(clinical$Chain=="alpha"))
betaclinical = subset(clinical, subset=(clinical$Chain=="beta"))

#########################################

#gini vs cd4 count by patient
ggplot(clinical, aes(x=clinical$CD4.count, y=clinical$Gini.coefficient, color=clinical$Individual, fill=clinical$Individual)) +
  geom_point(size=16) + labs(title = " ", x = "CD4 count", y = "Gini coeffient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("CD4vsGiniPatient.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#gini vs cd4 count by month
ggplot(clinical, aes(x=clinical$CD4.count, y=clinical$Gini.coefficient, color=clinical$Month.after.transplant3, fill=clinical$Month.after.transplant3)) +
  geom_point(size=16) + labs(title = " ", x = "CD4 count", y = "Gini coeffient", color = "Month after transplant") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("CD4vsGiniMonth.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs vdj reads
ggplot(noUCBs, aes(x=noUCBs$Number.of.reads, y=noUCBs$Number.of.VJ.rearrangements, color=noUCBs$Individual, fill=noUCBs$Individual)) +
  geom_point(size=12) + labs(title = " ", x = "Raw read count", y = "Number of identified reads", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsReads.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs vdj reads in alpha
ggplot(noUCBalpha, aes(x=noUCBalpha$Number.of.reads, y=noUCBalpha$Number.of.VJ.rearrangements, color=noUCBalpha$Individual, fill=noUCBalpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha Chain", x = "Raw read count", y = "Number of identified reads", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsReadsALPHA.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs vdj reads beta
ggplot(noUCBbeta, aes(x=noUCBbeta$Number.of.reads, y=noUCBbeta$Number.of.VJ.rearrangements, color=noUCBbeta$Individual, fill=noUCBbeta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain ", x = "Raw read count", y = "Number of identified reads", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsReadsBETA.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#gini vs months alpha
ggplot(alpha, aes(x=alpha$Month.after.transplant2, y=alpha$Gini.coefficient, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Months after transplant", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("GinivsMonthsAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#gini vs months beta
ggplot(beta, aes(x=beta$Month.after.transplant2, y=beta$Gini.coefficient, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Months after transplant", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("GinivsMonthsBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#shannon vs months alpha
ggplot(alpha, aes(x=alpha$Month.after.transplant2, y=alpha$Shannon.entropy, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Months after transplant", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ShannonvsMonthsAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#shannon vs months beta
ggplot(beta, aes(x=beta$Month.after.transplant2, y=beta$Shannon.entropy, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Months after transplant", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ShannonvsMonthsBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#chao1 vs months alpha
ggplot(alpha, aes(x=alpha$Month.after.transplant2, y=alpha$Chao1, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Months after transplant", y = "Chao1 index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("Chao1vsMonthsAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#chai1 vs months beta
ggplot(beta, aes(x=beta$Month.after.transplant2, y=beta$Chao1, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Months after transplant", y = "Chao1 index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("Chao1vsMonthsBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#gini vs shannon alpha
ggplot(alpha, aes(x=alpha$Gini.coefficient, y=alpha$Shannon.entropy, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Gini coefficient", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("GinivsShannonAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#gini vs shannon beta
ggplot(beta, aes(x=beta$Gini.coefficient, y=beta$Shannon.entropy, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Gini coefficient", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("GinivsShannonBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs shannon alpha
ggplot(alpha, aes(x=alpha$Number.of.reads, y=alpha$Shannon.entropy, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsShannonAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs shannon beta
ggplot(beta, aes(x=beta$Number.of.reads, y=beta$Shannon.entropy, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsShannonBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################


#########################################

#reads vs shannon alpha
ggplot(alpha, aes(x=alpha$Number.of.reads, y=alpha$Shannon.entropy, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsShannonAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#reads vs shannon beta
ggplot(beta, aes(x=beta$Number.of.reads, y=beta$Shannon.entropy, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsShannonBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################
alphapatients = subset(alpha, subset=(alpha$Patient.or.control=="Patient"))
betapatients = subset(beta, subset=(beta$Patient.or.control=="Patient"))

#########################################

#reads vs shannon alpha justpatients
ggplot(alphapatients, aes(x=alphapatients$Number.of.reads, y=alphapatients$Shannon.entropy, color=alphapatients$Individual, fill=alphapatients$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)
  #geom_smooth(method = 'glm', se = FALSE)


ggsave("ReadsvsShannonAlphaPatients.pdf", width = 40, height = 30, units = c('cm'))
############################################

shanReadAlpha.lm <- lm(alphapatients$Shannon.entropy ~ alphapatients$Number.of.reads, data = alphapatients)
summary(shanReadAlpha.lm)
#########################################

#reads vs shannon beta justpatinets
ggplot(betapatients, aes(x=betapatients$Number.of.reads, y=betapatients$Shannon.entropy, color=betapatients$Individual, fill=betapatients$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsShannonBetaPatients.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#reads vs gini alpha
ggplot(alpha, aes(x=alpha$Number.of.reads, y=alpha$Gini.coefficient, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsGiniAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#reads vs gin i beta
ggplot(beta, aes(x=beta$Number.of.reads, y=beta$Gini.coefficient, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ReadsvsGiniBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#shannon vs ng rna alpha
ggplot(alpha, aes(x=alpha$ng.of.RNA, y=alpha$Shannon.entropy, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "ng of RNA", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ShannonvNGAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#shannon  vs ng rna0 beta
ggplot(beta, aes(x=beta$ng.of.RNA, y=beta$Shannon.entropy, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "ng of RNA", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ShannonvNGGiniBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#gini  vs spectratyping
ggplot(clinical, aes(x=clinical$Month.after.transplant3, y=clinical$Gini.coefficient, color=clinical$Spectratyping, fill=clinical$Spectratyping)) +
  geom_point(size=12, position ="jitter") + labs(title = "", x = "Month after transplant", y = "Gini coefficient", color = "Spectratyping result") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("GinivsSpectra.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#shannon  vs spectratyping
ggplot(clinical, aes(x=clinical$Month.after.transplant3, y=clinical$Shannon.entropy, color=clinical$Spectratyping, fill=clinical$Spectratyping)) +
  geom_point(size=12, position ="jitter") + labs(title = "", x = "Month after transplant", y = "Shannon index", color = "Spectratyping result") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ShannonvsSpectra.pdf", width = 40, height = 30, units = c('cm'))
############################################
#largest clonal  vs month
ggplot(tcrstats_new, aes(x=tcrstats_new$Month.after.transplant2, y=tcrstats_new$Mean.abundance.of.top.100.clonotypes.CDR3, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_point(size=12) + labs(title = "", x = "Month after transplant", y = "Mean abundance of the 100 largest clonotypes (CDR3)", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("Top100Clono.pdf", width = 40, height = 30, units = c('cm'))
############################################

############################################
#largest clonal  vs month
ggplot(tcrstats_new, aes(x=tcrstats_new$Mean.abundance.of.top.100.clonotypes.CDR3, y=tcrstats_new$Largest.clonal.expansion, color=tcrstats_new$Individual, fill=tcrstats_new$Individual)) +
  geom_point(size=12) + labs(title = "", x = "Mean abundance of the 100 largest clonotypes (CDR3)", y = "Largest clonal expansion", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ClonovsClono.pdf", width = 40, height = 30, units = c('cm'))
############################################

############################################
#largest trecs  vs month
ggplot(clinical, aes(x=clinical$Month.after.transplant2, y=clinical$TRECs, color=clinical$Individual, fill=clinical$Individual)) +
  geom_point(size=12) + labs(title = "", x = "Month after transplant", y = "TRECS", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("TRECSvsMonth.pdf", width = 40, height = 30, units = c('cm'))
############################################


############################################
#number of vj vs error correct
ggplot(noUCBs, aes(x=noUCBs$Number.of.error.corrected.reads, y=noUCBs$Number.of.VJ.rearrangements, color=noUCBs$Individual, fill=noUCBs$Individual)) +
  geom_point(size=12) + labs(title = "", x = "Number of reads (error-corrected)", y = "Number of reads (VJ)", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsErrorCorrectRead.pdf", width = 40, height = 30, units = c('cm'))
############################################
spearman <-cor.test(noUCBs$Number.of.VJ.rearrangements, noUCBs$Number.of.error.corrected.reads,  method = "spearman")
spearman
############################################
############################################
#number of vj vs error correct ALPHA
ggplot(noUCBalpha, aes(x=noUCBalpha$Number.of.error.corrected.reads, y=noUCBalpha$Number.of.VJ.rearrangements, color=noUCBalpha$Individual, fill=noUCBalpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads (error-corrected)", y = "Number of reads (VJ)", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsErrorCorrectReadALPHA.pdf", width = 40, height = 30, units = c('cm'))
############################################
spearman <-cor.test(noUCBalpha$Number.of.VJ.rearrangements, noUCBalpha$Number.of.error.corrected.reads,  method = "spearman")
spearman
############################################
############################################
#number of vj vs error correct BETA
ggplot(noUCBbeta, aes(x=noUCBbeta$Number.of.error.corrected.reads, y=noUCBbeta$Number.of.VJ.rearrangements, color=noUCBbeta$Individual, fill=noUCBbeta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads (error-corrected)", y = "Number of reads (VJ)", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("VJvsErrorCorrectReadBETA.pdf", width = 40, height = 30, units = c('cm'))
############################################
spearman <-cor.test(noUCBbeta$Number.of.VJ.rearrangements, noUCBbeta$Number.of.error.corrected.reads,  method = "spearman")
spearman
############################################
spearman <-cor.test(noUCBbeta$Number.of.error.corrected.reads, noUCBbeta$Shannon.entropy,  method = "spearman")
spearman

spearman <-cor.test(noUCBalpha$Number.of.error.corrected.reads, noUCBalpha$Shannon.entropy,  method = "spearman")
spearman

#########################################

#error-corrected reads vs shannon alpha
ggplot(alpha, aes(x=alpha$Number.of.error.corrected.reads, y=alpha$Shannon.entropy, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads (error corrected)", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsShannonAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

# error-corrected reads vs shannon beta
ggplot(beta, aes(x=beta$Number.of.error.corrected.reads, y=beta$Shannon.entropy, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads (error corrected)", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsShannonBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################
alphapatients = subset(alpha, subset=(alpha$Patient.or.control=="Patient"))
betapatients = subset(beta, subset=(beta$Patient.or.control=="Patient"))

#########################################

#error-corrected reads vs shannon alpha justpatients
ggplot(alphapatients, aes(x=alphapatients$Number.of.error.corrected.reads, y=alphapatients$Shannon.entropy, color=alphapatients$Individual, fill=alphapatients$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads (error corrected)", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)
#geom_smooth(method = 'glm', se = FALSE)


ggsave("ErrorCorrectReadsvsShannonAlphaPatients.pdf", width = 40, height = 30, units = c('cm'))
############################################

shanReadAlpha.lm <- lm(alphapatients$Shannon.entropy ~ alphapatients$Number.of.reads, data = alphapatients)
summary(shanReadAlpha.lm)
#########################################

#error-corrected reads vs shannon beta justpatinets
ggplot(betapatients, aes(x=betapatients$Number.of.error.corrected.reads, y=betapatients$Shannon.entropy, color=betapatients$Individual, fill=betapatients$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads (error corrected)", y = "Shannon index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsShannonBetaPatients.pdf", width = 40, height = 30, units = c('cm'))
############################################

#########################################

#error-corrected reads vs gini alpha
ggplot(alpha, aes(x=alpha$Number.of.error.corrected.reads, y=alpha$Gini.coefficient, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads (error corrected)", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsGiniAlpha.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#rerror-corrected eads vs gin i beta
ggplot(beta, aes(x=beta$Number.of.error.corrected.reads, y=beta$Gini.coefficient, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads (error corrected)", y = "Gini coefficient", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsGiniBeta.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#error-corrected reads vs chao alpha
ggplot(alpha, aes(x=alpha$Number.of.error.corrected.reads, y=alpha$Chao1, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(size=12) + labs(title = "Alpha chain", x = "Number of reads (error corrected)", y = "Chao1 index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsChao1Alpha.pdf", width = 40, height = 30, units = c('cm'))
############################################
#########################################

#rerror-corrected eads vs chao beta
ggplot(beta, aes(x=beta$Number.of.error.corrected.reads, y=beta$Chao1, color=beta$Individual, fill=beta$Individual)) +
  geom_point(size=12) + labs(title = "Beta chain", x = "Number of reads (error corrected)", y = "Chao1 index", color = "Patient") + theme_bw() + theme(
    plot.title = element_text(color="black", size=26, face="bold"),
    axis.title.x = element_text(color="black", size=26, face="bold"),
    axis.title.y = element_text(color="black", size=26, face="bold"),
    axis.text=element_text(size=22),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + theme(axis.text.x = element_text(angle = 90, hjust = 0.5,size=22)) + guides(fill=FALSE)

ggsave("ErrorCorrectReadsvsChao1Beta.pdf", width = 40, height = 30, units = c('cm'))
############################################

