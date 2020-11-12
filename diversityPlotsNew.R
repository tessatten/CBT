library("devtools")
library("ggplot2")
library("vegan")
library("ineq")
library("tidyverse")
library("praise")


#import main data file
tcrstats_new = read.csv("TCR_statistics_March.csv")  # read csv file
class(tcrstats_new$Month.after.transplant)
tcrstats_new$Month.after.transplant2 <- factor(tcrstats_new$Month.after.transplant, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

alpha <- subset(tcrstats_new, subset=(tcrstats_new$Chain=="alpha"))
beta <- subset(tcrstats_new, subset=(tcrstats_new$Chain=="beta"))
alpha$Month.after.transplant2 <- factor(alpha$Month.after.transplant, levels = c("0", "0.75", "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))
beta$Month.after.transplant2 <- factor(beta$Month.after.transplant, levels = c("0", "0.75",  "1", "2", "3", "4", "6", "8", "12", "17", "18", "19", "22", "30", "C", "UCB"))

clinical = subset(tcrstats_new, subset=(tcrstats_new$Patient.or.control =="Patient" & tcrstats_new$CD4.count > 0.01))
clinical$Month.after.transplant3 <- factor(clinical$Month.after.transplant, levels = c("0.75", "1", "2", "3", "6", "12"))

alphapatients = subset(alpha, subset=(alpha$Patient.or.control=="Patient"))
betapatients = subset(beta, subset=(beta$Patient.or.control=="Patient"))
allpatients = subset(tcrstats_new, subset=(tcrstats_new$Patient.or.control=="Patient"))
allpatients$Month.after.transplant.numerical <- allpatients$Month.after.transplant
allpatients$Month.after.transplant.numerical <- as.numeric(levels(allpatients$Month.after.transplant))[allpatients$Month.after.transplant]

#in this file I will make new diversity plots for the individual patients

#control a
#Control_A <- A
#Control_A_alpha_cdr3 <-A_alpha_dcr_cdr3_freq
#Control_A_beta_cdr3 <- A_beta_dcr_cdr3_freq
#Control_A_alpha_dcr <- A_alpha_dcr_freq_genenames
#Control_A_beta_dcr <- A_beta_dcr_freq_genenames

#patient a/aaf
A_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="A"))
#A_alpha <- AAF_alpha
#A_alpha$ID <- "A"
#A_beta <- AAF_beta
#A_beta$ID <- "A"

#A_alpha_cdr3 <-A_alpha_dcr_cdr3_freq
#A_beta_cdr3 <- A_beta_dcr_cdr3_freq
#A_alpha_dcr <- A_alpha_dcr_freq_genenames
#A_beta_dcr <- A_beta_dcr_freq_genenames

A_info$Month.after.transplant <- factor(A_info$Month.after.transplant,levels = c("1", "2", "3", "6", "12"))

#gini
p = ggplot(A_info, aes(x=A_info$Month.after.transplant, y=A_info$Gini.coefficient, color=A_info$Chain, fill=A_info$Chain)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient A: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )  +geom_line(aes(group=Chain), size =2) + scale_x_continuous(breaks = integer_breaks())+ coord_cartesian(ylim = c(0, 1))
p
ggsave("patient_A_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(A_info, aes(x=A_info$Month.after.transplant, y=A_info$Shannon.entropy, color=A_info$Chain, fill=A_info$Chain)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient A: Shannon entropy", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )+geom_line(aes(group=Chain), size =2) + scale_x_continuous(breaks = integer_breaks())+ coord_cartesian(ylim = c(3, 12))
p
ggsave("patient_A_shannon.pdf", width = 27, height = 20, units = c('cm'))

#patient b/aam
B_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="B"))
B_alpha <- AAM_alpha
B_alpha$ID <- "B"
B_beta <- AAM_beta
B_beta$ID <- "B"
B_info$Month.after.transplant <- factor(B_info$Month.after.transplant,levels = c("1", "2"))

#gini
p = ggplot(B_info, aes(x=B_info$Month.after.transplant, y=B_info$Gini.coefficient, color=B_info$Chain, fill=B_info$Chain)) +
  geom_point(size = 8) + theme(axis.title.x = element_text(size=5),
                               axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient B: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )  +geom_line(aes(group=Chain), size =2) + scale_x_continuous(breaks = integer_breaks())+ coord_cartesian(ylim = c(0, 1))
p
ggsave("patient_B_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(B_info, aes(x=B_info$Month.after.transplant, y=B_info$Shannon.entropy, color=B_info$Chain, fill=B_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient B: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_B_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################
#patient b/aam
C_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="C"))
C_alpha <- AKM_alpha
C_alpha$ID <- "C"
C_beta <- AKM_beta
C_beta$ID <- "C"
C_info$Month.after.transplant <- factor(C_info$Month.after.transplant,levels = c("0", "2"))

#gini
p = ggplot(C_info, aes(x=C_info$Month.after.transplant, y=C_info$Gini.coefficient, color=C_info$Chain, fill=C_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
   labs(title = "Patient C: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_C_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(C_info, aes(x=C_info$Month.after.transplant, y=C_info$Shannon.entropy, color=C_info$Chain, fill=C_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + labs(title = "Patient C: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_C_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################
#patient e/amauf
E_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="E"))
E_alpha <- AmAuF_alpha
E_alpha$ID <- "E"
E_beta <- AmAuF_beta
E_beta$ID <- "E"

E_info$Month.after.transplant <- factor(E_info$Month.after.transplant,levels = c("1", "2", "3"))

#gini
p = ggplot(E_info, aes(x=E_info$Month.after.transplant, y=E_info$Gini.coefficient, color=E_info$Chain, fill=E_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient E: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_E_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(E_info, aes(x=E_info$Month.after.transplant, y=E_info$Shannon.entropy, color=E_info$Chain, fill=E_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient E: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_E_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

#patient f/anf
F_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="F"))
F_alpha <- ANF_alpha
F_alpha$ID <- "F"
F_beta <- ANF_beta
F_beta$ID <- "F"

F_info$Month.after.transplant <- factor(F_info$Month.after.transplant,levels = c("1", "3", "6", "19"))

#gini
p = ggplot(F_info, aes(x=F_info$Month.after.transplant, y=F_info$Gini.coefficient, color=F_info$Chain, fill=F_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient F: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_F_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(F_info, aes(x=F_info$Month.after.transplant, y=F_info$Shannon.entropy, color=F_info$Chain, fill=F_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient F: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_F_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

#patient b/aam
G_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="G"))
#B_alpha <- rbind(AAM1_alpha_dcr_cdr3_freq, AAM2_alpha_dcr_cdr3_freq)
#B_alpha$ID <- "B"
#B_beta <- rbind(AAM1_beta_dcr_cdr3_freq, AAM2_beta_dcr_cdr3_freq)
#B_beta$ID <- "B"
#doesn't have chain or month info
G_info$Month.after.transplant <- factor(G_info$Month.after.transplant,levels = c("1", "6"))

#gini
p = ggplot(G_info, aes(x=G_info$Month.after.transplant, y=G_info$Gini.coefficient, color=G_info$Chain, fill=G_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient G: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_G_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(G_info, aes(x=G_info$Month.after.transplant, y=G_info$Shannon.entropy, color=G_info$Chain, fill=G_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient G: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_G_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient h/hbf
H_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="H"))
H_alpha <- HBF_alpha
H_alpha$ID <- "H"
H_beta <- HBF_beta
H_beta$ID <- "H"

H_info$Month.after.transplant <- factor(H_info$Month.after.transplant,levels = c("1", "2", "3", "6", "12", "18"))

#gini
p = ggplot(H_info, aes(x=H_info$Month.after.transplant, y=H_info$Gini.coefficient, color=H_info$Chain, fill=H_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient H: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_H_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(H_info, aes(x=H_info$Month.after.transplant, y=H_info$Shannon.entropy, color=H_info$Chain, fill=H_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient H: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_H_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################
##############################################

#patient hilkm
I_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="I"))
I_alpha <- LKM_alpha
I_alpha$ID <- "I"
I_beta <- LKM_beta
I_beta$ID <- "I"

I_info$Month.after.transplant <- factor(I_info$Month.after.transplant,levels = c("2", "4", "6", "12"))

#gini
p = ggplot(I_info, aes(x=I_info$Month.after.transplant, y=I_info$Gini.coefficient, color=I_info$Chain, fill=I_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient I: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_I_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(I_info, aes(x=I_info$Month.after.transplant, y=I_info$Shannon.entropy, color=I_info$Chain, fill=I_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient I: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_I_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

#patient b/aam
J_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="J"))
J_alpha <- MLM_alpha
J_alpha$ID <- "J"
J_beta <- MLM_beta
J_beta$ID <- "J"
#doesn't have chain or month info
J_info$Month.after.transplant <- factor(J_info$Month.after.transplant,levels = c("1", "2", "3"))

#gini
p = ggplot(J_info, aes(x=J_info$Month.after.transplant, y=J_info$Gini.coefficient, color=J_info$Chain, fill=J_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient J: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_J_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(J_info, aes(x=J_info$Month.after.transplant, y=J_info$Shannon.entropy, color=J_info$Chain, fill=J_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient J: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_J_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient k/aam
K_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="K"))
K_alpha <- RBM_alpha
K_alpha$ID <- "K"
K_beta <- RBM_beta
K_beta$ID <- "K"
K_info$Month.after.transplant <- factor(K_info$Month.after.transplant,levels = c("1", "2", "3"))

#gini
p = ggplot(K_info, aes(x=K_info$Month.after.transplant, y=K_info$Gini.coefficient, color=K_info$Chain, fill=K_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient K: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_K_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(K_info, aes(x=K_info$Month.after.transplant, y=K_info$Shannon.entropy, color=K_info$Chain, fill=K_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient K: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_K_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient l/rpm
L_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="L"))
L_alpha <- RPM_alpha
L_alpha$ID <- "L"
L_beta <- RPM_beta
L_beta$ID <- "L"
L_info$Month.after.transplant <- factor(L_info$Month.after.transplant,levels = c("1", "2", "3", "12", "22"))

#gini
p = ggplot(L_info, aes(x=L_info$Month.after.transplant, y=L_info$Gini.coefficient, color=L_info$Chain, fill=L_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient L: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_L_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(L_info, aes(x=L_info$Month.after.transplant, y=L_info$Shannon.entropy, color=L_info$Chain, fill=L_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient L: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_L_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient m/ryf
M_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="M"))
M_alpha <- RYF_alpha
M_alpha$ID <- "M"
M_beta <- RYF_beta
M_beta$ID <- "M"
M_info$Month.after.transplant <- factor(M_info$Month.after.transplant,levels = c("1", "2", "6", "8", "12", "22"))

#gini
p = ggplot(M_info, aes(x=M_info$Month.after.transplant, y=M_info$Gini.coefficient, color=M_info$Chain, fill=M_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient M: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_M_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(M_info, aes(x=M_info$Month.after.transplant, y=M_info$Shannon.entropy, color=M_info$Chain, fill=M_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient M: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_M_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient n/WGM
N_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="N"))
N_alpha <- WGM_alpha
N_alpha$ID <- "N"
N_beta <- WGM_beta
N_beta$ID <- "N"
N_info$Month.after.transplant <- factor(N_info$Month.after.transplant,levels = c("1", "3", "12"))

#gini
p = ggplot(N_info, aes(x=N_info$Month.after.transplant, y=N_info$Gini.coefficient, color=N_info$Chain, fill=N_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient N: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_N_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(N_info, aes(x=N_info$Month.after.transplant, y=N_info$Shannon.entropy, color=N_info$Chain, fill=N_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient N: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_N_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient o/yom
O_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="O"))
O_alpha <- YOM_alpha
O_alpha$ID <- "O"
O_beta <- YOM_beta
O_beta$ID <- "O"
O_info$Month.after.transplant <- factor(O_info$Month.after.transplant,levels = c("2", "12"))

#gini
p = ggplot(O_info, aes(x=O_info$Month.after.transplant, y=O_info$Gini.coefficient, color=O_info$Chain, fill=O_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient O: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_O_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(O_info, aes(x=O_info$Month.after.transplant, y=O_info$Shannon.entropy, color=O_info$Chain, fill=O_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient O: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_O_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient o/zrm
P_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="P"))
P_alpha <- ZRM_alpha
P_alpha$ID <- "P"
P_beta <- ZRM_beta
P_beta$ID <- "P"
P_info$Month.after.transplant <- factor(P_info$Month.after.transplant,levels = c("1", "2", "3","6", "12", "17", "30"))

#gini
p = ggplot(P_info, aes(x=P_info$Month.after.transplant, y=P_info$Gini.coefficient, color=P_info$Chain, fill=P_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient P: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_P_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(P_info, aes(x=P_info$Month.after.transplant, y=P_info$Shannon.entropy, color=P_info$Chain, fill=P_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient P: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_P_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################


##############################################

#patient q/ekf
Q_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="Q"))
Q_alpha <- EKF_alpha
Q_alpha$ID <- "Q"
Q_beta <- EKF_beta
Q_beta$ID <- "Q"
Q_info$Month.after.transplant <- factor(Q_info$Month.after.transplant,levels = c("0", "0.75"))

#gini
p = ggplot(Q_info, aes(x=Q_info$Month.after.transplant, y=Q_info$Gini.coefficient, color=Q_info$Chain, fill=Q_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient Q: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_Q_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(Q_info, aes(x=Q_info$Month.after.transplant, y=Q_info$Shannon.entropy, color=Q_info$Chain, fill=Q_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient Q: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_Q_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################

##############################################

#patient q/ekf
R_info<- subset(tcrstats_new, subset=(tcrstats_new$Individual=="R"))
R_alpha <- YAM_alpha
R_alpha$ID <- "R"
R_beta <- YAM_beta
R_beta$ID <- "R"
R_info$Month.after.transplant <- factor(R_info$Month.after.transplant,levels = c("0", "1"))

#gini
p = ggplot(R_info, aes(x=R_info$Month.after.transplant, y=R_info$Gini.coefficient, color=R_info$Chain, fill=R_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient R: Gini coefficient", x = "Month after transplant", y = "Gini coefficient", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#E77EA1", "#5FC27D")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_R_gini.pdf", width = 27, height = 20, units = c('cm'))

#shannon
p = ggplot(R_info, aes(x=R_info$Month.after.transplant, y=R_info$Shannon.entropy, color=R_info$Chain, fill=R_info$Chain)) +
  geom_point(size = 14) + theme(axis.title.x = element_text(size=5),
                                axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) + 
  geom_line() + labs(title = "Patient R: Shannon index", x = "Month after transplant", y = "Shannon index", color = "Chain") +
  scale_color_manual(labels = c("Alpha", "Beta"), values = c("#F0E442", "#56B4E9")) +
  theme_bw() + guides(fill=FALSE)+ theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  )
p
ggsave("patient_R_shannon.pdf", width = 27, height = 20, units = c('cm'))

##############################################