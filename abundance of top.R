tcrstats = read.csv("TCR statistics 14 - 11.csv")  # read csv file
tcrstats$Month.after.transplant2 <- factor(tcrstats$Month.after.transplant, levels = c("1", "2", "3", "6", "8", "12", "18", "22"))

sum(A_alpha_dcr_freq_genenames$"Number found in sample" == 1)
sum(A_alpha_dcr_freq_genenames$"Number found in sample")
range(A_alpha_dcr_freq_genenames$"Number found in sample")

# based on variable values
alpha = subset(tcrstats, subset=(tcrstats$Chain=="alpha" & tcrstats$Patient.or.control =="Patient"))
beta = subset(tcrstats, subset=(tcrstats$Chain=="beta"& tcrstats$Patient.or.control =="Patient"))

# based on variable values
alpha = subset(tcrstats, subset=(tcrstats$Chain=="alpha"))
beta = subset(tcrstats, subset=(tcrstats$Chain=="beta"))
#################################################################
#plotting average abundance of the top 100 clones in each sample (CDR3)

#barplot of abundance of top 100 CDR3s
ggplot(tcrstats, aes(x=tcrstats$Sample.ID, y=tcrstats$Mean.abundance.of.top.100.clonotypes.CDR3, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_bar(colour="black", stat="identity") + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                                                    axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Sample ID") + 
  ylab("Mean abundance of the 100 largest clonotypes (CDR3)") +
  guides(fill=FALSE)

##########################################
# Manual month levels
tcrstats$Month.after.transplant2 <- factor(tcrstats$Month.after.transplant, levels = c("1", "2", "3", "6", "8", "12", "18", "22"))
############################################

#average abundance of the top 100 clones in each sample (CDR3)
ggplot(tcrstats, aes(x=tcrstats$Month.after.transplant2, y=tcrstats$Mean.abundance.of.top.100.clonotypes.CDR3, color=tcrstats$Individual, fill=tcrstats$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Mean abundance of the 100 largest clonotypes (CDR3)") +
  guides(fill=FALSE)

############################################

#average abundance of the top 100 clones in each sample (CDR3) ALPHA
ggplot(alpha, aes(x=alpha$Month.after.transplant, y=alpha$Mean.abundance.of.top.100.clonotypes.CDR3, color=alpha$Individual, fill=alpha$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Mean abundance of the 100 largest alpha clonotypes (CDR3)") +
  guides(fill=FALSE)

############################################

#average abundance of the top 100 clones in each sample (CDR3) BETA
ggplot(beta, aes(x=beta$Month.after.transplant2, y=beta$Mean.abundance.top.100.CDR3, color=beta$Individual, fill=beta$Individual)) +
  geom_point(shape=1) + theme(axis.title.x = element_text(face="bold", colour="#990000", size=5),
                              axis.text.x  = element_text(angle=90, vjust=0.5, size=5)) +   # Remove x-axis label
  xlab("Month after transplant") + 
  ylab("Mean abundance of the 100 largest beta clonotypes (CDR3)") +
  guides(fill=FALSE)

#################################################################
#calculating average abundance of the top 100 clones in each sample (CDR3)

x = (sum(tail(sort(A_alpha_dcr_cdr3_freq$`Number found in sample`), 100)))/100

sum(tail(sort(YOM12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 1
###############################################
#A alpha

sum(tail(sort(A_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#A beta

sum(tail(sort(A_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C2 alpha

sum(tail(sort(C2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C2 beta

sum(tail(sort(C2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#E alpha

sum(tail(sort(E_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#E beta

sum(tail(sort(E_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#M alpha

sum(tail(sort(M_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#M beta

sum(tail(sort(M_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#E1 alpha

sum(tail(sort(E1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#E1 beta

sum(tail(sort(E1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C3 alpha

sum(tail(sort(C3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C3 beta

sum(tail(sort(C3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100


#run 2
###############################################
#YOM2 alpha

sum(tail(sort(YOM2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#YOM2 beta

sum(tail(sort(YOM2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#YOM12 alpha

sum(tail(sort(YOM12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#YOM12 beta

sum(tail(sort(YOM12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C300 alpha

sum(tail(sort(C300_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C300 beta

sum(tail(sort(C300_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C400 alpha

sum(tail(sort(C400_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#C400 beta

sum(tail(sort(C400_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#DSF1 alpha

sum(tail(sort(DSF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#DSF1 beta

sum(tail(sort(DSF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#DSF6 alpha

sum(tail(sort(DSF6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#DSF6 beta

sum(tail(sort(DSF6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 3
###############################################
#AAF1 alpha

sum(tail(sort(AAF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF1 beta

sum(tail(sort(AAF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF2 alpha

sum(tail(sort(AAF2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF2 beta

sum(tail(sort(AAF2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF3 alpha

sum(tail(sort(AAF3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF3 beta

sum(tail(sort(AAF3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF6 alpha

sum(tail(sort(AAF6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF6 beta

sum(tail(sort(AAF6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF12 alpha

sum(tail(sort(AAF12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF12 beta

sum(tail(sort(AAF12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#L200 alpha

sum(tail(sort(L200_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#L200 beta

sum(tail(sort(L200_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 4
###############################################
#HBF1 alpha

sum(tail(sort(HBF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF1 beta

sum(tail(sort(HBF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF2 alpha

sum(tail(sort(HBF2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF2 beta

sum(tail(sort(HBF2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF3 alpha

sum(tail(sort(HBF3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF3 beta

sum(tail(sort(HBF3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF6 alpha

sum(tail(sort(HBF6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF6 beta

sum(tail(sort(HBF6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF12 alpha

sum(tail(sort(HBF12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF12 beta

sum(tail(sort(HBF12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF18 alpha

sum(tail(sort(HBF18_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#HBF18 beta

sum(tail(sort(HBF18_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 5
###############################################
#RYF1 alpha

sum(tail(sort(RYF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF1 beta

sum(tail(sort(RYF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF2 alpha

sum(tail(sort(RYF2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF2 beta

sum(tail(sort(RYF2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF6 alpha

sum(tail(sort(RYF6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF6 beta

sum(tail(sort(RYF6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF8 alpha

sum(tail(sort(RYF8_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF8 beta

sum(tail(sort(RYF8_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF12 alpha

sum(tail(sort(RYF12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF12 beta

sum(tail(sort(RYF12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF22 alpha

sum(tail(sort(RYF22_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RYF22 beta

sum(tail(sort(RYF22_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 6
###############################################
#RPM1 alpha

sum(tail(sort(RPM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM1 beta

sum(tail(sort(RPM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM2 alpha

sum(tail(sort(RPM2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM2 beta

sum(tail(sort(RPM2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM3 alpha

sum(tail(sort(RPM3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM3 beta

sum(tail(sort(RPM3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM12 alpha

sum(tail(sort(RPM12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM12 beta

sum(tail(sort(RPM12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM22 alpha

sum(tail(sort(RPM22_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RPM22 beta

sum(tail(sort(RPM22_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF19 alpha

sum(tail(sort(ANF19_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF19 beta

sum(tail(sort(ANF19_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAF3 alpha

sum(tail(sort(AAF3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 7
###############################################
#ANF1 alpha

sum(tail(sort(ANF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF1 beta

sum(tail(sort(ANF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF3 alpha

sum(tail(sort(ANF3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF3 beta

sum(tail(sort(ANF3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF6 alpha

sum(tail(sort(ANF6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ANF6 beta

sum(tail(sort(ANF6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAuF1 alpha

sum(tail(sort(AmAuF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAuF1 beta

sum(tail(sort(AmAuF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAuF2 alpha

sum(tail(sort(AmAuF2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAuF2 beta

sum(tail(sort(AmAuF2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAuF3 alpha

sum(tail(sort(AmAuF3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAuF3 beta

sum(tail(sort(AmAuF3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 8
###############################################
#AAM1 alpha

sum(tail(sort(AAM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAM1 beta

sum(tail(sort(AAM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAM2 alpha

sum(tail(sort(AAM2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AAM2 beta

sum(tail(sort(AAM2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM1 alpha

sum(tail(sort(ZRM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM1 beta

sum(tail(sort(ZRM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#WGM1 alpha

sum(tail(sort(WGM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#WGM1 beta

sum(tail(sort(WGM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#WGM3 alpha

sum(tail(sort(WGM3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#WGM3 beta

sum(tail(sort(WGM3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#WGM12 alpha

sum(tail(sort(WGM12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#WGM12 beta

sum(tail(sort(WGM12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 9
###############################################
#LKM2 alpha

sum(tail(sort(LKM2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM2 beta

sum(tail(sort(LKM2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM4 alpha

sum(tail(sort(LKM4_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM4 beta

sum(tail(sort(LKM4_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM6 alpha

sum(tail(sort(LKM6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM6 beta

sum(tail(sort(LKM6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM12 alpha

sum(tail(sort(LKM12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#LKM12 beta

sum(tail(sort(LKM12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#MLM1 alpha

sum(tail(sort(MLM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#MLM1 beta

sum(tail(sort(MLM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#MLM3 alpha

sum(tail(sort(MLM3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#MLM3 beta

sum(tail(sort(MLM3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 10
###############################################
#RBM1 alpha

sum(tail(sort(RBM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RBM1 beta

sum(tail(sort(RBM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RBM2 alpha

sum(tail(sort(RBM2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RBM2 beta

sum(tail(sort(RBM2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RBM3 alpha

sum(tail(sort(RBM3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#RBM3 beta

sum(tail(sort(RBM3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM2 alpha

sum(tail(sort(ZRM2_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM2 beta

sum(tail(sort(ZRM2_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM3 alpha

sum(tail(sort(ZRM3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM3 beta

sum(tail(sort(ZRM3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM6 alpha

sum(tail(sort(ZRM6_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM6 beta

sum(tail(sort(ZRM6_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run 11
###############################################
#ZRM12 alpha

sum(tail(sort(ZRM12_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM12 beta

sum(tail(sort(ZRM12_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM17 alpha

sum(tail(sort(ZRM17_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM17 beta

sum(tail(sort(ZRM17_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM30 alpha

sum(tail(sort(ZRM30_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#ZRM30 beta

sum(tail(sort(ZRM30_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAlF1 alpha

sum(tail(sort(AmAlF1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAlF1 beta

sum(tail(sort(AmAlF1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAlF3 alpha

sum(tail(sort(AmAlF3_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAlF3 beta

sum(tail(sort(AmAlF3_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAlF6CD4 alpha

sum(tail(sort(AmAlF6CD4_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#AmAlF6CD4 beta

sum(tail(sort(AmAlF6CD4_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#cordblood
###############################################
#CB2A alpha

sum(tail(sort(CB2A_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2A beta

sum(tail(sort(CB2A_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2B alpha

sum(tail(sort(CB2B_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2B beta

sum(tail(sort(CB2B_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2C alpha

sum(tail(sort(CB2C_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2C beta

sum(tail(sort(CB2C_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2D alpha

sum(tail(sort(CB2D_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2D beta

sum(tail(sort(CB2D_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2E alpha

sum(tail(sort(CB2E_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2E beta

sum(tail(sort(CB2E_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2F alpha

sum(tail(sort(CB2F_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100

#CB2F beta

sum(tail(sort(CB2F_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

#run12
###############################################
#AmAlF6CD8 alpha

sum(tail(sort(AmAlF6CD8_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(AmAlF6CD8_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(AmAlF8CD4_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(AmAlF8CD4_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(AmAlF8CD8_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(AmAlF8CD8_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(AKMCORD_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(AKMCORD_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(AKM2x_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(AKM2x_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(AKM2y_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(AKM2y_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100


#final ones
sum(tail(sort(EKMCORD_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(EKMCORD_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(EKM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(EKM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100


sum(tail(sort(YAMCORD_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(YAMCORD_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100

sum(tail(sort(YAM1_alpha_dcr_cdr3_freq$`Number found in sample`), 100))/100
sum(tail(sort(YAM1_beta_dcr_cdr3_freq$`Number found in sample`), 100))/100
