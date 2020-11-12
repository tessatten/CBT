library(dplyr)
library(immunarch)
library("devtools")
library("ggplot2")
library("vegan")
library("ineq")
library("tidyverse")
library("praise")
library("scales")
library("patchwork")
library("wesanderson")

Relative abundance by clonotype size (Patient A)
Clonotype size category
Sample (patient and month)
Relative abundance of sample (%)

#clono size categottry

p1= repClonality(alphaVDJdataA$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataA$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataD$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataD$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataF$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataF$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataH$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataH$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataI$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataI$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataL$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataL$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataM$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataM$data, "homeo") %>% vis()
p2
fixVis(p2)

p1= repClonality(alphaVDJdataP$data, "homeo") %>% vis()
p1
fixVis(p1)
p2 = repClonality(betaVDJdataP$data, "homeo") %>% vis()
p2
fixVis(p2)

Relative abundance by clonotype size (alpha chain)
Samples grouped by immune reconstitution success of each patient
Immune reconstitution 
Clonotype size category
Relative abundance of sample (%)


reconstPalette <- c("#1C77A3", "#008E90", "#5E2D30", "#C5A387")

immalpha= repClonality(alphaVDJdataNoControls$data, "homeo")
immalpha
p1 = vis(immalpha, .by = c("Immune.reconstitution"), .meta = alphaVDJdataNoControls$meta)+scale_fill_manual(values=reconstPalette)
p1
fixVis(p1)

immabeta= repClonality(betaVDJdataNoControls$data, "homeo")
immabeta
p2 = vis(immabeta, .by = c("Immune.reconstitution"), .meta = betaVDJdataNoControls$meta)+scale_fill_manual(values=reconstPalette)
p2
fixVis(p2)

V gene usage (beta chain)
Samples grouped by immune reconstitution success of each patient
Immune reconstitution 
V gene segment
Usage (proportion of clonotypes)

imm_gu <- geneUsage(betaVDJdataNoControls$data, "hs.trbv", .norm = T, .quant = "count")
gene_stats()
p1 = vis(imm_gu, .by = "Immune.reconstitution", .meta = betaVDJdataNoControls$meta, .plot = "box", .test = "FALSE")+scale_fill_manual(values=reconstPalette)
p1
fixVis(p1)

#okay, in A

patientAPalette <- c("#C3DAEA", "#ECE28B", "#B1D5BB", "#9DAFC3", "#88988D")

V gene usage in Patient A (alpha chain)
Month
V gene segment
Usage (proportion of clonotypes)

imm_gu <- geneUsage(alphaVDJdataA$data[c(1, 2, 3, 4, 5)], "hs.trav", .norm = T, .quant = "count")
p1 = vis(imm_gu)
p1
fixVis(p1)

imm_gu <- geneUsage(alphaVDJdataA$data[c(1, 2, 3, 4, 5)], "hs.traj", .norm = T, .quant = "count")
p2 = vis(imm_gu)
p2
fixVis(p2)

imm_gu <- geneUsage(betaVDJdataA$data[c(1, 2, 3, 4, 5)], "hs.trbv", .norm = T, .quant = "count")
p1 = vis(imm_gu)
p1
fixVis(p1)

imm_gu <- geneUsage(betaVDJdataA$data[c(1, 2, 3, 4, 5)], "hs.trbj", .norm = T, .quant = "count")
p2 = vis(imm_gu)
p2
fixVis(p2)

patientMPalette <- c("#1C77A3", "#008E90", "#5E2D30", "#C5A387", "#C5A387", "#C5A387")

V gene usage in Patient M (alpha chain)
Month
V gene segment
Usage (proportion of clonotypes)

imm_gu <- geneUsage(alphaVDJdataM$data[c(1, 2, 3, 4, 5, 6)], "hs.trav", .norm = T, .quant = "count")
p1 = vis(imm_gu)
p1
fixVis(p1)

imm_gu <- geneUsage(alphaVDJdataM$data[c(1, 2, 3, 4, 5, 6)], "hs.traj", .norm = T, .quant = "count")
p2 = vis(imm_gu)
p2
fixVis(p2)

imm_gu <- geneUsage(betaVDJdataM$data[c(1, 2, 3, 4, 5, 6)], "hs.trbv", .norm = T, .quant = "count")
p1 = vis(imm_gu)
p1
fixVis(p1)

imm_gu <- geneUsage(betaVDJdataM$data[c(1, 2, 3, 4, 5, 6)], "hs.trbj", .norm = T, .quant = "count")
p2 = vis(imm_gu)
p2
fixVis(p2)

