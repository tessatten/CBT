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

# VDJDB - For determining TCR antigen specificity

vdjdb <- read.csv("vdjdb.csv", stringsAsFactors = F, 
                 header = TRUE)
tra <- subset(vdjdb, subset=(vdjdb$gene=="TRA"))
trb <- subset(vdjdb, subset=(vdjdb$gene=="TRB"))

#don't run this willy-nilly! will override good stuff
#################################
write.csv(E_alpha_top10, "E_alpha_top10_all.csv")
write.csv(E_beta_top10, "E_beta_top10_all.csv")

write.csv(F_alpha_top10, "F_alpha_top10_all.csv")
write.csv(F_beta_top10, "F_beta_top10_all.csv")

write.csv(G_alpha_top10, "G_alpha_top10_all.csv")
write.csv(G_beta_top10, "G_beta_top10_all.csv")

write.csv(H_alpha_top10, "H_alpha_top10_all.csv")
write.csv(H_beta_top10, "H_beta_top10_all.csv")

write.csv(I_alpha_top10, "I_alpha_top10_all.csv")
write.csv(I_beta_top10, "I_beta_top10_all.csv")

write.csv(J_alpha_top10, "J_alpha_top10_all.csv")
write.csv(J_beta_top10, "J_beta_top10_all.csv")

write.csv(K_alpha_top10, "K_alpha_top10_all.csv")
write.csv(K_beta_top10, "K_beta_top10_all.csv")

write.csv(L_alpha_top10, "L_alpha_top10_all.csv")
write.csv(L_beta_top10, "L_beta_top10_all.csv")

write.csv(M_alpha_top10, "M_alpha_top10_all.csv")
write.csv(M_beta_top10, "M_beta_top10_all.csv")

write.csv(N_alpha_top10, "N_alpha_top10_all.csv")
write.csv(N_beta_top10, "N_beta_top10_all.csv")

write.csv(O_alpha_top10, "O_alpha_top10_all.csv")
write.csv(O_beta_top10, "O_beta_top10_all.csv")

write.csv(P_alpha_top10, "P_alpha_top10_all.csv")
write.csv(P_beta_top10, "P_beta_top10_all.csv")

write.csv(R_alpha_top10, "R_alpha_top10_all.csv")
write.csv(R_beta_top10, "R_beta_top10_all.csv")

write.csv(Q_alpha_top10, "Q_alpha_top10_all.csv")
write.csv(Q_beta_top10, "Q_beta_top10_all.csv")


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
A_alpha_top10 = A_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

write.csv(A_alpha_top10, "A_alpha_top10_all.csv")
A_alpha_top_by_percent <- read.csv("A_alpha_top10_all.csv", stringsAsFactors = F, 
                  header = TRUE)

# Blast the CDR3 data against the VDJDB
A_alpha_top10$antigen <- NA
A_alpha_top10$score <- NA
A_alpha_top10$antigen[A_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(A_alpha_top10$CDR3s[A_alpha_top10$Chain=='alpha'],tra$cdr3)]
A_alpha_top10$score[A_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(A_alpha_top10$CDR3s[A_alpha_top10$Chain=='alpha'],tra$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(A_alpha_top10, aes(x=A_alpha_top10$Month, y=A_alpha_top10$Count, color=A_alpha_top10$CDR3s, fill=A_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient A: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
                                           limits=c("1","2","3", "4", "5", "6","7", "8", "9", "10","11", "12"))
p
ggsave("patient_A_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(A_alpha_top10$CDR3s)

write.csv(theseuniques, "A_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
A_beta_top10 = A_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

write.csv(A_beta_top10, "A_beta_top10_all.csv")
A_beta_top_by_percent <- read.csv("A_beta_top10_all.csv", stringsAsFactors = F, 
                                   header = TRUE)


# Blast the CDR3 data against the VDJDB
A_beta_top10$antigen <- NA
A_beta_top10$score <- NA
A_beta_top10$antigen[A_beta_top10$Chain=='beta'] <- trb$antigen.species[match(A_beta_top10$CDR3s[A_beta_top10$Chain=='beta'],trb$cdr3)]
A_beta_top10$score[A_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(A_beta_top10$CDR3s[A_beta_top10$Chain=='beta'],trb$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(A_beta_top10, aes(x=A_beta_top10$Month, y=A_beta_top10$Count, color=A_beta_top10$CDR3s, fill=A_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient A: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9", "10","11", "12"))
p
ggsave("patient_A_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(A_beta_top10$CDR3s)
theseuniques
write.csv(theseuniques, "A_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
B_alpha_top10 = B_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

write.csv(B_alpha_top10, "B_alpha_top10_all.csv")
B_alpha_top_by_percent <- read.csv("B_alpha_top10_all.csv", stringsAsFactors = F, 
                                   header = TRUE)

# Blast the CDR3 data against the VDJDB
B_alpha_top10$antigen <- NA
B_alpha_top10$score <- NA
B_alpha_top10$antigen[B_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(B_alpha_top10$CDR3s[B_alpha_top10$Chain=='alpha'],tra$cdr3)]
B_alpha_top10$score[B_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(B_alpha_top10$CDR3s[B_alpha_top10$Chain=='alpha'],tra$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#b_alpha top10
p = ggplot(B_alpha_top10, aes(x=B_alpha_top10$Month, y=B_alpha_top10$Count, color=B_alpha_top10$CDR3s, fill=B_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient B: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2"))
p
ggsave("patient_B_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(B_alpha_top10$CDR3s)

write.csv(theseuniques, "B_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
B_beta_top10 = B_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

write.csv(B_beta_top10, "B_beta_top10_all.csv")
B_beta_top_by_percent <- read.csv("B_beta_top10_all.csv", stringsAsFactors = F, 
                                  header = TRUE)


# Blast the CDR3 data against the VDJDB
B_beta_top10$antigen <- NA
B_beta_top10$score <- NA
B_beta_top10$antigen[B_beta_top10$Chain=='beta'] <- trb$antigen.species[match(B_beta_top10$CDR3s[B_beta_top10$Chain=='beta'],trb$cdr3)]
B_beta_top10$score[B_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(B_beta_top10$CDR3s[B_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#b_beta top10
p = ggplot(B_beta_top10, aes(x=B_beta_top10$Month, y=B_beta_top10$Count, color=B_beta_top10$CDR3s, fill=B_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient B: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2"))
p
ggsave("patient_B_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(B_beta_top10$CDR3s)
theseuniques
write.csv(theseuniques, "B_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
C_alpha_top10 = C_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)
#C_alpha_top10[C_alpha_top10=="0"]<-"1"

write.csv(C_alpha_top10, "C_alpha_top10_all.csv")
C_alpha_top_by_percent <- read.csv("C_alpha_top10_all.csv", stringsAsFactors = F, 
                                   header = TRUE)

# Blast the CDR3 data against the VDJDB
C_alpha_top10$antigen <- NA
C_alpha_top10$score <- NA
C_alpha_top10$antigen[C_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(C_alpha_top10$CDR3s[C_alpha_top10$Chain=='alpha'],tra$cdr3)]
C_alpha_top10$score[C_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(C_alpha_top10$CDR3s[C_alpha_top10$Chain=='alpha'],tra$cdr3)]


C_alpha_top102 <- C_alpha_top10 %>%
  group_by(CDR3s) %>%
  summarize(Count = sum(Count))

praise("Hello, you are ${adjective}! This code is ${adjective}!")

#C_alpha top10
p = ggplot(C_alpha_top10, aes(x=C_alpha_top10$Month, y=C_alpha_top10$Count, color=C_alpha_top10$CDR3s, fill=C_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  labs(title = "Patient C: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE) + scale_x_discrete(limits=c("0", "1","2"))
p
ggsave("patient_C_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(C_alpha_top10$CDR3s)

write.csv(theseuniques, "C_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
C_beta_top10 = C_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

write.csv(C_beta_top10, "C_beta_top10_all.csv")
C_beta_top_by_percent <- read.csv("C_beta_top10_all.csv", stringsAsFactors = F, 
                                  header = TRUE)

# Blast the CDR3 data against the VDJDB
C_beta_top10$antigen <- NA
C_beta_top10$score <- NA
C_beta_top10$antigen[C_beta_top10$Chain=='beta'] <- trb$antigen.species[match(C_beta_top10$CDR3s[C_beta_top10$Chain=='beta'],trb$cdr3)]
C_beta_top10$score[C_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(C_beta_top10$CDR3s[C_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#c_beta top10
p = ggplot(C_beta_top10, aes(x=C_beta_top10$Month, y=C_beta_top10$Count, color=C_beta_top10$CDR3s, fill=C_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
   labs(title = "Patient C: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("0","1", "2"))
p
ggsave("patient_C_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(C_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)
write.csv(theseuniques, "C_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
D_alpha_top10 = D_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)
write.csv(D_alpha_top10, "D_alpha_top10_all.csv")
# Blast the CDR3 data against the VDJDB
D_alpha_top10$antigen <- NA
D_alpha_top10$score <- NA
D_alpha_top10$antigen[D_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(D_alpha_top10$CDR3s[D_alpha_top10$Chain=='alpha'],tra$cdr3)]
D_alpha_top10$score[D_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(D_alpha_top10$CDR3s[D_alpha_top10$Chain=='alpha'],tra$cdr3)]

praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(D_alpha_top10, aes(x=D_alpha_top10$Month, y=D_alpha_top10$Count, color=D_alpha_top10$CDR3s, fill=D_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient D: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8"))
p
ggsave("patient_D_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(D_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "D_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
D_beta_top10 = D_beta %>%
  group_by(Month) %>%
  top_n(10, Count)
write.csv(D_beta_top10, "D_beta_top10_all.csv")
# Blast the CDR3 data against the VDJDB
D_beta_top10$antigen <- NA
D_beta_top10$score <- NA
D_beta_top10$antigen[D_beta_top10$Chain=='beta'] <- trb$antigen.species[match(D_beta_top10$CDR3s[D_beta_top10$Chain=='beta'],trb$cdr3)]
D_beta_top10$score[D_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(D_beta_top10$CDR3s[D_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(D_beta_top10, aes(x=D_beta_top10$Month, y=D_beta_top10$Count, color=D_beta_top10$CDR3s, fill=D_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient D: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8"))
p
ggsave("patient_D_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(D_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "D_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
E_alpha_top10 = E_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)
write.csv(E_alpha_top10, "E_alpha_top10_all.csv")

# Blast the CDR3 data against the VDJDB
E_alpha_top10$antigen <- NA
E_alpha_top10$score <- NA
E_alpha_top10$antigen[E_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(E_alpha_top10$CDR3s[E_alpha_top10$Chain=='alpha'],tra$cdr3)]
E_alpha_top10$score[E_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(E_alpha_top10$CDR3s[E_alpha_top10$Chain=='alpha'],tra$cdr3)]

praise("Hello, you are ${adjective}! This code is ${adjective}!")

praise("Hello, you are ${adjective}! This code is ${adjective}!")

#alpha top10
p = ggplot(E_alpha_top10, aes(x=E_alpha_top10$Month, y=E_alpha_top10$Count, color=E_alpha_top10$CDR3s, fill=E_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient E: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3"))
p
ggsave("patient_E_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(E_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "E_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
E_beta_top10 = E_beta %>%
  group_by(Month) %>%
  top_n(10, Count)
write.csv(E_beta_top10, "E_beta_top10_all.csv")

# Blast the CDR3 data against the VDJDB
E_beta_top10$antigen <- NA
E_beta_top10$score <- NA
E_beta_top10$antigen[E_beta_top10$Chain=='beta'] <- trb$antigen.species[match(E_beta_top10$CDR3s[E_beta_top10$Chain=='beta'],trb$cdr3)]
E_beta_top10$score[E_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(E_beta_top10$CDR3s[E_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(E_beta_top10, aes(x=E_beta_top10$Month, y=E_beta_top10$Count, color=E_beta_top10$CDR3s, fill=E_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient E: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3"))
p
ggsave("patient_E_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(E_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "E_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
F_alpha_top10 = F_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
F_alpha_top10$antigen <- NA
F_alpha_top10$score <- NA
F_alpha_top10$antigen[F_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(F_alpha_top10$CDR3s[F_alpha_top10$Chain=='alpha'],tra$cdr3)]
F_alpha_top10$score[F_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(F_alpha_top10$CDR3s[F_alpha_top10$Chain=='alpha'],tra$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(F_alpha_top10, aes(x=F_alpha_top10$Month, y=F_alpha_top10$Count, color=F_alpha_top10$CDR3s, fill=F_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient F: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19"))
p
ggsave("patient_F_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(F_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "F_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
F_beta_top10 = F_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
F_beta_top10$antigen <- NA
F_beta_top10$score <- NA
F_beta_top10$antigen[F_beta_top10$Chain=='beta'] <- trb$antigen.species[match(F_beta_top10$CDR3s[F_beta_top10$Chain=='beta'],trb$cdr3)]
F_beta_top10$score[F_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(F_beta_top10$CDR3s[F_beta_top10$Chain=='beta'],trb$cdr3)]

praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(F_beta_top10, aes(x=F_beta_top10$Month, y=F_beta_top10$Count, color=F_beta_top10$CDR3s, fill=F_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient F: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19"))
p
ggsave("patient_F_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(F_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "F_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
G_alpha_top10 = G_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
G_alpha_top10$antigen <- NA
G_alpha_top10$score <- NA
G_alpha_top10$antigen[G_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(G_alpha_top10$CDR3s[G_alpha_top10$Chain=='alpha'],tra$cdr3)]
G_alpha_top10$score[G_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(G_alpha_top10$CDR3s[G_alpha_top10$Chain=='alpha'],tra$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(G_alpha_top10, aes(x=G_alpha_top10$Month, y=G_alpha_top10$Count, color=G_alpha_top10$CDR3s, fill=G_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient G: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6"))
p
ggsave("patient_G_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(G_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "G_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
G_beta_top10 = G_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
G_beta_top10$antigen <- NA
G_beta_top10$score <- NA
G_beta_top10$antigen[G_beta_top10$Chain=='beta'] <- trb$antigen.species[match(G_beta_top10$CDR3s[G_beta_top10$Chain=='beta'],trb$cdr3)]
G_beta_top10$score[G_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(G_beta_top10$CDR3s[G_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(G_beta_top10, aes(x=G_beta_top10$Month, y=G_beta_top10$Count, color=G_beta_top10$CDR3s, fill=G_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient G: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6"))
p
ggsave("patient_G_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(G_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "G_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
H_alpha_top10 = H_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
H_alpha_top10$antigen <- NA
H_alpha_top10$score <- NA
H_alpha_top10$antigen[H_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(H_alpha_top10$CDR3s[H_alpha_top10$Chain=='alpha'],tra$cdr3)]
H_alpha_top10$score[H_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(H_alpha_top10$CDR3s[H_alpha_top10$Chain=='alpha'],tra$cdr3)]

praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(H_alpha_top10, aes(x=H_alpha_top10$Month, y=H_alpha_top10$Count, color=H_alpha_top10$CDR3s, fill=H_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient H: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18"))
p
ggsave("patient_H_alpha_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(H_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "H_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
H_beta_top10 = H_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
H_beta_top10$antigen <- NA
H_beta_top10$score <- NA
H_beta_top10$antigen[H_beta_top10$Chain=='beta'] <- trb$antigen.species[match(H_beta_top10$CDR3s[H_beta_top10$Chain=='beta'],trb$cdr3)]
H_beta_top10$score[H_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(H_beta_top10$CDR3s[H_beta_top10$Chain=='beta'],trb$cdr3)]

praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(H_beta_top10, aes(x=H_beta_top10$Month, y=H_beta_top10$Count, color=H_beta_top10$CDR3s, fill=H_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient H: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18"))
p
ggsave("patient_H_beta_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(H_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "H_beta_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
I_alpha_top10 = I_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
I_alpha_top10$antigen <- NA
I_alpha_top10$score <- NA
I_alpha_top10$antigen[I_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(I_alpha_top10$CDR3s[I_alpha_top10$Chain=='alpha'],tra$cdr3)]
I_alpha_top10$score[I_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(I_alpha_top10$CDR3s[I_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(I_alpha_top10, aes(x=I_alpha_top10$Month, y=I_alpha_top10$Count, color=I_alpha_top10$CDR3s, fill=I_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient I: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12"))
p
ggsave("patient_I_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(I_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(I_alpha_top10, "I_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
I_beta_top10 = I_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
I_beta_top10$antigen <- NA
I_beta_top10$score <- NA
I_beta_top10$antigen[I_beta_top10$Chain=='beta'] <- trb$antigen.species[match(I_beta_top10$CDR3s[I_beta_top10$Chain=='beta'],trb$cdr3)]
I_beta_top10$score[I_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(I_beta_top10$CDR3s[I_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(I_beta_top10, aes(x=I_beta_top10$Month, y=I_beta_top10$Count, color=I_beta_top10$CDR3s, fill=I_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient I: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12"))
p
ggsave("patient_I_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(I_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(I_beta_top10, "I_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
J_alpha_top10 = J_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
J_alpha_top10$antigen <- NA
J_alpha_top10$score <- NA
J_alpha_top10$antigen[J_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(J_alpha_top10$CDR3s[J_alpha_top10$Chain=='alpha'],tra$cdr3)]
J_alpha_top10$score[J_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(J_alpha_top10$CDR3s[J_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(J_alpha_top10, aes(x=J_alpha_top10$Month, y=J_alpha_top10$Count, color=J_alpha_top10$CDR3s, fill=J_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient J: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3"))
p
ggsave("patient_J_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(J_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "J_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
J_beta_top10 = J_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
J_beta_top10$antigen <- NA
J_beta_top10$score <- NA
J_beta_top10$antigen[J_beta_top10$Chain=='beta'] <- trb$antigen.species[match(J_beta_top10$CDR3s[J_beta_top10$Chain=='beta'],trb$cdr3)]
J_beta_top10$score[J_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(J_beta_top10$CDR3s[J_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(J_beta_top10, aes(x=J_beta_top10$Month, y=J_beta_top10$Count, color=J_beta_top10$CDR3s, fill=J_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient J: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3"))
p
ggsave("patient_J_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(J_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "J_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
K_alpha_top10 = K_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
K_alpha_top10$antigen <- NA
K_alpha_top10$score <- NA
K_alpha_top10$antigen[K_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(K_alpha_top10$CDR3s[K_alpha_top10$Chain=='alpha'],tra$cdr3)]
K_alpha_top10$score[K_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(K_alpha_top10$CDR3s[K_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(K_alpha_top10, aes(x=K_alpha_top10$Month, y=K_alpha_top10$Count, color=K_alpha_top10$CDR3s, fill=K_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient K: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3"))
p
ggsave("patient_K_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(K_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "K_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
K_beta_top10 = K_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
K_beta_top10$antigen <- NA
K_beta_top10$score <- NA
K_beta_top10$antigen[K_beta_top10$Chain=='beta'] <- trb$antigen.species[match(K_beta_top10$CDR3s[K_beta_top10$Chain=='beta'],trb$cdr3)]
K_beta_top10$score[K_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(K_beta_top10$CDR3s[K_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(K_beta_top10, aes(x=K_beta_top10$Month, y=K_beta_top10$Count, color=K_beta_top10$CDR3s, fill=K_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient K: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3"))
p
ggsave("patient_K_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(K_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "K_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
L_alpha_top10 = L_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
L_alpha_top10$antigen <- NA
L_alpha_top10$score <- NA
L_alpha_top10$antigen[L_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(L_alpha_top10$CDR3s[L_alpha_top10$Chain=='alpha'],tra$cdr3)]
L_alpha_top10$score[L_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(L_alpha_top10$CDR3s[L_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(L_alpha_top10, aes(x=L_alpha_top10$Month, y=L_alpha_top10$Count, color=L_alpha_top10$CDR3s, fill=L_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient L: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19","20", "21","22"))
p
ggsave("patient_L_alpha_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(L_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "L_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
L_beta_top10 = L_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
L_beta_top10$antigen <- NA
L_beta_top10$score <- NA
L_beta_top10$antigen[L_beta_top10$Chain=='beta'] <- trb$antigen.species[match(L_beta_top10$CDR3s[L_beta_top10$Chain=='beta'],trb$cdr3)]
L_beta_top10$score[L_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(L_beta_top10$CDR3s[L_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a beta top10
p = ggplot(L_beta_top10, aes(x=L_beta_top10$Month, y=L_beta_top10$Count, color=L_beta_top10$CDR3s, fill=L_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient L: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19","20", "21","22"))
p
ggsave("patient_L_beta_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(L_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "L_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
M_alpha_top10 = M_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
M_alpha_top10$antigen <- NA
M_alpha_top10$score <- NA
M_alpha_top10$antigen[M_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(M_alpha_top10$CDR3s[M_alpha_top10$Chain=='alpha'],tra$cdr3)]
M_alpha_top10$score[M_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(M_alpha_top10$CDR3s[M_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}!")

#a alpha top10
p = ggplot(M_alpha_top10, aes(x=M_alpha_top10$Month, y=M_alpha_top10$Count, color=M_alpha_top10$CDR3s, fill=M_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient M: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19","20", "21","22"))
p
ggsave("patient_M_alpha_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(M_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "M_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
M_beta_top10 = M_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
M_beta_top10$antigen <- NA
M_beta_top10$score <- NA
M_beta_top10$antigen[M_beta_top10$Chain=='beta'] <- trb$antigen.species[match(M_beta_top10$CDR3s[M_beta_top10$Chain=='beta'],trb$cdr3)]
M_beta_top10$score[M_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(M_beta_top10$CDR3s[M_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a beta top10
p = ggplot(M_beta_top10, aes(x=M_beta_top10$Month, y=M_beta_top10$Count, color=M_beta_top10$CDR3s, fill=M_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient M: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19","20", "21","22"))
p
ggsave("patient_M_beta_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(M_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "M_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
N_alpha_top10 = N_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
N_alpha_top10$antigen <- NA
N_alpha_top10$score <- NA
N_alpha_top10$antigen[N_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(N_alpha_top10$CDR3s[N_alpha_top10$Chain=='alpha'],tra$cdr3)]
N_alpha_top10$score[N_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(N_alpha_top10$CDR3s[N_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a alpha top10
p = ggplot(N_alpha_top10, aes(x=N_alpha_top10$Month, y=N_alpha_top10$Count, color=N_alpha_top10$CDR3s, fill=N_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient N: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12"))
p
ggsave("patient_N_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(N_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "N_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
N_beta_top10 = N_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
N_beta_top10$antigen <- NA
N_beta_top10$score <- NA
N_beta_top10$antigen[N_beta_top10$Chain=='beta'] <- trb$antigen.species[match(N_beta_top10$CDR3s[N_beta_top10$Chain=='beta'],trb$cdr3)]
N_beta_top10$score[N_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(N_beta_top10$CDR3s[N_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a beta top10
p = ggplot(N_beta_top10, aes(x=N_beta_top10$Month, y=N_beta_top10$Count, color=N_beta_top10$CDR3s, fill=N_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient N: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12"))
p
ggsave("patient_N_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(N_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "N_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
O_alpha_top10 = O_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
O_alpha_top10$antigen <- NA
O_alpha_top10$score <- NA
O_alpha_top10$antigen[O_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(O_alpha_top10$CDR3s[O_alpha_top10$Chain=='alpha'],tra$cdr3)]
O_alpha_top10$score[O_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(O_alpha_top10$CDR3s[O_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a alpha top10
p = ggplot(O_alpha_top10, aes(x=O_alpha_top10$Month, y=O_alpha_top10$Count, color=O_alpha_top10$CDR3s, fill=O_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient O: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12"))
p
ggsave("patient_O_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(O_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "O_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
O_beta_top10 = O_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
O_beta_top10$antigen <- NA
O_beta_top10$score <- NA
O_beta_top10$antigen[O_beta_top10$Chain=='beta'] <- trb$antigen.species[match(O_beta_top10$CDR3s[O_beta_top10$Chain=='beta'],trb$cdr3)]
O_beta_top10$score[O_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(O_beta_top10$CDR3s[O_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a beta top10
p = ggplot(O_beta_top10, aes(x=O_beta_top10$Month, y=O_beta_top10$Count, color=O_beta_top10$CDR3s, fill=O_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient O: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12"))
p
ggsave("patient_O_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(O_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "O_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
P_alpha_top10 = P_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
P_alpha_top10$antigen <- NA
P_alpha_top10$score <- NA
P_alpha_top10$antigen[P_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(P_alpha_top10$CDR3s[P_alpha_top10$Chain=='alpha'],tra$cdr3)]
P_alpha_top10$score[P_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(P_alpha_top10$CDR3s[P_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a alpha top10
p = ggplot(P_alpha_top10, aes(x=P_alpha_top10$Month, y=P_alpha_top10$Count, color=P_alpha_top10$CDR3s, fill=P_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient P: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19","20", "21","22","23", "24", "25", "26","27", "28", "29","30"))
p
ggsave("patient_P_alpha_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(P_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "P_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
P_beta_top10 = P_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
P_beta_top10$antigen <- NA
P_beta_top10$score <- NA
P_beta_top10$antigen[P_beta_top10$Chain=='beta'] <- trb$antigen.species[match(P_beta_top10$CDR3s[P_beta_top10$Chain=='beta'],trb$cdr3)]
P_beta_top10$score[P_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(P_beta_top10$CDR3s[P_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a beta top10
p = ggplot(P_beta_top10, aes(x=P_beta_top10$Month, y=P_beta_top10$Count, color=P_beta_top10$CDR3s, fill=P_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient P: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("1","2","3", "4", "5", "6","7", "8", "9","10", "11","12","13", "14", "15", "16","17", "18", "19","20", "21","22","23", "24", "25", "26","27", "28", "29","30"))
p
ggsave("patient_P_beta_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(P_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "P_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
Q_alpha_top10 = Q_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
Q_alpha_top10$antigen <- NA
Q_alpha_top10$score <- NA
Q_alpha_top10$antigen[Q_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(Q_alpha_top10$CDR3s[Q_alpha_top10$Chain=='alpha'],tra$cdr3)]
Q_alpha_top10$score[Q_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(Q_alpha_top10$CDR3s[Q_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a alpha top10
p = ggplot(Q_alpha_top10, aes(x=Q_alpha_top10$Month, y=Q_alpha_top10$Count, color=Q_alpha_top10$CDR3s, fill=Q_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient Q: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("0","0.75"))
p
ggsave("patient_Q_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(Q_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "Q_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
Q_beta_top10 = Q_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
Q_beta_top10$antigen <- NA
Q_beta_top10$score <- NA
Q_beta_top10$antigen[Q_beta_top10$Chain=='beta'] <- trb$antigen.species[match(Q_beta_top10$CDR3s[Q_beta_top10$Chain=='beta'],trb$cdr3)]
Q_beta_top10$score[Q_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(Q_beta_top10$CDR3s[Q_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a beta top10
p = ggplot(Q_beta_top10, aes(x=Q_beta_top10$Month, y=Q_beta_top10$Count, color=Q_beta_top10$CDR3s, fill=Q_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient Q: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("0","0.75"))
p
ggsave("patient_Q_beta_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(Q_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "Q_beta_top10.csv")

#################################################


#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
R_alpha_top10 = R_alpha %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
R_alpha_top10$antigen <- NA
R_alpha_top10$score <- NA
R_alpha_top10$antigen[R_alpha_top10$Chain=='alpha'] <- tra$antigen.species[match(R_alpha_top10$CDR3s[R_alpha_top10$Chain=='alpha'],tra$cdr3)]
R_alpha_top10$score[R_alpha_top10$Chain=='alpha'] <- tra$vdjdb.score[match(R_alpha_top10$CDR3s[R_alpha_top10$Chain=='alpha'],tra$cdr3)]



praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a alpha top10
p = ggplot(R_alpha_top10, aes(x=R_alpha_top10$Month, y=R_alpha_top10$Count, color=R_alpha_top10$CDR3s, fill=R_alpha_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient R: top 10 CDR3s (alpha chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 26),
    legend.text = element_text(color = "black", size = 24)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("0", "1"))
p
ggsave("patient_R_alpha_top10.pdf", width = 47, height = 20, units = c('cm'))

theseuniques <- unique(R_alpha_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "R_alpha_top10.csv")

#################################################

#sekect the top ten clonotypes (abundance over time) and plot - so some are found in multiple points and some are are not
R_beta_top10 = R_beta %>%
  group_by(Month) %>%
  top_n(10, Count)

# Blast the CDR3 data against the VDJDB
R_beta_top10$antigen <- NA
R_beta_top10$score <- NA
R_beta_top10$antigen[R_beta_top10$Chain=='beta'] <- trb$antigen.species[match(R_beta_top10$CDR3s[R_beta_top10$Chain=='beta'],trb$cdr3)]
R_beta_top10$score[R_beta_top10$Chain=='beta'] <- trb$vdjdb.score[match(R_beta_top10$CDR3s[R_beta_top10$Chain=='beta'],trb$cdr3)]


praise("Hello, you are ${adjective}! This code is ${adjective}! Your sister Ros is ${adjective}!")

#a beta top10
p = ggplot(R_beta_top10, aes(x=R_beta_top10$Month, y=R_beta_top10$Count, color=R_beta_top10$CDR3s, fill=R_beta_top10$CDR3s)) +
  geom_point(size = 10)+ theme(legend.position="none")+
  geom_line(aes(color = CDR3s), size = 3) + labs(title = "Patient R: top 10 CDR3s (beta chain)", x = "Month after transplant", y = "Clonal expansion size", color = "CDR3") +
  theme_bw() + theme(
    plot.title = element_text(color="black", size=28, face="bold"),
    axis.title.x = element_text(color="black", size=28, face="bold"),
    axis.title.y = element_text(color="black", size=28, face="bold"),
    axis.text=element_text(size=26),
    legend.title = element_text(color = "black", size = 16),
    legend.text = element_text(color = "black", size = 14)
  ) + guides(fill=FALSE)+ scale_x_discrete( 
    limits=c("0", "1"))
p
ggsave("patient_R_beta_top10.pdf", width = 52, height = 20, units = c('cm'))

theseuniques <- unique(R_beta_top10$CDR3s)
theseuniques <-sort(theseuniques)

write.csv(theseuniques, "R_beta_top10.csv")

#################################################

most_abundant_alpha <- read.csv("most_abundant_taxa_alpha.csv", stringsAsFactors = F, 
                  header = TRUE)
most_abundant_alpha$antigen <- NA
most_abundant_alpha$score <- NA
most_abundant_alpha$antigen[most_abundant_alpha$Chain=='alpha'] <- trb$antigen.species[match(most_abundant_alpha$CDR3s[most_abundant_alpha$Chain=='alpha'],trb$cdr3)]
most_abundant_alpha$score[most_abundant_alpha$Chain=='alpha'] <- trb$vdjdb.score[match(most_abundant_alpha$CDR3s[most_abundant_alpha$Chain=='alpha'],trb$cdr3)]

most_abundant_beta <- read.csv("most_abundant_taxa_beta.csv", stringsAsFactors = F, 
                               header = TRUE)
most_abundant_beta$antigen <- NA
most_abundant_beta$score <- NA
most_abundant_beta$antigen[most_abundant_beta$Chain=='beta'] <- trb$antigen.species[match(most_abundant_beta$CDR3s[most_abundant_beta$Chain=='beta'],trb$cdr3)]
most_abundant_beta$score[most_abundant_beta$Chain=='beta'] <- trb$vdjdb.score[match(most_abundant_beta$CDR3s[most_abundant_beta$Chain=='beta'],trb$cdr3)]