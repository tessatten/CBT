# Load libraries

library(data.table)
library(zoo)
library(plyr)
library(scales)
library(ggplot2)
library(ggthemes)
library(dplyr)
library(plotly)
library(pracma)
library(Hmisc)
library(rpart)
library(survminer)
library(knitr)
library(grid)
library(stringr)
library(ggrepel)
library(ggjoy)
library(RColorBrewer)
library(stringdist)
library(ape)
library(reshape2)
library(gplots)
library(ineq)
library(vegan)
library(Rtsne)
library(praise)

tigar = P_alpha
#ID, month, chain, CDR3, freq
ID <- "P"
weekPost <- tigar$Month
chain <- "a"
V1 <- tigar$CDR3s #CDR3s
V2 <- tigar$Count #Freqs

dat <- data.frame(ID, weekPost, chain, V1, V2, stringsAsFactors = F)

# % shared clonotypes time

#  y | x1 | x2 | shared |
# ---|----|----|--------|
#  1 |  1 |  2 |   50   |
#  1 |  2 |  3 |   50   |
#  1 |  1 |  3 |   25   |


shared <- data.frame()

chains <- c('a', 'b')
for (chain in chains){
  for (i in unique(dat$ID)){
    weekPosts <- unique(dat$weekPost[dat$ID == i & dat$chain == chain])
    weekPosts <- weekPosts[order(weekPosts)]
    
    # loop of t1s
    for (j in seq(2, length(weekPosts))){
      
      # loop of t0s to t1
      for (k in seq(1, j-1)){
        
        if (i == 'P'){ID <- 1}
        if (i == 'L'){ID <- 2}
        if (i == 'Q'){ID <- 3}
        if (i == 'T'){ID <- 4}
        
        temp1 <- dat[dat$ID == i & dat$weekPost == weekPosts[k] & dat$chain == chain,] #t0
        temp2 <- dat[dat$ID == i & dat$weekPost == weekPosts[j] & dat$chain == chain,] #t1
        temp3 <- temp1[temp1$V1 %in% temp2$V1,] # overlap
        
        if (chain == 'a'){tempChain <- 1}
        if (chain == 'b'){tempChain <- 2}
        
        row <- c(ID, weekPosts[k], weekPosts[j], ((100/length(temp1$V1)) * length(temp1$V1[temp1$V1 %in% temp2$V1])), ((100/length(temp2$V1)) * length(temp2$V1[temp2$V1 %in% temp1$V1])), length(temp1$V1), length(temp2$V1), tempChain)
        
        shared <- rbind(shared, row)
        names(shared) <- c('ID', 'x1', 'x2', 'shared', 'sharedTo', 'x1clonotypes', 'x2clonotypes', 'chain')
        
        ### Shared clonotype plot ###
        # creates a vector identifying shared sequences
        temp1$shared <- NA
        temp2$shared <- NA
        temp1$shared[temp1$V1 %in% temp2$V1] <- T
        temp1$shared[!temp1$V1 %in% temp2$V1] <- F
        temp2$shared[temp2$V1 %in% temp1$V1] <- T
        temp2$shared[!temp2$V1 %in% temp1$V1] <- F
        
        # orders the dataset by the shared sequences, in order of clonal expansion size
        temp1 <- temp1[order(-temp1$shared, -temp1$V2),]
        
        # vector to order and colour sequence
        fillseq <- seq(1, length(temp1$ID[temp1$shared == T]))
        temp1$fillseq <- NA
        temp2$fillseq <- NA
        temp1$fillseq[temp1$shared == T] <- fillseq
        temp1$fillseq[temp1$shared == F] <- length(temp1$fillseq[temp1$shared == T])+1
        temp2$fillseq <- temp1$fillseq[match(temp2$V1, temp1$V1)]
        temp2$fillseq[temp2$shared == F] <- length(temp2$fillseq[temp2$shared == T])+1
        
        # orders the dataset by the shared sequences, in order of clonal expansion size
        temp2 <- temp2[order(temp2$fillseq),]
        
        # create 1 entry for all non-shared sequences
        temp1.mod <- temp1
        temp2.mod <- temp2
        temp1.mod$V2[length(temp1.mod$ID)] <- sum(temp1.mod$V2[temp1.mod$shared == F])
        temp2.mod$V2[length(temp2.mod$ID)] <- sum(temp2.mod$V2[temp2.mod$shared == F])
        temp1.mod$shared[length(temp1.mod$ID)] <- T
        temp2.mod$shared[length(temp2.mod$ID)] <- T
        temp1.mod <- temp1.mod[!temp1.mod$shared == F,]
        temp2.mod <- temp2.mod[!temp2.mod$shared == F,]
        temp1.mod$shared[length(temp1.mod$ID)] <- F
        temp2.mod$shared[length(temp2.mod$ID)] <- F
        
        # normalise sequence frequency
        temp1.mod$V2 <- (100/sum(temp1.mod$V2))*temp1.mod$V2 
        temp2.mod$V2 <- (100/sum(temp2.mod$V2))*temp2.mod$V2
        
        #temp1.mod.expanded <- temp1.mod[rep(row.names(temp1.mod), temp1.mod$V2),]
        #temp2.mod.expanded <- temp2.mod[rep(row.names(temp2.mod), temp2.mod$V2),]
        
        maxPcntShared <- max(c(100 - sum(temp1.mod$V2[temp1.mod$shared==T]), 100 - sum(temp2.mod$V2[temp2.mod$shared==T])))
        minPcntShared <- min(c(100 - sum(temp1.mod$V2[temp1.mod$shared==T]), 100 - sum(temp2.mod$V2[temp2.mod$shared==T])))
        
        # create x coordinates for line plotting
        temp1.mod$x <- 1.45
        temp2.mod$x <- 2.55
        
        # create the y point for the plot
        temp1.mod$y <- 100 - cumsum(temp1.mod$V2) + temp1.mod$V2
        temp2.mod$y <- 100 - cumsum(temp2.mod$V2) + temp2.mod$V2
        
        # merge the two datasets for plotting
        temp3.mod <- rbind(temp1.mod, temp2.mod)
        
        yadjst <- 8
        sampleName1 <- paste(unique(temp1.mod$ID), unique(temp1.mod$weekPost), unique(temp1.mod$chain), sep = '-')
        sampleName2 <- paste(unique(temp2.mod$ID), unique(temp2.mod$weekPost), unique(temp2.mod$chain), sep = '-')
        
        ggp <- ggplot()+
          geom_line(aes(x = temp3.mod$x, y = temp3.mod$y, group = temp3.mod$fillseq), colour = 'grey')+
          geom_bar(aes(x = 1, temp1.mod$V2, fill = as.factor(temp1.mod$fillseq), alpha = as.factor(temp1.mod$shared)), stat = 'identity', colour = 'black') +
          geom_errorbar(aes(x = 1, ymin = round(minPcntShared-yadjst), ymax = round(minPcntShared-yadjst)))+
          geom_bar(aes(x = 3, temp2.mod$V2, fill = as.factor(temp2.mod$fillseq), alpha = as.factor(temp2.mod$shared)), stat = 'identity', colour = 'black') + 
          scale_alpha_manual(values=c(0, 1)) +
          geom_errorbar(aes(x = 3, ymin = round(minPcntShared-yadjst), ymax = round(minPcntShared-yadjst)))+
          scale_y_continuous(breaks = c(round(minPcntShared-yadjst), seq(round(minPcntShared, -1), 100, by = 10)), labels = c(0,  seq(round(minPcntShared, -1), 100, by = 10)), expand = c(0,0))+
          scale_x_continuous(breaks = c(1, 3), labels = c(sampleName1, sampleName2))+
          geom_text(aes(x = 1, y = round(minPcntShared-(yadjst/2))), label = 'Not Shared')+
          geom_text(aes(x = 3, y = round(minPcntShared-(yadjst/2))), label = 'Not Shared')+
          annotate(geom = "segment", x = -Inf, xend = -Inf, y = -Inf, yend = Inf) +
          annotate(geom = "segment", x = -Inf, xend = -Inf, y =  round(minPcntShared-yadjst), yend = round(minPcntShared, -1),
                   linetype = "dashed", color = "white")+
          coord_cartesian(ylim=c(round(minPcntShared-yadjst), 100), expand = T) +
          labs(title = 'Shared Clonotypes', x = ' ', y = 'Cumulative % of Sample')+
          theme(legend.position="none", panel.background = element_blank()) + coord_flip()
        
        ggsave(paste(paste(sampleName1, sampleName2, sep = '_'), '.pdf', sep = ''), 
               ggp, width = 17, height = 10, units = c('cm'))
        
      }
    }
  }
}
shared$delta <- shared$x2 - shared$x1
shared$chain[shared$chain == 1] <- 'a'
shared$chain[shared$chain == 2] <- 'b'




# Shared clonotypes over time


cutoff1 <- 6
cutoff2 <- 20
shared.mod <- shared
shared.mod$group[shared.mod$delta <= cutoff1 | shared.mod$delta >= cutoff2] <- 1
shared.mod$group[shared.mod$delta > cutoff1 & shared.mod$delta < cutoff2] <- 2


for (i in unique(shared.mod$ID)){
  for (j in unique(shared.mod$chain)){
    temp <- shared.mod[shared.mod$ID == i & shared.mod$chain == j,] 
    
    if (i == 'P'){ID <- 1}
    if (i == 'L'){ID <- 2}
    if (i == 'Q'){ID <- 3}
    if (i == 'T'){ID <- 4}
    
    tempx <- rep(unique(c(temp$x1, temp$x2))[order(unique(c(temp$x1, temp$x2)))], length(unique(c(temp$x1, temp$x2))[order(unique(c(temp$x1, temp$x2)))]))
    tempy <- as.factor(rep(unique(c(temp$x1, temp$x2))[order(unique(c(temp$x1, temp$x2)))], each = length(unique(c(temp$x1, temp$x2))[order(unique(c(temp$x1, temp$x2)))])))
    
    
    ggp <- ggplot()+
      
      geom_point(aes(x = tempx, y = tempy), size = 5)+
      geom_line(aes(x = tempx, y = tempy))+
      
      geom_curve(data = temp, aes(x = x1, y = as.factor(x1), xend = x2, yend = as.factor(x1), 
                                  group = ID, colour = shared), curvature = -0.2, arrow = arrow(length = unit(0.01, "npc")))+
      geom_curve(data = temp, aes(x = x1, y = as.factor(x2), xend = x2, yend = as.factor(x2), 
                                  group = ID, colour = shared), curvature = 0.2, arrow = arrow(length = unit(0.01, "npc")))+
      
      geom_text(aes(x = unique(tempx), y = as.factor(unique(c(temp$x1, temp$x2))), 
                    label = paste(c(unique(temp$x1clonotypes), 
                                    unique(temp$x2clonotypes[temp$x2 == max(temp$x2)])), 'Clonotypes')), 
                nudge_y = 0.15)+
      
      geom_text(aes(x = unique(tempx), y = as.factor(unique(c(temp$x1, temp$x2))), label = '|'), nudge_y = 0.09)+
      geom_point(aes(x = unique(tempx), y = unique(tempy)), size = 5, shape = 21, colour = 'black', fill = 'white')+
      scale_colour_gradientn("% of Shared Clonotypes", colours=c("red","violet","blue"))+
      scale_alpha("% of Shared S1 Clonotypes", range = c(0.1, 1))+
      
      scale_x_continuous(expand = c(0.1,0.1))+
      scale_y_discrete(expand = c(0.05,0.05))+
      
      labs(x = 'Time (months)', y = 'Sample of interest (months)', nudge_y = 0.09,title = paste('PATIENT', 'P', 'alpha chain'))+
      theme_classic()+
      theme(legend.position = 'bottom', axis.title=element_text(size=15), axis.text=element_text(size=12))
    
    ggsave(paste(paste(i, j, sep = '_'), '.pdf', sep = ''), 
           ggp, width = 20, height = 25, units = c('cm'))
    
  }
}
