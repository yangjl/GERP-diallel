### Jinliang Yang
### jan 29th, 2016

library(ggplot2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")



res1 <- read.csv("largedata/var_explained_shuffling.csv")

r0 <- subset(res1, stime == 0)
r1 <- subset(res1, stime != 0)
r0 <- r0[, c("trait", "h2")]
names(r0)[2] <- "h20"
r1 <- merge(r1, r0, by="trait")
r1$trait <- toupper(r1$trait)

p1 <- ggplot(r1, aes(x=h2, y=..density..)) + 
    geom_histogram(position="identity") +
    #geom_density(aes(x=h2, y=..density..)) +
    geom_vline(aes(xintercept=h20), color="red") + 
    facet_wrap( ~ trait, scales="free") +
    ylab("Density") +
    xlab("Proportion of variance explained") +
    ggtitle("") + theme_bw() 
    #theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

pdf("graphs/SFig_var_explianed.pdf", width=8, height=8)
p1
dev.off()




