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



##########
res1 <- read.csv("largedata/var_explained_shuffling.csv")
r0 <- subset(res1, stime == 0)
r100k <- read.csv("largedata/var_explained_100kb.csv")

r2 <- rbind(r0, r100k)
r2$trait <- toupper(r2$trait)
p2 <- ggplot(r2, aes(x=factor(trait), y=h2, 
                       fill=factor(stime, levels=c("ws", "0"), labels=c("whole", "subset")))) + 
    geom_bar(position=position_dodge(), stat="identity") +
    xlab("") +
    ylim(c(0,1)) +
    ylab("Proportion of Variance Explained") +
    ggtitle("") + theme_bw() +
    labs(fill="SNP Set") +
    theme(axis.text.x = element_text(angle = 0, size=12))

#theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

pdf("graphs/SFig_var_sub_100k.pdf", width=5, height=5)
p2
dev.off()


