### Jinliang Yang
### Sept 19th, 2015
### update: 03-23-2017

############


library(ggplot2)
library(reshape2)
library("data.table")

### for Figure S3
t <- fread("largedata/lcache/noB73_kval_perse_5x.csv", data.table=FALSE)
t <- subset(t, !(trait %in% "Sample"))
#res2 <- nx_flt(res=dat, X=5)
res2 <- subset(t, abs(k) < 2)

med <- read.csv("cache/loh_pMPH_median.csv")
#bymed2 <- with(trait, reorder(trait, pBPHmax, median))
bymed <- med[order(med$h),]

#######
#res2 <- fread("cache/noB73_kval_perse_5x.csv", data.table=FALSE)
res2$trait <- toupper(res2$trait)
res2$trait <- factor(res2$trait, levels=bymed$trait)
p2 <- ggplot(data=res2) +
    geom_density(aes(x= k, y=-..scaled.., fill= as.factor(trait)) ) +
    #guides(fill=FALSE) +
    labs(y=NULL, fill="Trait") + theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    coord_flip() + xlab("Degree of Dominance (k)") + ylab("") + facet_grid(~ trait) 

pdf("graphs/SFig_k_dis_excludeB73.pdf", height=5, width=8)
p2
dev.off()








