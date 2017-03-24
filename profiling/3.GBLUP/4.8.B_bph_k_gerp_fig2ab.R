### Jinliang Yang
### June 14th, 2016
### updated BPH


############
getvar <- function(res){
  traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
  out <- data.frame()
  for(i in 1:7){
    sub <- subset(res, trait==traits[i])
    tem <- data.frame(trait=traits[i], A=sum(sub$h2_mrk_A), D=sum(sub$h2_mrk_D))
    out <- rbind(out, tem)
  }
  out$trait <- toupper(out$trait)
  return(out)
}

nx_flt <- function(res=dat, x=5){
  traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
  out <- data.frame()
  for(i in 1:7){
    sub <- subset(res, trait==traits[i])
    tot <- sum(sub$H2_mrk)
    tem <- subset(sub, H2_mrk > x*tot/nrow(sub))
    out <- rbind(out, tem)
  }
  out$trait <- toupper(out$trait)
  return(out)
}
##########################################

dat <- read.csv("largedata/lcache/kval_bph_0x.csv")
res1 <- getvar(res=dat)
res2 <- nx_flt(res=dat, x=1)
res2 <- subset(res2, abs(k) < 2)
dim(res2) #107346      8
#############################################
library(ggplot2)
library(reshape2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

med <- read.csv("cache/loh_pBPHmax_median.csv")
#bymed2 <- with(trait, reorder(trait, pBPHmax, median))
bymed <- med[order(med$h),]

#######
theme_set(theme_grey(base_size = 18)) 

out1 <- melt(res1, id.var="trait")
p1 <- ggplot(out1, aes(x=factor(trait, levels=bymed$trait), y=value, 
                       fill=factor(variable, levels =c("A", "D"), labels=c("A", "D")))) + 
  geom_bar(position=position_dodge(), stat="identity") +
  xlab("") +
  ylab("Cumulative Variance") +
  ggtitle("") + theme_bw() +
  labs(fill="Effect") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

res2$trait <- toupper(res2$trait)
res2$trait <- factor(res2$trait, levels=bymed$trait)
p2 <- ggplot(data=res2) +
  geom_density(aes(x= k, y=-..scaled.., fill= as.factor(trait)) ) +
  #guides(fill=FALSE) +
  labs(y=NULL, fill="Trait") + theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  coord_flip() + xlab("Degree of Dominance (k)") + ylab("") + facet_grid(~ trait) 

########
pdf("graphs/SFig2_eff_var.pdf", width=13, height=5)
multiplot(p1, p2, cols=2)
dev.off()










########################

for(i in 1:10){
  res2 <- getk(filepath="largedata/snpeff/perse/", H2_cutoff=i )
  res2$trait <- toupper(res2$trait)
  res2$trait <- factor(res2$trait, levels=bymed$trait)
  p2 <- ggplot(data=res2) +
    geom_density(aes(x= k, y=-..scaled.., fill= as.factor(trait)) ) +
    #guides(fill=FALSE) +
    labs(y=NULL, fill="Trait") + theme_bw() +
    theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
    coord_flip() + xlab("Degree of Dominance (k)") + ylab("") + facet_grid(~ trait) 
  
  
  pdf(paste0("largedata/Test_eff_var_", i, ".pdf"), width=13, height=5)
  multiplot(p2, p2, cols=2)
  dev.off()
  
}


###### calculate the porportion of positive
dat <- read.csv("largedata/lcache/kval_perse_0x.csv")
res1 <- getvar(res=dat)
res2 <- nx_flt(res=dat, x=5)
res2 <- subset(res2, abs(k) < 2)

excess_pos <- function(dat, myt="TW", cutoff=0){
  out2 <- dat
  out2tw <- subset(out2, trait == myt)
  a <- nrow(subset(out2tw, k >= cutoff))/nrow(out2tw) - nrow(subset(out2tw, k < -cutoff))/nrow(out2tw)
  print(a)
}

excess_pos(dat=res2, myt="GY", cutoff=0) #0.4598517
excess_pos(dat=res2, myt="GY", cutoff=1) #0.1057723




