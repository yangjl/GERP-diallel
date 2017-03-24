### Jinliang Yang
### 11/02/2016
### move Fig1a to Fig S3

############
getvar <- function(res, x=5){
  traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
  out <- data.frame()
  for(i in 1:7){
    sub <- subset(res, trait==traits[i])
    tot <- sum(sub$H2_mrk)
    sub2 <- subset(sub, H2_mrk > x*tot/nrow(sub))
    tem <- data.frame(trait=traits[i], A=sum(sub2$h2_mrk_A), D=sum(sub2$h2_mrk_D))
    out <- rbind(out, tem)
  }
  out$trait <- toupper(out$trait)
  return(out)
}

##########################################

library(data.table)
dat <- fread("largedata/lcache/kval_perse_0x.csv", data.table=FALSE)
res1 <- getvar(res=dat, x=1)
write.table(res1, "cache/gblup_var_1x.csv", sep=",", row.names=FALSE, quote=FALSE)

res1 <- read.csv("cache/gblup_var_1x.csv")
#############################################
library(ggplot2)
library(reshape2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

med <- read.csv("cache/loh_pBPHmax_median.csv")
#bymed2 <- with(trait, reorder(trait, pBPHmax, median))
bymed <- med[order(med$h),]

temp <- merge(med, res1, by="trait")
cor.test(temp$h, temp$D)

#######
theme_set(theme_grey(base_size = 18)) 

out1 <- melt(res1, id.var="trait")
p1 <- ggplot(out1, aes(x=factor(trait, levels=bymed$trait), y=value, 
                       fill=factor(variable, levels =c("A", "D"), labels=c("A", "D")))) + 
  geom_bar(position=position_dodge(), stat="identity") +
  xlab("") +
  ylab("Variance Explained") +
  ggtitle("") + theme_bw() +
  labs(fill="Effect") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))
p1
