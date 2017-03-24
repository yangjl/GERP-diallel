### Jinliang Yang
### Sept 19th, 2015
### update: 03-23-2017

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

nx_flt <- function(res=dat, X=5){
  traits <- tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"))
  out <- data.frame()
  for(i in 1:7){
    sub <- subset(res, trait==traits[i])
    tot <- sum(sub$H2_mrk)
    tem <- subset(sub, H2_mrk > X*tot/nrow(sub))
    out <- rbind(out, tem)
  }
  out$trait <- toupper(out$trait)
  return(out)
}
##########################################

files <- list.files(path="largedata/lcache", pattern="perse", full.names = TRUE)

dat <- data.frame()
for(i in 1:length(files)){
    t <- fread(files[i], data.table=FALSE)
    t$trait <- gsub(".*\\/|_perse.*", "", files[i])
    dat <- rbind(dat, t)
}

### for Figure S3
res1 <- getvar(res=dat)
write.table(res1, "cache/gblup_var_updated.csv", sep=",", row.names=FALSE, quote=FALSE)

res2 <- nx_flt(res=dat, X=5)
res2 <- subset(res2, abs(k) < 2)
dim(res2) #107346      8
write.table(res2, "cache/kval_perse_5x.csv", sep=",", row.names=FALSE, quote=FALSE)

#############################################
library(ggplot2)
library(reshape2)
library("data.table")
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

med <- read.csv("cache/loh_pMPH_median.csv")
#bymed2 <- with(trait, reorder(trait, pBPHmax, median))
bymed <- med[order(med$h),]

#######
res1 <- read.csv("cache/gblup_var_updated.csv")
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

res2 <- fread("cache/kval_perse_5x.csv", data.table=FALSE)
res2$trait <- toupper(res2$trait)
res2$trait <- factor(res2$trait, levels=bymed$trait)
p2 <- ggplot(data=res2) +
  geom_density(aes(x= k, y=-..scaled.., fill= as.factor(trait)) ) +
  #guides(fill=FALSE) +
  labs(y=NULL, fill="Trait") + theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  coord_flip() + xlab("Degree of Dominance (k)") + ylab("") + facet_grid(~ trait) 
p2


########
pdf("graphs/Fig_S3_var.pdf", width=6, height=5)
p1
dev.off()




