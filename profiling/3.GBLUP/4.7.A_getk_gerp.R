### Jinliang Yang
### Sept 19th, 2015

##########################################
getk <- function(filepath="largedata/snpeff/pBPH/", deff, q=0.9, method="q"){
    
  files <- list.files(path=filepath, pattern="snpe$", full.names=TRUE)
  
  output <- data.frame()
  for(i in 1:length(files)){
    h1 <- fread(files[i], header=TRUE, data.table=FALSE)
    #h1 <- as.data.frame(h1)
    names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                   "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
    #h1 <- subset(h1, Effect_A !=0)
    #tot <- sum(h1$H2_mrk)
    #out1 <- h1[h1$snpid %in% deff$snpid, ]
    #out1 <- h1
    #out1$Effect_A <- -out1$Effect_A
    #out1$Effect_D <- -out1$Effect_D
    
    #out2 <- h1[!(h1$snpid %in% deff$snpid), ]
    #h1 <- rbind(out1, out2)
    
    if(method == "q"){
      qv <- quantile(h1$H2_mrk, probs=q)
      message(sprintf("###>>> quantile: [ %s ]", names(qv)))
      h2 <- subset(h1, H2_mrk > qv & Effect_A !=0)
    }else if(method == "var"){
      h1 <- subset(h1, Effect_A !=0 )
      qv <- q*sum(h1$H2_mrk)/nrow(h1)
      message(sprintf("###>>> var/#snp >: [ %s ], time [ %s ]", qv, q))
      h2 <- subset(h1, H2_mrk > qv)
    }else{
      message(sprintf("###>>> No filteration (except Effect_A !=0 )"))
      h2 <- subset(h1, Effect_A !=0 )
    }

    h2$k <- h2$Effect_D/h2$Effect_A
    out <- h2[, c("snpid", "k", "Effect_A", "Effect_D", "h2_mrk_A", "h2_mrk_D", "H2_mrk")]
    if(sum(out$k > 1) > 0){
      if(sum(out$k > 2) > 0){
        out[out$k > 2, ]$k <- 2
      }
      #out[out$k > 1, ]$k <- rescale(out[out$k > 1, ]$k, c(1, 2))
    }
    if(sum(out$k < -1) > 0){
      if(sum(out$k < -2) > 0){
        out[out$k < -2, ]$k <- -2
      }
      #out[out$k < -1, ]$k <- rescale(out[out$k < -1, ]$k, c(-2, -1))
    }
    myt <- gsub(".*/|_.*", "", files[i])
    message(sprintf("###>>> trait [ %s ], snp # [ %s ], k ranged [ %s - %s ]", myt, nrow(out), max(out$k), min(out$k)))
    #write.table(out, paste0(outpwd, "/", trait, "_k.txt"), sep="\t", row.names=FALSE, quote=FALSE)
    out$trait <- myt
    output <- rbind(output, out)
  }
  return(output)
}

###############################
### determine direction of effects
deff <- read.csv("largedata/Alignment/conserved_alleles_AGPv2.csv")
deff$major <- as.character(deff$major)
deff$Zea <- as.character(deff$Zea)
deff <- subset(deff, major != Zea)



### var/#snp filteration
for(i in 0:10){
  res2 <- getk(filepath="largedata/GBLUP", deff, q=i, method="var")
  write.table(res2, paste0("largedata/lcache/kval_perse_", i, "x.csv"), sep=",", row.names=FALSE, quote=FALSE)
}


### quantile filtering
for(i in c(1:10, 2.5, 7.5)){
  res2 <- getk(filepath="largedata/snpeff/perse/", q=0.1*i, method="q")
  write.table(res2, paste0("largedata/lcache/kval_perse_q", i, ".csv"), sep=",", row.names=FALSE, quote=FALSE)
}

i =7.5
res2 <- getk(filepath="largedata/snpeff/perse/", q=0.1*i)
write.table(res2, paste0("largedata/lcache/kval_perse_q", i, ".csv"), sep=",", row.names=FALSE, quote=FALSE)

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

res1 <- getvar(res=res2)

#############################################
library(ggplot2)
library(reshape2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

med <- read.csv("cache/loh_pBPHmax_median.csv")
#bymed2 <- with(trait, reorder(trait, pBPHmax, median))
bymed <- med[order(med$h),]

###### calculate the porportion of positive
excess_pos <- function(myt="TW"){
  out2tw <- subset(out2, trait == myt)
  a <- nrow(subset(out2tw, h >= 0))/nrow(out2tw) - nrow(subset(out2tw, h < 0))/nrow(out2tw)
  print(a)
}

excess_pos(myt="TW") #0.2505068
excess_pos(myt="PHT") #0.1352641
excess_pos(myt="EHT") #0.17809
excess_pos(myt="GY") #0.1656072

###### calculate the porportion of positive
per_overd <- function(myt="TW"){
  out2tw <- subset(out2, trait == myt)
  a <- nrow(subset(out2tw, h >= 1))/nrow(out2tw)
  print(a)
}

per_overd(myt="TW") #0.1368988
per_overd(myt="PHT") #0.05725433
per_overd(myt="EHT") #0.04802514
per_overd(myt="GY") #0.2914191





#########################################
out1 <- melt(res1, id.var="trait")
theme_set(theme_grey(base_size = 18)) 

p1 <- ggplot(out1, aes(x=factor(trait, levels=bymed2$trait), y=value, 
                       fill=factor(variable, levels =c("A", "D"), labels=c("A", "D")))) + 
  geom_bar(position=position_dodge(), stat="identity") +
  xlab("") +
  ylab("Accumulative Variance") +
  ggtitle("") + theme_bw() +
  labs(fill="Effect") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

res2$trait <- toupper(res2$trait)
res2$trait <- factor(res2$trait, levels=bymed2$trait)
p2 <- ggplot(data=res2) +
  geom_density(aes(x= k, y=-..scaled.., fill= as.factor(trait)) ) +
  #guides(fill=FALSE) +
  labs(y=NULL, fill="Trait") + theme_bw() +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank()) +
  coord_flip() + xlab("Degree of Dominance (k)") + ylab("") + facet_grid(~ trait) 


for(i in 1:10){
  res2 <- getk(filepath="largedata/snpeff/perse/", H2_cutoff=i )
  res2$trait <- toupper(res2$trait)
  res2$trait <- factor(res2$trait, levels=bymed2$trait)
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

########
pdf("graphs/Fig_eff_var.pdf", width=13, height=5)
multiplot(p1, p2, cols=2)
dev.off()



