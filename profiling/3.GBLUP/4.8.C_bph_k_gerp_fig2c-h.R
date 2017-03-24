### Jinliang Yang
### June 14th, 2016

##############################################
library(wesanderson)
library(ggplot2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

plot_k_gerp <- function(dat,med2, out, outfile="largedata/lgraphs/gerp_k7x_others.pdf"){
  
  #cols <- wes_palette(7, name = "Zissou", type = "continuous")
  cols <- c("#f6546a", "#daa520", "#00ff00", "#66cdaa", "#3b5998", "#8a2be2", "#ff00ff")
  theme_set(theme_grey(base_size = 18)) 
  
  lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
  p1 <- ggplot(dat, aes(x=RS, y=Effect_A, colour=factor(trait, levels=med2$trait),
                        linetype=factor(trait, levels=med2$trait))) +
    labs(colour="Traits") +
    theme_bw() +
    xlab("GERP Score") +
    ylab("Additive Effect") +
    #  (by default includes 95% confidence region)
    #scale_fill_manual(values=c("#008080", "#003366", "#40e0d0", "#ffa500", "#f6546a", "#ff00ff", "#800000")) +
    #http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=lty1) +
    guides(colour=FALSE, linetype=FALSE) +
    
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
  
  
  lty2 <- getlty(df=out, eff="effd", cutoff=0.05)$l
  p2 <- ggplot(dat, aes(x=RS, y=Effect_D, colour=factor(trait, levels=med2$trait),
                        linetype=factor(trait, levels=med2$trait))) +
    labs(colour="Traits") +
    theme_bw() +
    xlab("GERP Score") +
    ylab("Dominant Effect") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=lty2) +
    guides(colour=FALSE, linetype=FALSE) +
    
    theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
    geom_smooth(method="gam", size=1.3)   # Add linear regression line 
  
  lty3 <- getlty(df=out, eff="effk", cutoff=0.05)$l
  p3 <- ggplot(dat, aes(x=RS, y=k, colour=factor(trait, levels=med2$trait),
                        linetype=factor(trait, levels=med2$trait))) +
    #geom_point(shape=1) +    # Use hollow circles
    labs(colour="Traits") +
    theme_bw() +
    xlab("GERP Score") +
    ylab("Degree of Domiance (k)") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=lty3) +
    guides(colour=FALSE, linetype=FALSE) +
    
    theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
    geom_smooth(method="gam", size=1.3)   # Add linear regression line 
  
  lty4 <- getlty(df=out, eff="h2a", cutoff=0.05)$l
  p4 <- ggplot(dat, aes(x=RS, y=h2_mrk_A, colour=factor(trait, levels=med2$trait),
                        linetype=factor(trait, levels=med2$trait))) +
    #geom_point(shape=1) +    # Use hollow circles
    labs(colour="Traits") +
    theme_bw() +
    xlab("GERP Score") +
    ylab("Additive Variance") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=lty4) +
    guides(colour=FALSE, linetype=FALSE) +
    
    theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
    geom_smooth(method="gam", size=1.3) 
  
  lty5 <- getlty(df=out, eff="h2d", cutoff=0.05)$l
  p5 <- ggplot(dat, aes(x=RS, y=h2_mrk_D, colour=factor(trait, levels=med2$trait),
                        linetype=factor(trait, levels=med2$trait))) +
    #geom_point(shape=1) +    # Use hollow circles
    labs(colour="Traits") +
    theme_bw() +
    xlab("GERP Score") +
    ylab("Dominant Variance") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=lty5) +
    guides(colour=FALSE, linetype=FALSE) +
    
    theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
    geom_smooth(method="gam", size=1.3) 
  
  lty6 <- getlty(df=out, eff="h2k", cutoff=0.05)$l
  p6 <- ggplot(dat, aes(x=RS, y=H2_mrk, colour=factor(trait, levels=med2$trait),
                        linetype=factor(trait, levels=med2$trait))) +
    #geom_point(shape=1) +    # Use hollow circles
    labs(colour="Traits") +
    theme_bw() +
    xlab("GERP Score") +
    ylab("Total Variance") +
    scale_color_manual(values=cols) +
    scale_linetype_manual(values=lty6) +
    guides(colour=FALSE, linetype=FALSE) +
    
    theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
    geom_smooth(method="gam", size=1.3) 
  
  #multiplot(p1, p4, p2, p5, p3, p6, cols=3)
  pdf(outfile, width=13, height=8)
  multiplot(p1, p4, p2, p5, p3, p6, cols=3)
  dev.off()
}

getlty <- function(df, eff, cutoff=0.05){
  df$l <- 2
  if(nrow(df[df[, eff] < cutoff, ]) >0) df[df[, eff] < cutoff, ]$l <- 1
  return(df)
}



##########################################################
#####
geno <- read.csv("largedata/GERPv2/gerpsnp_506898.csv")
geno <- geno[, 1:5]

for(i in 0:0){
  kval <- read.csv(paste0("largedata/lcache/kval_bph_", i, "x.csv"))
  dat <- merge(kval, geno, by.x="snpid", by.y="marker")
  #dat$Effect_A <- dat$Effect_A + mean(dat$Effect_A)
  #dat$Effect_D <- -dat$Effect_D
  #dat$k <- -dat$Effect_D/abs(dat$Effect_A)
  #dat$Effect_D <- -dat$Effect_D
  dat$Effect_A <- -dat$Effect_A
  dat$k <- dat$Effect_D/abs(dat$Effect_A)
  
  if(sum(dat$k > 1) > 0){
    if(sum(dat$k > 2) > 0){
      dat[dat$k > 2, ]$k <- 2
    }
    #out[out$k > 1, ]$k <- rescale(out[out$k > 1, ]$k, c(1, 2))
  }
  if(sum(dat$k < -1) > 0){
    if(sum(dat$k < -2) > 0){
      dat[dat$k < -2, ]$k <- -2
    }
    #out[out$k < -1, ]$k <- rescale(out[out$k < -1, ]$k, c(-2, -1))
  }
  med2 <- data.frame(trait=tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW")), 
                     phph=c(-0.24725, -0.08345, -0.10605,  0.29005,  1.24175,  0.25460, -0.00955 ))
  med2 <- med2[order(med2$phph),]
  
  out <- data.frame()
  #med2$trait <- as.character(med2$trait)
  for(j in 1:7){
    sub <- subset(dat, trait == med2$trait[j])
    t1 <- cor.test(sub$Effect_A, sub$RS)
    t2 <- cor.test(sub$Effect_D, sub$RS)
    t3 <- cor.test(sub$k, sub$RS)
    t4 <- cor.test(sub$h2_mrk_A, sub$RS)
    t5 <- cor.test(sub$h2_mrk_D, sub$RS)
    t6 <- cor.test(sub$H2_mrk, sub$RS)
    
    tem <- data.frame(trait=med2$trait[j], effa=t1$p.value, effar=t1$estimate,
                      effd=t2$p.value, effdr=t2$estimate,
                      effk=t3$p.value, effkr=t3$estimate,
                      h2a=t4$p.value, h2ar=t4$estimate,
                      h2d=t5$p.value, h2dr=t5$estimate,
                      h2k=t6$p.value, h2kr=t6$estimate)
    out <- rbind(out, tem)
  }
  print(i)
  print(out)
  #### start to plot:
  write.table(dat, "largedata/gerp_dat.csv", sep=",", row.names=FALSE, quote=FALSE)
  #plot_k_gerp(dat, med2, out, outfile=paste0("graphs/SFig_gerp_k", i, "x_bph.pdf"))
  
}















###########################
med2$l <- getlty(df=out, eff="effa", cutoff=0.05)$l

cols <- c("#f6546a", "#daa520", "#00ff00", "#66cdaa", "#3b5998", "#8a2be2", "#ff00ff")
theme_set(theme_grey(base_size = 18)) 
le <- ggplot(dat, aes(x=RS, y=Effect_A, colour=factor(trait, levels=med2$trait),
                      linetype=factor(trait, levels=med2$trait))) +
  #geom_point(shape=1) +    # Use hollow circles
  labs(colour="Traits", linetype="type") +
  theme_bw() +
  xlab("GERP Score") +
  ylab("Total Variance") +
  #guides(colour=FALSE) +
  scale_colour_manual(values=cols) +
  scale_linetype_manual(values=med2$l) +
  theme(axis.text.y = element_text(angle = 90, hjust = 1)) +
  geom_smooth(method="gam", size=1.3) 

pdf("graphs/Fig4_k_others_legend.pdf", width=4, height=4)
le
dev.off()


