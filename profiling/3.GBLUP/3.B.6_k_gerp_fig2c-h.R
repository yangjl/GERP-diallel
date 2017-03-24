### Jinliang Yang
### 9/24/2015

##############################################
library(wesanderson)
library(ggplot2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

plot_k_gerp <- function(dat,med2, out, outfile="largedata/lgraphs/gerp_k7x_others.pdf"){
  
  #cols <- wes_palette(7, name = "Zissou", type = "continuous")
  cols <- c( "#daa520", "#3b5998", "#ff00ff")
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
geno <- fread("largedata/gerpsnp_v3_410183.csv", data.table=FALSE)
geno <- geno[, 1:4]

res2 <- fread("cache/kval_perse_5x.csv", data.table=FALSE)

out <- merge(res2, geno,  by.x="snpid", by.y="marker")
write.table(out, "largedata/lcache/kval_perse_5x_gerpv3.csv", sep=",", row.names=FALSE, quote=FALSE)



med2 <- data.frame(trait=c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW"), 
                   phph=c(-0.24725, -0.08345, -0.10605,  0.29005,  1.24175,  0.25460, -0.00955) )
med2 <- med2[order(med2$phph), ]


kval <- fread("largedata/lcache/kval_perse_5x_gerpv3.csv", data.table=FALSE)
kval$Effect_A <- abs(kval$Effect_A)
kval$Effect_D <- abs(kval$Effect_D)


plot_adk_gerp <- function(outfile, getpdf){
    
    dat <- fread("largedata/lcache/kval_perse_5x_gerpv3.csv", data.table=FALSE)
    dat$Effect_A <- abs(dat$Effect_A)
    dat$Effect_D <- abs(dat$Effect_D)
    
    dat$trait <- toupper(dat$trait)
    med2 <- read.csv("cache/loh_pMPH_median.csv")
    #out <- read.csv("cache/eff_adk_1x.csv")
    #cols <- wes_palette(7, name = "Zissou", type = "continuous")
    cols <- c("#f6546a", "#daa520", "#00ff00", "#66cdaa", "#3b5998", "#8a2be2", "#ff00ff")
    theme_set(theme_grey(base_size = 18)) 
    
    getlty <- function(df, eff, cutoff=0.05){
        df$l <- 2
        if(nrow(df[df[, eff] < cutoff, ]) >0) df[df[, eff] < cutoff, ]$l <- 1
        return(df)
    }
    
    lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
    p1 <- ggplot(dat, aes(x=RS, y=Effect_A, colour=factor(trait, levels=med2$trait),
                          linetype=factor(trait, levels=med2$trait))) +
        labs(colour="Traits") +
        theme_bw() +
        xlab("GERP Score") +
        ylab("Additive Effect") +
        scale_color_manual(values=cols) +
        scale_linetype_manual(values=lty1) +
        guides(colour=FALSE, linetype=FALSE) +
        geom_smooth(method="gam", size=1.3) +
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
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
        geom_smooth(method="gam", size=1.3) +  # Add linear regression line 
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    lty3 <- getlty(df=out, eff="effk", cutoff=0.05)$l
    p3 <- ggplot(dat, aes(x=RS, y=k, colour=factor(trait, levels=med2$trait),
                          linetype=factor(trait, levels=med2$trait))) +
        labs(colour="Traits") +
        theme_bw() +
        xlab("GERP Score") +
        ylab("Degree of Domiance (k)") +
        scale_color_manual(values=cols) +
        scale_linetype_manual(values=lty3) +
        guides(colour=FALSE, linetype=FALSE) +
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize)) +
        geom_smooth(method="gam", size=1.3)   # Add linear regression line 
    
    multiplot(p1, p2, p3, cols=3)
    
    
    if(getpdf == TRUE){
        pdf(outfile, width=13, height=4)
        multiplot(p1, p2, p3, cols=3)
        dev.off()
    }
    
}

####
plot_adk_gerp(outfile="graphs/Fig2c-e_adk.pdf", getpdf)


