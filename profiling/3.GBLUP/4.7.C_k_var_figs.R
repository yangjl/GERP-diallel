### Jinliang Yang
### 8/24/2016
### answer howie's question 

dat <- read.csv("largedata/lcache/kval_perse_0x.csv")
res1 <- getvar(res=dat)
res2 <- nx_flt(res=dat, x=5)
res2 <- subset(res2, abs(k) < 2)
dim(res2) #107346      8

pdf("graphs/explore1.pdf")
for(mytrait in c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW")){
  c1 <- cor(abs(subset(res2, trait==mytrait)$Effect_D), abs(subset(res2, trait==mytrait)$Effect_A))
  c2 <- cor(abs(subset(res2, trait==mytrait)$Effect_D)/abs(subset(res2, trait==mytrait)$Effect_A), 
            abs(subset(res2, trait==mytrait)$Effect_A))
  par(mfrow=c(2,2))
  hist(abs(subset(res2, trait==mytrait)$Effect_A), xlab="Additive Effect", main=mytrait)
  hist(abs(subset(res2, trait==mytrait)$Effect_D), xlab="Dominant Effect", main=mytrait)
  plot(abs(subset(res2, trait==mytrait)$Effect_D), abs(subset(res2, trait==mytrait)$Effect_A),
       xlab="Additive Effect", ylab="Dominant Effect", main=paste("correlation =", round(c1,3)))
  plot(abs(subset(res2, trait==mytrait)$Effect_D)/abs(subset(res2, trait==mytrait)$Effect_A), 
       abs(subset(res2, trait==mytrait)$Effect_A),
       xlab="Additive Effect", ylab="Dominant/Additive", main=paste("correlation =", round(c2,3)))
  
}
dev.off()



i <- seq(from=0.6, to=0.75, by=0.05)

for(myi in i){
  mytrait <- "PHT"
  sub <- subset(res2, trait==mytrait)
  sub <- sort(abs(sub$Effect_A), decreasing=T)
  l <- length(sub)
  sub2 <- sub[round(l*myi):l]
  hist(sub2, xlab="Additive Effect", main= paste("removing", myi*100, "% large effect loci") )
  
}



c1 <- cor(abs(), abs(subset(res2, trait==mytrait)$Effect_A))
c2 <- cor(abs(subset(res2, trait==mytrait)$Effect_D)/abs(subset(res2, trait==mytrait)$Effect_A), 
          abs(subset(res2, trait==mytrait)$Effect_A))
par(mfrow=c(2,2))

hist(abs(subset(res2, trait==mytrait)$Effect_D), xlab="Dominant Effect", main=mytrait)
plot(abs(subset(res2, trait==mytrait)$Effect_D), abs(subset(res2, trait==mytrait)$Effect_A),
     xlab="Additive Effect", ylab="Dominant Effect", main=paste("correlation =", round(c1,3)))
plot(abs(subset(res2, trait==mytrait)$Effect_D)/abs(subset(res2, trait==mytrait)$Effect_A), 
     abs(subset(res2, trait==mytrait)$Effect_A),
     xlab="Additive Effect", ylab="Dominant/Additive", main=paste("correlation =", round(c2,3)))












##############################################
library(wesanderson)
library(ggplot2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")

plot_k_gerp <- function(dat,med2, out, outfile="largedata/lgraphs/gerp_k7x_others.pdf"){
  
  #cols <- wes_palette(7, name = "Zissou", type = "continuous")
  cols <- c("#f6546a", "#daa520", "#00ff00", "#66cdaa", "#3b5998", "#8a2be2", "#ff00ff")
  theme_set(theme_grey(base_size = 18)) 
  
  #lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
  p1 <- ggplot(dat$data, aes(x=abs(Effect_A), y=k, colour=factor(trait, levels=dat$med$trait),
                        linetype=factor(trait, levels=dat$med$trait))) +
    facet_grid(~trait, scales = "free") +
    labs(colour="Traits") +
    theme_bw() +
    xlab("Additive Effect") +
    ylab("Degree of Dominance") +
    #  (by default includes 95% confidence region)
    #scale_fill_manual(values=cols) +
    #http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
    scale_color_manual(values=cols) +
    #scale_linetype_manual(values=lty1) +
    guides(colour=FALSE, linetype=FALSE) +
    
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
  
  p2 <- ggplot(dat$data, aes(x=abs(Effect_D), y=k, colour=factor(trait, levels=dat$med$trait),
                             linetype=factor(trait, levels=dat$med$trait))) +
    labs(colour="Traits") +
    theme_bw() +
    facet_grid(~trait, scales = "free") +
    xlab("Dominace Effect") +
    ylab("Degree of Dominance") +
    #  (by default includes 95% confidence region)
    #scale_fill_manual(values=cols) +
    #http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
    #scale_color_manual(values=cols) +
    scale_colour_manual(values=cols) +
    labs(colour="Traits") +
    guides(linetype=FALSE) +
    geom_smooth(method="gam", size=1.3) +
    theme(axis.text.y = element_text(angle = 90, hjust = 1))
  
  
  
  #multiplot(p1, p4, p2, p5, p3, p6, cols=3)
  pdf(outfile, width=13, height=8)
  multiplot(p1, p2, cols=2)
  dev.off()
}

getlty <- function(df, eff, cutoff=0.05){
  df$l <- 2
  if(nrow(df[df[, eff] < cutoff, ]) >0) df[df[, eff] < cutoff, ]$l <- 1
  return(df)
}



##########################################################
##### prepare data for plotting!!!

getdata <- function(toremove=FALSE){
  geno <- read.csv("largedata/GERPv2/gerpsnp_506898.csv")
  geno <- geno[, 1:5]
  
  kval <- read.csv("largedata/lcache/kval_perse_1x.csv")
  dat <- merge(kval, geno, by.x="snpid", by.y="marker")
  #dat$Effect_A <- dat$Effect_A + mean(dat$Effect_A)
  #dat$Effect_D <- -dat$Effect_D
  #dat$k <- -dat$Effect_D/abs(dat$Effect_A)
  #dat$Effect_D <- -dat$Effect_D
  dat$Effect_A <- -dat$Effect_A
  dat$k <- dat$Effect_D/abs(dat$Effect_A)
  
  dat <- subset(dat, trait %in% tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW")))
  
  ### whether to remove k > or k < -2 sites
  if(toremove){
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
  }
  
  
  med2 <- data.frame(trait=tolower(c("ASI", "DTP", "DTS", "EHT", "GY", "PHT", "TW")), 
                     phph=c(-0.24725, -0.08345, -0.10605,  0.29005,  1.24175,  0.25460, -0.00955) )
  med2 <- med2[order(med2$phph), ]
  
  out <- data.frame()
  med2$trait <- as.character(med2$trait)
  for(j in 1:7){
    sub <- subset(dat, trait == med2$trait[j])
    fit1 <- lm(sub$k ~ sub$Effect_A)
    fit2 <- lm(sub$k ~ sub$Effect_D)
    
    tem1 <- data.frame(trait=med2$trait[j], a=coef(fit1)["(Intercept)"], beta=coef(fit1)["sub$Effect_A"], 
                       pval=summary(fit1)$coefficients[2,4])
    tem1$effect <- "add"
    tem2 <- data.frame(trait=med2$trait[j], a=coef(fit2)["(Intercept)"], beta=coef(fit2)["sub$Effect_D"],
                       pval=summary(fit2)$coefficients[2,4])
    tem2$effect <- "dom"
    out <- rbind(out, tem1, tem2)
  }
  
  return(list(data=dat, test=out, med=med2))
}

########
dat <- getdata(toremove=FALSE)  
save(list="dat", file="largedata/dat_08152016.RData")

ob <- load("largedata/dat_08152016.RData")

#####
dat2 <- getdata(toremove=TRUE)  
save(list="dat2", file="largedata/dat2_08152016.RData")

ob <- load("largedata/dat2_08152016.RData")
dat <- dat2

d <- dat$data
hist(abs(subset(d, trait == "gy")$Effect_A ))

#lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
p1 <- ggplot(d, aes(x=abs(Effect_A), y=k, colour=factor(trait, levels=dat$med$trait),
                           linetype=factor(trait, levels=dat$med$trait))) +
  facet_grid(~trait, scales = "free") +
  labs(colour="Traits") +
  theme_bw() +
  xlab("Additive Effect") +
  ylab("Degree of Dominance") +
  #  (by default includes 95% confidence region)
  #scale_fill_manual(values=cols) +
  #http://www.sthda.com/english/wiki/ggplot2-colors-how-to-change-colors-automatically-and-manually
  scale_color_manual(values=cols) +
  #scale_linetype_manual(values=lty1) +
  guides(colour=FALSE, linetype=FALSE) +
  
  geom_smooth(method="gam", size=1.3) +
  theme(axis.text.y = element_text(angle = 90, hjust = 1))

#### start to plot:
plot_k_gerp(dat, med2, out, outfile=paste0("largedata/lgraphs/gerp_k", i, "x_tem.pdf"))
#plot_k_gerp(dat, med2, out, outfile=paste0("largedata/lgraphs/gerp_k_q", i, ".pdf"))  
  

write.table(out, "cache/s_estimation.csv", sep=",", row.names=FALSE, quote=FALSE)







