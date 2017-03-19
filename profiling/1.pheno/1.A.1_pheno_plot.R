## Jinliang Yang
## 03-18-2017
## phenotypic data of diallel

#setwd("~/Documents/Github/pvpDiallel/")
trait <- read.csv("data/hyb_heterosis.csv")
trait$trait <- toupper(trait$trait)

out <- data.frame()
for(t in c("GY", "ASI", "DTP", "DTS", "EHT",  "PHT",  "TW")){
  sub <- subset(trait, trait == t)
  tem <- data.frame(trait=t, ybar=mean(sub$valHyb))
  out <- rbind(out, tem)
}
write.table(out, "cache/ybar.csv", sep=",", row.names=FALSE)



#######
plot_trait_per_se <- function(trait=trait, ...){
  par(mar=c(3,3,4,1))
  layout(matrix(c(1,1,2,3,4,1,1,5,6,7), 2, 5, byrow = TRUE))
  ti <- c("GY", "ASI", "DTP", "DTS", "EHT",  "PHT",  "TW")
  for(i in 1:length(ti)){
    myp <- subset(trait, trait==ti[i])
    myp$norv <- scale(myp$valHyb)
    d <- density(myp$norv)
    plot(d, main=ti[i], xlab="", cex.axis=1, bty="n", ...)
    polygon(d, col="antiquewhite3", border="antiquewhite3", lwd=3)
  }
}
#####
pdf("graphs/Fig_S1b.pdf", width=12, height=6)
plot_trait_per_se(trait=trait)
dev.off()


normality_test <- function(trait=trait){
  ti <- c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW")
  
  res <- vector()
  for(i in 1:length(ti)){
    myp <- subset(trait, trait==ti[i])
    out <- shapiro.test(myp$valHyb)
    res <- c(res, out$p.value)
  }
  return(res)
}

#######
res <- normality_test(trait=trait)
#0.05244055 0.18683006 0.19546192 0.86986939 0.59271848 0.82229468 0.48353513


