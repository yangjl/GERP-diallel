# Jinliang Yang
# 03-18-2017
# purpose: levels of heterosis (MPH)


#######
plot_loh <- function(trait,  ...){
  bymed2 <- with(trait, reorder(trait, pMPH, median))
  boxplot(pMPH ~ bymed2, data=trait,
          xlab = "", ylab= "MPH (100%)", col="antiquewhite3", 
          ...)
}
##### note: change to abs value

trait <- read.csv("data/hyb_heterosis.csv")
trait$pMPH <- abs(trait$pMPH*100)

pdf("graphs/Fig1a.pdf", width=5, height=5)
bymed <- plot_loh(trait, main="")
dev.off()
write.table(trait, "table/STable_heterosis.csv", sep=",", row.names=FALSE, quote=FALSE)


##########
library(plyr)
loh <- ddply(trait, .(trait), summarise,
             h = median(pMPH))
loh <- loh[order(loh$h),]
write.table(loh, "cache/loh_pMPH_median.csv", sep=",", row.names=FALSE, quote=FALSE)

