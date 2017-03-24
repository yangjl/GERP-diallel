### Jinliang Yang
### April 29th, 2015
### plot the variance explained by genome-wide markers

mc1 <- read.csv("cache/model_comparison_perse.csv")
mc2 <- read.csv("cache/model_comparison_bph.csv")



library(ggplot2)
source("~/Documents/Github/zmSNPtools/Rcodes/multiplot.R")
library(reshape2)
res1 <- melt(mc1, id.vars = c("trait", "mode"))
res1$qval <- p.adjust(res1$value, method="fdr")

res2 <- melt(mc2, id.vars = c("trait", "mode"))
res2$qval <- p.adjust(res2$value, method="fdr")

h <- read.csv("cache/loh_pBPHmax_median.csv")
med2 <- h
med2$traitlw <- tolower(med2$trait)
#bymed2 <- with(trait, reorder(trait, pBPHmax, median))
bymed2 <- med2[order(med2$h),]

res1 <- subset(res1, variable %in% c("p21", "p43"))
res1$variable <- gsub("p21", "m2 vs m1", res1$variable)
res1$variable <- gsub("p43", "m2 vs m3", res1$variable)
names(res1)[3] <- "Comparison"
p1 <- ggplot(subset(res1, mode=="a2"), aes(x=factor(trait, levels=bymed2$trait), y= -log10(qval) )) + 
  geom_point(aes(color = Comparison), size = 5) +
  xlab("") +
  ylab("-log10(Adjusted P-value)") +
  ggtitle("Additive") + theme_bw() +
  ylim(c(0, 30)) +
  #guides(fill=guide_legend(title="Comparison")) + 
  geom_hline(yintercept=-log10(0.01), colour = "red", lty=2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12), legend.position=c(0.8, 0.8))

p2 <- ggplot(subset(res1, mode=="d2"), aes(x=factor(trait, levels=bymed2$trait), y= -log10(qval) )) + 
  geom_point(aes(color = Comparison), size = 5) +
  xlab("") +
  ylab("-log10(Adjusted P-value)") +
  ggtitle("Dominance") + theme_bw() +
  ylim(c(0, 30)) +
  guides(colour=FALSE) + 
  geom_hline(yintercept=-log10(0.01), colour = "red", lty=2) +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

p3 <- ggplot(subset(res1, mode=="h2"), aes(x=factor(trait, levels=bymed2$trait), y= -log10(qval) )) + 
  geom_point(aes(color = Comparison), size = 5) +
  xlab("") +
  ylab("-log10(Adjusted P-value)") +
  ggtitle("Incomplete Dominance") + theme_bw() +
  geom_hline(yintercept=-log10(0.01), colour = "red", lty=2) +
  ylim(c(0, 30)) +
  guides(colour=FALSE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

########
pdf("graphs/Fig_model_comp.pdf", width=12, height=4)
multiplot(p1, p2, p3, cols=3)
dev.off()





p4 <- ggplot(subset(res2, mode=="a2"), aes(x=factor(trait, levels=bymed2$trait), y= -log10(value) )) + 
  geom_point(aes(color = variable), size = 5) +
  xlab("") +
  ylab("Posterior Variance Explained") +
  ggtitle("Additive") + theme_bw() +
  labs(fill="Models") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

p5 <- ggplot(subset(res2, mode=="d2"), aes(x=factor(trait, levels=bymed2$trait), y= -log10(value) )) + 
  geom_point(aes(color = variable), size = 5) +
  xlab("") +
  ylab("Posterior Variance Explained") +
  ggtitle("Dominance") + theme_bw() +
  labs(fill="Models") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))

p6 <- ggplot(subset(res2, mode=="h2"), aes(x=factor(trait, levels=bymed2$trait), y= -log10(value) )) + 
  geom_point(aes(color = variable), size = 5) +
  xlab("") +
  ylab("Posterior Variance Explained") +
  ggtitle("Incomplete Dominance") + theme_bw() +
  labs(fill="Models") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, size=12))







