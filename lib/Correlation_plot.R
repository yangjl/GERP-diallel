# Jinliang Yang
# 12.21.2010
# Correlation plot

# read in the data
#pheno <- read.csv("/Users/yangjl/Documents/workingSpace/SAM/pheno.ibm.ril_25.csv")

## put (absolute) correlations on the upper panels,
## with size proportional to the correlations.
panel.cor <- function(x, y, digits=2, prefix="*", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- cor(x, y, method="pearson", use="complete.obs")
  
  
  txt <- format(c(r, 0.123456789), digits=digits)[1]
  pv <- cor.test(x,y, method="pearson", alternative="two.sided")$p.value
  if (pv <=0.05) txt <- paste(txt, prefix, sep="")
  if(missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  m = max(r, 0.4)
  text(0.5, 0.5, txt, cex = cex.cor * m)
}

diag <- function(x, y, labels, cex=2, font=2){
  a <- unlist(strsplit(labels, split="_"))
  x=0.5
  for (i in 1:length(a)){
    text(x, 0.5, a[i], cex=2, font=font)}
}

#pairs(krntem[, 2:3], text.panel = diag, upper.panel=panel.smooth, 
#      lower.panel=panel.cor, gap=0, main="KRN Correlation Plots of three obs")
