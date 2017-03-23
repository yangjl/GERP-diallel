### Jinliang Yang
### 9/3/2015
### purpose: find the deleterious per line

gerpsnp <- read.csv("largedata/GERPv2/gerpsnp_506898.csv")
gerpsnp <- subset(gerpsnp, B73 != "N")

#### compute deleterious snps carried per line
del_perline <- function(gerpsnp, cutoff=0){
  
  gerpsnp <- subset(gerpsnp, RS > cutoff)
  message(sprintf("###>>> gerp > %s, total [ %s ]", cutoff, nrow(gerpsnp)))
  out <- data.frame()
  for(i in 11:ncol(gerpsnp)){
    tem <- subset(gerpsnp, B73 != gerpsnp[, i] & gerpsnp[, i]!= "N")
    temout <- data.frame(line=names(gerpsnp)[i], no=nrow(tem), score=mean(tem$RS))
    out <- rbind(out, temout)
  }
  return(out[order(out$no, decreasing=TRUE),])
}

out <- del_perline(gerpsnp, cutoff=0)
write.table(out, "data/deleterious_perline.csv", sep=",", row.names=FALSE, quote=FALSE)

out1 <- del_perline(gerpsnp, cutoff=1)
#write.table(out, "data/deleterious_perline.csv", sep=",", row.names=FALSE, quote=FALSE)

out2 <- del_perline(gerpsnp, cutoff=2)
#write.table(out, "data/deleterious_perline.csv", sep=",", row.names=FALSE, quote=FALSE)


line <- read.csv("data/deleterious_perline.csv")
###############################################################

cal_complementation <- function(gerpsnp, cutoff=0){
  #### compute complementary snps of two lines
  message(sprintf("###>>> converting trait mxtrix"))
  SCA <- read.csv("data/SCA_all_traits.csv")
  SCA$trait <- as.character(SCA$trait)
  traits <- unique(SCA$trait)
  tmx <- subset(SCA, trait== traits[1] & P1 != "B73")
  names(tmx)[c(1,4)] <- c("F1", traits[1])
  tmx$F1 <- paste(tmx$P1, tmx$P2, sep="x")
  for(i in 2:length(traits)){
    tem <- subset(SCA, trait== traits[i] & P1 != "B73")
    names(tem)[c(1,4)] <- c("F1", traits[i])
    tem$F1 <- paste(tem$P1, tem$P2, sep="x")
    
    tmx <- merge(tmx, tem[, c(1,4)], by="F1")
  }
  
  tmx$P1 <- as.character(tmx$P1)
  tmx$P2 <- as.character(tmx$P2)
  
  #### compute complementary snps of two lines
  message(sprintf("###>>> calculating # of complementation!"))
  gerpsnp <- subset(gerpsnp, B73 != "N")
  gerpsnp <- subset(gerpsnp, RS > cutoff)
  tmx$compno <- 0
  tmx$delscore <- 0
  gerpsnp$marker <- as.character(gerpsnp$marker)
  for(i in 1:nrow(tmx)){
    #### non-deleterious of the two parents
    tem1 <- subset(gerpsnp, B73 == gerpsnp[, tmx$P1[i]] )
    tem2 <- subset(gerpsnp, B73 == gerpsnp[, tmx$P2[i]] )
    unq <- unique(c(tem1$marker, tem2$marker))
    sub <- gerpsnp[!(gerpsnp$marker %in% unq), ]
    tmx$compno[i] <- length(unq)
    tmx$delscore[i] <- sum(sub$RS)
  }
  return(tmx)
}


tmx <- cal_complementation(gerpsnp, cutoff=0)
cor.test(tmx$GY, tmx$compno)
cor.test(tmx$GY, tmx$delscore)
write.table(tmx, "data/tmx_complemenation.csv", sep=",", row.names=FALSE, quote=FALSE)
idx <- which.max(tmx$compno)

tmx1 <- cal_complementation(gerpsnp, cutoff=1)
cor.test(tmx$GY, tmx1$compno)
cor.test(tmx$GY, tmx1$delscore)

tmx2 <- cal_complementation(gerpsnp, cutoff=2)
for(i in 4:7){
  print(cor.test(tmx2[, i], tmx$compno))
  cor.test(tmx2[, i], tmx$delscore)  
}

###############################################################

tmx <- read.csv("data/tmx_complemenation.csv")
pdf("graphs/Fig2_complementation.pdf", width=6, height=6)
## add extra space to right margin of plot within frame
par(mar=c(5, 4, 4, 5) + 0.1)

## Plot first set of data and draw its axis
plot(tmx$GY, tmx$compno, pch=16, col="red", yaxt="n", xlab="SCA of GY", ylab="", 
     type="p", main="")
abline(lm(compno ~ GY, data=tmx), col="red", lwd=3)
axis(2, ylim=range(tmx$compno), col.axis="red")  ## las=1 makes horizontal labels
mtext("No. of Complementations",side=2,line=2.5, col="red")

## Allow a second plot on the same graph
par(new=TRUE)

## Plot the second plot and put axis scale on right
plot(tmx$GY, tmx$delscore, pch=16,  xlab="", ylab="",
     axes=FALSE, type="p", col="blue")
abline(lm(delscore ~ GY, data=tmx), col="blue", lwd=3)
## a little farther out (line=4) to make room for labels
mtext("Accumulative GERP Scores",col="blue", side=4,line=2.5) 
axis(4, ylim= range(tmx$delscore), col.axis="blue",las=0)
dev.off()

####### confident interval
n<-50
x<-sample(40:70,n,rep=T)
y<-.7*x+rnorm(n,sd=5)
plot(x,y,xlim=c(20,90),ylim=c(0,80))
mylm<-lm(y~x)
abline(mylm,col="red")
newx<-seq(20,90)
prd<-predict(mylm,newdata=data.frame(x=newx),interval = c("confidence"), 
             level = 0.90,type="response")
lines(newx,prd[,2],col="red",lty=2)
lines(newx,prd[,3],col="red",lty=2)








