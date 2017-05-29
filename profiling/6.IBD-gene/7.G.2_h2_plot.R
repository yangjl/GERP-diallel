### Jinliang Yang
### April 29th, 2015
### plot the variance explained by genome-wide markers

format_h2 <- function(){
  gs1 <- read.csv("cache/allsnp_wholeset_h2.csv")
  gs1$type <- "allsnp"
  gs1$transf <- gsub("_..$", "", gs1$mode)
  gs1$transf <- gsub(".*_", "", gs1$transf)
  gs1$mode <- gsub(".*_", "", gs1$mode)
  
  gs2 <- read.csv("cache/gs_wholeset_h2.csv")
  gs2$type <- "genicsnp"
  
  h2 <- rbind(gs1, gs2)
  
  #idx <- grep("MPH", h2$transf)
  #h2 <- h2[-idx, ]
  
  ###>>> ASI part
  asi <- subset(h2, trait == "asi" )
  idx2 <- grep("BPHmax", asi$transf)
  asi <- asi[-idx2, ]
  asi$transf <- gsub("min", "max", asi$transf)
  
  ###>>> non ASI
  nonasi <- subset(h2, trait != "asi")
  nonasi <- subset(nonasi, transf != "BPHmin" & transf != "pBPHmin")
  
  h2 <- rbind(asi, nonasi)
  
  h2mx <- data.frame()
  for(typei in c("allsnp", "genicsnp")){
    mytype <- subset(h2, type == typei)
    for(modei in c("a2", "d2")){
      mymode <- subset(mytype, mode == modei)
      for(tsfi in c("perse", "BPHmax")){
        mytsf <- subset(mymode, transf == tsfi)
        tem <- data.frame(type=typei, mode=modei, tsf=tsfi, 
                          ASI=mytsf[mytsf$trait=="asi", ]$h2,
                          DTP=mytsf[mytsf$trait=="dtp", ]$h2,
                          DTS=mytsf[mytsf$trait=="dts", ]$h2,
                          EHT=mytsf[mytsf$trait=="eht", ]$h2,
                          GY=mytsf[mytsf$trait=="gy", ]$h2,
                          PHT=mytsf[mytsf$trait=="pht", ]$h2,
                          TW=mytsf[mytsf$trait=="tw", ]$h2)
        h2mx <- rbind(h2mx, tem)
      }
    }
  }
  return(h2mx)
}

###>>>>>>
h2 <- format_h2()
counts <- table(mtcars$vs, mtcars$gear)

pdf("manuscript/Figure_Table/Figure_h2.pdf", width=10, height=5)
par(mfrow=c(1,2))
#h2mx_all_a2 <- as.matrix(subset(h2, type == "allsnp" & mode == "a2")[, 4:10])
#barplot(h2mx_all_a2, main="Additive model using SNPs with GERP > 0", ylim=c(0, 1), ylab="Phenotypic variance explained",
#        xlab="Phenotypic traits", beside=TRUE)
#h2mx_all_d2 <- as.matrix(subset(h2, type == "allsnp" & mode == "d2")[, 4:10])
#barplot(h2mx_all_d2, main="Dominant model using SNPs with GERP > 0", ylim=c(0, 1), ylab="Phenotypic variance explained",
#        xlab="Phenotypic traits", beside=TRUE)

h2mx_all_a2 <- as.matrix(subset(h2, type == "genicsnp" & mode == "a2")[, 4:10])
barplot(h2mx_all_a2, main="Additive model", ylim=c(0, 1), ylab="Phenotypic variance explained",
        xlab="Phenotypic traits", beside=TRUE)
h2mx_all_d2 <- as.matrix(subset(h2, type == "genicsnp" & mode == "d2")[, 4:10])
barplot(h2mx_all_d2, main="Dominant model", ylim=c(0, 1), ylab="Phenotypic variance explained",
        xlab="Phenotypic traits", beside=TRUE)
dev.off()

