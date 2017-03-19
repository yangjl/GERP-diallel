## Jinliang Yang
## Oct. 13th, 2014
## phenotypic data of offpvp diallel

#setwd("~/Documents/Github/pvpDiallel/")
ca <- read.csv("data/trait_CA.csv")

#######
get_GCA_matrix <- function(trait=ca){
  ti <- c("ASI", "DTP", "DTS", "EHT",  "GY", "PHT",  "TW")
  out <- trait[, c("P1", "GCA1.all")]
  out <- out[!duplicated(out$P1), ]
  for(i in 1:length(ti)){
    myp <- subset(trait, trait==ti[i])
    myp <- myp[, c("P1", "GCA1.all")]
    myp <- myp[!duplicated(myp$P1), ]
    names(myp)[2] <- ti[i]
    
    out <- merge(out, myp)
  }
  row.names(out) <- out$P1
  return(out[, -1:-2])
}


get_SCA_matrix <- function(trait=ca, pheno="GY"){
  
  sub <- subset(trait, trait == pheno)
  sub$P1 <- as.character(sub$P1)
  sub$P2 <- as.character(sub$P2)
  line <- sort(unique(c(sub$P1, sub$P2) ))
  
  mx <- as.data.frame(matrix(data=NA, nrow=12, ncol=12))
  names(mx) <- line
  row.names(mx) <- line
  for(i in 1:nrow(sub)){
    mx[sub$P1[i], sub$P2[i]] <- sub$SCA.all[i]
  }

  return(mx[-12, -1])
}

#######
GCA <- get_GCA_matrix(trait=ca)
SCA <- ca[, c("trait", "P1", "P2", "SCA.all")]
write.table(SCA, "data/SCA_all_traits.csv", sep=",", row.names=FALSE, quote=FALSE)

SCA_mx <- get_SCA_matrix(trait=ca, pheno="GY")

pdf("graphs/Fig1d.pdf", width=5, height=5)
heatmap(as.matrix(SCA_mx), Rowv = NA, Colv=NA, revC=TRUE)
dev.off()

pdf("graphs/Fig1d.pdf", width=5, height=5)
ggplot(data = subset(SCA, trait=="GY"), aes(P2, P1, fill = SCA.all))+
  geom_tile(color = "white") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "yellow", name="") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 12, hjust = 1),
        axis.text.y = element_text(size = 12, hjust = 1)) +
  xlab("") +
  ylab("") +
  coord_fixed()
dev.off()
