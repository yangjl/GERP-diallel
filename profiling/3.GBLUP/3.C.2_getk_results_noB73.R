### Jinliang Yang
### 09/08/2015
### updated: 03-23-2017


##########################################
getk <- function(filepath="largedata/snpeff/pBPH/", deff, q=0.9, method="q"){
    files <- list.files(path=filepath, pattern="perse", full.names=TRUE)
    
    output <- data.frame()
    for(i in 1:length(files)){
        h1 <- fread(files[i], header=TRUE, data.table=FALSE)
        #h1 <- as.data.frame(h1)
        names(h1) <- c("snpid","chr","pos","Effect_A","Effect_D","Effect_A2","Effect_D2","h2_mrk_A", 
                       "h2_mrk_D","H2_mrk","h2_mrk_A_p","h2_mrk_D_p","H2_mrk_p","log10_h2_mrk_A","log10_h2_mrk_D","log10_H2_mrk")
        #h1 <- subset(h1, Effect_A !=0)
        #tot <- sum(h1$H2_mrk)
        
        # del == beneficial alleles!!!
        ba <- subset(deff, del == B73)$marker
        out1 <- h1[h1$snpid %in% ba, ]
        
        
        da <- subset(deff, del != B73)$marker
        out2 <- h1[h1$snpid %in% da, ]
        out2$Effect_A <- -out2$Effect_A
        out2$Effect_D <- -out2$Effect_D
        
        h1 <- rbind(out1, out2)
        
        if(method == "q"){
            qv <- quantile(h1$H2_mrk, probs=q)
            message(sprintf("###>>> quantile: [ %s ]", names(qv)))
            h2 <- subset(h1, H2_mrk > qv & Effect_A !=0)
        }else if(method == "var"){
            h1 <- subset(h1, Effect_A !=0 )
            qv <- q*sum(h1$H2_mrk)/nrow(h1)
            message(sprintf("###>>> var/#snp >: [ %s ], time [ %s ]", qv, q))
            h2 <- subset(h1, H2_mrk > qv)
        }else{
            message(sprintf("###>>> No filteration (except Effect_A !=0 )"))
            h2 <- subset(h1, Effect_A !=0 )
        }
        
        h2$k <- h2$Effect_D/h2$Effect_A
        out <- h2[, c("snpid", "k", "Effect_A", "Effect_D", "h2_mrk_A", "h2_mrk_D", "H2_mrk")]
        if(sum(out$k > 1) > 0){
            if(sum(out$k > 2) > 0){
                out[out$k > 2, ]$k <- 2
            }
            #out[out$k > 1, ]$k <- rescale(out[out$k > 1, ]$k, c(1, 2))
        }
        if(sum(out$k < -1) > 0){
            if(sum(out$k < -2) > 0){
                out[out$k < -2, ]$k <- -2
            }
            #out[out$k < -1, ]$k <- rescale(out[out$k < -1, ]$k, c(-2, -1))
        }
        myt <- gsub(".*/|_.*", "", files[i])
        message(sprintf("###>>> trait [ %s ], snp # [ %s ], k ranged [ %s - %s ]", myt, nrow(out), max(out$k), min(out$k)))
        #write.table(out, paste0(outpwd, "/", trait, "_k.txt"), sep="\t", row.names=FALSE, quote=FALSE)
        out$trait <- myt
        output <- rbind(output, out)
    }
    return(output)
}


###############################
### determine direction of effects
library(data.table)
deff <- fread("largedata/gerpsnp_v3_345176_del.csv", data.table=FALSE)

### var/#snp filteration
for(i in c(0,5)){
    res2 <- getk(filepath="largedata/snpeff/rsnp0", deff, q=i, method="var")
    write.table(res2, paste0("largedata/lcache/noB73_kval_perse_", i, "x.csv"), sep=",", row.names=FALSE, quote=FALSE)
}
