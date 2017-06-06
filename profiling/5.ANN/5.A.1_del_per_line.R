### Jinliang Yang
### 04/13/2017
### cal deleterious per line and complemantation


library(data.table)
g <- fread("largedata/gerpsnp_v3_345176_del.csv", data.table=FALSE)
## Note col "del" is the major from the alignment, should be the beneficial alleles

get_per_line <- function(g){
    g$del <- as.character(g$del)
    out <- data.frame()
    for(i in 10:ncol(g)){
        g[, i] <- as.character(g[, i])
        sub <- subset(g, !(g[, i] %in% "N"))
        b <- subset(sub, g[, "del"] == g[, i])
        b <- b[!is.na(b$snpid), ]
        
        d <- subset(sub, g[, "del"] != g[, i])
        d <- d[!is.na(d$snpid), ]
        lod <- sum(d$RS)
        tem <- data.frame(pid=names(g)[i], tot=nrow(sub), ben=nrow(b), del=nrow(d), sumRS=lod)
        
        out <- rbind(out, tem)
    }
    return(out)
}

### deleterious alleles carried by each line
res <- get_per_line(g)
write.table(res, "cache/del_per_line.csv", sep=",", row.names=FALSE, quote=FALSE)

perline <- read.csv("cache/del_per_line.csv")
perline$percent <- with(perline, round(del/tot, 3))
with(perline, mean(del))
with(perline, range(del))

comp_two_lines <- function(g){
    
    out <- data.frame()
    for(i in 10:(ncol(g)-1)){
        message(sprintf("###>>> working on [ %s ]", names(g)[i]))
        for(j in (i+1):ncol(g)){
            
            ### complement
            c <- subset(g, g[, "del"] == g[, i] | g[, "del"] == g[, j])
            
            ### additive genetic load
            # ad <- subset(sub, (g[, "del"] != g[, i] & g[, i] != "N")  | (g[, "del"] == g[, j] & g[, j] != "N"))
            
            ### non-complement
            dd <- subset(g, (g[, "del"] != g[, i] & g[, i] != "N")  & (g[, "del"] != g[, j] & g[, j] != "N"))
            
            tem <- data.frame(pid1=names(g)[i], pid2=names(g)[j], comp=nrow(c), load=nrow(dd))
            out <- rbind(out, tem)
        }
    }
    return(out)
}

### deleterious alleles carried by each line
ids <- sort(names(g)[10:ncol(g)])
newg <- g[, c(names(g)[1:9], ids)]

res2 <- comp_two_lines(g=newg)
write.table(res2, "cache/del_complemenation.csv", sep=",", row.names=FALSE, quote=FALSE)

res22 <- subset(res2, pid1 != "B73")
write.table(res22, "cache/del_complemenation_noB73.csv", sep=",", row.names=FALSE, quote=FALSE)
range(res22$load)
#[1] 47219 77210
res22[which.min(res22$load),]
res22[which.max(res22$load),]


c1 <- read.csv("cache/del_complemenation.csv")
c2 <- read.csv("cache/del_complemenation_noB73.csv")


##########
comp_trait <- function(res2, pheno, excludeB73=TRUE){
    
    out <- data.frame()
    t <- as.character(unique(pheno$trait))
    if(excludeB73){
        res2 <- subset(res2, !(pid1 %in% "B73"))
    }
    
    for(i in 1:length(t)){
        subp <- subset(pheno, trait %in% t[i])
        subp <- merge(res2[, 3:5], subp, by="Hyb")
        if(nrow(subp) > 50){
            fit1 <- lm(valHyb ~ comp, data=subp)
            fit2 <- lm(valHyb ~ load, data=subp)
            fit3 <- lm(BPHmax ~ comp, data=subp)
            fit4 <- lm(BPHmax ~ load, data=subp)
            fit5 <- lm(pBPHmax ~ comp, data=subp)
            fit6 <- lm(pBPHmax ~ load, data=subp)
            fit7 <- lm(MPH ~ comp, data=subp)
            fit8 <- lm(MPH ~ load, data=subp)
            fit9 <- lm(pMPH ~ comp, data=subp)
            fit10 <- lm(pMPH ~ load, data=subp)
            
            
            tem <- data.frame(trait=t[i], pheno=rep(c("perse", "BPHmax", "pBPHmax", "MPH", "pMPH"), each=2),
                              geno=rep(c("comp", "load"), times=5))
            pval <- c(summary(fit1)$coefficients[2,4],
                      summary(fit2)$coefficients[2,4],
                      summary(fit3)$coefficients[2,4],
                      summary(fit4)$coefficients[2,4],
                      summary(fit5)$coefficients[2,4],
                      summary(fit6)$coefficients[2,4],
                      summary(fit7)$coefficients[2,4],
                      summary(fit8)$coefficients[2,4],
                      summary(fit9)$coefficients[2,4],
                      summary(fit10)$coefficients[2,4])
            
            est <- c(summary(fit1)$coefficients[2,1],
                      summary(fit2)$coefficients[2,1],
                      summary(fit3)$coefficients[2,1],
                      summary(fit4)$coefficients[2,1],
                      summary(fit5)$coefficients[2,1],
                      summary(fit6)$coefficients[2,1],
                      summary(fit7)$coefficients[2,1],
                      summary(fit8)$coefficients[2,1],
                      summary(fit9)$coefficients[2,1],
                      summary(fit10)$coefficients[2,1])
            tem$pval <- pval  
            tem$est <- est
  
        }else{
            stop("err! less < 50 ids!")
        }
        out <- rbind(out, tem)
    }
    return(out)
}


res2 <- read.csv("cache/del_complemenation.csv")
res2$Hyb <- paste(res2$pid1, res2$pid2, sep="_")

mean(res2$load)
range(res2$load)

idx1 <- which.min(res2$load)
idx2 <- which.max(res2$load)
res2[idx1,]
res2[idx2,]

pheno <- read.csv("data/hyb_heterosis.csv")
pheno$Hyb <- paste(pheno$Par1, pheno$Par2, sep="_")


res3 <- comp_trait(res2, pheno, excludeB73 = TRUE)
write.table(res3, "cache/complementation_rest_noB73.csv", sep=",", row.names=FALSE, quote=FALSE)

res3 <- read.csv("cache/complementation_rest_noB73.csv")
subset(res3, pheno %in% "MPH" & geno %in% "load" & pval < 0.05)

stable <- subset(res3, pheno %in% c("perse", "MPH") & geno %in% "load")
write.table(stable[, -3], "table/Table_hybload_pheno.csv", sep=",", row.names=FALSE, quote=FALSE)


res4 <- comp_trait(res2, pheno, excludeB73 = FALSE)
write.table(res3, "cache/complementation_rest.csv", sep=",", row.names=FALSE, quote=FALSE)

res3 <- read.csv("cache/complementation_rest.csv")
subset(res3, pheno %in% "perse" & geno %in% "load" & pval < 0.05)
stable2 <- subset(res3, pheno %in% c("perse", "MPH") & geno %in% "load")

stable$withB73 <- "no"
stable2$withB73 <- "yes"
st <- rbind(stable, stable)
write.table(st[, -3], "table/Table_hyb_load_pheno.csv", sep=",", row.names=FALSE, quote=FALSE)






