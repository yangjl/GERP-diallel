### Jinliang Yang
### 04/13/2017
### cal deleterious per line and complemantation


library(data.table)
g <- fread("largedata/gerpsnp_v3_345176_del.csv", data.table=FALSE)
## Note col "del" is the major from the alignment, should be the beneficial alleles

get_per_line <- function(g){
    out <- data.frame()
    for(i in 10:ncol(g)){
        sub <- subset(g, g[, i] != "N")
        b <- subset(sub, g[, "del"] == g[, i])
        d <- subset(sub, g[, "del"] != g[, i])
        tem <- data.frame(pid=names(g)[i], tot=nrow(sub), ben=nrow(b), del=nrow(d))
        
        out <- rbind(out, tem)
    }
    return(out)
}

### deleterious alleles carried by each line
res <- get_per_line(g)
write.table(res, "cache/del_per_line.csv", sep=",", row.names=FALSE, quote=FALSE)

comp_two_lines <- function(g){
    out <- data.frame()
    for(i in 10:(ncol(g)-1)){
        message(sprintf("###>>> working on [ %s ]", names(g)[i]))
        for(j in (i+1):ncol(g)){
            c <- subset(sub, g[, "del"] == g[, i] | g[, "del"] == g[, j])
            tem <- data.frame(pid1=names(g)[i], pid2=names(g)[j], comp=nrow(c))
            out <- rbind(out, tem)
        }
    }
    return(out)
}

### deleterious alleles carried by each line
res2 <- comp_two_lines(g)
write.table(res2, "cache/del_complemenation.csv", sep=",", row.names=FALSE, quote=FALSE)


