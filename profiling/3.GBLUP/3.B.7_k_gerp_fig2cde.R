### Jinliang Yang
### 9/24/2015


library(data.table)
source("lib/multiplot.R")

plot_adk_gerp <- function(outfile, getpdf){
    
    dat <- fread("largedata/gerpdat0x.csv", data.table=FALSE)
    
    
    dat$trait <- toupper(dat$trait)
    med2 <- read.csv("cache/loh_pMPH_median.csv")
    out <- read.csv("cache/eff_adk_0x.csv")
    #cols <- wes_palette(7, name = "Zissou", type = "continuous")
    cols <- c("#f6546a", "#daa520", "#00ff00", "#66cdaa", "#3b5998", "#8a2be2", "#ff00ff")
    theme_set(theme_grey(base_size = 18)) 
    
    getlty <- function(df, eff, cutoff=0.05){
        df$l <- 2
        if(nrow(df[df[, eff] < cutoff, ]) >0) df[df[, eff] < cutoff, ]$l <- 1
        return(df)
    }
    
    lty1 <- getlty(df=out, eff="effa", cutoff=0.05)$l
    p1 <- ggplot(dat, aes(x=RS, y=Effect_A, colour=factor(trait, levels=med2$trait),
                          linetype=factor(trait, levels=med2$trait))) +
        labs(colour="Traits") +
        theme_bw() +
        xlab("GERP Score") +
        ylab("Additive Effect") +
        scale_color_manual(values=cols) +
        scale_linetype_manual(values=lty1) +
        guides(colour=FALSE, linetype=FALSE) +
        geom_smooth(method="gam", size=1.3) +
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    lty2 <- getlty(df=out, eff="effd", cutoff=0.05)$l
    p2 <- ggplot(dat, aes(x=RS, y=Effect_D, colour=factor(trait, levels=med2$trait),
                          linetype=factor(trait, levels=med2$trait))) +
        labs(colour="Traits") +
        theme_bw() +
        xlab("GERP Score") +
        ylab("Dominant Effect") +
        scale_color_manual(values=cols) +
        scale_linetype_manual(values=lty2) +
        guides(colour=FALSE, linetype=FALSE) +
        geom_smooth(method="gam", size=1.3) +  # Add linear regression line 
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize))
    
    lty3 <- getlty(df=out, eff="effk", cutoff=0.05)$l
    p3 <- ggplot(dat, aes(x=RS, y=k, colour=factor(trait, levels=med2$trait),
                          linetype=factor(trait, levels=med2$trait))) +
        labs(colour="Traits") +
        theme_bw() +
        xlab("GERP Score") +
        ylab("Degree of Domiance (k)") +
        scale_color_manual(values=cols) +
        scale_linetype_manual(values=lty3) +
        guides(colour=FALSE, linetype=FALSE) +
        theme(axis.text.y = element_text(angle = 90, hjust = 1),
              axis.text=element_text(size=fsize),
              axis.title=element_text(size=fsize, face="bold"),
              legend.title = element_text(size=fsize, face="bold"),
              legend.text = element_text(size=fsize)) +
        geom_smooth(method="gam", size=1.3)   # Add linear regression line 
    
    multiplot(p1, p2, p3, cols=3)
    if(getpdf == TRUE){
        pdf(outfile, width=13, height=4)
        multiplot(p1, p2, p3, cols=3)
        dev.off()
    }
    
}

####
plot_adk_gerp(outfile="graphs/Fig2c-e_adk.pdf", getpdf=TRUE)


