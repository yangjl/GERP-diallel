## Define parameters, functions:

numinds=1000
numloci=100
kneecap<-function(x){ if(x>1){ return(1) } else(return(x))} # turn loci from additive > dominant (two 1's == one 1)
bumpup<-function(x){ if(x!=1){ return(x) } else(return(2))} # turn loci from additive > dominant (two 1's == one 1)

library(ggplot2)
hetplot=function(a,b){  
  bob=data.frame(a,b)
  namex=deparse(substitute(a));
  namey=deparse(substitute(b));
  the_plot=ggplot(bob)+geom_point(aes(y=b,x=a)) + 
    geom_smooth(aes(y=b,x=a),method="lm") +   
    geom_text(aes(-Inf,Inf,label=paste("r2=",round(summary(lm(b~a))$r.squared,2))),vjust=2,hjust=-0.5) +
    xlab(namex)+ylab(namey);
  print(the_plot)
}

a_het=numeric()
a_phet=numeric()
d_het=numeric()
d_phet=numeric()
ad_het=numeric()
ad_phet=numeric()

for(run in 1:100){
  ## Make parents, F1
  
  #### Inbreds and phenotypes
  #There's no dominance in inbreds, so we just sum to get phenotypes, and multiply by two since homozygous for each locus.
  
  p1 <- lapply(1:numinds, function(X) sample(c(0,1),numloci,replace=T))
  p1_phenotype=2*unlist(lapply(p1,sum))
  p2 <- lapply(1:numinds, function(X) sample(c(0,1),numloci,replace=T))
  p2_phenotype=2*unlist(lapply(p2,sum))
  
  
  #Identify best parent
  best_p<-sapply(1:1000, function(i) max(p1_phenotype[i],p2_phenotype[i]))
  
  #### Make F1
  
  f1 <- lapply(1:numinds, function(X) p1[[X]]+p2[[X]])
  
  #purely dominant phenotype for the F1
  d1=f1
  d1<-lapply(1:numinds, function(X) sapply(1:numloci, function(Y) bumpup(d1[[X]][Y]) ))
  f1_d_phenotype<-unlist(lapply(d1,sum))
  
  ### Calculate heterosis for all three F1 types
  
  #heterosis dominant
  HPH_d<-f1_d_phenotype-best_p
  pHPH_d<-HPH_d/best_p
  
  
  a_het[run]=HPH_d
  a_phet[run]=pHPH_d
  ad_het[run]=HPH_d
  ad_phet[run]=pHPH_d
  d_het[run]=HPH_d
  d_phet[run]=pHPH_d
}


#a_het[run]=cor(f1_a_phenotype,HPH_d)
#a_phet[run]=cor(f1_a_phenotype,pHPH_d)
ad_het[run]=cor(f1_ad_phenotype,HPH_d)
ad_phet[run]=cor(f1_ad_phenotype,pHPH_d)
d_het[run]=cor(f1_d_phenotype,HPH_d)
d_phet[run]=cor(f1_d_phenotype,pHPH_d)

hetplot(a_het, f1_d_phenotype)
