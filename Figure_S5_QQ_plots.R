library(tidyverse)

setwd('~/Genetica/')

f2 <- paste0('Tanner/paper enero 2024/output/qqplot_tanner_v2.pdf')

set.seed(1155)
quants <-seq(0,1,length=5001)[2:5000]

pdf(file=f2, width=5, height=10)

par(mfrow=c(4,2),oma=c(2,2,1,1), mar=c(3,3,3,3))
for(i in c(12,23,24,34)){
  
  f1 <- NULL
  df1 <- NULL
  xlab=''
  ylab=''
  
  f1 <- paste0('Tanner/Tanner',i,'/output/Tanner',i,'_SNP_P.csv')
  df1 <- read_delim(f1,delim = ';') %>% 
    filter(name=='geno_snp',Seg!='All') %>% 
    mutate(y=-log10(p))
  
  t1 <- str_sub(i,1,1)
  t2 <- str_sub(i,2,2)
  
  for(j in c('boys','girls')){
    
    titulo <- NULL
    p <- NULL
    y <- NULL
    df2 <- NULL
    p_quants <- NULL
    fit_quants <- NULL
    data_quants <- NULL
    
    
    title <- paste0('T',t1,'-T',t2,': ',  j)
    
    
    if(j=='girls'){
      df2 <- df1 %>% filter(Seg=='Female') %>% 
        filter(FREQ_female>0.01) 
    }
    if(j=='boys'){
      df2 <- df1 %>% filter(Seg=='Male') %>% 
        filter(FREQ_male>0.01)
    }
    
    p <- df2$p
    
    u <- qunif(quants)
    fit_quants <- -log10(u)
    
    p_quants <- quantile(p,quants)
    data_quants <- -log10(p_quants)
    
    plot(x=fit_quants,y=data_quants,xlab='',ylab='',main=title)
    abline(0,1,col='red')
    
  }}
mtext('Theoretical Quantile',side=1,line=0,outer=TRUE,cex=1.3)
mtext('Sample Quantile',side=2,line=0,outer=TRUE,cex=1.3,las=0)
dev.off()
