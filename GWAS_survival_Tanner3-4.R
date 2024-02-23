library(data.table)
library(dplyr)
library(icenReg)

args <- commandArgs(TRUE)
chr <- args[1]

link <- '~/proyectos/'
setwd(link)

f1 <- 'tanner34/data/Tanner_34.tsv'
f2 <- paste0('bmi/detail/chr',chr,'_SNP_selected.tsv')
f3 <- paste0('/data/genetica/tomas_localancestry/gocs_pel_ibs_pel95_phen904/rfmix15.chr',chr,'.stats.phen.pel_ibs.txt')
f4 <- paste0('tanner34/output/chr',chr,'_tanner34.txt')
f5 <- paste0('tanner34/output/chr',chr,'_tanner34failed.txt')

df1 <- fread(f1, sep='\t', header=T)

df2 <- fread(f2,header = T,sep='\t')
k <- nrow(df2)

df3 <- fread(f3,
             sep='\t', 
             select =c(1,6, df2$POSITION, df2$POSITION+1),
             header=T)

colnames(df3)[c(1,2)] <- c('CODE','ga')

df3 <- df1 %>% inner_join(df3)
rm(df1)

myformula <- function(i, seg){
  
  p1 = paste0('cbind(Age_LastT3 , Age_FirstT4) ~ ga + ', names(df3)[5+i+k])
  p2 = ifelse(seg == 'All',
              paste0('Sex + ','Sex:',names(df3)[5+i]),
              names(df3)[5+i])
  p = paste0(c(p1, p2), collapse = '+') %>% as.formula()
  return(p)
}

myoutput <- function(m,i, seg){
  
  param = summary(m)$summaryParameters %>% 
    as.data.frame() %>% 
    add_rownames('name') %>% 
    mutate(name = case_when(
      stringr::str_detect(name, "^N_") ~ 'locanc_snp',
      stringr::str_detect(name, "^GT_") ~ 'geno_snp',
      stringr::str_detect(name, "^SexMale:GT_") ~ 'SexMale:geno_snp',
      stringr::str_detect(name, "^SexFemale:GT_") ~ 'SexFemale:geno_snp',
      TRUE ~ name))
  
  param = param[,c(1,2,4,5,6)]
  colnames(param) <- c('name','estimate','std_error','z','p')
  
  tabla <- bind_cols(
    CHR = chr,
    SNP = df2$SNP[i],
    Seg = seg,
    param)  
  return(tabla)
}

mydf <- function(seg){
  
  if(seg == 'All') return(df3)
  if(seg == 'Male') return(dfm)
  if(seg == 'Female') return(dff)
}

dfm <- df3 %>% filter(Sex == 'Male')
dff <- df3 %>% filter(Sex == 'Female')

for(i in 1:k){
  for(j in c('Male','Female','All')){
    
    model <- NULL
    fmla = myformula(i, seg = j) 
    
    tryCatch({
      
      model <- ic_par(fmla, data = mydf(j))
      
      write.table(myoutput(model,i, seg = j), file = f4, append=TRUE, col.names = ifelse(i==1 & j=='Male',TRUE,FALSE), row.names = FALSE, quote = FALSE)},
      error=function(e){ write.table(tibble(SNP=df2$SNP[i], POSTITION =df2$POSITION[i],Seg = j), file = f5, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE) })
  }
}


