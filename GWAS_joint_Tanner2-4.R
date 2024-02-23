library(data.table)
library(readr)
library(dplyr)
library(dtplyr)
library(nlme)
library(lme4)
library(survival)
library(parallel)

args <- commandArgs(TRUE)
chr <- args[1]
chr <- as.numeric(chr)

setwd('~/proyectos/')

df1f <- read_csv('bmi/two_stage_weibull/data/girls_bmi_T24.csv')
df1m <- read_csv('bmi/two_stage_weibull/data/boys_bmi_T24.csv')

f2 <- paste0('bmi/detail/chr',chr,'_SNP_selected.tsv')
f22 <- paste0('tanner24/output/chr',chr,'_tanner24.txt')
f23 <- 'data/freq_loc.tsv'

f3 <- paste0('/data/genetica/tomas_localancestry/gocs_pel_ibs_pel95_phen904/rfmix15/rfmix15.chr',chr,'.stats.phen.pel_ibs.txt')
f4 <- paste0('tanner24/output_conjunto/chr',chr,'_conjunto.txt')
f5 <- paste0('tanner24/output_conjunto/chr',chr,'_conjuntofailed.txt')

myformula <- function(i,df,k){
  
  p1 = 'Surv(time, status) ~  ga +'
  p2 = names(df)[c(9+i,9+k+i)] %>% paste0(., collapse = '+')
  p3 = '+ fasoc'
  p = paste0(p1,p2,p3) %>% as.formula()
  return(p)
}

loglike <- function(data,i=i,k,theta,...){
  
  n <- data %>% nrow
  p <- length(theta)
  bi <- data %>% select(b0,b1) %>% as.data.frame()
  li <- data$Age_LastT2 
  ri <- data$Age_FirstT4
  
  data <- data[,c(8,9,9+i,9+k+i)] %>% as.matrix()
  N <- dim(bi)[1]
  phi <- theta[1]	
  gammas <- theta[2:(p-1)]
  alpha <- theta[p]
  
  S <- function(tt,i,...){
    h <- function(s,...){
      eta <- as.vector(matrix(data[i,],nrow=1)%*%gammas)
      fval <- bi[i,2]*s 
      return(exp(log(phi)+(phi-1)*log(s)+eta+alpha*fval))
    }
    out <- integrate(h,lower=0,upper=tt[i])$value
    return(exp(-out))
  }
  S.li <- S.ri <- c()
  for(i in 1:n){
    S.li[i] <- S(li,i)
    S.ri[i] <- S(ri,i)
  }
  
  aux <- sum(log(S.li-S.ri))
  
  return(-aux)
}

myk <- function(seg){
  if(seg=='Male') return(km)
  if(seg=='Female') return(kf)
}

mydf <- function(seg){
  if(seg=='Male') return(dfm)
  if(seg=='Female') return(dff)
}

mydf2 <- function(seg){
  if(seg=='Male') return(df2m)
  if(seg=='Female') return(df2f)
}

myoutput <- function(param, i, s, n){
  
  tabla <- tibble(
    CHR = chr,
    SNP = mydf2(s)$SNP[i],
    Sex = s,
    name = c('inv_scale','intercept','ga','geno_snp','locanc_snp','fasoc'),
    value = param$par,
    std_error = sqrt(diag(solve(param$hessian))),
    z = param$par/std_error,
    p = 2*(1-pnorm(abs(z))) )
  
  return(tabla)
}

write.table(tibble(t(c('CHR','SNP','Sex','name','value','std_error','z','p'))), file = f4, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

df2 <- fread(f2,header = T,sep='\t')
df22 <- fread(f22,header = T,sep=' ')
df23 <- fread(f23,header = T,sep='\t')

df22 <- df22[(p<=1e-1) & (name=='geno_snp'),] %>% as_tibble()
df24 <- df22 %>% inner_join(df23 %>% as_tibble())
df24 <- df24 %>% 
  filter((Seg=='Male' & FREQ_male>0.01)|(Seg=='Female' & FREQ_female>0.01))

df2 <- df2 %>% as_tibble() %>% inner_join(df24) 

df2m <- df2 %>% filter(Seg == 'Male')
km <- nrow(df2m)

df2f <- df2 %>% filter(Seg == 'Female')
kf <- nrow(df2f)

df3m <- fread(f3,
              sep='\t', 
              select =c(1,6, df2m$POSITION, df2m$POSITION+1),
              header=T)
colnames(df3m)[c(1,2)] <- c('CODE','ga')

df3f <- fread(f3,
              sep='\t', 
              select =c(1,6, df2f$POSITION, df2f$POSITION+1),
              header=T)
colnames(df3f)[c(1,2)] <- c('CODE','ga')

dfm <- df1m %>% left_join(df3m)
dff <- df1f %>% left_join(df3f)
rm(df1m, df1f, df3m,df3f)

f <- function(seg){
  
  df <- mydf(seg)
  n <- df %>% nrow
  k <- myk(seg)
  
  for(i in 1:k){
    
    model.surv <- NULL
    result <- NULL
    fmla = myformula(i,df,k)
    
    tryCatch({
      
      model.surv <- survreg(fmla, dist = 'weibull', data = df)
      result <- optim(c(1/model.surv$scale,
                        -model.surv$coefficients[1]/model.surv$scale,
                        -model.surv$coefficients[2]/model.surv$scale,
                        -model.surv$coefficients[3]/model.surv$scale,
                        -model.surv$coefficients[4]/model.surv$scale,
                        -model.surv$coefficients[5]/model.surv$scale),
                      loglike, hessian=TRUE, data = df, i=i, k=k)
      
      write.table(myoutput(result,i, seg, n = n), file = f4, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)},
      error=function(e){ write.table(tibble(SNP=mydf2(seg)$SNP[i], POSTITION =mydf2(seg)$POSITION[i] ,Sex = seg), file = f5, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE) })
  }
}

r <- mclapply(c('Male','Female'), f, mc.cores = 2)

