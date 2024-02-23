library(data.table)
library(dplyr)
library(nlme)
library(lme4)
library(survival)
library(parallel)

args <- commandArgs(TRUE)
chr <- args[1]

link <- '~/proyectos/bmi/'
setwd(link)

df1 <- readr::read_rds('two_stage_weibull/data/df1.rds')

f2 <- paste0('detail/chr',chr,'_SNP_selected.tsv')
f3 <- paste0('/data/genetica/tomas_localancestry/gocs_pel_ibs_pel95_phen904/rfmix15.chr',chr,'.stats.phen.pel_ibs.txt')
f4 <- paste0('two_stage_weibull/output/chr',chr,'_2stage_weibull.txt')
f5 <- paste0('two_stage_weibull/output/chr',chr,'_2stage_weibullfailed.txt')

myformula <- function(i){
  
  p1 = 'Surv(time, status) ~'
  p2 = names(df1)[c(11,11+i,11+k+i)] %>% paste0(., collapse = '+')
  p3 <- '+ fasoc'
  p = paste0(p1,p2,p3) %>% as.formula()
  return(p)
}

loglike <- function(data,i=i,theta,...){
  n <- data %>% nrow
  p <- length(theta)
  bi <- data %>% select(b0,b1,b2) %>% as.data.frame()
  
  li <- data$Age_LasttT1 
  ri <- data$Age_FirstT2
  
  data <- data[,c(10,11,11+i,11+k+i)] %>% as.matrix()
  N <- dim(bi)[1]
  phi <- theta[1]	
  gammas <- theta[2:(p-1)]
  alpha <- theta[p]
  
  S <- function(tt,i,...){
    h <- function(s,...){
      eta <- as.vector(matrix(data[i,],nrow=1)%*%gammas)
      fval <- bi[i,1] + bi[i,2]*s + bi[i,3]*s^2 
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

myoutput <- function(param, i, s, n){
  
  tabla <- tibble(
    CHR = chr,
    SNP = df2$SNP[i],
    Sex = s,
    name = c('inv_scale','intercept','ga','geno_snp','locanc_snp','fasoc'),
    value = param$par,
    std_error = sqrt(diag(param$hessian)/n),
    z = param$par/sqrt(diag(param$hessian)/n),
    p = 2*(1-pnorm(abs(z))) )
  
  return(tabla)
}

write.table(tibble(t(c('CHR','SNP','Sex','name','value','std_error','z','p'))), file = f4, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)

df2 <- fread(f2,header = T,sep='\t')
k <- nrow(df2)


df3 <- fread(f3,
             sep='\t', 
             select =c(1,6, df2$POSITION, df2$POSITION+1),
             header=T)

colnames(df3)[c(1,2)] <- c('CODE','ga')

df1 <- df1 %>% left_join(df3)

df <- df1 %>% group_split(Sex)

f <- function(j){
  
  sex <- df[[j]]$Sex %>% unique()
  n <- df[[j]] %>% nrow
  
  for(i in 1:k){
    
    model.surv <- NULL
    result <- NULL
    fmla = myformula(i)
    
    tryCatch({
      
      model.surv <- survreg(fmla, dist = 'weibull', data = df[[j]])
      result <- optim(c(1/model.surv$scale,
                        -model.surv$coefficients[1]/model.surv$scale,
                        -model.surv$coefficients[2]/model.surv$scale,
                        -model.surv$coefficients[3]/model.surv$scale,
                        -model.surv$coefficients[4]/model.surv$scale,
                        -model.surv$coefficients[5]/model.surv$scale),
                      loglike, hessian=TRUE, data = df[[j]], i=i)
      
      write.table(myoutput(result,i, sex, n = n), file = f4, append=TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE)},
      error=function(e){ write.table(tibble(SNP=df2$SNP[i], POSTITION =df2$POSITION[i] ,Sex = sex), file = f5, append = TRUE, col.names = FALSE, row.names = FALSE, quote = FALSE) })
  }
}

r <- mclapply(1:2, f, mc.cores = 2)

