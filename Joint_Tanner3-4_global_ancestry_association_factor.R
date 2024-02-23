library(data.table)
library(readr)
library(dplyr)
library(dtplyr)
library(nlme)
library(lme4)
library(survival)

setwd('~/proyectos/')
chr <- 22


df1f <- read_csv('bmi/two_stage_weibull/data/girls_bmi_T34.csv')
df1m <- read_csv('bmi/two_stage_weibull/data/boys_bmi_T34.csv')


f3 <- paste0('/data/genetica/tomas_localancestry/gocs_pel_ibs_pel95_phen904/rfmix15/rfmix15.chr',chr,'.stats.phen.pel_ibs.txt')
f4 <- 'tanner34/output_conjunto/summary_ga_fasoc.csv'

myformula <- function(){
  
  p1 = 'Surv(time, status) ~  ga + fasoc'
  p = p1 %>% as.formula()
  return(p)
}

loglike <- function(data, theta,...){
  
  n <- data %>% nrow
  p <- length(theta)
  bi <- data %>% select(b0,b1) %>% as.data.frame()
  li <- data$Age_LastT3 
  ri <- data$Age_FirstT4
  
  data <- data[,c(8,9)] %>% as.matrix()
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

myoutput <- function(param, s){
  
  tabla <- tibble(
    Sex = s,
    name = c('inv_scale','intercept','ga','fasoc'),
    value = param$par,
    std_error = sqrt(diag(solve(param$hessian))),
    z = param$par/std_error,
    p = 2*(1-pnorm(abs(z))) )
  
  return(tabla)
}


df3 <- fread(f3,
             sep='\t', 
             select =c(1,6),
             header=T)
colnames(df3)[c(1,2)] <- c('CODE','ga')

fmla = myformula()

dfm <- df1m %>% left_join(df3) %>% as_tibble
dff <- df1f %>% left_join(df3) %>% as_tibble

df=dfm
seg='boys'

df=dff
seg='girls'

n <- df %>% nrow

model.surv <- survreg(fmla, dist = 'weibull', data = df)

result <- optim(c(1/model.surv$scale,
                  -model.surv$coefficients[1]/model.surv$scale,
                  -model.surv$coefficients[2]/model.surv$scale,
                  -model.surv$coefficients[3]/model.surv$scale),
                loglike, hessian=TRUE, data = df)

# output_boys=myoutput(result,seg)
# output_girls=myoutput(result,seg)

# bind_rows(output_girls, output_boys) %>% 
#   write_csv(f4)

# result.girls=result
# result.boys=result

ll.boys=result.boys$value
ll.girls=result.girls$value
k=2
AIC.boys = 2*k -2*(-ll.boys)
AIC.girls = 2*k -2*(-ll.girls)
