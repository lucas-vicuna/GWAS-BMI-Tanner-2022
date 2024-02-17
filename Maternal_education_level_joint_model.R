rm(list=ls(all=TRUE))

library(icenReg)
library(survival)

educ_madre <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/education_data_gocs",h=T,sep=",")
dim(educ_madre)
table(educ_madre$edu)
educ_madre$edu2 <- educ_madre$edu
educ_madre$edu2[educ_madre$edu2==0] <- "S/I"
educ_madre$edu2[educ_madre$edu2==4] <- "Media incompleta"
educ_madre$edu2[educ_madre$edu2==123] <- "Sin estudios, Básica incompleta y completa"
educ_madre$edu2[educ_madre$edu2==567] <- "Media completa, Instituto profesional incompleto y completo"
educ_madre$edu2[educ_madre$edu2==8910] <- "Universitaria incompleta y completa, Posgrado"
head(educ_madre)
str(educ_madre)

ancestry <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/global_ancestry.txt",h=T)

loglike <- function(data,theta,base,random,n,...){
p <- length(theta)
bi <- random 
li <- base[,4]
ri <- base[,5]
N <- length(bi)
phi <- theta[1]	
gammas <- theta[2:(p-1)]
alpha <- theta[p]
S <- function(tt,i,...){
	h <- function(s,...){
    eta <- as.vector(data[i,]%*%gammas)
    fval <- bi[i]# random slope
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

## BOYS

### T1-T2
boys_T12 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/boys_bmi_T12.csv",h=T,sep=",")
ind_T12 <- match(boys_T12$CODE,educ_madre$CODE)
ind_T12_gan <- match(boys_T12$CODE,ancestry$ID)
data_T12 <- data.frame(boys_T12[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T12], 
'gan' = ancestry$GAN[ind_T12_gan])
head(data_T12)

## puntos iniciales

time_T12 <- boys_T12$time
fasoc_T12 <- boys_T12$fasoc
status_T12 <- rep(1,length(time_T12))
edu_t12 <- data_T12$edu
gan_t12 <- data_T12$gan
model.surv_T12 <- survreg(Surv(time_T12,status_T12)~edu_t12 + gan_t12+fasoc_T12,
dist="weibull")
model.surv_T12$coefficients
model.surv_T12$scale

X_T12 <- model.matrix(model.surv_T12)[,-7]
n_T12 <- dim(X_T12)[1]

param_T12 <- optim(c(1/model.surv_T12$scale,
-model.surv_T12$coefficients[1]/model.surv_T12$scale,
-model.surv_T12$coefficients[2]/model.surv_T12$scale,
-model.surv_T12$coefficients[3]/model.surv_T12$scale,
-model.surv_T12$coefficients[4]/model.surv_T12$scale,
-model.surv_T12$coefficients[5]/model.surv_T12$scale,
-model.surv_T12$coefficients[6]/model.surv_T12$scale,
-model.surv_T12$coefficients[7]/model.surv_T12$scale),
loglike,hessian=TRUE, data=X_T12,base=boys_T12,
random=boys_T12[,3], n = n_T12)

Tabla_T12 <- cbind('value' = param_T12$par,
'z' = param_T12$par/sqrt(diag(solve(param_T12$hessian))),
'p'=2*(1-pnorm(abs(param_T12$par/sqrt(diag(solve(param_T12$hessian)))))))
Tabla_T12

### T2-T3
boys_T23 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/boys_bmi_T23.csv",h=T,sep=",")
head(boys_T23)
ind_T23 <- match(boys_T23$CODE,educ_madre$CODE)
ind_T23_gan <- match(boys_T23$CODE,ancestry$ID)
data_T23 <- data.frame(boys_T23[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T23],
'gan' = ancestry$GAN[ind_T23_gan])
head(data_T23)
## puntos iniciales

time_T23 <- boys_T23$time
fasoc_T23 <- boys_T23$fasoc
status_T23 <- rep(1,length(time_T23))
edu_t23 <- data_T23$edu
gan_t23 <- data_T23$gan

model.surv_T23 <- survreg(Surv(time_T23,status_T23)~edu_t23 + gan_t23+fasoc_T23,
dist="weibull")
model.surv_T23$coefficients
model.surv_T23$scale

X_T23 <- model.matrix(model.surv_T23)[,-7]
n_T23 <- dim(X_T23)[1]

param_T23 <- optim(c(1/model.surv_T23$scale,
-model.surv_T23$coefficients[1]/model.surv_T23$scale,
-model.surv_T23$coefficients[2]/model.surv_T23$scale,
-model.surv_T23$coefficients[3]/model.surv_T23$scale,
-model.surv_T23$coefficients[4]/model.surv_T23$scale,
-model.surv_T23$coefficients[5]/model.surv_T23$scale,
-model.surv_T23$coefficients[6]/model.surv_T23$scale,
-model.surv_T23$coefficients[7]/model.surv_T23$scale),
loglike,hessian=TRUE, data=X_T23,base=boys_T23,
random=boys_T23[,3], n = n_T23)

Tabla_T23 <- cbind('value' = param_T23$par,
'z' = param_T23$par/sqrt(diag(solve(param_T23$hessian))),
'p'=2*(1-pnorm(abs(param_T23$par/sqrt(diag(solve(param_T23$hessian)))))))
Tabla_T23

### T2-T4
boys_T24 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/boys_bmi_T24.csv",h=T,sep=",")
head(boys_T24)
ind_T24 <- match(boys_T24$CODE,educ_madre$CODE)
ind_T24_gan <- match(boys_T24$CODE,ancestry$ID)
data_T24<- data.frame(boys_T24[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T24],
'gan' = ancestry$GAN[ind_T24_gan])
head(data_T24)
## puntos iniciales

time_T24 <- boys_T24$time
fasoc_T24 <- boys_T24$fasoc
status_T24 <- rep(1,length(time_T24))
edu_t24 <- data_T24$edu
gan_t24 <- data_T24$gan
model.surv_T24 <- survreg(Surv(time_T24,status_T24)~edu_t24 + gan_t24+fasoc_T24,
dist="weibull")
model.surv_T24$coefficients
model.surv_T24$scale

X_T24 <- model.matrix(model.surv_T24)[,-7]
n_T24 <- dim(X_T24)[1]

param_T24 <- optim(c(1/model.surv_T24$scale,
-model.surv_T24$coefficients[1]/model.surv_T24$scale,
-model.surv_T24$coefficients[2]/model.surv_T24$scale,
-model.surv_T24$coefficients[3]/model.surv_T24$scale,
-model.surv_T24$coefficients[4]/model.surv_T24$scale,
-model.surv_T24$coefficients[5]/model.surv_T24$scale,
-model.surv_T24$coefficients[6]/model.surv_T24$scale,
-model.surv_T24$coefficients[7]/model.surv_T24$scale),
loglike,hessian=TRUE, data=X_T24,base=boys_T24,
random=boys_T24[,3], n = n_T24)

Tabla_T24 <- cbind('value' = param_T24$par,
'z' = param_T24$par/sqrt(diag(solve(param_T24$hessian))),
'p'=2*(1-pnorm(abs(param_T24$par/sqrt(diag(solve(param_T24$hessian)))))))
Tabla_T24

### T3-T4
boys_T34 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/boys_bmi_T34.csv",h=T,sep=",")
head(boys_T34)
ind_T34 <- match(boys_T34$CODE,educ_madre$CODE)
ind_T34_gan <- match(boys_T34$CODE,ancestry$ID)
data_T34<- data.frame(boys_T34[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T34],
'gan' = ancestry$GAN[ind_T34_gan])
head(data_T34)
## puntos iniciales

time_T34 <- boys_T34$time
fasoc_T34 <- boys_T34$fasoc
status_T34 <- rep(1,length(time_T34))
edu_t34 <- data_T34$edu
gan_t34 <- data_T34$gan
model.surv_T34 <- survreg(Surv(time_T34,status_T34)~edu_t34 +gan_t34+ fasoc_T34,
dist="weibull")
model.surv_T34$coefficients
model.surv_T34$scale

X_T34 <- model.matrix(model.surv_T34)[,-7]
n_T34 <- dim(X_T34)[1]

param_T34 <- optim(c(1/model.surv_T34$scale,
-model.surv_T34$coefficients[1]/model.surv_T34$scale,
-model.surv_T34$coefficients[2]/model.surv_T34$scale,
-model.surv_T34$coefficients[3]/model.surv_T34$scale,
-model.surv_T34$coefficients[4]/model.surv_T34$scale,
-model.surv_T34$coefficients[5]/model.surv_T34$scale,
-model.surv_T34$coefficients[6]/model.surv_T34$scale,
-model.surv_T34$coefficients[7]/model.surv_T34$scale),
loglike,hessian=TRUE, data=X_T34,base=boys_T34,
random=boys_T34[,3], n = n_T34)

Tabla_T34 <- cbind('value' = param_T34$par,
'z' = param_T34$par/sqrt(diag(solve(param_T34$hessian))),
'p'=2*(1-pnorm(abs(param_T34$par/sqrt(diag(solve(param_T34$hessian)))))))
Tabla_T34

## GIRLS

### T1-T2
girls_T12 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/girls_bmi_T12.csv",h=T,sep=",")
ind_T12G <- match(girls_T12$CODE,educ_madre$CODE)
ind_T12G_gan <- match(girls_T12$CODE,ancestry$ID)
data_T12G <- data.frame(girls_T12[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T12G],
'gan' = ancestry$GAN[ind_T12G_gan])
head(data_T12G)
## puntos iniciales

time_T12G <- girls_T12$time
fasoc_T12G <- girls_T12$fasoc
status_T12G <- rep(1,length(time_T12G))
edu_t12G <- data_T12G$edu
gan_t12G <- data_T12G$gan
model.surv_T12G <- survreg(Surv(time_T12G,status_T12G)~edu_t12G + gan_t12G +fasoc_T12G,
dist="weibull")
model.surv_T12G$coefficients
model.surv_T12G$scale

X_T12G <- model.matrix(model.surv_T12G)[,-7]
n_T12G <- dim(X_T12G)[1]

param_T12G <- optim(c(1/model.surv_T12G$scale,
-model.surv_T12G$coefficients[1]/model.surv_T12G$scale,
-model.surv_T12G$coefficients[2]/model.surv_T12G$scale,
-model.surv_T12G$coefficients[3]/model.surv_T12G$scale,
-model.surv_T12G$coefficients[4]/model.surv_T12G$scale,
-model.surv_T12G$coefficients[5]/model.surv_T12G$scale,
-model.surv_T12G$coefficients[6]/model.surv_T12G$scale,
-model.surv_T12G$coefficients[7]/model.surv_T12G$scale),
loglike,hessian=TRUE, data=X_T12G,base=girls_T12,
random=girls_T12[,3], n = n_T12G)

Tabla_T12G <- cbind('value' = param_T12G$par,
'z' = param_T12G$par/sqrt(diag(solve(param_T12G$hessian))),
'p'=2*(1-pnorm(abs(param_T12G$par/sqrt(diag(solve(param_T12G$hessian)))))))
Tabla_T12G

### T2-T3
girls_T23 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/girls_bmi_T23.csv",h=T,sep=",")
head(girls_T23)
ind_T23G <- match(girls_T23$CODE,educ_madre$CODE)
ind_T23G_gan <- match(girls_T23$CODE,ancestry$ID)
data_T23G <- data.frame(girls_T23[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T23G],
'gan' = ancestry$GAN[ind_T23G_gan])
head(data_T23G)
## puntos iniciales

time_T23G <- girls_T23$time
fasoc_T23G <- girls_T23$fasoc
status_T23G <- rep(1,length(time_T23G))
edu_t23G <- data_T23G$edu
gan_t23G <- data_T23G$gan
model.surv_T23G <- survreg(Surv(time_T23G,status_T23G)~edu_t23G + gan_t23G +fasoc_T23G,
dist="weibull")
model.surv_T23G$coefficients
model.surv_T23G$scale

X_T23G <- model.matrix(model.surv_T23G)[,-7]
n_T23G <- dim(X_T23G)[1]

param_T23G <- optim(c(1/model.surv_T23G$scale,
-model.surv_T23G$coefficients[1]/model.surv_T23G$scale,
-model.surv_T23G$coefficients[2]/model.surv_T23G$scale,
-model.surv_T23G$coefficients[3]/model.surv_T23G$scale,
-model.surv_T23G$coefficients[4]/model.surv_T23G$scale,
-model.surv_T23G$coefficients[5]/model.surv_T23G$scale,
-model.surv_T23G$coefficients[6]/model.surv_T23G$scale,
-model.surv_T23G$coefficients[7]/model.surv_T23G$scale),
loglike,hessian=TRUE, data=X_T23G,base=girls_T23,
random=girls_T23[,3], n = n_T23G)

Tabla_T23G <- cbind('value' = param_T23G$par,
'z' = param_T23G$par/sqrt(diag(solve(param_T23G$hessian))),
'p'=2*(1-pnorm(abs(param_T23G$par/sqrt(diag(solve(param_T23G$hessian)))))))
Tabla_T23G

### T2-T4
girls_T24 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/girls_bmi_T24.csv",h=T,sep=",")
head(girls_T24)
ind_T24G <- match(girls_T24$CODE,educ_madre$CODE)
ind_T24G_gan <- match(girls_T24$CODE,ancestry$ID)
data_T24G<- data.frame(girls_T24[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T24G],
'gan' = ancestry$GAN[ind_T24G_gan])
head(data_T24G)
## puntos iniciales

time_T24G <- girls_T24$time
fasoc_T24G <- girls_T24$fasoc
status_T24G <- rep(1,length(time_T24G))
edu_t24G <- data_T24G$edu
gan_t24G <- data_T24G$gan
model.surv_T24G <- survreg(Surv(time_T24G,status_T24G)~edu_t24G + gan_t24G+fasoc_T24G,
dist="weibull")
model.surv_T24G$coefficients
model.surv_T24G$scale

X_T24G <- model.matrix(model.surv_T24G)[,-7]
n_T24G <- dim(X_T24G)[1]

param_T24G <- optim(c(1/model.surv_T24G$scale,
-model.surv_T24G$coefficients[1]/model.surv_T24G$scale,
-model.surv_T24G$coefficients[2]/model.surv_T24G$scale,
-model.surv_T24G$coefficients[3]/model.surv_T24G$scale,
-model.surv_T24G$coefficients[4]/model.surv_T24G$scale,
-model.surv_T24G$coefficients[5]/model.surv_T24G$scale,
-model.surv_T24G$coefficients[6]/model.surv_T24G$scale,
-model.surv_T24G$coefficients[7]/model.surv_T24G$scale),
loglike,hessian=TRUE, data=X_T24G,base=girls_T24,
random=girls_T24[,3], n = n_T24G)

Tabla_T24G <- cbind('value' = param_T24G$par,
'z' = param_T24G$par/sqrt(diag(solve(param_T24G$hessian))),
'p'=2*(1-pnorm(abs(param_T24G$par/sqrt(diag(solve(param_T24G$hessian)))))))
Tabla_T24G

### T3-T4
girls_T34 <- read.table("/Users/vjleiva/Simulaciones/Scripts/Simulacion_per/GOCS/Graficos/data/girls_bmi_T34.csv",h=T,sep=",")
head(girls_T34)
ind_T34G <- match(girls_T34$CODE,educ_madre$CODE)
ind_T34G_gan <- match(girls_T34$CODE,ancestry$ID)
data_T34G<- data.frame(girls_T34[,c(1,4,5)],'edu'=educ_madre$edu2[ind_T34G],
'gan' = ancestry$GAN[ind_T34G_gan])
head(data_T34G)
## puntos iniciales

time_T34G <- girls_T34$time
fasoc_T34G <- girls_T34$fasoc
status_T34G <- rep(1,length(time_T34G))
edu_t34G <- data_T34G$edu
gan_t34G <- data_T34G$gan

model.surv_T34G <- survreg(Surv(time_T34G,status_T34G)~edu_t34G + gan_t34G+fasoc_T34G,
dist="weibull")
model.surv_T34G$coefficients
model.surv_T34G$scale

X_T34G <- model.matrix(model.surv_T34G)[,-6]
n_T34G <- dim(X_T34G)[1]

param_T34G <- optim(c(1/model.surv_T34G$scale,
-model.surv_T34G$coefficients[1]/model.surv_T34G$scale,
-model.surv_T34G$coefficients[2]/model.surv_T34G$scale,
-model.surv_T34G$coefficients[3]/model.surv_T34G$scale,
-model.surv_T34G$coefficients[4]/model.surv_T34G$scale,
-model.surv_T34G$coefficients[5]/model.surv_T34G$scale,
-model.surv_T34G$coefficients[6]/model.surv_T34G$scale,
-model.surv_T34G$coefficients[7]/model.surv_T34G$scale),
loglike,hessian=TRUE, data=X_T34G,base=girls_T34,
random=girls_T34[,3], n = n_T34G)

Tabla_T34G <- cbind('value' = param_T34G$par,
'z' = param_T34G$par/sqrt(diag(solve(param_T34G$hessian))),
'p'=2*(1-pnorm(abs(param_T34G$par/sqrt(diag(solve(param_T34G$hessian)))))))
Tabla_T34G
