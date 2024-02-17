library(readr)
library(dplyr)
library(ggplot2)
library(survival)

rm(list=ls(all=TRUE))

setwd("/data/genetica/Tanner/data_curves/")
getwd()

## Data de niños
boysT12 <- read.table("boys_bmi_T12.csv",h=T,sep=",")
boysT13 <- read.table("boys_bmi_T13.csv",h=T,sep=",")
boysT23 <- read.table("boys_bmi_T23.csv",h=T,sep=",")
boysT34 <- read.table("boys_bmi_T34.csv",h=T,sep=",")
boysT24 <- read.table("boys_bmi_T24.csv",h=T,sep=",")

## Data de niñas
girlsT12 <- read.table("girls_bmi_T12.csv",h=T,sep=",")
girlsT13 <- read.table("girls_bmi_T13.csv",h=T,sep=",")
girlsT23 <- read.table("girls_bmi_T23.csv",h=T,sep=",")
girlsT34 <- read.table("girls_bmi_T34.csv",h=T,sep=",")
girlsT24 <- read.table("girls_bmi_T24.csv",h=T,sep=",")

## ancestría
ga <- read.table("datos_con_areaPesoAncestria2.txt",h=T)

library(survival)
library(icenReg)
library(ggplot2)

Curvas <- function(x,y,stage,...){
  ### análisis de la data	
  ## x : data _bmi_
  ## y: data ga
  ## stage : vector con transición de tanner según ancestría
  y_unique <- unique(y[,c("id","genero","GA_N")])	
  id <- match(x$CODE,y_unique$id)
  base <- x
  base$id <- y_unique[id,"id"]
  base$GA_N <- y_unique[id,"GA_N"]
  Lower <- base[,4]
  Upper <- base[,5]
  fit <- ic_par(cbind(Lower,Upper) ~ GA_N, data=base, model="ph")
  SE <- -solve(fit$hessian)
  
  lambda <- exp(coef(fit)[1])  # forma
  lambda_L <- exp(coef(fit)[1] - 1.96 * sqrt(SE[1,1]))
  lambda_U <- exp(coef(fit)[1] + 1.96 * sqrt(SE[1,1]))
  
  gamma <- exp(coef(fit)[2])   # escala
  gamma_L <- exp(coef(fit)[2] - 1.96 * sqrt(SE[2,2]))
  gamma_U <- exp(coef(fit)[2] + 1.96 * sqrt(SE[2,2]))
  
  beta <- coef(fit)[3]
  beta_L <- coef(fit)[3] - 1.96 * sqrt(SE[3,3])
  beta_U <- coef(fit)[3] + 1.96 * sqrt(SE[3,3])

  Survcurve <- function(tt,i,lambda,gamma,beta,group,...){
    h <- function(s,...){
      eta <- as.vector(group*beta)
      return(  exp(log(lambda)  + (lambda-1)*log(s) - lambda*log(gamma) + eta) )
    }
    out <- integrate(h,lower=0,upper=tt[i])$value
    return(exp(-out))
  }
  
  S_M <- S_M_L <- S_M_U <- c()
  S_E <- S_E_L <- S_E_U <- c()
  S_05_M <- S_05_L <- S_05_U <- c()
  
  tau <- seq(0,16,len=100)
  for(i in 1:length(tau)){
    S_M[i] <- Survcurve(tau,i,lambda,gamma,beta,1)
    S_M_L[i] <- Survcurve(tau,i,lambda_L,gamma_L,beta_L,1)
    S_M_U[i] <- Survcurve(tau,i,lambda_U,gamma_U,beta_U,1)
    S_E[i] <- Survcurve(tau,i,lambda,gamma,beta,0)
    S_E_L[i] <- Survcurve(tau,i,lambda_L,gamma_L,beta_L,0)
    S_E_U[i] <- Survcurve(tau,i,lambda_U,gamma_U,beta_U,0) 
  }
  
  # Extrapolate some points. To get median S(t), replace with 0.5. :
  print ("Mapuche")
  print (tau[which(S_M_L<0.75)[1]]) # 100% Mapuche
  print (tau[which(S_M<0.75)[1]]) # 100% Mapuche
  print (tau[which(S_M_U<0.75)[1]]) # 100% Mapuche
  
  print ("Europeo")
  print (tau[which(S_E_L<0.75)[1]]) # 100% European
  print (tau[which(S_E<0.75)[1]]) # 100% European
  print (tau[which(S_E_U<0.75)[1]]) # 100% European

  df <- data.frame( time=tau,S_M=S_M,S_M_L=S_M_L,
                    S_M_U = S_M_U, S_E=S_E, S_E_L=S_E_L, S_E_U=S_E_U
  )
  
  ggplot(df,aes(x=time)) +
    scale_color_manual(values=c("orangered3","royalblue3"))+
    scale_x_continuous(breaks = seq(5, 16, by = 1)) +
    scale_y_continuous(breaks = seq(0, 1, by = 0.1)) +
    geom_line(aes(y=S_M,color=stage[1])) +
    geom_line(aes(y=S_E,color=stage[2]))+
    labs(x="Age (years)",y ="S(t)",color="",title="")+
    ylim(0,1)  +
    xlim(5.5,16) +
    theme(axis.text=element_text(size=10),
          axis.title=element_text(size=11),
          legend.position=c(0.75,1.15),
          legend.text=element_text(size=8),
          legend.key= element_rect(fill = "transparent"),
          legend.key.size = unit(0.2, "cm"),
          axis.line.x = element_line(color="black", size = 0.4),
          axis.line.y = element_line(color="black", size = 0.4)) +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    theme(panel.background = element_blank())
  
c1=Curvas(boysT12,ga, c("Mapuche boys", "European boys"))
c2=Curvas(girlsT12,ga, c("Mapuche girls", "European girls"))
c3=Curvas(boysT23,ga, c("T2-T3 Mapuche", "T2-T3 European"))
c4=Curvas(girlsT23,ga, c("T2-T3 Mapuche", "T2-T3 European"))
c5=Curvas(boysT24,ga, c("T2-T4 Mapuche", "T2-T4 European"))
c6=Curvas(girlsT24,ga, c("T2-T4 Mapuche", "T2-T4 European"))
c7=Curvas(boysT34,ga, c("T3-T4 Mapuche", "T3-T4 European"))
c8=Curvas(girlsT34,ga, c("T3-T4 Mapuche", "T3-T4 European"))
c9=Curvas(boysT13,ga, c("T1-T3 Mapuche", "T1-T3 European"))
c10=Curvas(girlsT13,ga, c("T1-T3 Mapuche", "T1-T3 European"))

# Create plot:

pdf("/data/lucas/tanner/figures/surv_curves_ancestry.pdf", 
    width=7, height=2.5)
grid.arrange(c1,c2,nrow=1,ncol=2)

dev.off()
