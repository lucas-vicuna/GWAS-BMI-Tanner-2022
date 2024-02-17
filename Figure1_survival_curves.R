### Survival curves

rm(list=ls(all=TRUE))

setwd("/data/genetica/Tanner/data_curves/")

## Data de niños
boysT12 <- read.table("boys_bmi_T12.csv",h=T,sep=",")
boysT23 <- read.table("boys_bmi_T23.csv",h=T,sep=",")
boysT34 <- read.table("boys_bmi_T34.csv",h=T,sep=",")
boysT24 <- read.table("boys_bmi_T24.csv",h=T,sep=",")

## Data de niñas
girlsT12 <- read.table("girls_bmi_T12.csv",h=T,sep=",")
girlsT23 <- read.table("girls_bmi_T23.csv",h=T,sep=",")
girlsT34 <- read.table("girls_bmi_T34.csv",h=T,sep=",")
girlsT24 <- read.table("girls_bmi_T24.csv",h=T,sep=",")

## ancestría
ga <- read.table("datos_con_areaPesoAncestria2.txt",h=T)

library(survival)
library(icenReg)
library(ggplot2)

Curvas <- function(x,y,...){
  ### análisis de la data	
  ## x : data _bmi_
  ## y: data ga
  y_unique <- unique(y[,c("id","genero","GA_N")])	
  id <- match(x$CODE,y_unique$id)
  base <- x
  base$id <- y_unique[id,"id"]
  Lower <- base[,4]
  Upper <- base[,5]
  fit <- ic_par(cbind(Lower,Upper) ~ 1, data=base, model="ph")
  SE <- -solve(fit$hessian)
  
  lambda <- exp(coef(fit)[1])  # forma
  lambda_L <- exp(coef(fit)[1] - 1.96 * sqrt(SE[1,1]))
  lambda_U <- exp(coef(fit)[1] + 1.96 * sqrt(SE[1,1]))
  
  gamma <- exp(coef(fit)[2])   # escala
  gamma_L <- exp(coef(fit)[2] - 1.96 * sqrt(SE[2,2]))
  gamma_U <- exp(coef(fit)[2] + 1.96 * sqrt(SE[2,2]))
  
  
  Survcurve <- function(tt,i,lambda,gamma,...){
    h <- function(s,...){
      return(  exp(log(lambda) + (lambda-1)*log(s) - lambda*log(gamma)) )
    }
    out <- integrate(h,lower=0,upper=tt[i])$value
    #out <- (tt[i]/gamma)^(lambda)
    return(exp(-out))
  }
  
  S_0 <- S_0_L <- S_0_U <- c()
  tau <- seq(0,16,len=100)
  for(i in 1:length(tau)){
    S_0[i] <- Survcurve(tau,i,lambda,gamma)
    S_0_L[i] <- Survcurve(tau,i,lambda_L,gamma_L)
    S_0_U[i] <- Survcurve(tau,i,lambda_U,gamma_U)
  }
  
  return(df <- data.frame(time=tau, Surv=S_0, Surv_L=S_0_L, Surv_U=S_0_U))
  
}


df12boys <- Curvas(boysT12, ga)
colnames(df12boys)=c("time_boysT12","Surv_boysT12","Surv_L_boysT12","Surv_U_boysT12")

df23boys <- Curvas(boysT23, ga)
colnames(df23boys)=c("time_boysT23","Surv_boysT23","Surv_L_boysT23","Surv_U_boysT23")

df34boys <- Curvas(boysT34, ga)
colnames(df34boys)=c("time_boysT34","Surv_boysT34","Surv_L_boysT34","Surv_U_boysT34")

df24boys <- Curvas(boysT24, ga)
colnames(df24boys)=c("time_boysT24","Surv_boysT24","Surv_L_boysT24","Surv_U_boysT24")

df12girls <- Curvas(girlsT12, ga)
colnames(df12girls)=c("time_girlsT12","Surv_girlsT12","Surv_L_girlsT12","Surv_U_girlsT12")

df23girls <- Curvas(girlsT23, ga)
colnames(df23girls)=c("time_girlsT23","Surv_girlsT23","Surv_L_girlsT23","Surv_U_girlsT23")

df34girls <- Curvas(girlsT34, ga)
colnames(df34girls)=c("time_girlsT34","Surv_girlsT34","Surv_L_girlsT34","Surv_U_girlsT34")

df24girls <- Curvas(girlsT24, ga)
colnames(df24girls)=c("time_girlsT24","Surv_girlsT24","Surv_L_girlsT24", "Surv_U_girlsT24")

df_f= cbind(df12boys,df23boys,df34boys,df24boys,
            df12girls,df23girls,df34girls,df24girls)

# Now get confidence intervals:
df34girls$Surv_U_girlsT34[which(df34girls$Surv_U_girlsT34<0.75)[1]]

# Make plots:
p1 = print(ggplot(df_f,aes(x=time_girlsT12)) +
             scale_color_manual(values=c("orangered3","green4","orange1","royalblue3"))+
             scale_x_continuous(breaks = seq(2, 16, by = 1)) +
             scale_y_continuous(breaks = seq(0, 1.1, by = 0.1)) +
             geom_line(aes(y=Surv_girlsT12,color="T1-T2")) +
             geom_line(aes(y=Surv_L_girlsT12),col="orangered3",lty=2)+
             geom_line(aes(y=Surv_U_girlsT12),col="orangered3",lty=2)+
             geom_line(aes(y=Surv_girlsT23,color="T2-T3"))+
             geom_line(aes(y=Surv_L_girlsT23),col="green4",lty=2)+
             geom_line(aes(y=Surv_U_girlsT23),col="green4",lty=2)+
             geom_line(aes(y=Surv_girlsT24,color="T2-T4"))+
             geom_line(aes(y=Surv_L_girlsT24),col="orange1",lty=2)+
             geom_line(aes(y=Surv_U_girlsT24),col="orange1",lty=2)+
             geom_line(aes(y=Surv_girlsT34,color="T3-T4"))+
             geom_line(aes(y=Surv_L_girlsT34),col="royalblue3",lty=2)+
             geom_line(aes(y=Surv_U_girlsT34),col="royalblue3",lty=2)+
             ylim(0,1)+
             xlim(5.5,16) +
             theme(axis.text=element_text(size=13),
                   axis.title=element_text(size=15),
                   legend.position=c(0.84,0.9),
                   legend.text=element_text(size=8),
                   legend.key= element_rect(fill = "transparent"),
                   legend.key.size = unit(0.2, "cm"),
                   axis.line.x = element_line(color="black", size = 0.7),
                   axis.line.y = element_line(color="black", size = 0.7)) +
             labs(x="Age girls (years)",y ="S(t)",color="",title="")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank())

p2 = print(ggplot(df_f,aes(x=time_boysT12)) +
             scale_color_manual(values=c("orangered3","green4","orange1","royalblue3"))+
             scale_x_continuous(breaks = seq(2, 16, by = 1)) +
             scale_y_continuous(breaks = seq(0, 1.1, by = 0.1)) +
             geom_line(aes(y=Surv_boysT12,color="T1-T2")) +
             geom_line(aes(y=Surv_L_boysT12),col="orangered3",lty=2)+
             geom_line(aes(y=Surv_U_boysT12),col="orangered3",lty=2)+
             geom_line(aes(y=Surv_boysT23,color="T2-T3"))+
             geom_line(aes(y=Surv_L_boysT23),col="green4",lty=2)+
             geom_line(aes(y=Surv_U_boysT23),col="green4",lty=2)+
             geom_line(aes(y=Surv_boysT24,color="T2-T4"))+
             geom_line(aes(y=Surv_L_boysT24),col="orange1",lty=2)+
             geom_line(aes(y=Surv_U_boysT24),col="orange1",lty=2)+
             geom_line(aes(y=Surv_boysT34,color="T3-T4"))+
             geom_line(aes(y=Surv_L_boysT34),col="royalblue3",lty=2)+
             geom_line(aes(y=Surv_U_boysT34),col="royalblue3",lty=2)+
             ylim(0,1)+
             xlim(5.5,16) +
             theme(axis.text=element_text(size=13),
                   axis.title=element_text(size=15),
                   legend.position=c(0.84,0.9),
                   legend.text=element_text(size=8),
                   legend.key= element_rect(fill = "transparent"),
                   legend.key.size = unit(0.2, "cm"),
                   axis.line.x = element_line(color="black", size = 0.7),
                   axis.line.y = element_line(color="black", size = 0.7)) +
             labs(x="Age boys (years)",y ="S(t)",color="",title="")) +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  theme(panel.background = element_blank())

# Create plot:

pdf("/data/lucas/tanner/figures/curves_T1toT4.pdf", 
    width=7, height=3)
grid.arrange(p2,p1,nrow=1,ncol=2)
dev.off()

