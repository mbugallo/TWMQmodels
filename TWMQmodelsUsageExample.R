
#####################################################
###     Temporal M-quantile models and robust     ###
###     bias-corrected small area predictors      ###
###                                               ###
### -- Both continuos and discrete covariates --  ###
##         (we work with dummy variables)         ###
#####################################################

# Model: y_dtj = 5*x1_dtj + 10*I(x2_dtj==0) + 15*I(x2_dtj==1) 
#              + 20*I(x2_dtj==2) + u_d + u_t + e_dtj
# x1_dtj ~ LogNormal(1, 0.5)
# x2_dtj ~ Binom(2,0.3)
# u_d ~ N(0,sigma.u), sigma.u=3  
# u_t ~ N_T(0,SIGMA), rho=0.2, sigma.t=1
# e_dtj ~ N(0,sigma.e), sigma.e=6 

# x1_dtj is a positive continuous variable. It could be the household 
# expenditure on housing, food, etc.

# x2_dtj is a categorical variable with 3 factors. It could be the 
# age group, the educational level, etc.

# y_dtj is a continuous variable. It could be the equivalent disposable
# income per unit of consumption.

# d: index for the regioncode (d=1,...,D). D = 40
# t: index for the timecode (t=1,...,T). T = 10 

# Remove all stored variables from memory
rm(list=ls())

# Set a seed
set.seed(1234)

# Load auxiliary functions stored in our own scripts
source("QRLM.R")
source("QRLMweights.R")

### HUBER FUNCTION

hub.psi <- function(x, k){ ifelse(abs(x) <= k, x, sign(x) * k) }
der.hub.psi <- function(x, k){ ifelse(abs(x) <= k, 1, 0) }

### MEDIAN ABSOLUTE DEVIATION

fun.MAD<-function(x){ median(abs( x-median(x)) )/0.6745 }

### LIBRARIES

library(MASS)
library(fpp2) 
library(dplyr)
library(fastDummies)
library(agricolae)

# DESIGN
m=40 # number of areas
T=10 # number of time periods

# Population sizes
Nit=matrix(100,nrow=m,ncol=T)
Ni=rowSums(Nit)
N=sum(Ni)

# Samples sizes
nit=matrix(5,nrow=m,ncol=T)
ni=rowSums(nit)
n=sum(ni)

# Area-level random effects
ui<-rnorm(m,0,sqrt(3))

# Time-level random effects
COV.AR <- function(p, rho) {
  i <- matrix(1:p, p, p, byrow = TRUE)
  (1 / (1 - rho^2)) * rho^abs(i - t(i))}

S1<-COV.AR(p=T,rho=0.2)

# Covariates
xdtj.discr<-rbinom(N,2,0.3)
xdtj.cont<-rlnorm(N,1,0.5)

# Model errors
edtj <- rnorm(N,0,sqrt(6))

# Target variables
ydtj<-xdtj.cont*5+as.matrix(dummy_cols(xdtj.discr))[,-1]%*%c(10,15,20)+rep(ui,Ni)+
		rep(mnormt::rmnorm(1, rep(0,dim(S1)[2]), S1),each=N/T)+edtj

# Region and timecodes
regioncode<-rep(1:m,Ni)
timecode<-rep(1:T,each=N/(m*T),times=m)

# Final data.frame
df<-data.frame(merge(data.frame('id'=1:N, 'regioncode'=regioncode),
	data.frame('id'=1:N,'timecode'=timecode),by='id'), ydtj, xdtj.cont, xdtj.discr)
head(df)

# Real values
truet<-aggregate(ydtj, by=list(timecode,regioncode),mean)[,3]

# Sample and non-sample sets
s<-NULL
for (j in 1:m){ 
   for (t in 1:T){
      s<-sort(c(s,sample(df$id[df$regioncode==j & df$timecode==t],nit[j,t])))  
   }
}
df.s<-df[s,]
df.r<-df[-s,]

# Sample data  
x.s<-cbind(xdtj.cont[s],as.matrix(dummy_cols(xdtj.discr[s]))[,-1])
y.s<-ydtj[s]
regioncode.s<-regioncode[s]
timecode.s<-timecode[s]
n.t<-colSums(nit)
n.ts<-c(0,cumsum(n.t))
  
# Non-sample data
x.r<-cbind(xdtj.cont[-s],as.matrix(dummy_cols(xdtj.discr[-s]))[,-1])
regioncode.r<-regioncode[-s]
timecode.r<-timecode[-s]


#############################
###     Area-level MQ     ###
#############################  

# 1º Quantile regression
tau=sort(c(seq(0.006,0.99,0.045),0.5,0.994,0.01,0.02,0.96,0.98))
mod<-QRLM(x=x.s, y=y.s, q=tau, maxit=25, k = 1.345)

## Linear Interpolation
qo<-matrix(c(gridfitinter(y.s,mod$fitted.values,mod$q.values)),nrow=n,ncol=1)
qmat<-matrix(c(qo,regioncode.s),nrow=n,ncol=2)
mqo<-aggregate(qmat[,1],by=list(qmat[,2]),mean)[,2]

# 2º Quantile regression
mod.SAE<-QRLM(x=x.s, y=y.s, q=mqo, maxit=25, k = 1.345)


########################################
###     T-Weighted Area-level MQ     ###
########################################

# Borrowing strength from time
wt.results <- compute.wt(x.s, y.s, regioncode.s, timecode.s, mod.SAE, x.r)

W <- wt.results$W
P <- wt.results$P

# Time-weighted M-quantile models
mod.SAE2t<-mod.50t<-list()

for (t in 1:T){
  if ((t == 1)|(P==0)){ incl.t = t } else { incl.t = t-P+1; if (incl.t<1){ incl.t=1 } }
  cond.t<-(timecode.s<=t & timecode.s>=incl.t)
  
  # 2º Quantile regression
  # SAE models
  mod.SAE2t[[t]]<-QRLMweights(x=x.s[cond.t,], y=y.s[cond.t], q=mqo, maxit=25, k = 1.345,
                              w1=W[(n.ts[incl.t]+1):(n.ts[t+1]), (n.ts[incl.t]+1):(n.ts[t+1])])
  
  # Median MQ regression needed to calculate the MSE of the predictors
  mod.50t[[t]]<-QRLMweights(x=x.s[cond.t,], y=y.s[cond.t], q=0.5, maxit=25, k = 1.345,
                            w1=W[(n.ts[incl.t]+1):(n.ts[t+1]), (n.ts[incl.t]+1):(n.ts[t+1])])
}

# Regression parameters
beta<-list()
for(t in 1:T){ beta[[t]]<-t(mod.SAE2t[[t]]$coef) }
df.beta<-data.frame(List=unlist(beta),Coef=rep(1:4,each=m, times=T), 
                    Time=rep(1:T,each=m*4), Region=rep(1:m, times=T))

boxplot(df.beta$List[df.beta$Coef==1]~df.beta$Time[df.beta$Coef==1])
boxplot(df.beta$List[df.beta$Coef==2]~df.beta$Time[df.beta$Coef==2])
boxplot(df.beta$List[df.beta$Coef==3]~df.beta$Time[df.beta$Coef==3])
boxplot(df.beta$List[df.beta$Coef==4]~df.beta$Time[df.beta$Coef==4])


################################################################# 
###     Prediction and MSE estimation for TWMQ predictors     ###
################################################################# 

# TWMQ.mse: External function that has all calculations implemented
TWMQ.results<-TWMQ.mse(x.s, y.s, x.r, mod.SAE2t, mod.50t, Nit, nit, regioncode.s, timecode.s,
                       regioncode.r, timecode.r, P, seq(0,10,by=0.1))

# Contents of the object TWMQ.results
## TWMQ: Plug-type predictor
## TWMQ.BC: Optimal robust bias-corrected plug-type predictor
## mse1, mse2, mse3, mse4: MSE estimation for the plug-type predictor
## bias.sq: Bias estimation
## mse5.BC, mse6.BC: MSE estimation for the optimal robust bias-corrected plug-type predictor
## k.dt: optimal area-time specific robustness parameters

str(TWMQ.results)


# Predicted values versus real values
plot(TWMQ.results$TWMQ.BC, truet, xlab='MQB', ylab='True')
abline(a=0, b=1, col='red', lty=1)

# Coefficients of variation
plot(100*sqrt(TWMQ.results$mse6.BC)/TWMQ.results$TWMQ.BC, ylab='CV (%)')


#############################
###   Outlier detection   ###
#############################

friedman<-data.frame('respuesta'=matrix(TWMQ.results$k.dt, ncol=1), 
                     'area'=rep(1:m,times=T), 'time'=rep(1:T,each=m))

# Area-level outliers         
friedman.test(y=friedman$respuesta, groups=friedman$area, blocks=friedman$time) 
Tukey<-HSD.test(aov(friedman$respuesta~friedman$time+friedman$area),"friedman$area", 
                group=TRUE,alpha=0.05); Tukey

# Time-level outliers 
friedman.test(y=friedman$respuesta, groups=friedman$time, blocks=friedman$area) 
Tukey<-HSD.test(aov(friedman$respuesta~friedman$time+friedman$area),"friedman$time", 
                group=TRUE,alpha=0.05); Tukey

# Optimal area-time specific robustness parameters
boxplot(TWMQ.results$k.dt, ylab=expression(hat(c)[phi][dt]), xlab='Time period')
abline(h=3, col='red')

boxplot(t(TWMQ.results$k.dt), ylab=expression(hat(c)[phi][dt]), xlab='Area')
abline(h=3, col='red')

