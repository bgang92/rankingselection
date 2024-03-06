#setwd("~/Desktop/ranking selection/code/new code")
source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsorgen.R')
library(CVXR)

#par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
d=read.csv("mean_se_5yrs.csv")
x=d$alpha_firsthalf
se=d$se_alpha_firsthalf


  qunt=quantile(se,c(0.001,0.999))
  x=x[which(se>qunt[1]&se<qunt[2])]
  se=se[which(se>qunt[1]&se<qunt[2])]


s=se

mu0=0
z=(x-mu0)/se
pv=1-pnorm(z, 0, 1)
sig=s

par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
#rrval=rvalue2.func(x,se,q=0.1,gd=50,mod=2,jacknife = T)




rrval1=rvalue.func(x,sig,mu0=0,gd=50,mod=2,jacknife = T)

rrval2=rvalue2.func(x,sig,q=0.1,gd=50,mod=2,jacknife = T)

mm=20
indr1=order(rrval1,decreasing = F)[1:mm]
indr2=order(rrval2,decreasing = T)[1:mm]


indboth=intersect(indr1,indr2)
indu=union(indr2,indr1)



ptcolor=rep('gray',length(indu))
for (i in 1:length(indu)) {
  if(indu[i] %in% indr2){
    ptcolor[i]='red'
  }
  if(indu[i] %in% indr1){
    ptcolor[i]='blue'
  }
  if((indu[i] %in% indr1)&& (indu[i] %in% indr2) ){
    ptcolor[i]='green'
  }
}

ptype=rep(1,length(indu))
plot(x,s,col='gray',pch=1,cex=2,ylab='se',xlab = 'x',main='(a) Top 20 mutual funds')
points(x[indu],s[indu],col=ptcolor,pch=18,cex=2)





d=read.csv("AYP_05.csv")
x=d$n_seapass/d$n_seatest-d$n_sedpass/d$n_sedtest
d$total=d$n_seatest+d$n_seapass+d$n_sedtest+d$n_sedpass


sea_var=d$n_seapass/d$n_seatest*(1-d$n_seapass/d$n_seatest)/d$n_seatest
sed_var=d$n_sedpass/d$n_sedtest*(1-d$n_sedpass/d$n_sedtest)/d$n_sedtest
se=sqrt(sea_var+sed_var)


if(T){
  qunt=quantile(se,c(0.01,0.99))
  x=x[which(se>qunt[1]&se<qunt[2])]
  se=se[which(se>qunt[1]&se<qunt[2])]
}


s=se
mu0=0.2
tot=d$total
z=(x-mu0)/se

#par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
#rrval=rvalue2.func(x,se,q=0.1,gd=50,mod=2,jacknife = T)

sig=s


rrval1=rvalue.func(x,sig,mu0=mu0,gd=50,mod=2,jacknife = T)

rrval2=rvalue2.func(x,sig,q=0.01,gd=50,mod=2,jacknife = T)

mm=20
indr1=order(rrval1,decreasing = F)[1:mm]
indr2=order(rrval2,decreasing = T)[1:mm]


indboth=intersect(indr1,indr2)
indu=union(indr2,indr1)



ptcolor=rep('gray',length(indu))
for (i in 1:length(indu)) {
  if(indu[i] %in% indr2){
    ptcolor[i]='red'
  }
  if(indu[i] %in% indr1){
    ptcolor[i]='blue'
  }
  if((indu[i] %in% indr1)&& (indu[i] %in% indr2) ){
    ptcolor[i]='green'
  }
}


ptype=rep(1,length(indu))
plot(x,s,col='gray',pch=1,cex=2,ylab='se',xlab = 'x',main='(b) Top 20 schools')
points(x[indu],s[indu],col=ptcolor,pch=18,cex=2)





