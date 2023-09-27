#setwd("~/Desktop/ranking selection/code/new code")
source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsorgen.R')
library(CVXR)
library(invgamma)
library(rvalues)
library(truncnorm)
par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
d=read.csv("AYP_05.csv")
x=d$n_seapass/d$n_seatest-d$n_sedpass/d$n_sedtest
d$total=d$n_seatest+d$n_seapass+d$n_sedtest+d$n_sedpass


sea_var=d$n_seapass/d$n_seatest*(1-d$n_seapass/d$n_seatest)/d$n_seatest
sed_var=d$n_sedpass/d$n_sedtest*(1-d$n_sedpass/d$n_sedtest)/d$n_sedtest
se=sqrt(sea_var+sed_var)
par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

if(1==1){
qunt=quantile(se,c(0.01,0.99))
x=x[which(se>qunt[1]&se<qunt[2])]
se=se[which(se>qunt[1]&se<qunt[2])]
}

s=se
plot(x,s,ylab = 'se')
hist(x,main = 'Distribution of x')
hist(s,main = 'Distribution of se',xlab = 'se')
mu0=0.2
z=(x-mu0)/se

q=0.01

pv=1-pnorm(z, 0, 1)



  d.dd=rih.func(x,se,mu0,q,gd=50,mod=2,jacknife = T)
  d.hart=sc.func(d.dd$clfdr,q)
  d.bh=bh.func(pv,q)
  se.dd=se[which(d.dd$de==1)]
  se.hart=se[which(d.hart$de==1)]
  se.bh=se[which(d.bh$de==1)]



ind=which(d.dd$de==1&d.hart$de==0)
ind2=which(d.dd$de==0&d.hart$de==1)
ind3=which(d.dd$de==1&d.hart$de==1)


ptcolor=rep('gray',length(x))
ptcolor[ind3]='green'

ptcolor[ind2]='blue'
ptcolor[ind]='red'

ptype=rep(1,length(x))

par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
ptype[which(ptcolor!='gray')]=3
plot(x,se,col='gray',pch=1,cex=2,main="(a) Rejection by DD and Clfdr ")
points(x[c(ind3,ind2,ind)],se[c(ind3,ind2,ind)],col=ptcolor[c(ind3,ind2,ind)],pch=18,cex=2)



rrval=rvalue2.func(x,se,q=0.01,gd=50,mod=2,jacknife = T)




mm=20 #ranking the top 20



indp=order(pv)[1:mm]
indr=order(rrval,decreasing = T)[1:mm]
indboth=intersect(indp,indr)
indu=union(indp,indr)

ptcolor=rep('gray',length(indu))
for (i in 1:length(indu)) {
  if(indu[i] %in% indr){
    ptcolor[i]='red'
  }
  if(indu[i] %in% indp){
    ptcolor[i]='blue'
  }
  if((indu[i] %in% indp)&& (indu[i] %in% indr) ){
    ptcolor[i]='green'
  }
}

ptype=rep(1,length(indu))


plot(x,s,col='gray',pch=1,cex=2,ylab='se',xlab = 'x',main='(b)Top 20 schools')
points(x[indu],s[indu],col=ptcolor,pch=18,cex=2)




