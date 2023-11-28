#setwd("~/Desktop/ranking selection/code/new code")
source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsorgen.R')
library(CVXR)

par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
d=read.csv("AYP_05.csv")
x=d$n_seapass/d$n_seatest-d$n_sedpass/d$n_sedtest
ranknumber=20

sea_var=d$n_seapass/d$n_seatest*(1-d$n_seapass/d$n_seatest)/d$n_seatest
sed_var=d$n_sedpass/d$n_sedtest*(1-d$n_sedpass/d$n_sedtest)/d$n_sedtest
se=sqrt(sea_var+sed_var)
par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

qunt=quantile(se,c(0.01,0.99))
x=x[which(se>qunt[1]&se<qunt[2])]
se=se[which(se>qunt[1]&se<qunt[2])]

s=se
par(mfrow=c(1, 4), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)

mu0=0.2
z=(x-mu0)/se

pv=1-pnorm(z, 0, 1)

d.dd=rih.func(x,se,mu0,0.1,gd=50,mod=1,jacknife = T)

bb=seq(0,max(s),length.out=10)

hist(s,breaks = bb,xlim = c(0,max(s)),ylab = 'count',xlab = 'se',main='(a) Distribution of se',col ='white')
hist(s[order(pv)[1:ranknumber]],breaks = bb,xlim = c(0,max(s)),ylim = c(0,12),main = ' (b) p-value' ,ylab = 'count',xlab = 'se',col = 'white')#pvalue
hist(s[which(rank(-x)<=ranknumber)],breaks=bb, xlim = c(0,max(s)),ylim = c(0,12),main = '(c) raw observation',ylab = 'count',xlab = 'se',col = 'white')
hist(s[which(rank(-d.dd$post)<=ranknumber)],breaks=bb, xlim = c(0,max(s)),ylim = c(0,12),main = '(d) posterior mean',ylab = 'count',xlab = 'se',col = 'white')




