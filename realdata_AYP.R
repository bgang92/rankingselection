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
d$total=d$n_seatest+d$n_sedtest
x=d$n_seapass/d$n_seatest-d$n_sedpass/d$n_sedtest
#d$total=d$n_seatest+d$n_seapass+d$n_sedtest+d$n_sedpass


sea_var=d$n_seapass/d$n_seatest*(1-d$n_seapass/d$n_seatest)/d$n_seatest
sed_var=d$n_sedpass/d$n_sedtest*(1-d$n_sedpass/d$n_sedtest)/d$n_sedtest
se=sqrt(sea_var+sed_var)
par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
d$x=x
d$se=se

if(1==1){
qunt=quantile(se,c(0.01,0.99))
x.rd=x[which(se>qunt[1]&se<qunt[2])]
se.rd=se[which(se>qunt[1]&se<qunt[2])]
d.rd=d[which(se>qunt[1]&se<qunt[2]),]
}


plot(x.rd,se.rd,ylab = 'se')
hist(x.rd,main = 'Distribution of x',col='white',breaks=12)
hist(se.rd,main = 'Distribution of se',xlab = 'se',col='white',breaks=12)
mu0=0.2
z=(x.rd-mu0)/se.rd

q=0.01

pv=1-pnorm(z, 0, 1)



  d.dd=rih.func(x.rd,se.rd,mu0,q,gd=50,mod=2,jacknife = T)
  d.hart=sc.func(d.dd$clfdr,q)
  d.bh=bh.func(pv,q)
  se.dd=se[which(d.dd$de==1)]
  se.hart=se[which(d.hart$de==1)]
  se.bh=se[which(d.bh$de==1)]

sum(x.rd[which(d.dd$de==1)]-0.2)

ind=which(d.dd$de==1&d.hart$de==0)
ind2=which(d.dd$de==0&d.hart$de==1)
ind3=which(d.dd$de==1&d.hart$de==1)


ptcolor=rep('gray',length(x.rd))
ptcolor[ind3]='green'

ptcolor[ind2]='blue'
ptcolor[ind]='red'

ptype=rep(1,length(x))

par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
ptype[which(ptcolor!='gray')]=3
plot(x.rd,se.rd,col='gray',pch=1,cex=2,main="(a) Rejection by DD and Clfdr ")
points(x.rd[c(ind3,ind2,ind)],se.rd[c(ind3,ind2,ind)],col=ptcolor[c(ind3,ind2,ind)],pch=18,cex=2)



rrval=rvalue2.func(x.rd,se.rd,q=0.01,gd=50,mod=2,jacknife = T)
rihres=rih.func(x.rd,se.rd,mu0=0.2)
posterior.x=rihres$post
length(posterior.x)


mm=20 #ranking the top 20


indp=order(pv)[1:mm]
raw_rank_r=order(rrval,decreasing = T)
rankr=data.frame(x=x.rd[raw_rank_r],s=se.rd[raw_rank_r])

#rank_postx=data.frame(x=x.rd[order(posterior.x,decreasing = T)],s=se.rd[order(posterior.x,decreasing = T)])


rankp=data.frame(x=x.rd[indp],s=se.rd[indp])

sorted_df <- sort_dataframe(rankr)

#sorted_df_postx=sort_dataframe(rank_postx)

# View the sorted data frame
print(sorted_df)
######## tyding up rank




indr_full=rep(0,length(x.rd))
for (i in 1:length(x.rd)) {
  print(which(x.rd==sorted_df$x[i] & se.rd==sorted_df$s[i]))
  indr_full[i]=which(x.rd==sorted_df$x[i] & se.rd==sorted_df$s[i])
}


order(pv)

d.rd$total[order(pv)[1:10]]
d.rd$total[indr_full[1:10]]
d.rd[order(indr_full)[1:10],]

ecdf_func <- ecdf(d.rd$total)
ecdf_func(mean(d.rd$total[order(pv)[1:20]]))
ecdf_func(mean(d.rd$total[indr_full[1:20]]))
## want to find the p-value rank in r-value
#top p in r
top_p_r=rep(0,length(pv))
for (i in 1:length(pv)) {
  top_p_r[i]=which(indr_full==order(pv)[i])
}

top_r_p=rep(0,length(pv))
for (i in 1:length(pv)) {
  top_r_p[i]=which(order(pv)==indr_full[i])
}


indr=indr_full[1:mm]
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


plot(x.rd,se.rd,col='gray',pch=1,cex=2,ylab='se',xlab = 'x',main='(b)Top 20 schools')
points(x.rd[indu],se.rd[indu],col=ptcolor,pch=18,cex=2)

y_up=5
k=20
common_breaks <- seq(0,max(se.rd),length.out = 10)
par(mfrow=c(1, 5), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
hist(se.rd,col='white',breaks=common_breaks,xlim=c(0,max(se.rd)),main='(a) Distribution of se',xlab='se',cex.lab=1.5,ylab='')
hist(se.rd[order(pv)][1:k],col='white',breaks=common_breaks,xlim=c(0,max(se.rd)),ylim=c(0,y_up),xlab='se',main='(b) p-value',cex.lab=1.5,ylab='')
hist(se.rd[order(x.rd,decreasing = T)][1:k],col='white',breaks=common_breaks,xlim=c(0,max(se.rd)),ylim=c(0,y_up),xlab='se',main='(c) raw observation',cex.lab=1.5,ylab='')
hist(se.rd[order(posterior.x,decreasing = T)][1:k],col='white',breaks=common_breaks,xlim=c(0,max(se.rd)),ylim=c(0,y_up),xlab='se',main='(d) posterior mean',cex.lab=1.5,ylab='')
hist(se.rd[indr_full][1:k],col='white',breaks=common_breaks,xlim=c(0,max(se.rd)),ylim=c(0,y_up),xlab='se',main='(e) r-value',cex.lab=1.5,ylab='')


par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
common_breaks <- seq(0,3500,length.out = 14)
common_breaks1<- seq(0,3500,length.out = 14)
hist(d.rd$total,col='white',breaks=common_breaks1,xlim=c(0,3500),main='(a) Distribution of size',xlab='size',cex.lab=1.5,ylab='')
hist(d.rd$total[order(pv)][1:k],col='white',breaks=common_breaks,xlim=c(0,max(common_breaks)),ylim=c(0,y_up),xlab='size',main='(b) p-value',cex.lab=1.5,ylab='')
#hist(d.rd$total[order(x.rd,decreasing = T)][1:k],col='white',breaks=common_breaks,xlim=c(0,max(common_breaks)),ylim=c(0,y_up),xlab='size',main='(c) raw observation',cex.lab=1.5,ylab='')
#hist(d.rd$total[order(posterior.x,decreasing = T)][1:k],col='white',breaks=common_breaks,xlim=c(0,max(common_breaks)),ylim=c(0,y_up),xlab='size',main='(d) posterior mean',cex.lab=1.5,ylab='')
hist(d.rd$total[indr_full][1:k],col='white',breaks=common_breaks,xlim=c(0,max(common_breaks)),ylim=c(0,y_up),xlab='size',main='(c) r-value',cex.lab=1.5,ylab='')


k=10
d.rd[order(pv)[1:k],]

d.rd[order(x.rd,decreasing = T)[1:k],]
d.rd[order(posterior.x,decreasing = T)[1:k],]
d.rd[indr_full[1:k],]

