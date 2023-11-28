#setwd("~/Desktop/ranking selection/code/new code")
source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsor.R')
source('rsorgen.R')
source('cb.colors.R')
library(CVXR)
#library(invgamma)
#library(rvalues)
#library(truncnorm)


m<-10000   
nrep<-100
q=0.1

epsilon=10^(-5)

rs.fdp<-rep(0, nrep)
rs2.fdp=rep(0,nrep)


rs.etp<-rep(0, nrep)
rs2.etp<-rep(0, nrep)

rs.etp.new<-rep(0, nrep)
rs2.etp.new<-rep(0, nrep)

rs.th=rep(NA,nrep)
clfdr.th=rep(NA,nrep)


for (i in 1:nrep){
 
  a1=-3
  b1=-1
  a2=1
  b2=2
  prob=0.2
  sig=runif(m,0.5,3)
  d=eg.data(m,a1,b1,a2,b2,prob,sig)
  x=d$x
  theta=d$theta
  fnull<- function(z,s) integrate(function(x,z) dunif(x,a1,b1)*dnorm(z,x,s),a1,b1,z)$value
  falt<- function(z,s) integrate(function(x,z) dunif(x,a2,b2)*dnorm(z,x,s),a2,b2,z)$value
  
  num=mapply(fnull,x,sig)*(1-prob)
  denom=mapply(falt,x,sig)*prob+num
  clfdr=num/denom
  
  drs=rsorddgen.func(x-0,clfdr,q)
  drs2=sc.func(clfdr,q)
  rs.fdp[i]=sum((1-theta)*drs$de)/max(sum(drs$de), 1)
  rs2.fdp[i]=sum((1-theta)*drs2$de)/max(sum(drs2$de), 1) 
  rs.etp[i]=sum((x-0)*drs$de)
  rs2.etp[i]=sum((x-0)*drs2$de)
  rs.etp.new[i]=sum((theta)*drs$de)/max(sum(theta), 1)/m
  rs2.etp.new[i]=sum((theta)*drs2$de)/max(sum(theta), 1)/m
  rs.th[i]=drs$th
  clfdr.th[i]=drs2$th
  
  
  
}#end of i loop




findt=function(clfdr,s){
  prob=0.2
  q=0.1
  fnull<- function(z,s) integrate(function(x,z) dunif(x,a1,b1)*dnorm(z,x,s),a1,b1,z)$value
  falt<- function(z,s) integrate(function(x,z) dunif(x,a2,b2)*dnorm(z,x,s),a2,b2,z)$value
  potentialx=seq(-5,5,by=0.01)
  num=mapply(fnull,potentialx,s)*(1-prob)
  denom=mapply(falt,potentialx,s)*prob+num
  potential.clfdr=num/denom
  xx=potentialx[which.min(abs(potential.clfdr-clfdr))]
  tt=xx/clfdr-q
  return(tt)
}



findclfdr=function(tt,s){
  prob=0.2
  q=0.1
  fnull<- function(z,s) integrate(function(x,z) dunif(x,a1,b1)*dnorm(z,x,s),a1,b1,z)$value
  falt<- function(z,s) integrate(function(x,z) dunif(x,a2,b2)*dnorm(z,x,s),a2,b2,z)$value
  potentialx=seq(-5,5,by=0.01)
  num=mapply(fnull,potentialx,s)*(1-prob)
  denom=mapply(falt,potentialx,s)*prob+num
  potential.clfdr=num/denom
  potential.t=potentialx/potential.clfdr-q
  clfdr=potential.clfdr[which.min(abs(potential.t-tt))]
  return(clfdr)
}




findclfdrx=function(tt,x){
  clfdr=x/tt+q
 
  return(clfdr)
}



plot.regions=function(clfdr.cut, rs.cut,N=500,lower,upper,prob=.2,mu=3){
  source("cb.colors.R")
  q=0.1
  par(mfrow=c(1,2))
  par(mar=c(4,4,3.5,1))
  require(latex2exp)

  xgrid=seq(2,6,length=N)
 
  clfdr.grid=seq(0.15,0.7,length=N)
 
  clfdr.mat=matrix(rep(NA,N*N),nrow = N)
  t.mat=matrix(rep(NA,N*N),nrow = N)
  for (i in 1:N) {
    clfdr.mat[i,]=clfdr.grid[i]
   
    t.mat[i,]=xgrid/(clfdr.grid[i]-q)
  }
 
  clfdr.t=rep(clfdr.cut,N)
 
  image(xgrid,clfdr.grid,t(clfdr.mat),col=cm.colors(n=10,alpha=.7)[seq(0,10,1)],ylab="Clfdr",xlab=TeX("$x$"),cex.axis=.7,cex.lab=.8)
  lines(xgrid,clfdr.t,lwd=1.5)
 
   t2.mat=log(t.mat)
 
  t2.t=rep(0,N)
  for (i in 1:N){
   
    t2.t[i]=max(0.1,findclfdrx(rs.cut,xgrid[i]))
  }
  
 
  image(xgrid,clfdr.grid,t(t2.mat),col=cm.colors(n=10,alpha=1)[seq(10,0,-1)],ylab="Clfdr",xlab=TeX("$x$"),cex.axis=.7,cex.lab=.8)
  lines(xgrid,t2.t,lwd=1.5)

}


plot.regions(mean(clfdr.th),mean(rs.th),lower = 0.5, upper = 3)

