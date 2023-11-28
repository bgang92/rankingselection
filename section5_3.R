source("3M.funcs.R")
source("datagen.R")
source("nest.func.R")
source('rih.R')
source('rsorgen.R')
library(CVXR)
#library(invgamma)
#library(rvalues)
#library(truncnorm)



m<-10000   #number of hypothesis
nrep<-1#number of reptitions
q<-0.1 
mu0=1 #we are interested in mu>mu_0
nullmean=-0.5
nullsd=0.25
epsilon=10^(-5)


altmu1=1.5
altmu2=3
altsd=0.25
sig.vec=seq(1.5,2,by=0.1)
np<-length(sig.vec)


dd.fdr<-rep(0, np)
rs.fdr<-rep(0, np)
bh.fdr=rep(0,np)
rs2.fdr=rep(0,np)
z.fdr=rep(0,np)

dd.etp<-rep(0, np)
rs.etp<-rep(0, np)
rs2.etp<-rep(0, np)
bh.etp<-rep(0, np)
z.etp=rep(0,np)
ddave.etp=rep(0,np)
rs2ave.etp=rep(0,np)
orave.etp=rep(0,np)

diff1.etp=rep(0,np)
diff2.etp=rep(0,np)
diff3.etp=rep(0,np)
diff4.etp=rep(0,np)


dd.etp.old<-rep(0, np)
rs.etp.old<-rep(0, np)
rs2.etp.old<-rep(0, np)
bh.etp.old<-rep(0, np)
z.etp.old=rep(0,np)

diff1.etp.old<-rep(0, np)
diff2.etp.old<-rep(0, np)
diff3.etp.old<-rep(0, np)
diff4.etp.old<-rep(0, np)

dd.fdp<-matrix(rep(0, nrep*np), np, nrep)
rs.fdp<-matrix(rep(0, nrep*np), np, nrep)
rs2.fdp<-matrix(rep(0, nrep*np), np, nrep)
bh.fdp<-matrix(rep(0, nrep*np), np, nrep)
z.fdp<-matrix(rep(0, nrep*np), np, nrep)


dd.ntp<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp<-matrix(rep(0, nrep*np), np, nrep)
rs2.ntp<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp<-matrix(rep(0, nrep*np), np, nrep)
z.ntp<-matrix(rep(0, nrep*np), np, nrep)
ddave.ntp<-matrix(rep(0, nrep*np), np, nrep)
rs2ave.ntp<-matrix(rep(0, nrep*np), np, nrep)
orave.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff1.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff2.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff3.ntp<-matrix(rep(0, nrep*np), np, nrep)
diff4.ntp<-matrix(rep(0, nrep*np), np, nrep)

dd.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
rs.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
rs2.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
bh.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
z.ntp.old<-matrix(rep(0, nrep*np), np, nrep)

diff1.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff2.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff3.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
diff4.ntp.old<-matrix(rep(0, nrep*np), np, nrep)
for (i in 1:np){
  
 
  sigmult=c(0.25,1.25)*sig.vec[i]
  p=c(1/10,1/10,1/10)
  
  
 
  
 
  sig=c(rep(sigmult[1],m/2),rep(sigmult[2],m/2))
  
  
  

    
    fnull<- function(z,s) integrate(function(x,z) dnorm(x,nullmean,nullsd)*dnorm(z,x,s),-Inf,mu0,z)$value
    fnull1=function(z,s) integrate(function(x,z) dnorm(x,altmu1,altsd)*dnorm(z,x,s),-Inf,mu0,z)$value
    fnull2=function(z,s) integrate(function(x,z) dnorm(x,altmu2,altsd)*dnorm(z,x,s),-Inf,mu0,z)$value
   
     falt0<- function(z,s) integrate(function(x,z) dnorm(x,nullmean,nullsd)*dnorm(z,x,s),mu0,Inf,z)$value
    falt1<- function(z,s) integrate(function(x,z) dnorm(x,altmu1,altsd)*dnorm(z,x,s),mu0,Inf,z)$value
    falt2<- function(z,s) integrate(function(x,z) dnorm(x,altmu2,altsd)*dnorm(z,x,s),mu0,Inf,z)$value
  
  

  #sig=rep(1,m)
  for (j in 1:nrep) {
    oldtime=Sys.time()
    cat("j=",j,"i=",i)

    
    mu1=mudata.norm2(m/2,nullmean,nullsd,altmu1,altsd,p[1])
    mu2=mudata.norm2(m/2,nullmean,nullsd,altmu2,altsd,p[2])
    mu=c(mu1,mu2)
    
    d=rdata.rs.func(m,mu,sig,mu0)
    x=d$x
    z=(x-mu0)/sig
    theta=d$theta
    pv=1-pnorm(z, 0, 1)
    

  
    num0=(1-p[1])*mapply(fnull,x,sig)
    num1=(p[1])*mapply(fnull1,x[1:(m/2)],sig[1:(m/2)])
    num2=(p[1])*mapply(fnull2,x[(m/2+1):m],sig[(m/2+1):m])
    
    altpart0=(1-p[1])*mapply(falt0,x,sig)
 
    altpart1=p[1]*mapply(falt1,x[1:(m/2)],sig[1:(m/2)])
    altpart2=p[2]*mapply(falt2,x[(m/2+1):m],sig[(m/2+1):m])

    num=num0+c(num1,num2)
    denom=c(altpart1,altpart2)
    clfdr=num/(num+denom+altpart0)
    clfdr=pmin(clfdr,1)
    
    
  

     
      tst.dd1=rsstathomo.func(x[1:(m/2)],sig[1:(m/2)],mu0,gd=50)
      tst.dd2=rsstathomo.func(x[(m/2+1):m],sig[(m/2+1):m],mu0,gd=50)
   
     plot(tst.dd1$gd,tst.dd1$w)
     plot(tst.dd2$gd,tst.dd2$w)
       clfdr.dd=c(tst.dd1$clfdr,tst.dd2$clfdr)
     
      d.dd=rsorddgen.func(x-mu0,clfdr.dd,q)
      dd.fdp[i, j]<-sum((1-theta)*d.dd$de)/max(sum(d.dd$de), 1)
      dd.ntp[i, j]<-sum((x-mu0)*d.dd$de)
      dd.ntp.old[i,j]=sum((theta)*d.dd$de)/max(sum(theta), 1)
      ddave.ntp[i,j]=dd.ntp[i,j]/max(sum(d.dd$de), 1)
 

    drs=rsorddgen.func(x-mu0,clfdr,q)

    rs.fdp[i, j]<-sum((1-theta)*drs$de)/max(sum(drs$de), 1)
    rs.ntp[i, j]<-sum((x-mu0)*drs$de)
    rs.ntp.old[i,j]=sum((theta)*drs$de)/max(sum(theta), 1)
    orave.ntp[i,j]=rs.ntp[i,j]/max(sum(drs$de), 1)

    drs2=sc.func(clfdr,q)
    rs2.fdp[i, j]<-sum((1-theta)*drs2$de)/max(sum(drs2$de), 1)
    rs2.ntp[i, j]<-sum((x-mu0)*drs2$de)
    rs2.ntp.old[i,j]=sum((theta)*drs2$de)/max(sum(theta), 1)
    rs2ave.ntp[i,j]=rs2.ntp[i,j]/max(sum(drs2$de), 1)
  
    bh.res<-bh.func(pv, q)
    bh.de<-bh.res$de
    bh.fdp[i, j]<-sum((1-theta)*bh.de)/max(sum(bh.de), 1)
    bh.ntp[i, j]<-sum((x-mu0)*bh.de)
    bh.ntp.old[i,j]=sum((theta)*bh.de)/max(sum(theta), 1)
    
 

    ind=which(d.dd$de==1&drs2$de==0)
    ind2=which(d.dd$de==0&drs2$de==1)
    ind3=which(drs$de==1&drs2$de==0)
    ind4=which(drs$de==0&drs2$de==1)
    
    diff1.ntp[i,j]=sum(x[ind]-mu0)
    diff2.ntp[i,j]=sum(x[ind2]-mu0)
    diff3.ntp[i,j]=sum(x[ind3]-mu0)
    diff4.ntp[i,j]=sum(x[ind4]-mu0)
    
    diff1.ntp.old[i,j]=length(ind)
    diff2.ntp.old[i,j]=length(ind2)
    diff3.ntp.old[i,j]=length(ind3)
    diff4.ntp.old[i,j]=length(ind4)
    
    
    newtime=Sys.time()
    print(newtime-oldtime)
  }#end of j loop
  dd.fdr[i]<-mean(dd.fdp[i,])
  rs.fdr[i]<-mean(rs.fdp[i,])
  rs2.fdr[i]<-mean(rs2.fdp[i,])
  bh.fdr[i]<-mean(bh.fdp[i,])
  
  
  dd.etp[i]<-mean(dd.ntp[i,])
  rs.etp[i]<-mean(rs.ntp[i,])
  rs2.etp[i]<-mean(rs2.ntp[i,])
  bh.etp[i]<-mean(bh.ntp[i,])
  
  
  dd.etp.old[i]<-mean(dd.ntp.old[i,])
  rs.etp.old[i]<-mean(rs.ntp.old[i,])
  rs2.etp.old[i]<-mean(rs2.ntp.old[i,])
  bh.etp.old[i]<-mean(bh.ntp.old[i,])
  
  ddave.etp[i]=mean(ddave.ntp[i,])
  orave.etp[i]=mean(orave.ntp[i,])
  rs2ave.etp[i]=mean(rs2ave.ntp[i,])
  
  diff1.etp[i]=mean(diff1.ntp[i,])
  diff2.etp[i]=mean(diff2.ntp[i,])
  diff3.etp[i]=mean(diff3.ntp[i,])
  diff4.etp[i]=mean(diff4.ntp[i,])
  
  diff1.etp.old[i]=mean(diff1.ntp.old[i,])
  diff2.etp.old[i]=mean(diff2.ntp.old[i,])
  diff3.etp.old[i]=mean(diff3.ntp.old[i,])
  diff4.etp.old[i]=mean(diff4.ntp.old[i,])
  
}#end of i loop

fdr1<-cbind(dd.fdr, rs.fdr,bh.fdr,rs2.fdr)
etp1<-cbind(dd.etp, rs.etp,bh.etp,rs2.etp)
etp.old<-cbind(dd.etp.old,rs.etp.old,bh.etp.old,rs2.etp.old)
fdr1
etp1
etp.old
ave.etp=cbind(ddave.etp,rs2ave.etp,orave.etp)

par(mfrow=c(1, 3), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
if(T){ 
  matplot(sig.vec, fdr1, type="o", pch=1:5, lwd=2, main=expression(paste("a)FDR Comparison varing ", sigma)), xlab=expression(sigma), ylim=c(-0.01, 0.20), ylab="FDR")
 
  legend("top", c("DD", "OR", "BH","Clfdr"), pch=1:4, col=1:4, lwd=2)
  matplot(sig.vec, etp1, type="o", pch=1:5, lwd=2, main=expression(paste("b)" ,"ETP* Comparison varing ", sigma)), xlab=expression(sigma), ylab="ETP*",ylim = c(700,1200))
 
  matplot(sig.vec, etp.old*988.6249, type="o", pch=1:5, lwd=2, main=expression(paste("c)ETP Comparison varing ", sigma)), xlab=expression(sigma), ylab="ETP",ylim = c(0.3,0.7)*988.6249)
  
}


if(T){
  ind=which(d.dd$de==1&drs2$de==0)
  ind2=which(d.dd$de==0&drs2$de==1)
  ind3=which(d.dd$de==1&drs2$de==1)
  

  
  length(ind)
  length(ind2)
  sum(d.dd$de)
  
  ptcolor=rep('gray',length(x))
  ptcolor[ind]='red'
  ptcolor[ind2]='blue'
  ptcolor[ind3]='green'
  par(mfrow=c(1, 2), mgp=c(2, 0.5, 0), mar=c(3, 3, 2, 1)+0.1)
  ptype=rep(1,length(x))
  
  ptype[which(ptcolor!='gray')]=3
  plot(x,sig,col='gray',pch=1,cex=2,ylab=expression(sigma),main=expression(paste("a)" ,"Rejection by DD and Clfdr ", sigma,"=2")))
  points(x[c(ind3,ind,ind2)],sig[c(ind3,ind,ind2)],col=ptcolor[c(ind3,ind,ind2)],pch=18,cex=2)
  
 
  etp.diff=cbind(diff1.etp,diff2.etp)
  matplot(sig.vec, etp.diff, type="o", pch=1:5,col=c('red','blue'),  lwd=2, main=expression(paste("b)" ,"ETP* Comparison varing ", sigma)), xlab=expression(sigma), ylab="ETP*",ylim = c(0,200))
  legend("topleft", c("DD", "Clfdr"), pch=1:4, col=c('red','blue'), lwd=2)
}
