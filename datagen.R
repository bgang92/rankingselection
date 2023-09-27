rdata.rs.func<-function(m,x1mean,sigma,mu0,reverse=F)
{
  theta=rep(0,m)
  theta[which(x1mean>mu0)]=1
  if(reverse==T){
    theta=1-theta
  }
  
  x=rnorm(m,x1mean,sigma)
  
  y<-list(x=x, theta=theta)
  return(y)
}

mudata.sun<-function(m,p,mu){
  theta=rbinom(m,1,p)
  res=(1-theta)*rep(0,m)+theta*rep(mu,m)
}


rdata.mult.func<-function(m,x1mean,sigma,mu0,nn)
{
  theta=rep(0,m)
  theta[which(x1mean>mu0)]=1
  x=rep(NA,m)
  sigest=rep(NA,m)
  for (i in 1:m) {
    da=rnorm(nn,x1mean[i],sigma[i]*sqrt(nn))
    x[i]=mean(da)
    sigest[i]=sd(da)/sqrt(nn)
    
  }
 
  y<-list(x=x, theta=theta,sigest=sigest)
  return(y)
}


rdata.mu.func<-function(a1,la1,pa1,lo,up,m)#unif
{
  theta=rbinom(la1,1,pa1)
  mu1=(1-theta)*runif(la1,lo,up)+theta*a1
  mu2=runif(m-la1,lo,up)
  mu=c(mu1,mu2)
  return(mu)
}

rdata.mu.func2<-function(a1,la1,pa1,m.vec,s.vec,p.vec,lobound,upbound,m)
{
  theta1=rbinom(la1,1,pa1)
  theta2=rmultinom(m,1,p.vec)
  ind=rep(NA,m)
  for (i in 1:m) {
    ind[i]=which(theta2[,i]==1)
  }
  mualt=rep(NA,m)
  for (i in 1:m) {
    mualt[i]=rtruncnorm(1,a=lobound,b=upbound,m.vec[ind[i]],s.vec[ind[i]])
  }
  theta2=rbinom(la1,1,pa1)
  mu1=(1-theta1)*mualt[1:la1]+theta1*a1
  mu=c(mu1,mualt[(la1+1):m])
  return(mu)
}

rdata.mu.func3<-function(a1,la1,pa1,m.vec,p.vec)
{
  theta1=rbinom(la1,1,pa1)
  theta2=rmultinom(m,1,p.vec)
  ind=rep(NA,m)
  for (i in 1:m) {
    ind[i]=which(theta2[,i]==1)
  }
  mualt=rep(NA,m)
  for (i in 1:m) {
    mualt[i]=m.vec[ind[i]]
  }
  theta2=rbinom(la1,1,pa1)
  mu1=(1-theta1)*mualt[1:la1]+theta1*a1
  mu=c(mu1,mualt[(la1+1):m])
  return(mu)
}

mudata=function(m,lo,mu0,mualt,prob){
  theta=rbinom(m,1,prob)
  munif=runif(m,lo,mu0)
  x1mean=munif*(1-theta)+mualt*theta
  return(x1mean)
}

mudata2=function(m,lo,mu0,altlo,altup,prob){
  theta=rbinom(m,1,prob)
  munif=runif(m,lo,mu0)
  mualt=runif(m,altlo,altup)
  x1mean=munif*(1-theta)+mualt*theta
  return(x1mean)
}
mudata3=function(m,lo,up,altlo,altup,prob){
  theta=rbinom(m,1,prob)
  munif=runif(m,lo,up)
  mualt=runif(m,altlo,altup)
  x1mean=munif*(1-theta)+mualt*theta
  return(x1mean)
}


mudata.norm=function(m,nullmean,nullsd,mualt,prob){
  theta=rbinom(m,1,prob)
  munif=rnorm(m,nullmean,nullsd)
  x1mean=munif*(1-theta)+mualt*theta
  return(x1mean)
}

mudata.norm2=function(m,nullmean,nullsd,altmu,altsd,prob){
  theta=rbinom(m,1,prob)
  munif=rnorm(m,nullmean,nullsd)
  mualt=rnorm(m,altmu,altsd)
  x1mean=munif*(1-theta)+mualt*theta
  return(x1mean)
}


eg.data=function(m,a1,b1,a2,b2,prob,sig){
  theta=rbinom(m,1,prob)
  mu.null=runif(m,a1,b1)
  mu.alt=runif(m,a2,b2)
  x1mean=mu.null*(1-theta)+mu.alt*theta
  x1=x1mean+rnorm(m,0,sig)
  res=list(x=x1,theta=theta,mu=x1mean)
  return(res)
  
}

rdata.unif.func4<-function(up,lo,mu0,sig,m)
{
  mu=runif(m,lo,up)
  theta=rep(0,m)
  theta[which(mu>mu0)]=1
  x=mu+rnorm(m,0,sig)
  res=list(x=x,theta=theta,mu=mu)
  return(res)
}

rdata.unif.func5<-function(up,lo,mu0,sig,m,copies)
{
  mu=runif(m,lo,up)
  theta=rep(0,m)
  theta[which(mu>mu0)]=1
  x.mat=matrix(rep(NA,m*copies),nrow = copies,ncol = m)
  for (i in 1:copies) {
    x.mat[i,]=mu+rnorm(m,0,sig)
  }
  x=colMeans(x.mat)
  std=apply(x.mat, 2, sd)
  res=list(x=x,std=std,theta=theta,mu=mu)
  return(res)
}

rdata.unif.func4corr<-function(up,lo,mu0,sig,m,rho,bksize)
{
  mu=runif(m,lo,up)
  theta=rep(0,m)
  theta[which(mu>mu0)]=1
  bknum=m/bksize
  noisemat=matrix(rep(NA,m),nrow = bksize)
  sigmat=matrix(sig,nrow = bksize)
  for (j in 1:bknum) {
    covmat.tent=matrix(rep(NA,bksize*bksize),nrow = bksize)
    covmat=matrix(rep(NA,bksize*bksize),nrow = bksize)
    for (i in 1:bksize) {
      covmat.tent[i,]=rho^(abs(i-seq(1,bksize,1)))
    }
    for (i in 1:bksize) {
      covmat[i,]=covmat.tent[i,]*sigmat[,j]*sigmat[i,j]
    }
  nn=mvrnorm(1,mu=rep(0,bksize),covmat)
  noisemat[,j]=nn
  }
  noise=as.vector(noisemat)


  x=mu+noise
  res=list(x=x,theta=theta,mu=mu)
  return(res)
}


