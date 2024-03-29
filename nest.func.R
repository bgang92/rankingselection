falt.nest=function(z,x,sig,weight=rep(1,length(x)),jacknife=T){

  #This function compute the weighted bivariate density estimate at z (see section 3.1 of the paper for detail)
  #x is the vector of raw observations.
  #sig is the vector of corresponding standard deviations.
  #z is a 2-tuple,the first is the observation, the second is standard deviation
 
  xds=density(x,from = min(x)-IQR(x)/2,to=max(x)+IQR(x)/2,n=length(x)/5)   #this line and the next are used to find bandwidth
  sigds=density(sig,from = min(sig)-IQR(sig)/2,to=max(sig)+IQR(sig)/2,n=length(x)/5)
  h=xds$bw
  hsig=sigds$bw
  if(jacknife==T){
    if(is.element(z[1],x)){
      z.ind=which(x==z[1])
      x=x[-z.ind[1]]
      sig=sig[-z.ind[1]]
      weight=weight[-z.ind[1]]/sum(weight[-z.ind[1]])
    }
  }
  obs=rep(z[1],length(x))
  obssig=rep(z[2],length(x))
  k=dnorm((obs-x)/h,0,sig)
  
  nestw=dnorm(obssig-sig,0,hsig)
  w=weight*nestw
  
  y=sum(w*k)/(sum(w)*h)
  return(y)
  
}






nest.func=function(x,sig,jacknife=T){
  #This function calculate the weighted bivariate density at every point of (x,sig)
  fm=function(a,b){
    res=falt.nest(c(a,b),x,sig,rep(1/length(x),length(x)),jacknife=jacknife)
    return(res)
  }

  denom=mapply(fm, x,sig)

  result=list(dens=denom)
  return(result)
}







dens.func=function(w,gd,sig,x){
  #w is the weigt
  #gd is the grid
  #sig is the standard deviation (vector)
  #compute the density of convolution with N(0,sig) at x (vector)
  est=rep(0,length(x))
  for (i in 1:length(x)) {
    a=dnorm(gd-x[i],0,sig[i])
    est[i]=sum(a*w)
  }
 
  return(est)
}

west.func=function(dens.dd,x,sig,gdsize=500,mu=NULL,adjust=T){
  #This function solves the optimization problem (3.2) in the paper
  
  if(is.null(mu)){
  
  gd=seq(quantile(x,0.01),quantile(x,0.99),length.out = gdsize)
  }
  else{
    gd=seq(min(mu),max(mu),length.out = gdsize)
  }
 
  
  w_hat=Variable(gdsize)
  
  f1=function(a,b){
    res=dnorm(gd-a,0,b)
    return(res)
  }
  A=mapply(f1, x,sig)
  
  objective<-Minimize(sum((dens.dd-t(A)%*%w_hat)^2))
  cons1<- (sum(w_hat)==1)
  cons2<- w_hat>=0
  problem <- Problem(objective,constraints = list(cons1, cons2))
  result <- solve(problem,solver='SCS')
  w=try(result$getValue(w_hat))
  if('try-error'%in% class(w)){
    result <- solve(problem,solver='ECOS')
    w=try(result$getValue(w_hat))
  }
  if('try-error'%in% class(w)){
    result <- solve(problem,solver='OSQP')
    w=try(result$getValue(w_hat))
  }
  if('try-error'%in% class(w)){
    result <- solve(problem,solver='GLPK')
    w=try(result$getValue(w_hat))
  }
  if('try-error'%in% class(w)){
    result <- solve(problem,solver='GLPK_MI')
    w=try(result$getValue(w_hat))
  }
  
  if(adjust==T){
    w=pmax(w,0)
    w=w/sum(w)
  }
  
  res=list(w=w,gd=gd)
  return(res)
}
 




