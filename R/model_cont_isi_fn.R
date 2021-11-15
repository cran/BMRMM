# Supporting functions for model_cont_isi.R

############## likelihood function ##############
isi_logml <- function(z,y,M,pM,alpha,zmax) {   # z:z_isi (only for important predictors), y:z_tau
  df <- cbind(y,z)
  for(pp in 1:size(df,2))
    df <- sortrows(df,size(df,2)-pp+1)
  z0 <- unique(df)                   # z0 are the sorted unique combinations of (y,z_{pp,1},...,z_{pp,p0}), 
  m <- which(duplicated(df)==FALSE)  # m contains the positions (indices) when they first appear on {1,...,n}
  m <- c(diff(m),size(df,1)-m[length(m)]+1)
  C <- array(0,dim=c(zmax,M))        # d0=levels of the response y, M=number of clustered levels of x_{pp,1},...,x_{pp,p0}
  C[z0] <- C[z0]+m                   # add the differences in positions to cells of clT corresponding to the unique combinations -> gives the number of times (y,z_{pp,1},...,z_{pp,p0}) appears 
  Cdata <- matrix(C,nrow=zmax)       # matrix representation of the array C, with rows of the matrix corresponding to dimension 1 i.e. the levels of y, # col: # of unique combinations of (y,z_{pp,1},...,z_{pp,p0})
  
  JJ <- size(Cdata)
  if(length(alpha)==1)
    loglik <- sum(sum(lgamma(Cdata+alpha)))-sum(lgamma(colSums(Cdata)+zmax*alpha))-JJ[2]*(zmax*lgamma(alpha)-lgamma(zmax*alpha))
  else
    loglik <- sum(lgamma(Cdata+t(repmat(alpha,JJ[2],1))))-sum(lgamma(colSums(Cdata)+t(repmat(sum(alpha),JJ[2],1))))-JJ[2]*(sum(lgamma(alpha))-lgamma(sum(alpha)))
  
  # likelihood for k_pp, #clusters
  p <- size(pM,1) # number of important predictors
  if(p==1)
    pM <- as.matrix(t(pM))
  for(j in 1:p) {
    d <- sum(pM[j,]>0) # number of possible number of clusters
    loglik <- loglik+log(pM[j,M[j]]/as.numeric(Stirling2(d,M[j]))) # S(n,k): # of ways to partition n objects into k groups 
  }
  return(loglik)
}

############## approximate gamma shape function ############## 
approx_gamma_shape <- function(dat,mu,a0,b0,it=10,e=10^(-8)) {
  R <- sum(log(dat))
  S <- sum(dat)
  n <- length(dat)
  Q <- S/mu-R+n*log(mu)-n
  A <- a0+n/2
  B <- b0+Q
  for(j in 1:it){
    a <- A/B
    A <- a0-n*a+n*(a^2)*trigamma(a)
    B <- b0+(A-a0)/a-n*log(a)+n*digamma(a)+Q
    if(abs(a/(A/B)-1)<e){
      return(rgamma(1,shape=A,rate=B))
    }
  }
  return(rgamma(1,shape=A,rate=B))
}
