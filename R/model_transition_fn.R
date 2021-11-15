# Supporting functions for model_transition.R

############## likelihood function ##############
trans_logml <- function(zz,Yt,Ytminus1,v,lambda0_mat_expanded,dpreds) {
  d0 <- length(unique(Yt))
  v0 <- cbind(Yt[v==1],Ytminus1[v==1],zz[v==1,])
  v0 <- v0[do.call(order,as.data.frame(v0)),]
  v00 <- unique(v0)
  m00starts <- which(duplicated(v0)==FALSE) 
  Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1) 
  C <- array(0,dim=c(d0,d0,dpreds))
  C[as.matrix(v00)] <- C[as.matrix(v00)]+Ts
  Cdata <- matrix(C,nrow=d0)
  lambda0_mat_expanded[lambda0_mat_expanded==0] <- 10^-5
  loglik <- sum(gammaln(Cdata+lambda0_mat_expanded))-sum(gammaln(colSums(Cdata)+colSums(lambda0_mat_expanded)))-sum(gammaln(lambda0_mat_expanded))+sum(gammaln(colSums(lambda0_mat_expanded)))
  return(loglik)
}
