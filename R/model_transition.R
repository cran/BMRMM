################################################################################
### This function is used when ISI is ignored or modeled as a discrete state ###
################################################################################

model_transition <- function(data,switch_off_random,simsize,burnin) {

  ################## process data ##################
  # data columns should be: id, selected covariates, prev state, current state
  # construct isi data and covs for simulation
  num_row <- nrow(data)
  Xexgns <- as.matrix(data[,2:(ncol(data)-2)]) # load covariate values
  if(ncol(Xexgns)==1) {
    colnames(Xexgns) <- colnames(data)[2]
  }
  Id <- data[,1]
  Ytminus1 <- data[,ncol(data)-1] # previous state
  Yt <- data[,ncol(data)] # current state

  d0 <- length(unique(Yt)) # number of unique states
  sz <- dim(Xexgns)
  p <- sz[2] # number of covariates
  dpreds <- as.vector(apply(Xexgns,2,function(x) length(unique(x)))) # size of each covariate
  num_pairs <- rep(0,p)
  for(pp in 1:p) {
    num_pairs[pp] <- choose(dpreds[pp],2)
  }

  v0 <- cbind(Xexgns,Id)
  v0 <- v0[do.call(order,as.data.frame(v0)),]
  v00 <- unique(v0) # unique combination of (x_{s,1},x_{s,2},...,id)
  m00starts <- which(duplicated(v0)==FALSE) # start of each unique combination
  Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1) # frequency of each combination

  ################ assign priors ###################
  pialpha <- ones(1,p)
  dmax <- max(dpreds)
  lambdaalpha_exgns <- 1
  lambdaalpha_anmls <- 1
  lambdaalpha0 <- 1

  ################ MCMC sampler ####################
  N_MCMC <- simsize
  N_Thin <- 5
  N_Store <- 0
  if(is.null(burnin))
    burnin <- floor(N_MCMC/2)
  np <- p
  Xnew <- Xexgns
  M <- repmat(dpreds,N_MCMC+1,1)
  G <- zeros(p,dmax)
  pi <- zeros(p,dmax)
  logmarginalprobs <- zeros(p,dmax)
  Ntot <- sz[1]
  z <- ones(Ntot,p) # covariate values

  for(j in 1:p) {
    G[j,1:dpreds[j]] <- 1:dpreds[j]
    z[,j] <- G[j,Xnew[,j]]
    pi[j,1:dpreds[j]] <- 1/dpreds[j]
  }
  GG <- G
  log0 <- zeros(N_MCMC,1)
  ind_all <- unique(sortrows(as.matrix(cbind(Xnew,Id))))

  if (switch_off_random) {
    piv <- ones(max(Id),d0) # probability for population-level effect
  } else {
    piv <- 0.8*ones(max(Id),d0) # probability for population-level effect
  }
  v <- sample(0:1,Ntot,TRUE,prob=c(1-piv[1,1],piv[1,1])) # v=0:anmls; v=1:exgns

  ### estimate TP_All (transition probabilities for all)
  C <- array(0,dim=c(d0,d0,dpreds,max(Id)))
  v0 <- cbind(Yt,Ytminus1,Xnew,Id)
  v0 <- v0[do.call(order,as.data.frame(v0)),]
  v00 <- unique(v0)
  m00starts <- which(duplicated(v0)==FALSE)
  Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
  C[as.matrix(v00)] <- C[as.matrix(v00)]+Ts
  T_All <- C
  TP_All <- array(0,dim=c(d0,d0,dpreds,max(Id)))
  for(i in 1:nrow(ind_all)) {
    attr <- ind_all[i,]
    for(j in 1:d0) {
      index <- cbind(1:d0,repmat(j,d0,1),repmat(attr,d0,1))
      TP_All[index] <- C[index]/sum(C[index])
    }
  }

  ### estimate TP_Anmls
  T_Anmls <- apply(T_All,c(1,2,length(size(T_All))),sum)
  TP_Anmls <- array(0,dim=c(d0,d0,max(Id)))
  for(i in 1:max(Id)) {
    for(j in 1:d0) {
      index <- cbind(1:d0,repmat(j,d0,1),repmat(i,d0,1))
      TP_Anmls[index] <- T_Anmls[index]/sum(T_Anmls[index])
    }
  }

  ### estimate TP_Exgns
  T_Exgns <- apply(T_All,1:(length(size(T_All))-1),sum)
  TP_Exgns <- array(0,dim=c(d0,d0,dpreds))
  v0 <- unique(Xnew)
  for(i in 1:nrow(v0)) {
    attr <- as.numeric(v0[i,])
    for(j in 1:d0) {
      index <- cbind(1:d0,repmat(j,d0,1),repmat(attr,d0,1))
      TP_Exgns[index] <- T_Exgns[index]/sum(T_Exgns[index])
    }
  }

  ### initialize lambda00
  C <- as.matrix(table(Yt))
  lambda00 <- C/sum(C)
  lambda00_mat_expanded <- repmat(lambda00,1,d0)

  ### initialize lambda0
  C <- zeros(d0,d0)
  v0 <- cbind(Yt,Ytminus1)
  v0 <- v0[do.call(order,as.data.frame(v0)),]
  v00 <- unique(v0)
  m00starts <- which(duplicated(v0)==FALSE)
  Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
  C[as.matrix(v00)] <- C[as.matrix(v00)]+Ts
  lambda0 <- C/repmat(colSums(C),d0,1)
  lambda0_mat_expanded_exgns <- repmat(lambda0,1,prod(dpreds))
  lambda0_mat_expanded_anmls <- repmat(lambda0,1,max(Id))

  ### initialize lambda_exgns
  C <- array(0,dim=c(d0,d0,dpreds))
  v0 <- cbind(Yt,Ytminus1,z)
  v0 <- v0[do.call(order,as.data.frame(v0)),]
  v00 <- unique(v0)
  m00starts <- which(duplicated(v0)==FALSE)
  Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
  C[as.matrix(v00)] <- C[as.matrix(v00)]+Ts
  Cdata <- matrix(C,nrow=d0)
  lambda_exgns_mat <- zeros(d0,size(Cdata)[2])
  for(i in 1:ncol(Cdata)) {
    lambda_exgns_mat[,i] <- rgamma(d0,shape=Cdata[,i]+lambdaalpha_exgns*lambda0_mat_expanded_exgns[,i],scale=1)
    lambda_exgns_mat[,i] <- lambda_exgns_mat[,i]/sum(lambda_exgns_mat[,i])
  }
  lambda_exgns_mat[lambda_exgns_mat==0] <- min(min(lambda_exgns_mat[lambda_exgns_mat>0]),0.0001)
  lambda_exgns <- array(lambda_exgns_mat,dim=c(d0,d0,dpreds))

  ### initialize lambda_anmls
  C <- array(0,dim=c(d0,d0,max(Id)))
  v0 <- cbind(Yt,Ytminus1,Id)
  v0 <- v0[do.call(order,as.data.frame(v0)),]
  v00 <- unique(v0)
  m00starts <- which(duplicated(v0)==FALSE)
  Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
  C[as.matrix(v00)] <- C[as.matrix(v00)]+Ts
  Cdata <- matrix(C,nrow=d0)
  lambda_anmls_mat <- zeros(d0,size(Cdata)[2])
  for(i in 1:ncol(Cdata)) {
    lambda_anmls_mat[,i] <- rgamma(d0,shape=Cdata[,i]+lambdaalpha_anmls*lambda0_mat_expanded_anmls[,i],scale=1)
    lambda_anmls_mat[,i] <- lambda_anmls_mat[,i]/sum(lambda_anmls_mat[,i])
  }
  lambda_anmls_mat[lambda_anmls_mat==0] <- min(min(lambda_anmls_mat[lambda_anmls_mat>0]),0.0001)
  lambda_anmls <- array(lambda_anmls_mat,dim=c(d0,d0,max(Id)))

  ### MCMC storage
  lambda_anmls_mat_tmp <- array(0,dim=c(d0,d0,max(Id)))
  TP_All_Post_Mean <- array(0,dim=c(d0,d0,dpreds,max(Id)))
  TP_Exgns_Post_Mean <- array(0,dim=c(d0,d0,dpreds))
  TP_Exgns_Store <- array(0,dim=c(d0,d0,dpreds,floor((N_MCMC-burnin)/N_Thin)))
  TP_Exgns_Diffs_Store <- array(0,dim=c(p,max(num_pairs),d0,d0,rep(max(dpreds),length(dpreds)-1),floor((N_MCMC-burnin)/N_Thin)))
  TP_Anmls_Store <- array(0,dim=c(d0,d0,max(Id),floor((N_MCMC-burnin)/N_Thin)))
  TP_Exgns_Comp_Post_Mean <- array(0,dim=c(d0,d0,dpreds))
  TP_Anmls_Comp_Post_Mean <- array(0,dim=c(d0,d0,max(Id)))
  TP_Anmls_Comp_Store <- array(0,dim=c(d0,d0,max(Id),floor((N_MCMC-burnin)/N_Thin)))
  TP_Exgns_All_Itns <- array(0,dim=c(d0,d0,dpreds,N_MCMC))

  ### start sampler
  for(kkk in 1:N_MCMC) {

    ### updating z (#clusters)
    M00 <- M[kkk,]
    if(kkk>1) {
      for(j in 1:p) {
        for(k in 1:dpreds[j]) {
          for(l in 1:dpreds[j]) {
            GG[j,k] <- l
            zz <- z
            zz[,j] <- GG[j,Xnew[,j]]
            logmarginalprobs[j,l] <- trans_logml(zz,Yt,Ytminus1,v,lambdaalpha_exgns*lambda0_mat_expanded_exgns,dpreds)
          }
          logmarginalprobs[j,1:dpreds[j]] <- logmarginalprobs[j,1:dpreds[j]]-max(logmarginalprobs[j,1:dpreds[j]])
          logprobs <- log(pi[j,1:dpreds[j]])+logmarginalprobs[j,1:dpreds[j]]
          probs <- exp(logprobs)/sum(exp(logprobs))
          GG[j,k] <- sample(1:dpreds[j],1,TRUE,probs)
        }
        G[j,] <- GG[j,]
        z[,j] <- GG[j,Xnew[,j]]
        M00[j] <- length(unique(z[,j]))
      }
      M[kkk+1,] <- M00
    }

    ### updating pi
    for(j in 1:p) {
      uzj <- unique(z[,j])
      ztab2 <- zeros(1,dpreds[j])
      ztab2[1:max(uzj)] <- rep(1,max(uzj))
      pi[j,1:dpreds[j]] <- rgamma(dpreds[j],shape=pialpha[j]+ztab2,scale=pialpha[j]+1)
      pi[j,1:dpreds[j]] <- pi[j,1:dpreds[j]]/sum(pi[j,1:dpreds[j]])
    }

    if (!switch_off_random) {
      ### updating v (v=0: anmls; v=1: exgns)
      prob <- zeros(Ntot,2)
      piv_mat <- matrix(piv,nrow=max(Id))
      prob[,1] <- (1-piv_mat[cbind(Id,Ytminus1)])*lambda_anmls[cbind(Yt,Ytminus1,Id)]
      prob[,2] <- piv_mat[cbind(Id,Ytminus1)]*lambda_exgns[cbind(Yt,Ytminus1,z)]
      prob <- prob/rowSums(prob)
      v <- rbinom(Ntot,1,prob[,2])

      ### updating piv
      C <- array(0,dim=c(2,d0,max(Id)))
      v0 <- cbind(v+1,Ytminus1,Id)
      v0 <- v0[do.call(order,as.data.frame(v0)),]
      v00 <- unique(v0)
      m00starts <- which(duplicated(v0)==FALSE)
      Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
      C[as.matrix(v00)] <- C[as.matrix(v00)]+Ts
      for(j in 1:max(Id)) {
        Cdata <- C[,,j]
        piv[j,] <- 1-rbeta(d0,Cdata[1,]+1,Cdata[2,]+1)
      }
    }

    ### updating lambda_exgns
    C1 <- array(0,dim=c(d0,d0,dpreds))
    v0 <- cbind(Yt[v==1],Ytminus1[v==1],z[v==1,])
    v0 <- v0[do.call(order,as.data.frame(v0)),]
    v00 <- unique(v0)
    m00starts <- which(duplicated(v0)==FALSE)
    Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
    C1[as.matrix(v00)] <- C1[as.matrix(v00)]+Ts
    C1data <- matrix(C1,nrow=d0)
    sz1 <- size(C1data)
    lambda_exgns_mat <- zeros(d0,sz1[2])+NaN
    for(j in 1:sz1[2]) {
      while(sum(is.na(lambda_exgns_mat[,j]))>0) {
        lambda_exgns_mat[,j] <- rgamma(d0,shape=C1data[,j]+lambdaalpha_exgns*lambda0_mat_expanded_exgns[,j],scale=1)
        lambda_exgns_mat[,j] <- lambda_exgns_mat[,j]/sum(lambda_exgns_mat[,j])
      }
    }
    lambda_exgns_mat[lambda_exgns_mat==0] <- min(min(lambda_exgns_mat[lambda_exgns_mat>0]),0.0001)
    lambda_exgns <- array(lambda_exgns_mat,dim=c(d0,d0,dpreds))
    cov_combs <- unique(Xexgns)
    for(row in 1:nrow(cov_combs)) {
      v0 <- expand.grid(1:d0,1:d0)
      v0 <- as.matrix(cbind(v0,repmat(as.numeric(cov_combs[row,]),nrow(v0),1)))
      v00 <- cbind(v0,rep(kkk,nrow(v0)))
      TP_Exgns_All_Itns[v00] <- lambda_exgns[v0]
    }

    ### updating lambda_anmls
    C2 <- array(0,dim=c(d0,d0,max(Id)))
    v0 <- cbind(Yt[v==0],Ytminus1[v==0],Id[v==0])
    v0 <- v0[do.call(order,as.data.frame(v0)),]
    v00 <- unique(v0)
    m00starts <- which(duplicated(v0)==FALSE)
    Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
    C2[as.matrix(v00)] <- C2[as.matrix(v00)]+Ts
    C2data <- matrix(C2,nrow=d0)
    sz2 <- size(C2data)
    lambda_anmls_mat <- zeros(d0,sz2[2])+NaN
    for(j in 1:sz2[2]) {
      while(sum(is.na(lambda_anmls_mat[,j]))>0) {
        lambda_anmls_mat[,j] <- rgamma(d0,shape=C2data[,j]+lambdaalpha_anmls*lambda0_mat_expanded_anmls[,j],scale=1)
        lambda_anmls_mat[,j] <- lambda_anmls_mat[,j]/sum(lambda_anmls_mat[,j])
      }
    }
    lambda_anmls_mat[lambda_anmls_mat==0] <- min(min(lambda_anmls_mat[lambda_anmls_mat>0]),0.0001)
    lambda_anmls <- array(lambda_anmls_mat,dim=c(d0,d0,max(Id)))

    ### updating hyper-parameters
    if (kkk>N_MCMC/4) {

      ### updating lambdaalpha_exgns
      C1 <- array(0,dim=c(d0,d0,dpreds))
      v0 <- cbind(Yt[v==1],Ytminus1[v==1],z[v==1,])
      v0 <- v0[do.call(order,as.data.frame(v0)),]
      v00 <- unique(v0)
      m00starts <- which(duplicated(v0)==FALSE)
      Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
      C1[as.matrix(v00)] <- C1[as.matrix(v00)]+Ts
      C1data <- matrix(C1,nrow=d0)
      sz1 <- size(C1data)
      vmat_exgns <- zeros(sz1[1],sz1[2])
      for(j in 1:sz1[2]) {
        jtemp <- j%%d0
        jtemp[jtemp==0] <- d0
        for(i in 1:sz1[1]) {
          if(C1data[i,j]>0) {
            prob <- lambdaalpha_exgns*lambda0[i,jtemp]/((1:C1data[i,j])-1+lambdaalpha_exgns*lambda0[i,jtemp])
            prob[prob<0] <- 0
            vmat_exgns[i,j] <- vmat_exgns[i,j]+sum(rbinom(length(prob),1,prob))
          }
        }
      }
      vmat_exgns_colsums <- colSums(vmat_exgns)
      vmat_exgns <- array(vmat_exgns,dim=c(d0,d0,prod(dpreds)))
      vmat_exgns <- apply(vmat_exgns,c(1,2),sum)
      C1data_colsums <- colSums(C1data)
      rmat_exgns <- rbeta(sz1[2],lambdaalpha_exgns+1,C1data_colsums)
      smat_exgns <- rbinom(sz1[2],1,C1data_colsums/(C1data_colsums+lambdaalpha_exgns))
      lambdaalpha_exgns <- rgamma(1,shape=1+sum(vmat_exgns_colsums)-sum(smat_exgns),scale=1/(1-sum(log(rmat_exgns))))

      ### updating lambdaalpha_anmls
      C2 <- array(0,dim=c(d0,d0,max(Id)))
      v0 <- cbind(Yt[v==0],Ytminus1[v==0],Id[v==0])
      v0 <- v0[do.call(order,as.data.frame(v0)),]
      v00 <- unique(v0)
      m00starts <- which(duplicated(v0)==FALSE)
      Ts <- c(diff(m00starts),nrow(v0)-m00starts[length(m00starts)]+1)
      C2[as.matrix(v00)] <- C2[as.matrix(v00)]+Ts
      C2data <- matrix(C2,nrow=d0)
      sz2 <- size(C2data)
      vmat_anmls <- zeros(sz2[1],sz2[2])
      for(j in 1:sz2[2]) {
        jtemp <- j%%d0
        jtemp[jtemp==0] <- d0
        for(i in 1:sz2[1]) {
          if(C2data[i,j]>0) {
            prob <- lambdaalpha_anmls*lambda0[i,jtemp]/((1:C2data[i,j])-1+lambdaalpha_anmls*lambda0[i,jtemp])
            prob[prob<0] <- 0
            vmat_anmls[i,j] <- vmat_anmls[i,j]+sum(rbinom(length(prob),1,prob))
          }
        }
      }
      vmat_anmls_colsums <- colSums(vmat_anmls)
      vmat_anmls <- array(vmat_anmls,dim=c(d0,d0,max(Id)))
      vmat_anmls <- apply(vmat_anmls,c(1,2),sum)
      C2data_colsums <- colSums(C2data)
      rmat_anmls <- rbeta(sz2[2],lambdaalpha_anmls+1,C2data_colsums)
      smat_anmls <- rbinom(sz2[2],1,C2data_colsums/(C2data_colsums+lambdaalpha_anmls))
      lambdaalpha_anmls <- rgamma(1,shape=1+sum(vmat_anmls_colsums)-sum(smat_anmls),scale=1/(1-sum(log(rmat_anmls))))

      ### updating lambda0
      a <- vmat_exgns+vmat_anmls
      sz <- size(lambda0)
      for(j in 1:sz[2]) {
        lambda0[,j] <- rgamma(d0,shape=a[,j]+lambdaalpha0*lambda00_mat_expanded[,j],scale=1);
        lambda0[,j] <- lambda0[,j]/sum(lambda0[,j])
      }
      lambda0[lambda0==0] <- min(min(lambda0[lambda0>0]),0.0001)
      lambda0_mat_expanded_exgns <- repmat(lambda0,1,prod(dpreds))
      lambda0_mat_expanded_anmls <- repmat(lambda0,1,max(Id))
    }

    ### print progress of MCMC sampler
    if(kkk%%100==0) print(paste("Transition Probabilities: Iteration ",kkk))
    #ind1 <- which(M00>1)
    #np <- length(ind1)
    #cat(sprintf('k=%i, %i important predictors = {',kkk,np))
    #for(i in 1:length(ind1)) {
    #  cat(sprintf(' X(%i)(%i)',ind1[i],M00[ind1[i]]))
    #}
    #cat(sprintf(' }. mean piv=%f, lambdaalpha_exgns=%f, lambdaalpha_anmls=%f \n',mean(piv),lambdaalpha_exgns,lambdaalpha_anmls))

    ### storage after burn-in
    if(kkk>burnin & kkk%%N_Thin==0) {
      N_Store <- N_Store+1
      pistar <- array(0,dim=c(p,dmax,dmax))
      idx <- repmat(((1:dmax)-1)*dmax,p,1)+G
      for(j in 1:p) {
        V <- zeros(dmax,dmax)
        V[idx[j,]] <- 1
        pistar[j,1:dmax,1:dmax] <- t(V)
      }
      TP_Exgns_Comp <- lambda_exgns
      TP_Exgns_Comp_Mean <- apply(array(lambda_exgns,dim=c(d0,d0,prod(dpreds))),c(1,2),sum)/prod(dpreds)
      TP_Exgns_Comp_Mean_Anmls <- array(repmat(TP_Exgns_Comp_Mean,max(Id)),dim=c(d0,d0,max(Id)))
      all_trans <- cbind(rep(1:d0,each=d0),rep(1:d0,d0))

      for(i in 1:nrow(ind_all)) {
        attr <- ind_all[i,1:(ncol(ind_all)-1)]
        id <- ind_all[i,ncol(ind_all)]
        lambda_anmls_mat_tmp[,,id] <- lambda_anmls[,,id]
        all_index <- as.matrix(cbind(all_trans,repmat(c(attr,id),nrow(all_trans),1)))
        anmls_mat <- matrix(lambda_anmls[all_index[,c(1,2,ncol(all_index))]],nrow=d0)
        exgns_mat <- matrix(lambda_exgns[all_index[,c(1:(ncol(all_index)-1))]],nrow=d0)
        TP_All_Post_Mean_Temp <- repmat(1-piv[id,],d0,1)*anmls_mat+repmat(piv[id,],d0,1)*exgns_mat
        TP_All_Post_Mean[all_index] <- TP_All_Post_Mean[all_index]+TP_All_Post_Mean_Temp[all_trans]
      }
      TP_Anmls_Comp <- lambda_anmls
      TP_Anmls_Comp_Mean <- apply(lambda_anmls_mat_tmp,c(1,2),sum)/max(Id)
      TP_Anmls_Comp_Mean_Exgns <- array(repmat(TP_Anmls_Comp_Mean,1,prod(dpreds)),dim=c(d0,d0,dpreds))
      TP_Exgns_kkk <- (lambda0_mat_expanded_exgns+array(TP_Exgns_Comp,dim=c(d0,d0*prod(dpreds))))/2
      TP_Exgns_kkk <- array(TP_Exgns_kkk,dim=c(d0,d0,dpreds))
      TP_Exgns_Post_Mean <- TP_Exgns_Post_Mean + TP_Exgns_kkk
      TP_Exgns_Comp_Post_Mean <- TP_Exgns_Comp_Post_Mean + TP_Exgns_Comp
      TP_Anmls_Comp_Post_Mean <- TP_Anmls_Comp_Post_Mean + TP_Anmls_Comp
      v0 <- unique(cbind(Yt,Ytminus1,Xnew))
      v00 <- cbind(v0,rep(N_Store,nrow(v0)))
      TP_Exgns_Store[as.matrix(v00)] <- TP_Exgns_kkk[as.matrix(v0)]
      for(j in 1:max(Id)) {
        TP_Anmls_Comp_Store[,,j,N_Store] <- (1-piv[j,])*(lambda_anmls_mat_tmp[,,j]-lambda0)
      }
      for(pp in 1:p) {
        num_pair <- num_pairs[pp]
        pair <- combn(dpreds[pp],2)
        inds <- 1:p
        inds <- inds[inds!=pp]
        for(np in 1:num_pair) {
          lvl1 <- pair[1,np]
          lvl2 <- pair[2,np]
          v0 <- unique(cbind(Ytminus1,Yt,Xexgns[,inds]))
          v1 <- zeros(nrow(v0),2+p)
          v1[,1:2] <- as.matrix(v0[,1:2])
          if (length(inds)>0) {
            v1[,inds+2] <- as.matrix(v0[,3:ncol(v0)])
          }
          v1[,pp+2] <- lvl1
          v2 <- v1
          v2[,pp+2] <- lvl2
          v0 <- cbind(rep(pp,nrow(v0)),rep(np,nrow(v0)),v0,rep(N_Store,nrow(v0)))
          TP_Exgns_Diffs_Store[as.matrix(v0)] <- TP_Exgns_kkk[v1]-TP_Exgns_kkk[v2]
        }
      }
    }
  }

  TP_Exgns_Post_Mean <- TP_Exgns_Post_Mean/N_Store
  TP_Exgns_Post_Std <- apply(TP_Exgns_Store,1:(length(size(TP_Exgns_Store))-1),sd)

  TP_Exgns_Comp_Post_Mean <- TP_Exgns_Comp_Post_Mean/N_Store
  TP_Anmls_Comp_Post_Mean <- TP_Anmls_Comp_Post_Mean/N_Store
  TP_All_Post_Mean <- TP_All_Post_Mean/N_Store

  TP_Anmls_Comp_Post_Std <- apply(TP_Anmls_Comp_Store,c(1,2),sd)

  results <- list("Num_States"=d0,
                  "Xexgns"=Xexgns,
                  "dpreds"=dpreds,
                  "MCMCparams"=list("simsize"=simsize,"burnin"=burnin,"N_Thin"=N_Thin),
                  "TP_Exgns_Post_Mean"=TP_Exgns_Post_Mean,
                  "TP_Exgns_Post_Std"=TP_Exgns_Post_Std,
                  "TP_Anmls_Post_Mean"=TP_Anmls_Comp_Post_Mean,
                  "TP_All_Post_Mean"=TP_All_Post_Mean,
                  "TP_Anmls_Post_Std"=TP_Anmls_Comp_Post_Std,
                  "TP_Exgns_Diffs_Store"=TP_Exgns_Diffs_Store,
                  "TP_Exgns_All_Itns"=TP_Exgns_All_Itns,
                  "Clusters"=M,
                  "Type"="Transition Probabilities")
  return(results)
}
