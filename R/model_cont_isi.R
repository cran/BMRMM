###################################################################
### This function is used when ISI is modeled as mixture gammas ###
###################################################################

model_cont_isi <- function(data,switch_off_random,simsize,burnin,K,init_alpha,init_beta) {

  ################## process data ##################
  # construct isi data and covs for simulation
  num_row <- nrow(data)
  covs <- data[,2:(ncol(data)-1)] # load covariate values
  mouseId <- data[,1]
  isi <- data[,ncol(data)]


  ################## initialization ##################

  p <- ncol(covs) # number of covariates
  covs.max <- as.vector(apply(covs,2,function(x) length(unique(x)))) # size of each covariate
  all_combs <- unique(covs)
  num_combs <- nrow(all_combs)
  z.max <- K # number of mixture components

  M00 <- covs.max
  M <- matrix(0,nrow=simsize,ncol=p) # record # of clusters for each iteration

  clusts <- kmeans(isi,z.max)$cluster # get clusters via kmeans
  z.vec <- clusts # initialization for mixture component assignment
  lambda00 <- sapply(1:z.max,function(k) sum(clusts==k)/length(clusts)) # lambda_{isi,00}
  lambdaalpha0 <- 1 # alpha_{isi,00}
  lambda0 <- lambda00 # lambda_{isi,0}
  lambdaalpha <- 0.01 # alpha_{isi,0}
  lambda <- array(lambda00,dim=c(z.max,1)) # lambda_{g1,g2,g3}
  # individual probability and lambda
  if (switch_off_random) {
    piv <- matrix(1,nrow=max(mouseId),ncol=z.max) # probability for population-level effect
  } else {
    piv <- matrix(0.8,nrow=max(mouseId),ncol=z.max) # probability for population-level effect
  }
  v <- sample(0:1,num_row,prob=c(piv[1,1],1-piv[1,1]),replace=TRUE) # v=0: exgns, v=1: anmls
  lambda_mice <- matrix(rep(lambda00,each=max(mouseId)),nrow=max(mouseId))
  lambda_mice_alpha <- 0.01 # alpha_{isi}^{0}

  # record all marginal prob for gamma mixture of each iteration
  prob_mat <- matrix(0,nrow=simsize,ncol=z.max)

  # mixture gamma parameters
  Alpha <- matrix(1,nrow=simsize,ncol=z.max) # shape
  Beta <- matrix(1,nrow=simsize,ncol=z.max) # rate

  # cluster assignment for each covariate
  z.cov <- as.matrix(covs) # initially, each covariate forms its own cluster
  G <- matrix(0,nrow=p,ncol=max(covs.max))
  for(pp in 1:p){
    G[pp,1:covs.max[pp]] <- 1:covs.max[pp]
  }

  # pM: prior probabilities of # of clusters for each covariate
  phi <- 0.25
  pM <- matrix(0,p,max(covs.max))
  for(pp in 1:p){
    pM[pp,1:covs.max[pp]] <- exp(-(pp*(1:covs.max[pp])*phi-1)) # prior probability for k_{pp}, #cluster for covariate pp
    pM[pp,] <- pM[pp,]/sum(pM[pp,])
  }
  # cM: for each covariate, # of ways to partition each j elements into two empty groups, 1<=j<=(size of covariate)
  cM <- matrix(0,p,max(covs.max))                   # There are (2^r-2)/2 ways to split r objects into two non-empty groups.
  for(pp in 1:p)                                  # (Consider 2 cells and r objects, subtract 2 for the two cases in which
    cM[pp,1:covs.max[pp]] <- c(2^(1:covs.max[pp])-2)/2 # one cell receives all objects, divide by 2 for symmetry.)

  # MCMC storage
  lambda0_all <- matrix(0,nrow=simsize,ncol=z.max) # store lambda_{0}
  lambda_all <- array(0,dim=c(simsize,covs.max,z.max)) # store lambda_{g1,g2,g3}
  lambda_mice_all <- array(0,dim=c(simsize,max(mouseId),z.max)) # store lambda^{(i)}
  lambdaalpha_all <- rep(0,simsize) # store alpha_0
  lambdaalpha_mice_all <- rep(0,simsize) # store alpha^0

  ################## simulation ##################

  for(iii in 1:simsize) {
    M[iii,] <- M00
    prob_mat[iii,] <- as.numeric(table(factor(z.vec,levels=1:z.max))/num_row)
    if(iii%%100==0) print(paste("Duration Times: Iteration ",iii))

    # update # cluster for each covariate and cluster mapping
    if(iii > 1) {
      # update cluster indicator z
      M00 <- M[iii,]
      for(pp in 1:p){
        M0 <- M00[pp] # current # of clusters of covariate pp
        if(M0==1) {
          # propose to split one cluster into 2
          new <- rbinom(covs.max[pp]-1,size=1,prob=0.5)	# propose new mapping for (d_{pp}-1) values at 1
          while(sum(new)==0)
            new <- rbinom(covs.max[pp]-1,size=1,prob=0.5)
          GG <- G[pp,1:covs.max[pp]] + c(0,new)   # keep the first one at 1, propose new cluster mappings for the other (d_{pp}-1) levels of x_{pp}
          zz <- z.cov                # zz initiated at z
          zz[,pp] <- GG[covs[,pp]]   # proposed new zz by mapping to new cluster configurations for the observed values of x_{pp}
          ind2 <- which(M00>1)       # current set of relevant predictors
          if(length(ind2)==0)        # if no predictor is currently important
            ind2 <- 1
          MM <- M00                  # MM initiated at {k_{1},...,k_{p}}, current values, with k_{pp}=1 by the if condition
          MM[pp] <- 2                # proposed new value of k_{pp}=2
          ind1 <- which(MM>1)        # proposed set of important predictors, now includes x_{pp}, since MM[pp]=2
          logR <- isi_logml(zz[,ind1],z.vec,MM[ind1],pM[ind1,],lambdaalpha*lambda0,z.max)-isi_logml(z.cov[,ind2],z.vec,M00[ind2],pM[ind2,],lambdaalpha*lambda0,z.max)
          logR <- logR+log(0.5)+log(cM[pp,covs.max[pp]]) # computational reasons
          if(log(runif(1))<logR) {
            G[pp,1:covs.max[pp]] <- GG
            M00 <- MM
            z.cov <- zz
            np <- length(ind2)
          }
        }
        if((M0>1)&&(M0<covs.max[pp])){
          if(runif(1)<0.5) {  # with prob 0.5 split one mapped value into two
            df <- sort(G[pp,1:covs.max[pp]]) # cluster mapping for covariate pp
            z0 <- unique(df)   # z0 are unique cluster mappings, mm contains their positions
            mm <- which(duplicated(df)==FALSE)             # mm contains the positions when they first appear on {1,...,n}
            gn <- c(diff(mm),size(df,2)-mm[length(mm)]+1)  # frequencies of z0
            pgn <- cM[pp,gn]/sum(cM[pp,gn])                # see the definition of cM
            rr <- sum(rmultinom(1,size=1,prob=pgn)*(1:M0)) # rr is the state to split, gn[rr] is the frequency of rr
            new <- rmultinom(1,size=1,0.5*rep(1,gn[rr]-1)) # propose new mapping for (gn(rr)-1) values at rr
            while(sum(new)==0)
              new <- rmultinom(1,size=1,0.5*rep(1,gn[rr]-1))
            GG <- G[pp,1:covs.max[pp]]
            GG[GG==rr] <- rr+(M0+1-rr)*c(0,new)         # keep first value at rr, propose new mapping (M0+1) for the rest of (gn[rr]-1) values at rr
            zz <- z.cov        # zz initiated at z.cov
            zz[,pp] <- GG[covs[,pp]]   # proposed new zz by mapping to new cluster configurations for the observed values of x_{pp}
            MM <- M00            # MM initizted at current values {k_{1},...,k_{p}}
            MM[pp] <- M0+1       # proposed new value of k_{pp}, since one original mapped value is split into two
            ind1 <- which(MM>1)  # proposed set of important predictors, now includes x_{pp}, since MM[pp]=2
            ind2 <- which(M00>1) # current set of important predictors
            if(length(ind2)==0)   # if no predictor is currently important
              ind2 <- 1
            logR <- isi_logml(zz[,ind1],z.vec,MM[ind1],pM[ind1,],lambdaalpha*lambda0,z.max)-isi_logml(z.cov[,ind2],z.vec,M00[ind2],pM[ind2,],lambdaalpha*lambda0,z.max)
            if(M00[pp]<(covs.max[pp]-1)){
              logR = logR-log(M0*(M0+1)/2)+log(sum(cM[pp,gn]))
            } else {
              logR = logR-log(covs.max[pp]*(covs.max[pp]-1)/2)-log(0.5) # note that by replacing -log(0.5) with +log(2), we get the same reasoning as above
            }
            if(log(runif(1))<logR) {
              G[pp,1:covs.max[pp]] <- GG
              M00 <- MM
              z.cov <- zz
              np <- length(ind2)
            }
          } else {        # with prob 0.5 merge two mapped values into one
            znew <- sample(M0,2)
            lnew <- max(znew)
            snew <- min(znew)
            GG <- G[pp,1:covs.max[pp]]
            GG[GG==lnew] <- snew # replace all lnews by snews
            GG[GG==M0] <- lnew   # replace the largest cluster mapping by lnew (lnew itself may equal M0, in which case GG remains unchanged by this move)
            zz <- z.cov        # zz initiated at z.cov
            zz[,pp] <- GG[covs[,pp]]   # proposed new zz by mapping to new cluster configurations for the observed values of x_{pp}
            MM <- M00            # MM initizted at current values {k_{1},...,k_{p}}, with k_{pp}=d_{pp} by the if condition
            MM[pp] <- M00[pp]-1  # proposed new value of k_{pp}, since two mappings are merged
            ind1 <- which(MM>1)  # proposed set of important predictors, may not include x_{pp} if original k_(j) was at 2
            ind2 <- which(M00>1)  # current set of important predictors
            if(length(ind1)==0)
              ind1 <- 1
            if(length(ind2)==0)
              ind2 <- 1
            logR <- isi_logml(zz[,ind1],z.vec,MM[ind1],pM[ind1,],lambdaalpha*lambda0,z.max)-isi_logml(z.cov[,ind2],z.vec,M00[ind2],pM[ind2,],lambdaalpha*lambda0,z.max)
            if(M0>2) {
              df <- sort(GG)
              z0 <- unique(df)   # z0 are unique cluster mappings, mm contains their positions
              mm <- which(duplicated(df)==FALSE)             # mm contains the positions when they first appear on {1,...,n}
              gn <- c(diff(mm),size(df,2)-mm[length(mm)]+1)  # frequencies of z0
              logR <- logR-log(sum(cM[pp,gn]))+log(M00[pp]*(M00[pp]-1)/2)
            } else {
              logR <- logR-log(cM[pp,covs.max[pp]])-log(0.5)
            }
            if(log(runif(1))<logR) {
              G[pp,1:covs.max[pp]] <- GG
              M00 <- MM
              z.cov <- zz
              np <- length(ind2)
            }
          }
        }
        if(M0==covs.max[pp]){
          znew <- sample(covs.max[pp],2)
          lnew <- max(znew)
          snew <- min(znew)
          GG <- G[pp,1:covs.max[pp]]
          GG[GG==lnew] <- snew     # replace all lnews by snews
          GG[GG==M0] <- lnew       # replace the largest cluster mapping d_{pp} by lnew (lnew itself can be d_{pp}, in which case GG remains unchanged by this move)
          zz <- z.cov              # zz initiated at z.cov
          zz[,pp] <- GG[covs[,pp]] # proposed new z_{,pp} as per the proposed new cluster mappings of the levels of x_{pp}
          MM <- M00                # MM initizted at current values {k_{1},...,k_{p}}, with k_{pp}=d_{pp} by the if condition
          MM[pp] <- covs.max[pp]-1 # proposed new value of k_{pp}, since originally k_{pp}=d_{pp} and now two mappings are merged
          ind1 <- which(MM>1)      # proposed set of important predictors, does not include x_{pp} when d_{pp}=2
          if(length(ind1)==0)
            ind1 <- 1
          ind2 <- which(M00>1)     # current set of important predictors
          if(length(ind2)==0)
            ind2 <- 1
          logR <- isi_logml(zz[,ind1],z.vec,MM[ind1],pM[ind1,],lambdaalpha*lambda0,z.max)-isi_logml(z.cov[,ind2],z.vec,M00[ind2],pM[ind2,],lambdaalpha*lambda0,z.max)
          logR <- logR+log(0.5)+log(covs.max[pp]*(covs.max[pp]-1)/2)
          if(log(runif(1))<logR) {
            G[pp,1:covs.max[pp]] <- GG
            M00 <- MM
            z.cov <- zz
          }
        }
        # end of #cluster appointment
        # start of cluster mapping
        if(M00[pp]>1){ # propose a new cluster index mapping
          zz <- z.cov              # initiate zz at z.cov
          per <- sample(covs.max[pp],covs.max[pp])
          GG <- G[pp,per]            # proposed new cluster mappings of different levels of x_{pp}
          zz[,pp] <- GG[covs[,pp]] # proposed new z_{,pp} as per the proposed new cluster mappings of the levels of x_{pp}
          MM <- M00
          ind1 <- which(M00>1)
          ind2 <- which(M00>1)
          logR <- isi_logml(zz[,ind1],z.vec,MM[ind1],pM[ind1,],lambdaalpha*lambda0,z.max)-isi_logml(z.cov[,ind2],z.vec,M00[ind2],pM[ind2,],lambdaalpha*lambda0,z.max)
          if(log(runif(1))<logR){
            G[pp,1:covs.max[pp]] <- GG
            z.cov <- zz
            np <- length(ind2)
          }
        }
        if((M00[pp]==1)&&(length(unique(covs.max))==1)&&(iii<simsize/2)){
          ind2 <- which(M00>1)     # the current set of important predictors
          tempind <- sample(length(ind2),1)
          temp <- ind2[tempind]      # choose one, x_{temp}, from the current set of important predictors
          zz <- z.cov            # initiate zz at z.cov
          zz[,temp] <- rep(1,num_row)  # propose removal of x_{temp} by setting z_{i,temp}=1 for all i
          per <- sample(covs.max[pp],covs.max[pp]) # permutation of {1,...,d_{pp}
          GG <- G[temp,per]        # proposed new cluster mappings of different levels of x_{pp} obtained by permuting the cluster mappings of the levels of x_{temp}
          zz[,pp] <- GG[covs[,pp]]  # proposed new z_{i,pp} as per the proposed new cluster mappings of the levels of x_{pp}
          MM <- M00                # MM initiated at current values {k_{1},...,k_{p}}, with k_{pp}=1 and k_{temp}>1 by the conditions
          MM[temp] <- 1              # proposed new value of k_{temp}, since x_{temp} is removed from the set of important predictors
          MM[pp] <- M00[temp]      # propose k_{pp}=k_{temp}
          ind1 <- which(MM>1)        # proposed set of important predictors, now this set excludes x_{temp} but includes x_{pp}
          logR <- isi_logml(zz[,ind1],z.vec,MM[ind1],pM[ind1,],lambdaalpha*lambda0,z.max)-isi_logml(z.cov[,ind2],z.vec,M00[ind2],pM[ind2,],lambdaalpha*lambda0,z.max)
          if(log(runif(1))<logR) {
            G[pp,1:covs.max[pp]] <- GG
            G[temp,] <- rep(1,covs.max[pp])
            M00 <- MM
            z.cov <- zz
            np <- length(ind2)
          }
        }
      }
    }

    ind00 <- which(M00>1)  # selected predictors whose #clusters>1
    if(length(ind00)==0)
      ind00 <- 1
    K00 <- M00[ind00]    # k_{pp}'s for the selected predictors i.e. #clusters for selected predicators
    p00 <- length(ind00)   # number of selected predictors
    z00 <- as.matrix(z.cov[,ind00])

    clT <- array(0,dim=c(2,z.max,K00))   # dmax levels of the response y, K00={k_{1},..,k_{p00}} clustered levels of x_{1},...,x_{p00}
    df <- cbind(v+1,z.vec,z00)
    df <- df[do.call(order,as.data.frame(df)),]
    z0 <- unique(df)                     # z0 are the sorted unique combinations of (y,z_{pp,1},...,z_{pp,p00}),
    m <- which(duplicated(df)==FALSE)    # m contains the positions when they first appear on {1,...,n}
    m <- c(diff(m),size(df,1)-m[length(m)]+1) # m contains their frequencies
    clT[z0] <- clT[z0]+m                 # add the differences in positions to cells of clT corresponding to the unique combinations -> gives the number of times (y,z_{1},...,z_{p00}) appears
    clTdata <- array(clT,dim=c(2,z.max,prod(K00)))  # matrix representation of the tensor clT, with rows of the matrix corresponding to dimension 1 i.e. the levels of y
    clTdata_pop <- matrix(clTdata[1,,],nrow=z.max)
    clTdata_all <- matrix(clTdata[1,,]+clTdata[2,,],nrow=z.max)
    sz <- size(clTdata_all)

    # update lambda (lambda_{isi,g_1,g_2,g_3})
    lambdamat <- matrix(0,z.max,sz[2])
    for(jj in 1:sz[2]) {
      lambdamat[,jj] <- rgamma(z.max,clTdata_pop[,jj]+lambdaalpha*lambda0,1)
      lambdamat[,jj] <- lambdamat[,jj]/sum(lambdamat[,jj])   # normalization
    }
    lambda <- array(lambdamat,dim=c(z.max,K00))  # there is a probability vector of dim z.max for each uniq combination
    for(cc in 1:num_combs){ # store population effect
      exog_val <- as.numeric(all_combs[cc,])
      cls <- G[cbind(1:p,exog_val)][ind00]
      for(zm in 1:z.max){
        lambda_all[t(c(iii,exog_val,zm))] <- lambda[t(c(zm,cls))]
      }
    }

    if (!switch_off_random) {
      # update v (v=0: exgns; v=1: anmls)
      mouse_k_pairs <- cbind(mouseId,z.vec)
      prob_exgns <- piv[mouse_k_pairs]*lambda[cbind(z.vec,z00)]
      prob_anmls <- (1-piv[mouse_k_pairs])*lambda_mice[mouse_k_pairs]
      probs <- prob_exgns/(prob_exgns+prob_anmls)
      probs[is.na(probs)] <- 0
      v <- sapply(probs,function(x) sample(0:1,1,prob=c(x,1-x)))

      # update piv (pi_{isi,0}^{i})
      v0 <- cbind(v+1,z.vec,mouseId)
      v0 <- v0[do.call(order,as.data.frame(v0)),]
      v00 <- cbind(unique(v0),1) # unique combinations of {v,k,mouse_id}
      v0m <- which(duplicated(v0)==FALSE)
      v0m <- c(diff(v0m),size(v0,1)-v0m[length(v0m)]+1) # contain occurrences of each {v,k,id}
      vMat <- array(rep(0,2*z.max*max(mouseId)),dim=c(2,z.max,max(mouseId),1))
      vMat[v00] <- v0m
      for(id in 1:max(mouseId)){
        vMouseMat <- vMat[,,id,]
        piv[id,] <- 1-rbeta(z.max,vMouseMat[1,]+1,vMouseMat[2,]+1)
      }

      # update lambda_mice (lambda^(i))
      for(id in 1:max(mouseId)){
        vMouseCt <- vMat[,,id,][2,]
        lambda_mice[id,] <- rdirichlet(1,lambda_mice_alpha*lambda0+vMouseCt)
        lambda_mice_all[iii,id,] <- lambda_mice[id,]
      }
    }

    # update lambda0 (lambda_{isi,0})
    mmat <- matrix(0,z.max,sz[2])
    for(ii in 1:z.max) {
      for(mm in 1:sz[2]) {
        if(clTdata_all[ii,mm]>0) {
          prob <- lambdaalpha*lambdamat[ii,mm]/(c(1:clTdata_all[ii,mm])-1+lambdaalpha*lambdamat[ii,mm])
          prob[is.na(prob)] <- 10^(-5)    # For numerical reasons
          mmat[ii,mm] <- mmat[ii,mm] + sum(rbinom(length(prob),size=rep(1,clTdata_all[ii,mm]),prob=prob))
        }
      }
    }
    lambda0 <- rgamma(z.max,rowSums(mmat)+lambdaalpha0/z.max,1)
    lambda0[lambda0==0] <- 10^(-5)    # For numerical reasons
    lambda0[is.na(lambda0)] <- 10^(-5)    # For numerical reasons
    lambda0 <- lambda0/sum(lambda0)
    lambda0_all[iii,] <- lambda0

    # update alpha_{isi,0} and alpha_{isi}^{(0)}
    lambdaalpha_all[iii] <- lambdaalpha
    lambdaalpha_mice_all[iii] <- lambda_mice_alpha
    lambdaalpha <- rgamma(1,shape=1+z.max-1,scale=1-digamma(1)+log(num_row))
    lambda_mice_alpha <- rgamma(1,shape=1+z.max-1,scale=1-digamma(1)+log(num_row))

    # update mixed gamma parameters alpha & beta
    # idea: use approximation for conjugate gamma prior
    if (iii > 50) {
      for(kk in 1:z.max) {
        df <- isi[which(z.vec==kk)]
        a_cur <- Alpha[iii-1,kk]
        b_cur <- Beta[iii-1,kk]
        a_new <- approx_gamma_shape(df,a_cur/b_cur,1,1) # mu=shape/rate
        b_new <- rgamma(1,shape=1+length(df)*a_new,rate=1+sum(df))
        Alpha[iii,kk] <- a_new
        Beta[iii,kk] <- b_new
      }
    } else {
      Alpha[iii,] <- init_alpha
      Beta[iii,] <- init_beta
    }

    # update z_{tau,k}
    for(ii in 1:num_row) {
      z.covs.indices <- array(c(1:z.max,rep(z00[ii,],each=z.max)),dim=c(z.max,p00+1)) # ii's values of important covariates
      if(iii > burnin/2) {
        id <- mouseId[ii]
        prob <- (piv[id,]*lambda[z.covs.indices]+(1-piv[id,])*lambda_mice[id,])*dgamma(isi[ii],shape=Alpha[iii,],rate=Beta[iii,])
      } else {
        prob <- lambda0*dgamma(isi[ii],shape=Alpha[iii,],rate=Beta[iii,])
      }
      prob[is.nan(prob)] <- 0
      prob[is.infinite(prob)] <- max(prob[is.finite(prob)])   # Numerical Stability
      if (sum(prob==0)==z.max){
        prob <- rep(1/z.max,z.max)
      } else {
        prob <- prob/sum(prob)
      }
      z.vec[ii] <- sample(z.max,1,TRUE,prob)
    }
  }

  results <- list("K"=K,
                  "Xexgns"=covs,
                  "dpreds"=covs.max,
                  "MCMCparams"=list("simsize"=simsize,"burnin"=burnin,"N_Thin"=5),
                  "Duration_Times"=isi,
                  "Comp_Assign"=z.vec,
                  "Duration_Exgns_Store"=lambda_all,
                  "Marginal_Prob"=prob_mat,
                  "Shape_Samples"=Alpha,
                  "Rate_Samples"=Beta,
                  "Clusters"=M,
                  "Type"="Duration Times")
  return(results)

} # end of function
