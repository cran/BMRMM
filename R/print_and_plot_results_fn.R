# Supporting functions for print_and_plot_results.R

# global test 
global_test <- function(results) {
  simsize <- results$MCMCparams$simsize
  burnin <- results$MCMCparams$burnin
  N_Thin <- results$MCMCparams$N_Thin
  ind <- seq(burnin+1,simsize,N_Thin)
  M <- results$clusters
  cluster_data <- c()
  label_data <- c()
  for(pp in 1:ncol(M)) {
    cluster_data <- c(cluster_data,M[ind,pp])
    label_data <- c(label_data,rep(colnames(results$covs)[pp],length(ind)))
  }
  counts <- table(cluster_data,label_data)/length(ind)
  return(counts)
}

# mixture parameters
get_mixture_params <- function(results,digits) {
  simsize <- results$MCMCparams$simsize
  m <- cbind(shape.k=results$shape.samples[simsize,],rate.k=results$rate.samples[simsize,])
  m <- round(m,digits)
  rownames(m) <- paste0(rep('Comp ',nrow(m)),1:nrow(m))
  return(m)
}

# mixture probabilities
get_mixture_probs_by_cov <- function(results,digits) {
  mix_probs_res <- list()
  ind_arr <- list()
  for (i in 1:length(results$dpreds))
    ind_arr <- c(ind_arr,list(1:results$dpreds[i]))
  ind_mat <- do.call(expand.grid,ind_arr)
  ind_mat <- cbind(rep(results$MCMCparams$simsize,nrow(ind_mat)),ind_mat)
  for (pp in 1:ncol(results$covs)) {
    prob_to_plot <- matrix(0,nrow=size(results$shape.samples,2),ncol=results$dpreds[pp])
    for (cc in 1:size(results$shape.samples,2)) {
      select_ind <- as.matrix(cbind(ind_mat,rep(cc,nrow(ind_mat))))
      dat <- array(results$duration.exgns.store[select_ind],dim=results$dpreds)
      prob_to_plot[cc,] <- apply(dat,c(pp),mean)
    }
    colnames(prob_to_plot) <- results$cov.labels[[pp]]
    prob_to_plot <- round(prob_to_plot, digits)
    rownames(prob_to_plot) <- paste0(rep('Comp ',nrow(prob_to_plot)),1:nrow(prob_to_plot))
    mix_probs_res[[length(mix_probs_res)+1]] <- prob_to_plot
  }
  names(mix_probs_res) <- colnames(results$covs)
  return(mix_probs_res)
}

# get transition probabilities by covariate levels
get_trans_mat_by_level <- function(results,cov_combs,digits) {
  mat_mean <- list()
  mat_sd <- list()
  mat_lab <- list()
  for(row in cov_combs) {
    label <- paste(lapply(1:length(row), function(x) results$cov.labels[[x]][row[x]]),collapse=',')
    label <- paste('(',label,')',sep='')
    v0 <- expand.grid(1:length(results$state.labels),1:length(results$state.labels))
    v0 <- as.matrix(cbind(v0,repmat(row,nrow(v0),1)))
    trans_mat <- t(round(matrix(results$tp.exgns.post.mean[v0],nrow=length(results$state.labels)),digits))
    sd_mat <- t(round(matrix(results$tp.exgns.post.std[v0],nrow=length(results$state.labels)),digits))
    colnames(trans_mat) <- results$state.labels
    rownames(trans_mat) <- results$state.labels
    colnames(sd_mat) <- results$state.labels
    rownames(sd_mat) <- results$state.labels
    mat_mean[[length(mat_mean)+1]] <- trans_mat
    mat_sd[[length(mat_sd)+1]] <- sd_mat
    mat_lab[[length(mat_lab)+1]] <- label
  }
  names(mat_mean) <- mat_lab
  names(mat_sd) <- mat_lab
  return(list(mat_mean,mat_sd))
}

# local test: get probabilities for H0 for given delta
postprob <- function(x,delta) {
  return(length(x[x<delta])/length(x))
}

# local test: get local test results given covariates levels
get_tp_local_test_results_by_level <- function(results,cov,delta,levels,comp_pairs,digits) {
  tp.exgns.diffs.store_Mean <- apply(results$tp.exgns.diffs.store,1:(length(size(results$tp.exgns.diffs.store))-1),mean)
  tp.exgns.diffs.store_Prob <- apply(abs(results$tp.exgns.diffs.store),1:(length(size(results$tp.exgns.diffs.store))-1),function(x) postprob(x,delta))
  inds <- 1:length(results$dpreds)
  inds <- inds[inds!=cov]
  if (ncol(results$clusters)==1) {
    levels <- matrix(rep(0,2),1)
  }
  res_by_levels_mean.diff <- vector("list",nrow(levels))
  res_by_levels_null.test <- vector("list",nrow(levels))
  for(row in 1:nrow(levels)) {
    level <- as.numeric(levels[row,])
    res_by_comp.pairs_mean.diff <- vector("list",nrow(comp_pairs))
    res_by_comp.pairs_null.test <- vector("list",nrow(comp_pairs))
    for(np in 1:nrow(comp_pairs)) {
      covs_levels_1 <- rep(0,length(results$dpreds))
      covs_levels_1[inds] <- level
      covs_levels_1[cov] <- comp_pairs[np,1]
      covs_levels_2 <- covs_levels_1
      covs_levels_2[cov] <- comp_pairs[np,2]
      covs_levels_1 <- sapply(1:size(covs_levels_1)[2],function(i) results$cov.labels[[i]][covs_levels_1[i]])
      covs_levels_2 <- sapply(1:size(covs_levels_2)[2],function(i) results$cov.labels[[i]][covs_levels_2[i]])
      v0 <- expand.grid(1:length(results$state.labels),1:length(results$state.labels))
      pair_ind <- which(apply(comp_pairs,1,function(row) all.equal(row,comp_pairs[np,])=='TRUE'))
      if (ncol(results$clusters)==1) {
        v0 <- as.matrix(cbind(repmat(c(cov,pair_ind),nrow(v0),1),v0))
      } else {
        v0 <- as.matrix(cbind(repmat(c(cov,pair_ind),nrow(v0),1),v0,repmat(level,nrow(v0),1)))
      }
      mean_diff_mat <- abs(matrix(tp.exgns.diffs.store_Mean[v0],nrow=size(results$tp.exgns.post.mean)[1],dimnames=list(results$state.labels,results$state.labels)))
      prob_diff_mat <- matrix(tp.exgns.diffs.store_Prob[v0],nrow=size(results$tp.exgns.post.mean)[1],dimnames=list(results$state.labels,results$state.labels))
      mean_diff_mat <- t(round(mean_diff_mat,digits))
      prob_diff_mat <- t(round(prob_diff_mat,digits))
      res_by_comp.pairs_mean.diff[[np]] <- mean_diff_mat
      res_by_comp.pairs_null.test[[np]] <- prob_diff_mat
      comp_pair_name <- sapply(comp_pairs[np,],function(x) results$cov.labels[[cov]][x])
      comp_pair_name <- paste('Compare:',paste(comp_pair_name,collapse='&'),sep='')
      names(res_by_comp.pairs_mean.diff)[[np]] <- comp_pair_name
      names(res_by_comp.pairs_null.test)[[np]] <- comp_pair_name
    }
    res_by_levels_mean.diff[[row]] <- res_by_comp.pairs_mean.diff
    res_by_levels_null.test[[row]] <- res_by_comp.pairs_null.test
    level_name <- 'Fixing:'
    if (ncol(results$clusters) > 1) {
      other_level_name <- sapply(1:length(inds),function(x) results$cov.labels[[inds[x]]][levels[row,x]])
      level_name <- paste('Fixing:',paste(other_level_name,collapse=' '),sep='')
    }
    names(res_by_levels_mean.diff)[[row]] <- level_name
    names(res_by_levels_null.test)[[row]] <- level_name
  }
  return(list('mean.diff'=res_by_levels_mean.diff,'null.test'=res_by_levels_null.test))
}

# transition probabilities diagnostic plots by covariates levels
get_tp_diagnostic_plots_by_level <- function(results,transitions,cov_combs) {
  for(row in cov_combs) {
    label <- unlist(lapply(1:length(row), function(x) results$cov.labels[[x]][row[x]]))
    title <- paste(c("Covariates {",label,'}'),collapse=" ")
    d0 <- size(results$tp.exgns.all.itns)[1]
    simsize <- results$MCMCparams$simsize
    burnin <- results$MCMCparams$burnin
    thin_ind <- seq(burnin,simsize,results$MCMCparams$N_Thin)
    for(transition in transitions) {
      i <- transition[1]
      j <- transition[2]
      ind <- cbind(repmat(c(i,j,row),simsize,1),1:simsize)
      plot(results$tp.exgns.all.itns[ind],type='l',xlab='Iteration number',col='#00A600',
            ylab=paste(c('Prob(',results$state.labels[j],'->',results$state.labels[i],')'),collapse=''),main=paste('Traceplot for',title),ylim=c(0,1))
      abline(v=burnin,col='red',lwd=2)
      abline(h=mean(results$tp.exgns.all.itns[ind][burnin:simsize]),col='blue',lty=2,lwd=2)
      acf(results$tp.exgns.all.itns[ind][thin_ind],main=paste('Autocorrelation Plot for',title))
    }
  }
}

# mixture gamma distributions
mix_gam <- function(x,pi,alpha,beta) {
  ret <- 0
  for(i in 1:length(pi)) {
    ret <- ret + pi[i]*dgamma(x,shape=alpha[i],rate=beta[i])
  }
  return(ret)
} # mixture gamma distribution

# get discretized x-axis
get_x_axis <- function(x_range) {
  return(seq(x_range[1],x_range[2],(x_range[2]-x_range[1])/200))
}

# ISI: get 5%, 50% & 95% quantile for posterior distribution
quantile_calc <- function(results,prob_mat,Alpha,Beta,x_range) {
  x_axis <- get_x_axis(x_range)
  ind <- seq(results$MCMCparams$burnin+1,results$MCMCparams$simsize,5)
  a_s <- Alpha[ind,]
  pi_s <- prob_mat[ind,]
  pi_s[pi_s==0] <- 10^-5
  b_s <- Beta[ind,]
  percent5 <- rep(0,length(x_axis))
  percent50 <- rep(0,length(x_axis))
  percent95 <- rep(0,length(x_axis))
  for(i in 1:length(x_axis)) {
    val <- sapply(1:length(ind),function(x) mix_gam(x_axis[i],pi_s[x,],a_s[x,],b_s[x,]))
    res <- quantile(val,c(.05,.5,.95))
    percent5[i] <- res['5%']
    percent50[i] <- res['50%']
    percent95[i] <- res['95%']
  }
  return(list('p5'=percent5,'p50'=percent50,'p95'=percent95))
}

# heatmap function
plot.heatmap <- function(x,xlab=NULL,ylab=NULL,main=NULL,...) {
  if(!is.matrix(x)) {
    stop("'x' must be of class 'matrix'")
  }
  mat <- x
  image.plot(1:ncol(mat),1:nrow(mat),t(mat)[,c(nrow(mat):1)],col=terrain.colors(60),
             axes=FALSE,ylab=ylab,xlab=xlab,main=main,...)
  grid(nx = ncol(mat), ny = nrow(mat), lty = 'solid', col='black')
  axis(1, 1:ncol(mat), colnames(mat))
  axis(2, 1:nrow(mat), rownames(mat)[nrow(mat):1])
  for (x in 1:ncol(mat))
    for (y in 1:nrow(mat))
      text(x,y,mat[nrow(mat)+1-y,x])
}

