# Supporting functions for print_and_plot_results.R

# heat-map function
heatmap_of_matrix <- function(results,m,title) {
  if (results$Type=='Duration Times'){
    xlabel <- ''
    ylabel <- ''
  } else {
    xlabel <- 'y_{t}'
    ylabel <- 'y_{t-1}'
  }
  image.plot(1:ncol(m),1:nrow(m),t(m)[,c(nrow(m):1)],col=terrain.colors(60),
             axes=FALSE,ylab=ylabel,xlab=xlabel,main=title)
  grid(nx = ncol(m), ny = nrow(m), lty = 'solid', col='black')
  axis(1, 1:ncol(m), colnames(m))
  if (results$Type=='Duration Times')
    rownames(m) <- paste0(rep('Comp ',nrow(m)),1:nrow(m))
  axis(2, 1:nrow(m), rownames(m)[nrow(m):1])
  for (x in 1:ncol(m))
    for (y in 1:nrow(m))
      text(x,y,m[nrow(m)+1-y,x])
}

# get transition probabilities by covariate levels
get_trans_mat_by_level <- function(results,cov_levels,level_labels,state_labels,decimal_pts,include_plot) {
  for(row in 1:nrow(cov_levels)) {
    level <- cov_levels[row,]
    label <- level_labels[row,]
    v0 <- expand.grid(1:results$Num_States,1:results$Num_States)
    v0 <- as.matrix(cbind(v0,repmat(level,nrow(v0),1)))
    trans_mat <- t(round(matrix(results$TP_Exgns_Post_Mean[v0],nrow=results$Num_States),decimal_pts))
    sd_mat <- t(round(matrix(results$TP_Exgns_Post_Std[v0],nrow=results$Num_States),decimal_pts))
    colnames(trans_mat) <- state_labels
    rownames(trans_mat) <- state_labels
    colnames(sd_mat) <- state_labels
    rownames(sd_mat) <- state_labels
    title <- paste(c("Covariates {",label,'}'),collapse=" ")
    title_mean <- "Trans. Probs. Posterior Mean"
    title_sd <- "Trans. Probs. Posterior Std"
    print(title)
    print(title_mean)
    print(trans_mat)
    print(title_sd)
    print(sd_mat)
    if(include_plot) {
      #par(mfrow=c(1,2),mar=c(5,5,5,5))
      heatmap_of_matrix(results,trans_mat,paste(c(title,title_mean)))
      heatmap_of_matrix(results,sd_mat,paste(c(title,title_sd)))
    }
  }
}

# local test: get probabilities for H0 for given delta
postprob <- function(x,delta) {
  return(length(x[x<delta])/length(x))
}

# local test: get local test results given covariates levels
get_tp_local_test_results_by_level <- function(results,cov,delta,levels,comp_pairs,cov_labels,state_labels,decimal_pts,include_plot) {
  TP_Exgns_Diffs_Store_Mean <- apply(results$TP_Exgns_Diffs_Store,1:(length(size(results$TP_Exgns_Diffs_Store))-1),mean)
  TP_Exgns_Diffs_Store_Prob <- apply(abs(results$TP_Exgns_Diffs_Store),1:(length(size(results$TP_Exgns_Diffs_Store))-1),function(x) postprob(x,delta))
  all_pairs <- t(combn(results$dpreds[cov],2))
  if (is.null(comp_pairs)) {
    comp_pairs <- all_pairs
  } else if (!is.matrix(comp_pairs)) {
    stop("'comp_pairs' should be a matrix")
  }
  inds <- 1:length(results$dpreds)
  inds <- inds[inds!=cov]
  num_states <- size(results$TP_Exgns_Post_Mean)[1]
  if (ncol(results$Clusters)==1) {
    levels <-matrix(rep(0,2),1)
  }
  for(row in 1:nrow(levels)) {
    level <- as.numeric(levels[row,])
    for(np in 1:nrow(comp_pairs)) {
      covs_levels_1 <- rep(0,length(results$dpreds))
      covs_levels_1[inds] <- level
      covs_levels_1[cov] <- comp_pairs[np,1]
      covs_levels_2 <- covs_levels_1
      covs_levels_2[cov] <- comp_pairs[np,2]
      if(!is.null(cov_labels)) {
        covs_levels_1 <- sapply(1:size(covs_levels_1)[2],function(i) cov_labels[i,covs_levels_1[i]])
        covs_levels_2 <- sapply(1:size(covs_levels_2)[2],function(i) cov_labels[i,covs_levels_2[i]])
      }
      title <- paste(c('Compare covs {',covs_levels_1,'} & {',covs_levels_2,'}'),collapse=" ")
      title_mean <- "Mean Diff. in Trans. Probs."
      title_prob <- "Probability for H_0"
      v0 <- expand.grid(1:num_states,1:num_states)
      pair_ind <- which(apply(all_pairs,1,function(row) all.equal(row,comp_pairs[np,])=='TRUE'))
      if (ncol(results$Clusters)==1) {
        v0 <- as.matrix(cbind(repmat(c(cov,pair_ind),nrow(v0),1),v0))
      } else {
        v0 <- as.matrix(cbind(repmat(c(cov,pair_ind),nrow(v0),1),v0,repmat(level,nrow(v0),1)))
      }
      mean_diff_mat <- abs(matrix(TP_Exgns_Diffs_Store_Mean[v0],nrow=size(results$TP_Exgns_Post_Mean)[1],dimnames=list(state_labels,state_labels)))
      prob_diff_mat <- matrix(TP_Exgns_Diffs_Store_Prob[v0],nrow=size(results$TP_Exgns_Post_Mean)[1],dimnames=list(state_labels,state_labels))
      mean_diff_mat <- t(round(mean_diff_mat,decimal_pts))
      prob_diff_mat <- t(round(prob_diff_mat,decimal_pts))
      print(title)
      print(title_mean)
      print(mean_diff_mat)
      print(title_prob)
      print(prob_diff_mat)
      if(include_plot) {
        #par(mfrow=c(1,2),mar=c(5,5,5,5))
        heatmap_of_matrix(results,mean_diff_mat,paste(c(title,title_mean),collapse="\n"))
        heatmap_of_matrix(results,prob_diff_mat,paste(c(title,title_prob),collapse="\n"))
      }
    }
  }
}

# transition probabilities diagnostic plots by covariates levels
get_tp_diagnostic_plots_by_level <- function(results,from,to,cov_levels,level_labels,state_labels) {
  for(row in 1:nrow(cov_levels)) {
    level <- cov_levels[row,]
    label <- level_labels[row,]
    title <- paste(c("Covariates {",label,'}'),collapse=" ")
    d0 <- size(results$TP_Exgns_All_Itns)[1]
    simsize <- results$MCMCparams$simsize
    burnin <- results$MCMCparams$burnin
    thin_ind <- seq(burnin,simsize,results$MCMCparams$N_Thin)
    for(j in from) {
      for(i in to) {
        ind <- cbind(repmat(c(i,j,level),simsize,1),1:simsize)
        #par(mfrow=c(2,1),mar=c(4,4,4,1))
        plot(results$TP_Exgns_All_Itns[ind],type='l',xlab='Iteration number',col='#00A600',
             ylab=paste(c('Prob(',state_labels[j],'->',state_labels[i],')'),collapse=''),main=paste('Traceplot for',title),ylim=c(0,1))
        abline(v=burnin,col='red',lwd=2)
        abline(h=mean(results$TP_Exgns_All_Itns[ind][burnin:simsize]),col='blue',lty=2,lwd=2)
        acf(results$TP_Exgns_All_Itns[ind][thin_ind],main=paste('Autocorrelation Plot for',title))
      }
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
get_x_axis <- function(results,x_range) {
  return(seq(x_range[1],x_range[2],(x_range[2]-x_range[1])/200))
}

# ISI: get 5%, 50% & 95% quantile for posterior distribution
quantile_calc <- function(results,prob_mat,Alpha,Beta,x_range) {
  x_axis <- get_x_axis(results,x_range)
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

