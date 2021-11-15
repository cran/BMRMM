##################################################
################## Plot Results ##################
##################################################

##########################################
############## Global Tests ##############
##########################################

#' Global Test Results for Transition Probabilities and Duration Times
#'
#' For each covariate, print (and plot) the percentage of the number of clusters for each covariate in stored MCMC iterations.
#'
#' The height of the bar represents the percentage of the number of the clusters within the stored MCMC samples. \cr
#' Pr(# clusters for covariate j ==1) is the probability for the null hypothesis that the covairate j is not significant
#' for transition probabilities or duration times, depending on the result type.
#'
#' @param results either results$results_trans or results$results_duration
#' @param decimal_pts specify the number of decimal points of the results; default is 2
#' @param include_plot display plot if TRUE; default is TRUE
#' @return No return value, called for printing and plotting global test results.
#' 
#' @examples
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=50)
#' get_global_test_results(results$results_trans)
#' get_global_test_results(results$results_duration)
#' 
#' @export
get_global_test_results <- function(results,decimal_pts=2,include_plot=TRUE) {
  simsize <- results$MCMCparams$simsize
  burnin <- results$MCMCparams$burnin
  N_Thin <- results$MCMCparams$N_Thin
  ind <- seq(burnin+1,simsize,N_Thin)
  M <- results$Clusters
  cluster_data <- c()
  label_data <- c()
  for(pp in 1:ncol(M)) {
    cluster_data <- c(cluster_data,M[ind,pp])
    label_data <- c(label_data,rep(colnames(results$Xexgns)[pp],length(ind)))
  }
  counts <- round(table(cluster_data,label_data)/length(ind),decimal_pts)
  print(counts)
  if(include_plot) {
    #par(mfrow=c(1,1),mar=c(5,5,5,5))
    barplot(counts, main=paste("Global Test Results for",results$Type),ylab="Proportions",
            xlab="",col=terrain.colors(nrow(counts)),ylim=c(0,1),legend.text=TRUE,
            args.legend=list(x="topright",legend=rownames(counts),inset=c(-0.2,0),xpd=TRUE),
            beside=TRUE)
  }
}


##################################################
################ Transition Probs ################
##################################################


#' Transition Probabilities: Posterior Mean and Standard Deviation
#'
#' Print and plot the posterior mean and standard deviation for
#' transition probabilities from MCMC samples under given different
#' combinations of covariate levels.
#'
#' For each row of 'cov_levels', the function returns two matrices of size d0xd0 where d0 is the number of states:
#' (1) the posterior mean and (2) the posterior standard deviation
#' of transition probabilities, computed from MCMC samples after burn-ins and thinning.
#' The default for 'cov_levels' is all possible combinations of covariate levels.
#'
#' @param results results of transition probabilities, i.e., results$results_trans
#' @param cov_labels a matrix such that row i represents the labels for covariate i; default labels for covariate i is 1:i
#' @param state_labels a vector of strings that represent the state labels; default is 1:Num_States
#' @param cov_levels a matrix such that each row is a combination of covariate levels; default is all possible combinations of covariates
#' @param decimal_pts specify the number of decimal points of the results; default is 2
#' @param include_plot display plot if TRUE; default is TRUE
#' 
#' @return No return value, called for printing and plotting posterior distribution of transition probabilities.
#' 
#' @examples
#'
#' # Examples using the shortened version of the simulated Foxp2 data set, foxp2_sm
#'
#' # get results for all combinations of covariate levels
#' results <- BMRMM(foxp2_sm,num_cov=2,duration_type='None',simsize=50)
#' get_estimated_post_mean_and_sd(results$results_trans)
#'
#' # get results for covariate levels ("HET","U") and ("WT","U")
#' cov_labels <- matrix(c("HET","WT","","U","L","A"),nrow=2,byrow=TRUE)
#' cov_levels <- matrix(c(1,1,2,1),nrow=2,byrow=TRUE)
#' get_estimated_post_mean_and_sd(results$results_trans,cov_labels,cov_levels=cov_levels)
#'
#' @export
get_estimated_post_mean_and_sd <- function(results,cov_labels=NULL,state_labels=1:results$Num_States,cov_levels=NULL,decimal_pts=2,include_plot=TRUE) {
  if(is.null(cov_levels)) {
    cov_levels <- as.matrix(unique(results$Xexgns))
  } else if (!is.matrix(cov_levels)) {
    stop("'cov_levels' should be a matrix")
  }
  if(size(cov_levels)[1]==1)
    cov_levels <- matrix(cov_levels,nrow=1)
  if(is.null(cov_labels)) {
    level_labels <- cov_levels
  } else {
    level_labels <- sapply(1:ncol(cov_levels),function(i) cov_labels[i,cov_levels[,i]])
  }
  if(size(level_labels)[1]==1)
    level_labels <- matrix(level_labels,nrow=1)
  get_trans_mat_by_level(results,cov_levels,level_labels,state_labels,decimal_pts,include_plot)
}


#' Transition Probabilities: Local Tests for a Specific Covariate
#'
#' Given some covariate 'cov' and the levels of other covariates,
#' prints and plots the local test results for pairs of levels of 'cov',
#' including the absolute differences in transition probabilities,
#' and the probabilities for the null hypothesis.
#'
#' Fix a covariate i as the target of the local test. \cr
#' The function provides two matrices of size d0xd0 where d0 is the number of states: \cr
#' (1) the mean of the absolute difference in all transition types for each pair of levels of covariate i; \cr
#' (2) the probability of the null hypothesis of the local test for all transition types.
#'
#' Given a pair of levels of covariate i, say i_1 and i_2, and given the levels of other covariates,
#' the null hypothesis is that the difference between i_1 and i_2 is not significant for transition probabilities.
#' It is calculated as the percentage of the samples with absolute difference less than 'delta'.
#'
#' 'comp_pairs' are user-specified pairs of 'cov' for the local test.
#' Each row of 'comp_pairs' is a pair of indices of 'cov'.
#'
#' @param results results of transition probabilities, i.e., results$results_trans
#' @param cov the index of the covariate to perform local tests on
#' @param delta threshold for the null hypothesis H_0; see Details
#' @param all_cov_labels a matrix such that row i represents the labels for covariate i; default labels for covariate i is 1:i
#' @param state_labels a vector of strings that represent the state labels; default is 1:Num_States
#' @param other_cov_levels a matrix such that each row is a combination of covariate levels excluding 'cov'; default is all possible combinations of covariates excluding 'cov'
#' @param comp_pairs a matrix such that each row is a pair of levels of 'cov' to compare
#' @param decimal_pts specify the number of decimal points of the results; default is 2
#' @param include_plot display plot if TRUE; default is TRUE
#' 
#' @return No return value, called for printing and plotting local test results for transition probabilities.
#' 
#' @examples
#'
#' # Examples using the shortened built-in simulated Foxp2 data set, foxp2_sm
#' # recall the covariate labels are: ("HET","WT") and ("U","L","A")
#' 
#' results <- BMRMM(foxp2_sm,num_cov=2,duration_type='None',simsize=50)
#'
#' # results for genotypes (HET, WT) under all three social contexts
#' get_tp_local_test_results(results$results_trans,cov=1,delta=0.02)
#'
#' # results for genotypes (HET, WT) under social contexts U and A
#' get_tp_local_test_results(results$results_trans,cov=1,delta=0.02,
#'                           other_cov_levels=matrix(c(1,3),nrow=2))
#'
#' # results for contexts (U,L), (U,A) and (L,A) under two genotypes
#' get_tp_local_test_results(results$results_trans,cov=2,delta=0.02)
#'
#' # results for contexts (U,L), (U,A) under genotype HET (specify pairs)
#' get_tp_local_test_results(results$results_trans,cov=2,delta=0.02,
#'                           other_cov_levels=matrix(1,nrow=1),
#'                           comp_pairs=matrix(c(1,2,1,3),nrow=2,byrow=2))
#'
#' @export
get_tp_local_test_results <- function(results,cov,delta,all_cov_labels=NULL,state_labels=1:results$Num_States,other_cov_levels=NULL,comp_pairs=NULL,decimal_pts=2,include_plot=TRUE) {
  if(is.null(other_cov_levels)){
    inds <- 1:length(results$dpreds)
    inds <- inds[inds!=cov]
    other_cov_levels <- as.matrix(unique(results$Xexgns[,inds]))
  } else if (!is.matrix(other_cov_levels)) {
    stop("'other_cov_levels' should be a matrix")
  }
  if(!is.matrix(other_cov_levels) && size(other_cov_levels)[1]==1) {
    other_cov_levels <- matrix(other_cov_levels,nrow=size(other_cov_levels)[2])
  }
  get_tp_local_test_results_by_level(results,cov,delta,other_cov_levels,comp_pairs,all_cov_labels,state_labels,decimal_pts,include_plot)
}


#' Transition Probabilities: MCMC Diagnostic Plots
#'
#' Provides the traceplots and autocorrelation plots for each transition type under combinations of covariate levels
#'
#' @param results results of transition probabilities, i.e., results$results_trans
#' @param from list of "from" states in target state transitions; default is all states
#' @param to list of "to" states in target state transitions; default is all states
#' @param cov_levels a matrix such that each row is a combination of covariate levels; default is all possible combinations of covariates
#' @param cov_labels a matrix such that row i represents the labels for covariate i; default labels for covariate i is 1:i
#' @param state_labels a vector of strings that represent the state labels; default is 1:Num_States
#' 
#' @return No return value, called for plotting MCMC diagnostic plots for transition probabilities.
#' 
#' @examples
#'
#' # Examples using the shortened built-in simulated Foxp2 data set, foxp2_sm
#' # recall the covariate labels are: ("HET","WT") and ("U","L","A")
#' # recall the state labels are 'd', 'm', 's', 'u'
#' 
#' results <- BMRMM(foxp2_sm,num_cov=2,duration_type='None',simsize=50)
#'
#' # results for all transition types for every combinations of covariates
#' get_tp_diagnostic_plots(results$results_trans)
#'
#' # results for all transition types for covariate levels (HET,U) and (WT,A)
#' get_tp_diagnostic_plots(results$results_trans,cov_levels=matrix(c(1,1,2,3),nrow=2,byrow=TRUE))
#'
#' # results for all transition from (d,m) to (d,s) for covariate levels (HET,U) and (WT,A)
#' get_tp_diagnostic_plots(results$results_trans,from=c(1),to=c(2,3),
#'                         cov_levels=matrix(c(1,1,2,3),nrow=2,byrow=TRUE))
#'
#' @export
get_tp_diagnostic_plots <- function(results,from=1:results$Num_States,to=1:results$Num_States,cov_levels=NULL,cov_labels=NULL,state_labels=1:results$Num_States) {
  if(is.null(cov_levels)) {
    cov_levels <- sortrows(as.matrix(unique(results$Xexgns)))
  } else if (!is.matrix(cov_levels)) {
    stop("'cov_levels' should be a matrix")
  }
  if(!is.matrix(cov_levels) && size(cov_levels)[1]==1)
    cov_levels <- matrix(cov_levels,nrow=size(cov_levels)[2])
  if(is.null(cov_labels)) {
    level_labels <- cov_levels
  } else {
    level_labels <- sapply(1:size(cov_levels)[2],function(i) cov_labels[i,cov_levels[,i]])
  }
  if(size(level_labels)[1]==1)
    level_labels <- matrix(level_labels,nrow=1)
  get_tp_diagnostic_plots_by_level(results,from,to,cov_levels,level_labels,state_labels)
}


##################################################
############### Duration Times ###################
##################################################

#' Duration Times Plots: Histogram of duration times with Estimated Posterior Mean
#'
#' Plots the histogram of Duration Times superimposed the posterior mean mixture gamma distribution
#'
#' @param results results of duration times, i.e., results_duration
#' @param x_range custom range for x-axis; default is (min(durations),max(durations))
#' @param breaks custom breaks for the histogram; default breaks follows the Freedman-Diaconis rule, which would be (max(durations)-min(durations))*n^(1/3)/2/IQR(durations)
#' @examples
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=50)
#' get_histogram_with_post_mean(results$results_duration,x_range=c(0,2),breaks=50)
#' @return No return value, called for plotting histogram superimposed with the posterior distribution of mixture gamma.
#'
#' @export
get_histogram_with_post_mean <- function(results,
                                         x_range=c(min(results$Duration_Times),max(results$Duration_Times)),
                                         breaks=(max(results$Duration_Times)-min(results$Duration_Times))*length(results$Duration_Times)^(1/3)/2/IQR(results$Duration_Times)) {
  if(results$Type != 'Duration Times')
    stop("'results' must be duration times results")
  x_axis <- get_x_axis(results,x_range)
  quan <- quantile_calc(results,results$Marginal_Prob,results$Shape_Samples,results$Rate_Samples,x_range)
  #par(mfrow=c(1,1),mar=c(5,5,5,5))
  hist(results$Duration_Times,freq=FALSE,breaks=breaks,xlim=x_range,col='gray',
       main='Histogram with Posterior Mean',xlab='Duration times',ylab='Density',
       cex.lab=1.5,cex.axis=1.5,cex.main=1.2)
  polygon(c(x_axis,rev(x_axis)),c(quan$p95,rev(quan$p5)), # draw credible region
          col=rgb(190,190,190,130,maxColorValue=255), border=NA)
  lines(x=x_axis,y=quan$p50,col="red",lwd=2) # superimpose posterior mean
}


#' Duration Times Plots: Histogram of duration times by Component
#'
#' Plots the histogram of each mixture component superimposed the gamma distribution with shape and rate
#' parameters taken from the last MCMC iteration
#'
#' @param results results of duration times, i.e., results$results_duration
#' @param comp a vector of the components to plot; default is 1:isi_num_comp
#' 
#' @return No return value, called for plotting histogram for each mixture component.
#'
#' @examples
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=50)
#' get_histogram_by_component(results$results_duration)
#' 
#' @export
get_histogram_by_component <- function(results,comp=1:results$K) {
  if(results$Type != 'Duration Times')
    stop("'results' must be duration times results")
  #par(mfrow=c(ceil(length(comp)/2),2),mar=c(4,5,4,4))
  for (kk in comp) {
    dat <- results$Duration_Times[results$Comp_Assign==kk]
    if (length(dat)==0) {
      print(paste("Component",kk,'is empty'))
    } else {
      x_range <- c(min(dat),max(dat))
      x_axis <- get_x_axis(results,x_range)
      tle <- paste("Component",kk)
      hist(dat,freq=FALSE,xlim=x_range,
           main=tle,col='gray',border=FALSE,xlab='Duration times',ylab='Density',
           cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
      lines(x=x_axis,y=dgamma(x_axis,shape=results$Shape_Samples[nrow(results$Shape_Samples),kk],
                              rate=results$Rate_Samples[nrow(results$Shape_Samples),kk]),col="red",lwd=2)
    }
  }
}


#' Duration Times Plot: Heat-map of Mixture Gamma Shape & Rate Parameters
#'
#' Print and plot the heat-map for the gamma mixture parameters for each component
#'
#' Print a matrix of size K*2, where K is the number of components.
#' Each row i represents the shape and rate parameter of component i.
#'
#' @param results results of duration times, i.e., results$results_duration
#' @param decimal_pts specify the number of decimal points of the results; default is 2
#' @param include_plot display plot if TRUE; default is TRUE
#' 
#' @return No return value, called for printing and plotting heatmaps for mixture parameters.
#' 
#' @examples
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=50)
#' get_heatmap_mix_param(results$results_duration)
#' 
#' @export
get_heatmap_mix_param <- function(results,decimal_pts=2,include_plot=TRUE) {
  if(results$Type != 'Duration Times')
    stop("'results' must be duration times results")
  simsize <- results$MCMCparams$simsize
  m <- cbind(shape.k=results$Shape_Samples[simsize,],rate.k=results$Rate_Samples[simsize,])
  m <- round(m,decimal_pts)
  title <- 'Heatmap for Mixture Probabilities'
  cat(title,"\n")
  print(m)
  if(include_plot) {
    #par(mfrow=c(1,1),mar=c(4,4,4,4))
    heatmap_of_matrix(results,m,title)
  }
}


#' Duration Times Plot: Heat-map of Mixture Probabilities for Different Covariates
#'
#' Print and plot the heat-map of mixture probabilities for covariate levels
#'
#' Print a matrix for each covariate.
#' Each matrix has K rows, where K is the number of components.
#' The number of columns of a matrix equals to the total number of levels of the covariate.
#'
#' @param results results of duration times, i.e., results$results_duration
#' @param cov_labels a matrix such that row i represents the labels for covariate i; default labels for covariate i is 1:i
#' @param cov_index a vector of indices of covariate to plot; default is 1:Num_Covariates
#' @param decimal_pts specify the number of decimal points of the results; default is 2
#' @param include_plot display plot if TRUE; default is TRUE
#' @examples
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=50)
#' get_heatmap_by_cov(results$results_duration) 
#' get_heatmap_by_cov(results$results_duration,cov_index=c(1,2))
#' @return No return value, called for printing and plotting heatmaps for mixture probabilities for each covariate.
#'
#'
#' @export
get_heatmap_by_cov <- function(results,cov_labels=NULL,cov_index=1:ncol(results$Xexgns),decimal_pts=2,include_plot=TRUE) {
  if(results$Type != 'Duration Times')
    stop("'results' must be duration times results")
  ind_arr <- list()
  for (i in 1:length(results$dpreds))
    ind_arr <- c(ind_arr,list(1:results$dpreds[i]))
  ind_mat <- do.call(expand.grid,ind_arr)
  ind_mat <- cbind(rep(results$MCMCparams$simsize,nrow(ind_mat)),ind_mat)
  #par(mfrow=c(1,length(cov_index)),mar=c(5,5,5,5))
  for (pp in cov_index) {
    prob_to_plot <- matrix(0,nrow=results$K,ncol=results$dpreds[pp])
    for (cc in 1:results$K) {
      select_ind <- as.matrix(cbind(ind_mat,rep(cc,nrow(ind_mat))))
      dat <- array(results$Duration_Exgns_Store[select_ind],dim=results$dpreds)
      prob_to_plot[cc,] <- apply(dat,c(pp),mean)
    }
    colnames(prob_to_plot) <- cov_labels[pp,1:results$dpreds[pp]]
    title <- paste("Heatmap for Covariate",colnames(results$Xexgns)[pp])
    prob_to_plot <- round(prob_to_plot,decimal_pts)
    cat(title,"\n")
    print(prob_to_plot)
    if(include_plot) {
      heatmap_of_matrix(results,prob_to_plot,title)
      #count <- t(prob_to_plot)
      #colnames(count) <- paste0(rep('Comp ',ncol(count)),1:ncol(count))
      #barplot(count, main="",col=terrain.colors(nrow(count)),ylab="Mixture Probability",
      #        legend=rownames(count),beside=TRUE,ylim=c(0,1))
    }
  }
}


#' Duration Times Plots: MCMC Diagnostic Plots
#'
#' Provides the traceplots and autocorrelation plots for each gamma shape and rate parameter
#'
#' @param results results of duration times, i.e., results$results_duration
#' @param comps a vector of the components to plot; default is 1:duration_num_comp
#' 
#' @return No return value, called for plotting MCMC diagnostic plots for duration times.
#'
#' @examples
#'
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=80)
#' 
#' # get diagnostic plots for all components
#' get_duration_diagnostic_plots(results$results_duration)
#'
#' # get diagnostic plots for components 1,2
#' get_duration_diagnostic_plots(results$results_duration,comps=c(1,2))
#'
#' @export
get_duration_diagnostic_plots <- function(results,comps=1:results$K) {
  if(results$Type != 'Duration Times')
    stop("'results' must be duration times results")
  burnin <- results$MCMCparams$burnin
  simsize <- results$MCMCparams$simsize
  af.thin <- seq(burnin+1,simsize,5)
  #par(mfrow=c(2,2),mar=c(3,3,3,3))
  for(kk in comps) {
    plot(results$Shape_Samples[,kk],type='l',xlab='',ylab='',col='#00A600',main=paste("Traceplot for Shape, Comp",kk))
    abline(v=burnin,col='red',lwd=2)
    abline(h=mean(results$Shape_Samples[burnin:simsize,kk]),col='blue',lty=2,lwd=2)
    plot(results$Rate_Samples[,kk],type='l',xlab='',ylab='',col='#00A600',main=paste("Traceplot for Rate, Comp",kk))
    abline(v=burnin,col='red',lwd=2)
    abline(h=mean(results$Rate_Samples[burnin:simsize,kk]),col='blue',lty=2,lwd=2)
    acf(results$Shape_Samples[af.thin,kk],main=paste('Autocorrelation after thin. for Shape, Comp',kk))
    acf(results$Rate_Samples[af.thin,kk],main=paste('Autocorrelation after thin. for Rate, Comp',kk))
  }
}

#' Duration Times Result: Model Selection Scores for the Number of Components
#'
#' Provides the LPML (Geisser and Eddy, 1979) and WAIC (Watanabe, 2010) scores of the Bayesian Markov renewal mixture models
#'
#' The two scores can be used to compare different choices of isi_num_comp, i.e., the number
#' of the mixture gamma components. Larger values of LPML and smaller values of WAIC
#' indicate better model fits.
#' 
#' @examples
#'
#' results <- BMRMM(foxp2_sm,num_cov=2,simsize=50)
#' get_lpml_waic_scores(results$results_duration)
#'
#' @references
#' Geisser, S. and Eddy, W. F. (1979). A predictive approach to model selection. Journal of the American Statistical Association, 74, 153–160. \cr\cr
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11, 3571–3594.
#'
#' @param results results of duration times, i.e., results$results_duration
#' 
#' @return No return value, called for printing LPML and WAIC scores for gamma mixture models.
#' 
#' @export
get_lpml_waic_scores <- function(results) {
  if(results$Type != 'Duration Times')
    stop("'results' must be duration times results")
  isi <- results$Duration_Times
  pis <- results$Marginal_Prob
  alphas <- results$Shape_Samples
  betas <- results$Rate_Samples
  lpml <- 0
  waic <- 0
  p_waic <- 0
  for(k in 1:length(isi)) {
    expectation <- 0
    inv_expectation <- 0
    log_expectation <- 0
    for(j in 1:nrow(alphas)) {
      expectation <- expectation + mix_gam(isi[k],pis[j,],alphas[j,],betas[j,])
      inv_expectation <- inv_expectation + 1/mix_gam(isi[k],pis[j,],alphas[j,],betas[j,])
      log_expectation <- log_expectation + log(mix_gam(isi[k],pis[j,],alphas[j,],betas[j,]))
    }
    expectation <- expectation/nrow(alphas)
    inv_expectation <- inv_expectation/nrow(alphas)
    log_expectation <- log_expectation/nrow(alphas)
    lpml <- lpml + log(1/inv_expectation)
    waic <- waic - 2*log(expectation)
    p_waic <- p_waic + 2*(log(expectation)-log_expectation)
  }
  print(list('LPML'=lpml,'WAIC'=waic+2*p_waic))
}


