##################################################
################## Plot Results ##################
##################################################

############################################
############## Summary Method ##############
############################################

#' Summary Method for Objects of Class `BMRMM`
#'
#' Summarizing an object of class `BMRMM`, including results for transition probabilities and duration times, if applicable. 
#'
#' @param object an object of class `BMRMM`.
#' @param delta threshold for the null hypothesis for the local tests of transition probabilities (see Details). Default is 0.02.
#' @param digits integer used for number formatting. Default is 2.
#' @param ... further arguments for the summary function.
#' 
#' @return An object of class `BMRMMsummary` with the following elements: \tabular{ll}{
#' `trans.global` \tab global test results for transition probabilities (see Details). \cr
#' `trans.probs.mean` \tab mean for the posterior transition probabilities. \cr
#' `trans.probs.sd` \tab standard deviation for the posterior transition probabilities. \cr
#' `trans.local.mean.diff` \tab the absolute difference in transition probabilities for a pair of covariate levels (see Details). \cr
#' `trans.local.null.test` \tab probability for the null hypothesis that the difference between two covariate levels is not significant (see Details). \cr
#' `dur.global` \tab global test results for duration times (see Details). \cr
#' `dur.mix.params` \tab mixture parameters taken from the last MCMC iteration if duration times follow a mixture gamma distribution. \cr
#' `dur.mix.probs`\tab mixture probabilities for each covariate taken from the last MCMC iteration if duration times follow a mixture gamma distribution. \cr
#' }
#' 
#' @details We give more explanation for the global tests and local tests results.
#' \itemize{
#' \item{Global tests (for both transition probabilities and duration times)}
#' 
#' {Global tests are presented as a matrix, where the row denote the number of clusters and the column represents covariates.
#' For each row `i` and column `j`, the matrix entry is the percentage of the number of the clusters within the stored MCMC samples 
#' for this covariate, i.e., an estimation for `Pr(# clusters for covariate j == i)`. We note that the probability 
#' `Pr(# clusters for covariate j > 1)` would be the probability for the null hypothesis that the covariate `j` is significant.}
#' \item{Local tests (for transition probabilities only)}
#' 
#' {Local tests focus on a particular covariate and compare the influence among its levels when the other covariates values are fixed. \cr
#' Given a pair of levels of covariate `j`, say `j_1` and `j_2`, and given the levels of other covariates,
#' the null hypothesis is that the difference between `j_1` and `j_2` is not significant for transition probabilities.
#' It is calculated as the percentage of the samples with absolute difference less than `delta`.
#' 
#' The local tests provide two matrices of size `d0` x `d0` where `d0` is the number of states: \cr
#' \enumerate{
#' \item `mean.diff` -- the mean of the absolute difference in each transition type between levels `j_1` and `j_2`; 
#' \item `null.test` -- the probability of the null hypothesis that `j_1` and `j_2` have the same significance for each transition type.
#' }}
#' }
#' 
#' @seealso 
#' [BMRMM::plot.BMRMMsummary] for plotting the summary results.
#' 
#' @examples
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50, 
#'                  cov.labels = list(c("F", "W"), c("U", "L", "A")),
#'                  duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
#' sm <- summary(results)
#' sm
#' 
#' @export
summary.BMRMM <- function(object, delta=0.02, digits=2, ...) {
  if (!inherits(object, "BMRMM")) {
    stop("'object' must be of class 'BMRMM'")
  }
  summary_list <- list()
  if ("results.trans" %in% names(object)) {
    # transition probabilities posterior mean & sd
    summary_list$trans.global <- global_test(object$results.trans) 
    cov.combs <- sortrows((unique(object$results.trans$covs)))
    cov.combs <- as.list(data.frame(t(cov.combs)))
    trans.probs.res <- get_trans_mat_by_level(object$results.trans,cov.combs,digits)
    summary_list$trans.probs.mean <- trans.probs.res[[1]]
    summary_list$trans.probs.sd <- trans.probs.res[[2]]
    # transition probabilities local test - mean diff & null test
    trans.local.mean.diff <- vector("list",size(object$results.trans$covs,2))
    trans.local.null.test <- vector("list",size(object$results.trans$covs,2))
    names(trans.local.mean.diff) <- colnames(object$results.trans$covs)
    names(trans.local.null.test) <- colnames(object$results.trans$covs)
    for(pp in 1:length(trans.local.mean.diff)) {
      inds <- 1:length(object$results.trans$dpreds)
      inds <- inds[inds != pp]
      other_cov_levels <- as.matrix(unique(object$results.trans$covs[,inds]))
      if(size(other_cov_levels)[1]==1) {
        other_cov_levels <- matrix(other_cov_levels,nrow=size(other_cov_levels)[2])
      }
      all_pairs <- t(combn(object$results.trans$dpreds[pp],2))
      trans.local.res <- get_tp_local_test_results_by_level(object$results.trans,pp,delta,other_cov_levels,all_pairs,digits)
      trans.local.mean.diff[[pp]] <- trans.local.res$mean.diff
      trans.local.null.test[[pp]] <- trans.local.res$null.test
    }
    summary_list$trans.local.mean.diff <- trans.local.mean.diff
    summary_list$trans.local.null.test <- trans.local.null.test
  }
  if ("results.duration" %in% names(object)) {
    summary_list$dur.global <- global_test(object$results.duration) 
    summary_list$dur.mix.params <- get_mixture_params(object$results.duration,digits)
    summary_list$dur.mix.probs <- get_mixture_probs_by_cov(object$results.duration,digits)
  }
  class(summary_list) <- "BMRMMsummary"
  return(summary_list)
}

#############################################################
############## Plot Function for BMRMM Summary ##############
#############################################################
#' Plot Method for Visualizing BMRMM Summary
#'
#' Visualization of a specified field of a `BMRMMsummary` object.
#'
#' @param x an object of class `BMRMMsummary`.
#' @param type a string indicating the plot(s) to draw. Must be named after a field of `object`.
#' @param xlab x-axis label. Default is NULL.
#' @param ylab y-axis label. Default is NULL.
#' @param main main title. Default is NULL.
#' @param col color of the plot. Default is NULL.
#' @param ... further arguments for the plot function.
#' 
#' @return None
#' 
#' @examples
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50, 
#'                  cov.labels = list(c("F", "W"), c("U", "L", "A")),
#'                  duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
#' fit.summary <- summary(results)
#' plot(fit.summary, 'trans.probs.mean')
#' plot(fit.summary, 'dur.mix.probs')
#' 
#' @seealso [BMRMM::summary.BMRMM()]
#' 
#' @export
plot.BMRMMsummary <- function(x,type,xlab=NULL,ylab=NULL,main=NULL,col=NULL,...) {
  if (!inherits(x, "BMRMMsummary")) {
    stop("'x' must be of class 'BMRMMsummary'")
  } 
  if (!type %in% names(x)) {
    stop(paste("'x' does not have field", type), sep=' ')
  }
  object <- x
  if (type == 'trans.global') {
    if (missing(main)) main <- 'Global Test Results for Transition Probabilities'
    if (missing(xlab)) xlab <- ''
    if (missing(ylab)) ylab <- 'Proportions'
    if (missing(col)) col <- terrain.colors(nrow(object$trans.global))
    args.legend <- list(x="topright",legend=rownames(object$trans.global),xpd=TRUE)
    barplot(object$trans.global,main=main,ylab=ylab,col=col,ylim=c(0,1),
            legend.text=TRUE,beside=TRUE,xlab=xlab,args.legend=args.legend,...)
  }
  else if (type == 'trans.probs.mean') {
    if (missing(main)) main <- 'Trans. Probs. Posterior Mean'
    if (missing(xlab)) xlab <- 'y_{t-1}'
    if (missing(ylab)) ylab <- 'y_{t}'
    for (ii in 1:length(object$trans.probs.mean)) {
      cov.name <- paste('Covariates',names(object$trans.probs.mean)[ii],sep=' ')
      plot.heatmap(object$trans.probs.mean[[ii]],xlab=xlab,ylab=ylab,main=paste(cov.name,main,sep='\n'),...)
    }
  } 
  else if (type == 'trans.probs.sd') {
    if (missing(main)) main <- 'Trans. Probs. Posterior Std'
    if (missing(xlab)) xlab <- 'y_{t-1}'
    if (missing(ylab)) ylab <- 'y_{t}'
    for (ii in 1:length(object$trans.probs.sd)) {
      cov.name <- paste('Covariates',names(object$trans.probs.sd)[ii],sep=' ')
      plot.heatmap(object$trans.probs.sd[[ii]],xlab=xlab,ylab=ylab,main=paste(cov.name,main,sep='\n'),...)
    }
  } 
  else if (type == 'trans.local.mean.diff') {
    if (missing(main)) main <- 'Mean Diff. in Trans. Probs.'
    if (missing(xlab)) xlab <- 'y_{t-1}'
    if (missing(ylab)) ylab <- 'y_{t}'
    for(pp in 1:length(object$trans.local.mean.diff)) {
      for(rr in 1:length(object$trans.local.mean.diff[[pp]])) {
        fixed.cov.names <- names(object$trans.local.mean.diff[[pp]])[rr]
        for(kk in 1:length(object$trans.local.mean.diff[[pp]][[rr]])) {
          comp.cov.names <- names(object$trans.local.mean.diff[[pp]][[rr]])[kk]
          cov.names <- paste(fixed.cov.names,comp.cov.names,sep='\n')
          plot.heatmap(object$trans.local.mean.diff[[pp]][[rr]][[kk]],xlab=xlab,ylab=ylab,main=paste(cov.names,main,sep='\n'),...)
        }
      }
    }
  } 
  else if (type == 'trans.local.null.test') {
    if (missing(main)) main <- 'Probability for H_0'
    if (missing(xlab)) xlab <- 'y_{t-1}'
    if (missing(ylab)) ylab <- 'y_{t}'
    for(pp in 1:length(object$trans.local.null.test)) {
      for(rr in 1:length(object$trans.local.null.test[[pp]])) {
        fixed.cov.names <- names(object$trans.local.null.test[[pp]])[rr]
        for(kk in 1:length(object$trans.local.null.test[[pp]][[rr]])) {
          comp.cov.names <- names(object$trans.local.null.test[[pp]][[rr]])[kk]
          cov.names <- paste(fixed.cov.names,comp.cov.names,sep='\n')
          plot.heatmap(object$trans.local.null.test[[pp]][[rr]][[kk]],xlab=xlab,ylab=ylab,main=paste(cov.names,main,sep='\n'),...)
        }
      }
    }
  } 
  else if (type == 'dur.global') {
    if (missing(main)) main <- 'Global Test Results for Duration Times'
    if (missing(xlab)) xlab <- ''
    if (missing(ylab)) ylab <- 'Proportions'
    if (missing(col)) col <- terrain.colors(nrow(object$dur.global))
    args.legend <- list(x="topright",legend=rownames(object$dur.global),xpd=TRUE)
    barplot(object$dur.global,main=main,ylab=ylab,col=col,ylim=c(0,1),
            legend.text=TRUE,beside=TRUE,xlab=xlab,args.legend=args.legend,...)
  }
  else if (type == 'dur.mix.params') {
    if (missing(main)) {
      main <- 'Mixture Parameters'
    }
    plot.heatmap(object$dur.mix.params,xlab='',ylab='',main=main,...)
  } 
  else if (type == 'dur.mix.probs') {
    if(missing(main)) {
      main <- 'Mixture Probabilities for'
    }
    for(ii in 1:length(object$dur.mix.probs)) {
      cov.name <- paste("Covariate", names(object$dur.mix.probs)[ii], sep=' ')
      plot.heatmap(object$dur.mix.probs[[ii]],xlab='',ylab='',main=paste(main,cov.name,sep=' '),...)
    }
  } 
  else {
    stop("'type' must be one of the following: 'trans.global', trans.probs.mean', 
          'trans.probs.sd', 'trans.local.mean.diff', 'trans.local.null.test', 
          'dur.global', 'dur.mix.params', 'dur.mix.probs'")
  }
}

###########################################################
############## Histograms for Duration Times ##############
###########################################################

#' Histogram of Duration Times
#'
#' Plots the histogram of duration times in two ways as the users desire:
#' \enumerate{
#' \item Histogram of all duration times superimposed the posterior mean mixture gamma distribution;
#' \item Histogram of a specified mixture component superimposed the gamma distribution with shape and rate
#' parameters taken from the last MCMC iteration.
#' }
#'
#' @param x an object of class `BMRMM`.
#' @param comp one of 
#' \itemize{
#' \item `NULL`, which means the histogram for all duration times is plotted with the posterior mean mixture gamma distribution. Default option.
#' \item an integer specifying the mixture component for which the corresponding histogram is plotted with mixture gamma parameters taken from the last MCMC iteration.
#' }
#' @param xlim a range of x values with sensible defaults. Default is `NULL`, which is to use `c(min(duration), max(duration))`.
#' @param breaks an integer giving the number of cells for the histogram. Default is `NULL`, which is to use the Freedman-Diaconis rule, i.e., `(max(duration)-min(duration))*n^(1/3)/2/IQR(duration)`.
#' @param main main title. Default is `NULL`, which is to use `"Histogram with Posterior Mean"` when `comp` is `NULL` and `"Component X"` if `comp` is specified.
#' @param col color of the histogram bars. Default is `gray`.
#' @param xlab x-axis label. Default is `"Duration times"`.
#' @param ylab y-axis label. Default is `"Density"`.
#' @param ... further arguments for the hist function.
#' 
#' @examples
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50, 
#'                  duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
#' 
#' # plot the histogram of all duration times superimposed with 
#' # the posterior mixture gamma distribution
#' hist(results, xlim = c(0, 1), breaks = 50)
#' 
#' # plot the histogram for components 1 superimposed with 
#' # the mixture gamma distribution of the last MCMC iteration
#' hist(results, components = 1)
#' 
#' @return An object of class `histogram`.
#'
#' @export
hist.BMRMM <- function(x,comp=NULL,xlim=NULL,breaks=NULL,main=NULL,
                       col='gray',xlab='Duration times',ylab='Density',...) {
  if (!inherits(x, "BMRMM")) stop("'x' must be of class 'BMRMM'")
  if (!"results.duration" %in% names(x)) stop("'x' does not contain results for duration times")
  object <- x$results.duration
  if (missing(comp)) {
    dat <- object$duration.times
  } else if (!is.numeric(comp) || comp < 1 || comp > size(object$shape.samples,2)) {
    stop(sprintf("'comp' must be an integer within range [%s, %s]",1,size(object$shape.samples,2)))
  } else {
    dat <- object$duration.times[object$comp.assignment==comp]
    if (length(dat)==0) stop(paste("Component ",comp,'is empty'))
  }
  if(missing(xlim)) xlim <- c(min(dat),max(dat))
  if(missing(breaks)) breaks <- (max(dat)-min(dat))*length(dat)^(1/3)/2/IQR(dat)
  if(missing(main) || is.null(comp)) main <- 'Histogram with Posterior Mean'
  if(missing(main) || !is.null(comp)) main <- paste('Component',comp)
  x_axis <- get_x_axis(xlim)
  if (missing(comp)) {
    quan <- quantile_calc(object,object$marginal.prob,object$shape.samples,object$rate.samples,xlim)
    res <- hist(dat,freq=FALSE,breaks=breaks,xlim=xlim,col=col,main=main,xlab=xlab,ylab=ylab,...)
    polygon(c(x_axis,rev(x_axis)),c(quan$p95,rev(quan$p5)), # draw credible region
            col=rgb(190,190,190,130,maxColorValue=255), border=NA)
    lines(x=x_axis,y=quan$p50,col="red",lwd=2) # superimpose posterior mean
  } else {
    res <- hist(dat,freq=FALSE,xlim=xlim,main=main,breaks=breaks,col=col,border=FALSE,xlab=xlab,ylab=ylab,...)
    lines(x=x_axis,y=dgamma(x_axis,shape=object$shape.samples[nrow(object$shape.samples),comp],
                            rate=object$rate.samples[nrow(object$shape.samples),comp]),col="red",lwd=2)
  }
  return(res)
}

##############################################
############## Diagnostic Plots ##############
##############################################

#' MCMC Diagnostic Plots for Transition Probabilities and Duration Times
#'
#' Provides the traceplots and autocorrelation plots for (i) transition probabilities and (ii) mixture gamma shape and rate parameters.
#'
#' @param object an object of class `BMRMM`
#' @param cov.combs a list of covariate level combinations. Default is `NULL`, which is all possible combination of covariate levels.
#' @param transitions a list of pairs denoting state transitions. Default is `NULL`, which is all possible state transitions.
#' @param components a numeric vector denoting the mixture components of interest. Default is `NULL`, which is a list of all mixture components.
#' 
#' @return None
#'
#' @examples
#'
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 80, 
#'                  duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
#' diag.BMRMM(results)
#' diag.BMRMM(results, cov.combs = list(c(1,1),c(1,2)), 
#'            transitions = list(c(1,1)), components = c(3))
#'
#' @export
diag.BMRMM <- function(object,cov.combs=NULL,transitions=NULL,components=NULL) {
  if (!inherits(object, "BMRMM")) {
    stop("'object' must be of class 'BMRMM'")
  }
  if ("results.trans" %in% names(object)) {
    ind <- match("results.trans", names(object))
    result <- object[[ind]]
    if(missing(cov.combs)) {
      cov.combs <- sortrows((unique(result$covs)))
      rownames(cov.combs) <- 1:nrow(cov.combs)
      cov.combs <- as.list(data.frame(t(cov.combs)))
    }
    if(missing(transitions)) {
      transitions <- expand.grid(1:length(result$state.labels),1:length(result$state.labels))
      transitions <- as.list(data.frame(t(transitions)))
    }
    for(row in cov.combs) {
      label <- unlist(lapply(1:length(row), function(x) result$cov.labels[[x]][row[x]]))
      title <- paste(c("Covariates {",label,'}'),collapse=" ")
      d0 <- size(result$tp.exgns.all.itns)[1]
      simsize <- result$MCMCparams$simsize
      burnin <- result$MCMCparams$burnin
      thin_ind <- seq(burnin,simsize,result$MCMCparams$N_Thin)
      for(transition in transitions) {
        i <- transition[1]
        j <- transition[2]
        ind <- cbind(repmat(c(i,j,row),simsize,1),1:simsize)
        plot(result$tp.exgns.all.itns[ind],type='l',xlab='Iteration number',col='#00A600',
             ylab=paste(c('Prob(',result$state.labels[j],'->',result$state.labels[i],')'),collapse=''),
             main=paste('Traceplot for',title),ylim=c(0,1))
        abline(v=burnin,col='red',lwd=2)
        abline(h=mean(result$tp.exgns.all.itns[ind][burnin:simsize]),col='blue',lty=2,lwd=2)
        acf(result$tp.exgns.all.itns[ind][thin_ind],main=paste('Autocorrelation Plot for',title))
      }
    }
  } 
  if ("results.duration" %in% names(object)) {
    ind <- match("results.duration", names(object))
    result <- object[[ind]]
    burnin <- result$MCMCparams$burnin
    simsize <- result$MCMCparams$simsize
    af.thin <- seq(burnin+1,simsize,5)
    if(missing(components))
      components <- 1:ncol(result$shape.samples)
    for(kk in components) {
      plot(result$shape.samples[,kk],type='l',xlab='',ylab='',col='#00A600',main=paste("Traceplot for Shape, Comp",kk))
      abline(v=burnin,col='red',lwd=2)
      abline(h=mean(result$shape.samples[burnin:simsize,kk]),col='blue',lty=2,lwd=2)
      plot(result$rate.samples[,kk],type='l',xlab='',ylab='',col='#00A600',main=paste("Traceplot for Rate, Comp",kk))
      abline(v=burnin,col='red',lwd=2)
      abline(h=mean(result$rate.samples[burnin:simsize,kk]),col='blue',lty=2,lwd=2)
      acf(result$shape.samples[af.thin,kk],main=paste('Autocorrelation after thin. for Shape, Comp',kk))
      acf(result$rate.samples[af.thin,kk],main=paste('Autocorrelation after thin. for Rate, Comp',kk))
    }
  } 
}

#############################################
############## Model Selection ##############
#############################################

#' Model Selection Scores for the Number of Components for Duration Times
#'
#' Provides the LPML (Geisser and Eddy, 1979) and WAIC (Watanabe, 2010) scores of the Bayesian Markov renewal mixture models
#'
#' The two scores can be used to compare different choices of isi_num_comp, i.e., the number
#' of the mixture gamma components. Larger values of LPML and smaller values of WAIC
#' indicate better model fits.
#' 
#' @examples
#'
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50, 
#'                  duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
#' model.selection.scores(results)
#'
#' @references
#' Geisser, S. and Eddy, W. F. (1979). A predictive approach to model selection. Journal of the American Statistical Association, 74, 153–160. \cr\cr
#' Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11, 3571–3594.
#'
#' @param object An object of class BMRMM.
#' 
#' @return a list consisting of LPML and WAIC scores for gamma mixture models.
#' 
#' @export
model.selection.scores <- function(object) {
  if (!inherits(object, "BMRMM")) {
    stop("'object' must be of class 'BMRMM'")
  }
  if (!"results.duration" %in% names(object)) {
    stop("'object' does not contain results for duration times")
  }
  ind <- match("results.duration",names(object))
  isi <- object[[ind]]$duration.times
  pis <- object[[ind]]$marginal.prob
  alphas <- object[[ind]]$shape.samples
  betas <- object[[ind]]$rate.samples
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
  return(list('LPML'=lpml,'WAIC'=waic+2*p_waic))
}

