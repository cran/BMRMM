#' @import pracma 
#' @importFrom fields image.plot
#' @importFrom logOfGamma gammaln
#' @importFrom multicool Stirling2
#' @importFrom MCMCpack rdirichlet
#' @importFrom grDevices rgb terrain.colors
#' @importFrom graphics abline axis barplot grid hist lines polygon text
#' @importFrom stats IQR acf dgamma kmeans quantile rbeta rbinom rgamma rmultinom runif sd
#' @importFrom utils combn
NULL

#' Bayesian Markov Renewal Mixed Models (BMRMMs)
#'
#' Provides inference results of both transition probabilities and duration times using BMRMMs.
#'
#' Users have the option to ignore duration times or model duration times as 
#' a discrete or continuous variable via defining `duration.distr`. 
#' 
#' `duration.distr` can be one of the following: \cr
#' \itemize{
#' \item{`NULL`}{ duration times are ignored. This is the default setting.}
#' \item{`list('mixgamma', shape, rate)`}{ duration times are modeled as a mixture gamma variable. `shape` and `rate`
#'       must be numeric vectors of the same length. The length indicates the number of mixture components.}
#' \item{`list('mixDirichlet', unit)`}{ duration times are modeled as a new state with discretization `unit`. The duration
#'       state is then analyzed along with the original states. For example, if an duration time entry is 20 and `unit` is 5, 
#'       then the model will add 4 consecutive new states. If an duration time entry is 23.33 and `unit` is 5, then the model 
#'       will still add 4 consecutive new states as the blocks are calculated with the floor operation.}
#' }
#'
#'@return An object of class `BMRMM` consisting of `results.trans` and `results.duration` if duration times are analyzed as a continuous variable. \cr
#'
#' The field `results.trans` is a data frame giving the inference results of transition probabilities. 
#'\tabular{ll}{
#' `covs` \tab covariates levels for each row of the data. \cr
#' `dpreds` \tab maximum level for each related covariate. \cr
#' `MCMCparams` \tab MCMC parameters including simsize, burnin and thinning factor. \cr
#' `tp.exgns.post.mean`\tab posterior mean of transition probabilities for different combinations of covariates. \cr
#' `tp.exgns.post.std` \tab posterior standard deviation of transition probabilities for different combinations of covariates.  \cr
#' `tp.anmls.post.mean` \tab  posterior mean of transition probabilities for different individuals. \cr
#' `tp.anmls.post.std` \tab  posterior standard deviation of transition probabilities for different individuals. \cr
#' `tp.all.post.mean` \tab  posterior mean of transition probabilities for different combinations of covariates AND different individuals. \cr
#' `tp.exgns.diffs.store` \tab  difference in posterior mean of transition probabilities for every pair of covariate levels given levels of the other covariates. \cr
#' `tp.exgns.all.itns` \tab population-level transition probabilities for every MCMC iteration. \cr
#' `clusters` \tab number of clusters for each covariate for each MCMC iteration. \cr
#' `type` \tab a string identifier for results, which is "Transition Probabilities". \cr
#' `cov.labels` \tab a list of string vectors giving labels of covariate levels. \cr
#' `state.labels` \tab a list of strings giving labels of states. \cr
#'}
#' The field `results.duration` is a data frame giving the inference results of duration times. 
#'\tabular{ll}{
#' `covs` \tab covariates related to duration times.\cr
#' `dpreds` \tab maximum level for each related covariate.\cr
#' `MCMCparams` \tab MCMC parameters: simsize, burnin and thinning factor.\cr
#' `duration.times` \tab  duration times from the data set.\cr
#' `comp.assignment` \tab mixture component assignment for each data point in the last MCMC iteration.\cr
#' `duration.exgns.store` \tab posterior mean of mixture probabilities for different combinations of covariates of each MCMC iteration.\cr
#' `marginal.prob` \tab estimated marginal mixture probabilities for each MCMC iteration.\cr
#' `shape.samples` \tab estimated shape parameters for gamma mixtures for each MCMC iteration.\cr
#' `rate.samples` \tab estimated rate parameters for gamma mixtures for each MCMC iteration.\cr
#' `clusters` \tab number of clusters for each covariate for each MCMC iteration.\cr
#' `type` \tab a string identifier for results, which is "Duration Times".\cr
#' `cov.labels` \tab a list of string vectors giving labels of covariate levels. \cr
#'}
#'
#'@examples
#'
#' # In the examples, we use a shorted version of the foxp2 dataset, foxp2sm
#'
#' # ignores duration times and only models transition probabilities using all three covariates
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50)
#'
#' # models duration times as a continuous variable with 3 gamma mixture components,
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50,
#'                  duration.distr = list('mixgamma', shape = rep(1,3), rate = rep(1,3)))
#'
#' # models duration times as a discrete state with discretization 0.025 and
#' results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50, 
#'                  duration.distr = list('mixDirichlet', unit = 0.025))
#'
#'
#' @author Yutong Wu, \email{yutong.wu@@utexas.edu}
#'
#' @param data a data frame containing -- individual ID, covariate values, previous state, current state, duration times (if applicable), in that order.
#' @param num.cov total number of covariates provided in `data`.
#' @param cov.labels a list of vectors giving names of the covariate levels. Default is a list of numerical vectors. 
#' @param state.labels a vector giving names of the states. Default is a numerical vector.
#' @param random.effect `TRUE` if population-level effects are considered. Default is `TRUE`.
#' @param fixed.effect `TRUE` if individual-level effects are considered. Default is `TRUE`.
#' @param trans.cov.index a numeric vector indicating the indices of covariates that are used for transition probabilities. Default is all of the covariates.
#' @param duration.cov.index a numeric vector indicating the indices of covariates that are used for duration times. Default is all of the covariates.
#' @param duration.distr a list of arguments indicating the distribution of duration times. Default is `NULL`, which is ignoring duration times.
#' @param duration.incl.prev.state `TRUE` if the previous state is included in the inference of duration times. Default is `TRUE`.
#' @param simsize total number of MCMC iterations. Default is 10000.
#' @param burnin number of burn-ins for the MCMC iterations. Default is `simsize`/2.
#' 
#' @export
BMRMM <- function(data,num.cov,cov.labels=NULL,state.labels=NULL,random.effect=TRUE,fixed.effect=TRUE,
                  trans.cov.index=1:num.cov,duration.cov.index=1:num.cov,duration.distr=NULL,
                  duration.incl.prev.state=TRUE,simsize=10000,burnin=simsize/2) {

  if(num.cov>5 || num.cov==0) {
    stop("The model takes 1 to 5 covariates.")
  }

  # data should be in format: id, cov 1, ..., cov p, prev state, cur state, (duration if applicable)
  num_col <- ncol(data)
  id <- data[,1]
  covs <- as.matrix(data[,2:(num.cov+1)])
  prev_state <- data[,1+num.cov+1]
  cur_state <- data[,1+num.cov+2]
  
  # assign covariate labels if not given
  if(missing(cov.labels)){
    cov.labels <- apply(covs, 2, function(x) 1:length(unique(x)))
  }
    
  # assign state labels if not given
  if(missing(state.labels)){
    state.labels <- sort(union(unique(prev_state),unique(cur_state)))
  }

  # all covs are taken if none specified
  if (ncol(covs)==1) {
    trans_covs <- as.matrix(data[,2])
    duration_covs <- as.matrix(data[,2])
    colnames(trans_covs) <- colnames(data)[2]
    colnames(duration_covs) <- colnames(data)[2]
  } else {
    trans_covs <- covs[,trans.cov.index]
    duration_covs <- covs[,duration.cov.index]
  }
  
  trans.cov.labels <- cov.labels[trans.cov.index]
  duration.cov.labels <- cov.labels[duration.cov.index]

  if(duration.incl.prev.state) {
    duration_covs <- cbind(duration_covs,prev_state)
    duration.cov.labels$prev_state <- state.labels
  }

  # get results
  res_duration <- NULL
  if(missing(duration.distr)) {
    res_trans <- model_transition(data.frame(cbind(id,trans_covs,prev_state,cur_state)),random.effect,fixed.effect,simsize,burnin)
  } else if('mixDirichlet' %in% duration.distr) {
    if(is.null(duration.distr$unit)) {
      stop("Must specify 'unit' if 'duration.distr' contains 'mixDirichlet'.")
    } else {
      state.labels <- append(state.labels,'dur.state')
      duration <- data[,num_col]
      aug_data <- add_isi_as_state(data.frame(cbind(id,trans_covs,prev_state,cur_state,duration)),duration.distr$unit)
      res_trans <- model_transition(aug_data,random.effect,fixed.effect,simsize,burnin)
    }
  } else if('mixgamma' %in% duration.distr) {
    if (is.null(duration.distr$shape) || is.null(duration.distr$rate)) {
      stop("Must specify both 'shape' and 'rate' if 'duration.distr' contains 'mixgamma'.")
    } else if (length(duration.distr$shape) != length(duration.distr$rate)) {
      stop("'shape' and 'rate' must have the same dimension.")
    } else {
      duration <- data[,num_col]
      data_for_trans <- data.frame(cbind(id,trans_covs,prev_state,cur_state))
      data_for_duration <- data.frame(cbind(id,duration_covs,duration))
      res_trans <- model_transition(data_for_trans,random.effect,fixed.effect,simsize,burnin)
      res_duration <- model_cont_isi(data_for_duration,random.effect,fixed.effect,simsize=simsize,burnin=burnin,K=length(duration.distr$shape),duration.distr$shape,duration.distr$rate)
    }
  } else {
    stop("The argument 'duration.distr' must be a list containing either 'mixDirichlet' or 'mixgamma'.")
  }

  # return results
  res_trans$state.labels <- state.labels
  res_trans$cov.labels <- trans.cov.labels
  if (!is.null(res_duration)) {
    res_duration$cov.labels <- duration.cov.labels
    results <- list("results.trans"=res_trans,"results.duration"=res_duration)
  } else {
    results <- list("results.trans"=res_trans)
  }
  class(results) <- "BMRMM"
  return(results)
}
