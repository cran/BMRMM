#' @import fields
#' @import stats
#' @import pracma
#' @import graphics
#' @import grDevices
#' @import gsubfn
#' @import logOfGamma
#' @import mhsmm
#' @import multicool
#' @importFrom MCMCpack rdirichlet
#' @importFrom utils combn
NULL

#' Bayesian Markov Renewal Mixed Models (BMRMMs)
#'
#' Provides inference results of both transition probabilities and duration times using BMRMMs.
#'
#' Returns list(trans_results, duration_results) \cr
#'
#' (1) Users have the option to ignore duration times or model duration times as a discrete or continuous variable via defining 'duration_type':\cr
#'       "None": ignores duration times \cr
#'       "Continuous": treat duration times as a continuous variable (DEFAULT)\cr
#'       "Discrete": treat duration times as a new state with discretization 'duration_unit'.
#'       If 'Discrete' is used, duration times becomes a new state type. 'duration_unit' must be specified if 'duration_type' is 'Discrete'.
#'       For example, if an duration time entry is 20 and 'duration_unit' is 5, then the model will add 4 consecutive new states.
#'       If an duration time entry is 23.33 and 'duration_unit' is 5, then the model will still add 4 consecutive new states
#'       as the blocks are calculated with the floor operation  \cr
#'
#' (2) 'trans_results' contains the results of the transition probabilities:\cr
#'                   - "Num_States": total number of states \cr
#'                   - "Xexgns": covariates related to transition probabilities \cr
#'                   - "dpreds": maximum level for each related covariate \cr
#'                   - "MCMCparams": MCMC parameters: simsize, burnin and thinning factor \cr
#'                   - "TP_Exgns_Post_Mean": posterior mean of transition probabilities for different combinations of exogenous predictors \cr
#'                   - "TP_Exgns_Post_Std": posterior standard deviation of transition probabilities for different combinations of exogenous predictors  \cr
#'                   - "TP_Anmls_Post_Mean": posterior mean of transition distribution components for different individuals \cr
#'                   - "TP_All_Post_Mean": posterior mean of transition distribution components for different combinations of exogenous predictors AND different individuals \cr
#'                   - "TP_Anmls_Post_Std": standard deviation of transition probabilities among different mice \cr
#'                   - "TP_Exgns_Diffs_Store": difference in posterior mean of transition probabilities for every pair of covariate levels given levels of the other covariates \cr
#'                   - "TP_Exgns_All_Itns": exogenous transition probabilities for every MCMC iteration \cr
#'                   - "Clusters": number of clusters for each covariate for each MCMC iteration \cr
#'                   - "Type": an identifier for results, which is "Transition Probabilities". \cr
#'
#' (3) 'duration_results' contains the results of the duration times ('duration_results' is NULL if 'duration_type' is 'None'):\cr
#'                 "K" <- number of gamma mixture components \cr
#'                 "Xexgns" <- covariates related to duration times\cr
#'                 "dpreds" <- maximum level for each related covariate\cr
#'                 "MCMCparams" <- MCMC parameters: simsize, burnin and thinning factor\cr
#'                 "Duration_Times" <- input duration times\cr
#'                 "Comp_Assign" <- mixture component assignment for each data point in the last MCMC iteration\cr
#'                 "Duration_Exgns_Store" <- posterior mean of mixture probabilities for different combinations of exogenous predictors of each MCMC iteration\cr
#'                 "Marginal_Prob" <- estimated marginal mixture probabilities for each iteration\cr
#'                 "Shape_Samples" <- estimated shape parameters for gamma mixtures for each iteration\cr
#'                 "Rate_Samples" <- estimated rate parameters for gamma mixtures for each iteration\cr
#'                 "Clusters" <- number of clusters for each covariate for each MCMC iteration\cr
#'                 "Type" <- an identifier for results, which is "Inter-Syllable Intervals".\cr
#'
#'@return List of results for transition probabilities and durations times, list(trans_results, duration_results). See details.
#'
#'@examples
#'
#' # In the examples, we use a shorted version of the foxp2 dataset, foxp2_sm
#'
#' # ignores duration times and only models transition probabilities using all three covariates
#' results <- BMRMM(foxp2_sm,num_cov=2,duration_type='None',simsize=50)
#'
#' # models duration times as a continuous variable with 5 gamma mixture components,
#' # using covariate 1 and the previous state
#' results <- BMRMM(foxp2_sm,num_cov=2,trans_cov_index=c(1),'duration_type'='Continuous',
#'                  duration_cov_index=c(1),duration_num_comp=5,simsize=50)
#'
#' # models duration times as a discrete state with discretization 0.25 and
#' # do not include the previous state as a covariate
#' results <- BMRMM(foxp2_sm,num_cov=2,duration_type='Discrete',duration_excl_prev_state=TRUE,
#'                  duration_unit=0.025,simsize=50)
#'
#'
#' @author Yutong Wu, \email{yutong.wu@@utexas.edu}
#'
#' @param data a data frame containing -- individual ID, covariate values, previous state, current state, duration times (if applicable), in that order
#' @param num_cov total number of covariates provided in 'data'
#' @param switch_off_random TRUE if only population-level effects are considered, default is FALSE
#' @param trans_cov_index indices of covariates that are used for transition probabilities, default is all of the covariates
#' @param duration_type one of 'None', 'Discrete', 'Continuous', default is 'Continuous'
#' @param duration_unit the discretization of duration times, only used when 'duration_type' is 'Discrete'
#' @param duration_excl_prev_state TRUE if the previous state is excluded from duration times inference, default is FALSE
#' @param duration_cov_index indices of covariates that are used for duration times, default is all of the covariates plus the previous state
#' @param duration_num_comp number of gamma mixture components for duration times, default is 4
#' @param duration_init_shape initialization of mixture gamma shape parameters, default is a vector of 1 of size 'duration_num_comp'
#' @param duration_init_rate initialization of mixture gamma rate parameters, default is a vector of 1 of size 'duration_num_comp'
#' @param simsize total number of MCMC iterations, default is 10000
#' @param burnin number of burn-ins for the MCMC iterations, default is simsize/2
#' 
#' @export
BMRMM <- function(data,num_cov,switch_off_random=FALSE,trans_cov_index=1:num_cov,duration_type='Continuous',
                  duration_cov_index=1:num_cov,duration_excl_prev_state=FALSE,duration_unit=NULL,
                  duration_num_comp=4,duration_init_shape=rep(1,duration_num_comp),
                  duration_init_rate=rep(1,duration_num_comp),simsize=10000,burnin=simsize/2) {

  if(num_cov>5 || num_cov==0) {
    stop("The model takes 1 to 5 covariates.")
  }

  # data should be in format: id, cov 1, ..., cov p, prev state, cur state, (duration if applicable)
  num_col <- ncol(data)
  id <- data[,1]
  covs <- as.matrix(data[,2:(num_cov+1)])
  prev_state <- data[,1+num_cov+1]
  cur_state <- data[,1+num_cov+2]

  # all covs are taken if none specified
  if (ncol(covs)==1) {
    trans_covs <- as.matrix(data[,2])
    duration_covs <- as.matrix(data[,2])
    colnames(trans_covs) <- colnames(data)[2]
    colnames(duration_covs) <- colnames(data)[2]
  } else {
    trans_covs <- covs[,trans_cov_index]
    duration_covs <- covs[,duration_cov_index]
  }

  if(!duration_excl_prev_state)
    duration_covs <- cbind(duration_covs,prev_state)

  # get results
  res_duration <- NULL
  if(duration_type=='None') {
    res_trans <- model_transition(data.frame(cbind(id,trans_covs,prev_state,cur_state)),switch_off_random,simsize,burnin)
  } else if(duration_type=='Discrete') {
    duration <- data[,num_col]
    if(is.null(duration_unit)) {
      stop("Must specify 'duration_unit' if treating duration time as a discrete variable.")
    } else {
      aug_data <- add_isi_as_state(data.frame(cbind(id,trans_covs,prev_state,cur_state,duration)),duration_unit)
      res_trans <- model_transition(aug_data,switch_off_random,simsize,burnin)
    }
  } else if(duration_type=='Continuous') {
    duration <- data[,num_col]
    data_for_trans <- data.frame(cbind(id,trans_covs,prev_state,cur_state))
    data_for_duration <- data.frame(cbind(id,duration_covs,duration))
    res_trans <- model_transition(data_for_trans,switch_off_random,simsize,burnin)
    res_duration <- model_cont_isi(data_for_duration,switch_off_random,simsize=simsize,burnin=burnin,K=duration_num_comp,duration_init_shape,duration_init_rate)
  } else {
    stop("'duration_type' must be one of 'None', 'Discrete' and 'Continuous'.")
  }

  # return results
  results <- list("results_trans"=res_trans,"results_duration"=res_duration)
  return(results)
}
