
<!-- README.md is generated from README.Rmd. Please edit that file -->

# BMRMM

<!-- badges: start -->
<!-- badges: end -->

The goal of BMRMM is to analyze sequential categorical data with
continuous duration times, being either state durations or inter-state
durations. BMRMMs can comprehensively analyze the stochastic dynamics of
both state transitions and duration times under the influence of
multiple exogenous factors and random individual effect. The BMRMM
package provides a suite of functions implementing such models. The
default setting flexibly models the transition probabilities using
mixtures of Dirichlet distributions and the duration times using
mixtures of gamma kernels. It also provides users the flexibility of
modeling the categorical sequences using Bayesian Markov mixed models
alone, either ignoring the duration times altogether or treating them as
an additional category in the sequence when their values are discretized
by a user-specified unit. The package allows extensive inference of the
state transition probabilities and the duration times as well as
relevant plots and graphs.

## Installation and loading the package

You can install the development version of BMRMM like so:

``` r
install.packages('BMRMM')
library(BMRMM)
```

## Main function

Here are some examples of running the main BMRMM function.

``` r
library(BMRMM)

# The included synthetic foxp2 data set contains 3 covariates

# ignores duration times and only models transition probabilities using all three covariates
results <- BMRMM(foxp2,num_cov=3,duration_type='None')

# models duration times as a continuous variable with 5 gamma mixture components,
# using covariates 2,3 and the previous state
results <- BMRMM(foxp2,num_cov=3,trans_cov_index=c(1,3),'duration_type'='Continuous',
                                duration_cov_index=c(2,3),duration_num_comp=5)

# models duration times as a discrete state with discretization 0.25 and
# do not include the previous state as a covariate
results <- BMRMM(foxp2,num_cov=3,duration_type='Discrete',
                 duration_excl_prev_state=TRUE,duration_unit=0.25)
```

## Plotting function (for transition probabilities)

Here are some examples of plotting results for transition probabilities.

``` r
# GLOBAL TESTS
get_global_test_results(results$results_trans,decimal_pts=3)

# POSTERIOR MEAN & SD
# get posterior results for all combinations of covariate levels
get_estimated_post_mean_and_sd(results$results_trans)

# The covariate labels for transition probabilities for foxp2 data set are
#        [,1]  [,2] [,3]
#   [1,] "HET" "WT" ""
#   [2,] "U"   "L"  "A"

# get posterior results for covariate levels ("HET","U") and ("WT","U")
cov_labels <- matrix(c("HET","WT","","U","L","A"),nrow=2,byrow=TRUE)
cov_levels <- matrix(c(1,1,2,1),nrow=2,byrow=TRUE)
get_estimated_post_mean_and_sd(results$results_trans,cov_levels,cov_labels)

# models duration times as a discrete state with discretization 0.25 and
# do not include the previous state as a covariate
results <- BMRMM(foxp2,num_cov=3,duration_type='Discrete',
                 duration_excl_prev_state=TRUE,duration_unit=0.25)

# LOCAL TESTS
# results for genotypes (HET, WT) under all three social contexts
get_tp_local_test_results(results$results_trans,cov=1)

# results for genotypes (HET, WT) under social contexts U and A
get_tp_local_test_results(results$results_trans,cov=1,other_cov_levels=c(1,3))

# results for contexts (U,L), (U,A) and (L,A) under two genotypes
get_tp_local_test_results(results$results_trans,cov=2)

# results for contexts (U,L), (U,A) under genotype HET (specify pairs)
get_tp_local_test_results(results$results_trans,cov=2,other_cov_levels=c(1),
                           comp_pairs=matrix(c(1,1,1,3),nrow=2,byrow=2))

# DIAGNOSTIC PLOTS
# results for all transition types for every combinations of covariates
get_tp_diagnostic_plots(results$results_trans)

# results for all transition types for covariate levels (HET,U) and (WT,A)
get_tp_diagnostic_plots(results$results_trans,cov_levels=matrix(c(1,1,2,3),nrow=2,byrow=TRUE))

# results for all transition from (d,s) to (s,u) for covariate levels (HET,U) and (WT,A)
get_tp_diagnostic_plots(results$results_trans,from=c(1,3),to=c(3,4),
                         cov_levels=matrix(c(1,1,2,3),nrow=2,byrow=TRUE))
```

## Plotting function (for duration times)

Here are some examples of plotting results for duration times.

``` r
# GLOBAL TESTS
get_global_test_results(results$results_duration,decimal_pts=3)

# POSTERIOR MEAN
# histogram for all data superimposed with posterior mean
get_histogram_with_post_mean(results$results_duration,x_range=c(0,1),breaks=400)
# histogram by mixture component superimposed with posterior mean
get_histogram_by_component(results$results_duration,comps=c(1,2,5))

# HEATMAP FOR MIXTURE PARAMETERS
get_heatmap_mix_param(results$results_duration,decimal_pts=3)

# The covariate labels for duration times for foxp2 data set are
#        [,1]  [,2] [,3] [,4]
#   [1,] "HET" "WT"  ""   ""
#   [2,] "U"   "L"  "A"   ""
#   [3,] "d"   "m"  "s"   "u"

# get heatmap for mixture probabilities for all covariates
get_heatmap_by_cov(results$results_duration)
# get heatmap for mixture probabilities for just genotype and context
get_heatmap_by_cov(results$results_duration,cov_index=c(1,2))

# DIAGNOSTIC PLOTS
# get diagnostic plots for all components
get_duration_diagnostic_plots(results$results_duration)
# get diagnostic plots for components 1,2,5
get_duration_diagnostic_plots(results$results_duration,comps=c(1,2,5))

# MODEL SELECTION SCORES (LPML & WAIC SCORES)
# larger values of LPML and smaller values of WAIC indicate better model fits.
get_lpml_waic_scores(results$results_duration)
```
