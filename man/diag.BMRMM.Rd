% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_and_plot_results.R
\name{diag.BMRMM}
\alias{diag.BMRMM}
\title{MCMC Diagnostic Plots for Transition Probabilities and Duration Times}
\usage{
diag.BMRMM(object, cov.combs = NULL, transitions = NULL, components = NULL)
}
\arguments{
\item{object}{an object of class \code{BMRMM}}

\item{cov.combs}{a list of covariate level combinations. Default is \code{NULL}, which is all possible combination of covariate levels.}

\item{transitions}{a list of pairs denoting state transitions. Default is \code{NULL}, which is all possible state transitions.}

\item{components}{a numeric vector denoting the mixture components of interest. Default is \code{NULL}, which is a list of all mixture components.}
}
\value{
None
}
\description{
Provides the traceplots and autocorrelation plots for (i) transition probabilities and (ii) mixture gamma shape and rate parameters.
}
\examples{

results <- BMRMM(foxp2sm, num.cov = 2, simsize = 80, 
                 duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
diag.BMRMM(results)
diag.BMRMM(results, cov.combs = list(c(1,1),c(1,2)), 
           transitions = list(c(1,1)), components = c(3))

}
