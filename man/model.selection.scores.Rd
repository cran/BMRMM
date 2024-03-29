% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/print_and_plot_results.R
\name{model.selection.scores}
\alias{model.selection.scores}
\title{Model Selection Scores for the Number of Components for Duration Times}
\usage{
model.selection.scores(object)
}
\arguments{
\item{object}{An object of class BMRMM.}
}
\value{
a list consisting of LPML and WAIC scores for gamma mixture models.
}
\description{
Provides the LPML (Geisser and Eddy, 1979) and WAIC (Watanabe, 2010) scores of the Bayesian Markov renewal mixture models
}
\details{
The two scores can be used to compare different choices of isi_num_comp, i.e., the number
of the mixture gamma components. Larger values of LPML and smaller values of WAIC
indicate better model fits.
}
\examples{

results <- BMRMM(foxp2sm, num.cov = 2, simsize = 50, 
                 duration.distr = list('mixgamma',shape=rep(1,3),rate=rep(1,3)))
model.selection.scores(results)

}
\references{
Geisser, S. and Eddy, W. F. (1979). A predictive approach to model selection. Journal of the American Statistical Association, 74, 153–160. \cr\cr
Watanabe, S. (2010). Asymptotic equivalence of Bayes cross validation and widely applicable information criterion in singular learning theory. Journal of Machine Learning Research, 11, 3571–3594.
}
