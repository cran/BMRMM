#' Simulated FoxP2 Data Set.
#'
#' A simulated data set of the original FoxP2 data set, which contains
#' the sequences of syllables sung by male mice of
#' different genotypes under various social contexts.
#'
#' @format A data frame with 70818 rows and 6 variables:
#' \describe{
#'   \item{Id}{Mouse Id}
#'   \item{Genotype}{Genotype of the mouse, 1 = FoxP2 knocked out, 2 = wild type}
#'   \item{Context}{Social context for the mouse, 1 = U (urine sample placed in the cage), 2 = L (living female mouse placed in the cage), 3 = A (an anesthetized female placed on the lid of the cage)}
#'   \item{Prev_State}{The previous syllable, \{1,2,3,4\}=\{d,m,s,u\}}
#'   \item{Cur_State}{The current syllable, \{1,2,3,4\}=\{d,m,s,u\}}
#'   \item{ISI}{Modified inter-syllable interval times, log(original ISI + 1)}
#' }
"foxp2"
