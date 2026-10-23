#' gpqyouden: Generalized Pivotal Quantity Inference for the Youden Index
#'
#' Monte Carlo simulation study of GPQ-based confidence intervals for the
#' Youden index and its optimal cut-point, following Lai, Tian and
#' Schisterman (2012).
#'
#' @keywords internal
#' @importFrom stats dnorm optimize pnorm quantile rbinom rchisq rnorm rt sd
#' @importFrom stats uniroot var
#' @importFrom utils globalVariables
"_PACKAGE"

# Columns referenced via NSE inside dplyr/tidyr pipes (tables.R), plus the
# foreach loop index (simulation.R): neither is a real global, but static
# analysis in R CMD check cannot see that they are supplied at call time.
utils::globalVariables(c(
  "n1", "n2", "J", "var1", "var2", "mu2", "AUC",
  "coverage_J", "coverage_c", "length_J", "length_c", "i"
))
