#' Generate a sample from the contaminated normal mixture
#'
#' Draws \code{n} observations from the mixture distribution used in the
#' robustness study of Lai et al. (2012):
#' \deqn{Y \sim (1 - B)\,N(\mu,\,\sigma_a^2) + B\,N(\mu,\,\sigma_b^2)}
#' where \eqn{B \sim \mathrm{Bernoulli}(p)},
#' \eqn{\sigma_a^2 = \sigma^2 / 1.1} and \eqn{\sigma_b^2 = 2\sigma^2 / 1.1}.
#'
#' @param n Integer. Number of observations to generate.
#' @param mu Numeric. Mean of the distribution.
#' @param sigma Numeric. Base standard deviation (before scaling).
#' @param p Numeric. Contamination probability (default: 0.1).
#'
#' @export
#' @return A numeric vector of length \code{n}.
#'
#' @seealso \code{\link{dmixnorm}}, \code{\link{pmixnorm}}
#' @examples
#' set.seed(1)
#' y <- generate_mixture_sample(n = 50, mu = 0, sigma = 1, p = 0.1)
#' summary(y)
generate_mixture_sample <- function(n, mu, sigma, p = 0.1) {
  B <- rbinom(n, size = 1, prob = p)

  sigma1_sq <- (1 / 1.1) * sigma^2
  sigma2_sq <- (2 / 1.1) * sigma^2

  y    <- numeric(n)
  idx0 <- which(B == 0)
  idx1 <- which(B == 1)

  if (length(idx0) > 0) y[idx0] <- rnorm(length(idx0), mu, sqrt(sigma1_sq))
  if (length(idx1) > 0) y[idx1] <- rnorm(length(idx1), mu, sqrt(sigma2_sq))

  y
}

#' Standard deviations of the two mixture components
#'
#' Returns the standard deviations of the non-contaminated and contaminated
#' components of the mixture as defined in Lai et al. (2012).
#'
#' @param sigma Numeric. Base standard deviation.
#'
#' @export
#' @return A named numeric vector with elements \code{sd_a}
#'   (\eqn{\sigma / \sqrt{1.1}}) and \code{sd_b}
#'   (\eqn{\sigma\sqrt{2 / 1.1}}).
#' @examples
#' mix_sd(sigma = 1)
mix_sd <- function(sigma) {
  c(
    sd_a = sqrt((1 / 1.1) * sigma^2),
    sd_b = sqrt((2 / 1.1) * sigma^2)
  )
}

#' Density of the contaminated normal mixture
#'
#' @param x Numeric vector. Points at which to evaluate the density.
#' @param mu Numeric. Mean of the distribution.
#' @param sigma Numeric. Base standard deviation.
#' @param p Numeric. Contamination probability (default: 0.1).
#'
#' @export
#' @return A numeric vector of density values.
#'
#' @seealso \code{\link{pmixnorm}}, \code{\link{generate_mixture_sample}}
#' @examples
#' dmixnorm(x = seq(-3, 3, by = 1), mu = 0, sigma = 1, p = 0.1)
dmixnorm <- function(x, mu, sigma, p = 0.1) {
  sds <- mix_sd(sigma)
  (1 - p) * dnorm(x, mu, sds["sd_a"]) +
    p      * dnorm(x, mu, sds["sd_b"])
}

#' CDF of the contaminated normal mixture
#'
#' @param x Numeric vector. Points at which to evaluate the CDF.
#' @param mu Numeric. Mean of the distribution.
#' @param sigma Numeric. Base standard deviation.
#' @param p Numeric. Contamination probability (default: 0.1).
#'
#' @export
#' @return A numeric vector of cumulative probabilities.
#'
#' @seealso \code{\link{dmixnorm}}, \code{\link{generate_mixture_sample}}
#' @examples
#' pmixnorm(x = 0, mu = 0, sigma = 1, p = 0.1)
pmixnorm <- function(x, mu, sigma, p = 0.1) {
  sds <- mix_sd(sigma)
  (1 - p) * pnorm(x, mu, sds["sd_a"]) +
    p      * pnorm(x, mu, sds["sd_b"])
}

#' True Youden index and optimal cut-point under the mixture model
#'
#' Numerically maximises the Youden index over all possible cut-points for
#' two contaminated normal mixture distributions.
#'
#' @param mu1 Numeric. Mean of the diseased population.
#' @param mu2 Numeric. Mean of the non-diseased population.
#' @param sigma1 Numeric. Base standard deviation of the diseased population.
#' @param sigma2 Numeric. Base standard deviation of the non-diseased
#'   population.
#' @param p Numeric. Contamination probability (default: 0.1).
#'
#' @export
#' @return A list with two elements:
#' \describe{
#'   \item{J}{The maximum Youden index.}
#'   \item{c}{The optimal cut-point that achieves it.}
#' }
#'
#' @seealso \code{\link{find_mu1_from_J_mixture}}
#' @examples
#' J_true_mixture(mu1 = 1.68, mu2 = 0, sigma1 = 1, sigma2 = 1, p = 0.1)
J_true_mixture <- function(mu1, mu2, sigma1, sigma2, p = 0.1) {
  J_of_c <- function(c) {
    pmixnorm(c, mu2, sigma2, p) - pmixnorm(c, mu1, sigma1, p)
  }

  lo  <- min(mu1, mu2) - 20 * max(sigma1, sigma2)
  hi  <- max(mu1, mu2) + 20 * max(sigma1, sigma2)
  opt <- optimize(J_of_c, interval = c(lo, hi), maximum = TRUE)

  list(J = opt$objective, c = opt$maximum)
}
