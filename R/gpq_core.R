#' GPQ for the mean
#'
#' Computes the generalized pivotal quantity for the mean of one population
#' (Equation 7 in Lai et al., 2012).
#'
#' @param ybar Numeric. Observed sample mean.
#' @param s Numeric. Observed sample standard deviation.
#' @param n Integer. Sample size.
#' @param t_value Numeric. Realization of a Student's t with \code{n - 1}
#'   degrees of freedom.
#'
#' @export
#' @return A numeric scalar.
#' @examples
#' gpq_mu(ybar = 5, s = 1.2, n = 30, t_value = 0.5)
gpq_mu <- function(ybar, s, n, t_value) {
  ybar - t_value * s / sqrt(n)
}

#' GPQ for the standard deviation
#'
#' Computes the generalized pivotal quantity for the standard deviation of one
#' population (Equation 8 in Lai et al., 2012).
#'
#' @param s2 Numeric. Observed sample variance.
#' @param v Numeric. Realization of a chi-squared with \code{n - 1} degrees
#'   of freedom.
#' @param n Integer. Sample size.
#'
#' @export
#' @return A numeric scalar (standard deviation scale).
#' @examples
#' gpq_sigma(s2 = 1.5, v = 28.3, n = 30)
gpq_sigma <- function(s2, v, n) {
  sqrt((n - 1) * s2 / v)
}

#' GPQ for the optimal cut-point
#'
#' Computes the generalized pivotal quantity for the optimal cut-point
#' (Equations 10-11 in Lai et al., 2012).
#'
#' @param R_mu1 Numeric. GPQ realization for the mean of the diseased group.
#' @param R_mu2 Numeric. GPQ realization for the mean of the non-diseased group.
#' @param R_sigma1 Numeric. GPQ realization for the SD of the diseased group.
#' @param R_sigma2 Numeric. GPQ realization for the SD of the non-diseased
#'   group.
#' @param tol Numeric. Tolerance for the equal-variance case (default: 1e-2).
#'
#' @export
#' @return A numeric scalar.
#' @examples
#' gpq_c(R_mu1 = 2, R_mu2 = 0, R_sigma1 = 1, R_sigma2 = 1.2)
gpq_c <- function(R_mu1, R_mu2, R_sigma1, R_sigma2, tol = 1e-2) {
  Ra <- R_mu1 - R_mu2
  Rb <- R_sigma1 / R_sigma2

  if (abs(Rb - 1) < tol) {
    return((R_mu1 + R_mu2) / 2)
  }

  numerator <- R_mu2 * (Rb^2 - 1) -
               Ra +
               Rb * sqrt(Ra^2 + (Rb^2 - 1) * R_sigma2^2 * log(Rb^2))

  denominator <- Rb^2 - 1

  numerator / denominator
}

#' GPQ for the Youden index
#'
#' Computes the generalized pivotal quantity for the Youden index
#' (Equation 12 in Lai et al., 2012).
#'
#' @param R_mu1 Numeric. GPQ realization for the mean of the diseased group.
#' @param R_mu2 Numeric. GPQ realization for the mean of the non-diseased group.
#' @param R_sigma1 Numeric. GPQ realization for the SD of the diseased group.
#' @param R_sigma2 Numeric. GPQ realization for the SD of the non-diseased
#'   group.
#' @param Rc Numeric. GPQ realization for the optimal cut-point.
#'
#' @export
#' @return A numeric scalar in \eqn{(-1, 1)}.
#' @examples
#' gpq_J(R_mu1 = 2, R_mu2 = 0, R_sigma1 = 1, R_sigma2 = 1.2, Rc = 1)
gpq_J <- function(R_mu1, R_mu2, R_sigma1, R_sigma2, Rc) {
  pnorm((R_mu1 - Rc) / R_sigma1) +
    pnorm((Rc - R_mu2) / R_sigma2) - 1
}

#' Generate generalized pivotal quantities for one simulated sample
#'
#' Implements the GPQ algorithm of Lai et al. (2012) for a single dataset
#' drawn from two normal populations, producing Monte Carlo realizations of
#' the GPQs for the Youden index and its optimal cut-point.
#'
#' @param y1 Numeric vector. Observed sample from the diseased population.
#' @param y2 Numeric vector. Observed sample from the non-diseased population.
#' @param K Integer. Number of GPQ iterations (default: 250).
#' @param tol Numeric. Tolerance for the equal-variance case in
#'   \code{\link{gpq_c}} (default: 1e-2).
#'
#' @export
#' @return A list with two numeric vectors of length \code{K}:
#' \describe{
#'   \item{RJ}{GPQ realizations for the Youden index.}
#'   \item{Rc}{GPQ realizations for the optimal cut-point.}
#' }
#'
#' @references Lai, C.-Y., Tian, L., & Schisterman, E. F. (2012).
#'   Exact confidence interval estimation for the Youden index and its
#'   corresponding optimal cut-point.
#'   \emph{Computational Statistics & Data Analysis}, 56(5), 1103-1114.
#' @examples
#' set.seed(1)
#' y1 <- rnorm(30, mean = 2, sd = 1)
#' y2 <- rnorm(30, mean = 0, sd = 1)
#' gpq <- gpq_from_sample(y1, y2, K = 100)
#' quantile(gpq$RJ, probs = c(0.025, 0.975))
gpq_from_sample <- function(y1, y2, K = 250, tol = 1e-2) {
  n1 <- length(y1)
  n2 <- length(y2)

  ybar1 <- mean(y1)
  ybar2 <- mean(y2)

  s1   <- sd(y1)
  s2   <- sd(y2)
  s2_1 <- var(y1)
  s2_2 <- var(y2)

  RJ <- numeric(K)
  Rc <- numeric(K)

  for (k in seq_len(K)) {
    t1 <- rt(1, df = n1 - 1)
    t2 <- rt(1, df = n2 - 1)
    v1 <- rchisq(1, df = n1 - 1)
    v2 <- rchisq(1, df = n2 - 1)

    R_mu1    <- gpq_mu(ybar1, s1, n1, t1)
    R_mu2    <- gpq_mu(ybar2, s2, n2, t2)
    R_sigma1 <- gpq_sigma(s2_1, v1, n1)
    R_sigma2 <- gpq_sigma(s2_2, v2, n2)

    Rc[k] <- gpq_c(R_mu1, R_mu2, R_sigma1, R_sigma2, tol = tol)
    RJ[k] <- gpq_J(R_mu1, R_mu2, R_sigma1, R_sigma2, Rc[k])
  }

  list(RJ = RJ, Rc = Rc)
}
