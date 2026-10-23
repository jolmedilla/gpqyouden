#' Theoretical Youden index under the normal model
#'
#' Computes the true Youden index for two normal populations given their
#' means and standard deviations.
#'
#' @param parameters Numeric vector or matrix with four columns in the order
#'   \eqn{(\mu_2, \mu_1, \sigma_2, \sigma_1)}.  If a vector is supplied it is
#'   coerced to a one-row matrix.
#' @param model Character. Currently only \code{"norm"} is supported.
#'
#' @export
#' @return A numeric vector of Youden index values, one per row of
#'   \code{parameters}.
#'
#' @seealso \code{\link{find_mu1_from_J2}}
#' @examples
#' # (mu2, mu1, sigma2, sigma1)
#' J2.theoretical(c(0, 1.68, 1, 1), model = "norm")
J2.theoretical <- function(parameters, model) {
  if (model != "norm") stop("Only 'norm' model implemented")

  if (is.null(dim(parameters))) {
    parameters <- matrix(parameters, nrow = 1)
  }

  mu2    <- parameters[, 1]
  mu1    <- parameters[, 2]
  sigma2 <- parameters[, 3]
  sigma1 <- parameters[, 4]

  same_sigma <- !is.na(sigma2) & !is.na(sigma1) &
    (abs(sigma1 - sigma2) < 1e-10)

  c_val <- ifelse(
    same_sigma,
    (mu1 + mu2) / 2,
    {
      a <- mu1 - mu2
      b <- sigma1 / sigma2
      (mu2 * (b^2 - 1) - a +
         b * sqrt(a^2 + (b^2 - 1) * sigma2^2 * log(b^2))) / (b^2 - 1)
    }
  )

  pnorm((mu1 - c_val) / sigma1) +
    pnorm((c_val - mu2) / sigma2) - 1
}

#' Equation whose root gives \eqn{\mu_1}
#'
#' Auxiliary function used by \code{\link{find_mu1_from_J2}}: evaluates
#' \eqn{J(\mu_1) - J_{\text{target}}}.
#'
#' @param x Numeric. Candidate value of \eqn{\mu_1}.
#' @param parameters Numeric vector \eqn{(\mu_2, \sigma_2, \sigma_1)}.
#' @param J2 Numeric. Target Youden index value.
#' @param model Character. Currently only \code{"norm"} is supported.
#'
#' @export
#' @return A numeric scalar.
#' @examples
#' # parameters = (mu2, sigma2, sigma1)
#' eq.J2(x = 1.68, parameters = c(0, 1, 1), J2 = 0.6, model = "norm")
eq.J2 <- function(x, parameters, J2, model) {
  if (model != "norm") stop("Only 'norm' model implemented")
  params <- c(parameters[1], x, parameters[2:3])
  J2.theoretical(matrix(params, nrow = 1), model) - J2
}

#' Find \eqn{\mu_1} from the Youden index under the normal model
#'
#' Numerically solves for the mean of the diseased population \eqn{\mu_1}
#' such that the theoretical Youden index equals \code{J2}, given the
#' remaining parameters.
#'
#' @param J2 Numeric. Target Youden index value.
#' @param mu2 Numeric. Mean of the non-diseased population.
#' @param sigma2 Numeric. Standard deviation of the non-diseased population.
#' @param sigma1 Numeric. Standard deviation of the diseased population.
#' @param model Character. Currently only \code{"norm"} is supported.
#'
#' @export
#' @return A numeric scalar, or \code{NA_real_} if the root cannot be found.
#'
#' @seealso \code{\link{J2.theoretical}}, \code{\link{find_mu1_from_J_mixture}}
#' @examples
#' find_mu1_from_J2(J2 = 0.6, mu2 = 0, sigma2 = 1, sigma1 = 1, model = "norm")
find_mu1_from_J2 <- function(J2, mu2, sigma2, sigma1, model) {
  if (anyNA(c(J2, mu2, sigma2, sigma1))) return(NA_real_)
  if (sigma2 <= 0 || sigma1 <= 0)        return(NA_real_)

  res <- tryCatch(
    uniroot(
      eq.J2,
      interval   = c(mu2 - 10 * sigma1, mu2 + 10 * sigma1),
      parameters = c(mu2, sigma2, sigma1),
      J2         = J2,
      model      = model,
      tol        = 1e-6
    ),
    error = function(e) NULL
  )

  if (is.null(res)) NA_real_ else res$root
}

#' Find \eqn{\mu_1} from the Youden index under the mixture model
#'
#' Numerically solves for \eqn{\mu_1} such that the true Youden index of two
#' contaminated normal mixture distributions equals \code{J_target}.
#'
#' @param J_target Numeric. Target Youden index value (must be in (0, 1)).
#' @param mu2 Numeric. Mean of the non-diseased population.
#' @param sigma1 Numeric. Base standard deviation of the diseased population.
#' @param sigma2 Numeric. Base standard deviation of the non-diseased
#'   population.
#' @param p Numeric. Contamination probability (default: 0.1).
#' @param tol Numeric. Root-finding tolerance (default: 1e-6).
#' @param max_expand Integer. Maximum number of bracket-expansion iterations
#'   (default: 60).
#'
#' @export
#' @return A numeric scalar, or \code{NA_real_} if the root cannot be found.
#'
#' @seealso \code{\link{find_mu1_from_J2}}, \code{\link{J_true_mixture}}
#' @examples
#' find_mu1_from_J_mixture(J_target = 0.6, mu2 = 0, sigma1 = 1, sigma2 = 1, p = 0.1)
find_mu1_from_J_mixture <- function(J_target, mu2, sigma1, sigma2, p = 0.1,
                                    tol = 1e-6, max_expand = 60) {
  if (anyNA(c(J_target, mu2, sigma1, sigma2))) return(NA_real_)
  if (sigma1 <= 0 || sigma2 <= 0)              return(NA_real_)
  if (J_target <= 0 || J_target >= 1)          return(NA_real_)

  g    <- function(mu1) J_true_mixture(mu1, mu2, sigma1, sigma2, p)$J - J_target
  lo   <- mu2 + 1e-6
  hi   <- mu2 + sigma1 + sigma2 + 1
  g_lo <- g(lo)
  g_hi <- g(hi)

  iter <- 0
  while (is.finite(g_hi) && (g_lo * g_hi > 0) && iter < max_expand) {
    hi   <- hi * 1.5 + 1
    g_hi <- g(hi)
    iter <- iter + 1
  }

  if (!is.finite(g_lo) || !is.finite(g_hi) || (g_lo * g_hi > 0)) {
    return(NA_real_)
  }

  res <- tryCatch(
    uniroot(function(x) g(x), interval = c(lo, hi), tol = tol),
    error = function(e) NULL
  )

  if (is.null(res)) NA_real_ else res$root
}

#' Compute \eqn{\mu_1} for each row of a parameter data frame
#'
#' Convenience wrapper that applies \code{\link{find_mu1_from_J2}} row-wise
#' to a data frame containing columns \code{mu2}, \code{var1}, \code{var2}
#' and \code{J}.
#'
#' @param df A data frame with columns \code{mu2}, \code{var1}, \code{var2}
#'   and \code{J}.
#'
#' @export
#' @return A data frame with columns \code{mu1}, \code{original_J} and
#'   \code{theoretical_J}.
#' @examples
#' df <- data.frame(mu2 = 0, var1 = 1, var2 = 1, J = c(0.4, 0.6))
#' calculate_mu1_for_dataframe(df)
calculate_mu1_for_dataframe <- function(df) {
  mu2    <- df$mu2
  sigma2 <- sqrt(df$var2)
  sigma1 <- sqrt(df$var1)
  J_orig <- df$J
  model  <- "norm"

  mu1 <- mapply(find_mu1_from_J2, J_orig, mu2, sigma2, sigma1,
                MoreArgs = list(model = model))

  J_theoretical <- mapply(
    function(mu2, mu1, sigma2, sigma1) {
      J2.theoretical(c(mu2, mu1, sigma2, sigma1), model)
    },
    mu2, mu1, sigma2, sigma1
  )

  data.frame(mu1 = mu1, original_J = J_orig, theoretical_J = J_theoretical)
}
