#' Build the default parameter grid used in Lai et al. (2012)
#'
#' Returns the full factorial combination of variances, Youden index values,
#' and sample-size pairs used in the simulation study of Lai et al. (2012):
#' \eqn{\sigma_1^2 \in \{0.5, 1, 3, 5\}}, \eqn{\sigma_2^2 = 1},
#' \eqn{\mu_2 = 0}, \eqn{J \in \{0.2, 0.4, 0.6, 0.8, 0.9\}}, and eight
#' \eqn{(n_1, n_2)} pairs.
#'
#' @export
#' @return A data frame with 160 rows and columns \code{mu2}, \code{sigma1},
#'   \code{sigma2}, \code{J}, \code{n1}, \code{n2}, \code{var1}, \code{var2}.
#'
#' @seealso \code{\link{generate_param_grid}}, \code{\link{run_parallel_grid}}
#' @examples
#' grid <- default_param_grid()
#' dim(grid)
#' head(grid)
default_param_grid <- function() {
  var1 <- c(0.5, 1, 3, 5)
  var2 <- c(1)
  mu2  <- c(0)
  J    <- c(0.2, 0.4, 0.6, 0.8, 0.9)

  n_pairs <- matrix(
    c(
      10,  10,
      20,  20,
      30,  30,
      10,  30,
      30,  20,
      50,  20,
      50,  50,
     100, 100
    ),
    ncol = 2, byrow = TRUE
  )

  generate_param_grid(var1, var2, mu2, J, n_pairs)
}

#' Build a custom parameter grid
#'
#' Constructs a data frame of all combinations of the supplied parameter
#' vectors and sample-size pairs, adding derived standard deviation columns.
#'
#' @param var1 Numeric vector. Variances for the diseased population.
#' @param var2 Numeric vector. Variances for the non-diseased population.
#' @param mu2 Numeric vector. Means for the non-diseased population.
#' @param J Numeric vector. Target Youden index values.
#' @param n_pairs Integer matrix with two columns (\code{n1}, \code{n2}).
#'
#' @export
#' @return A data frame with columns \code{mu2}, \code{sigma1}, \code{sigma2},
#'   \code{J}, \code{n1}, \code{n2}, \code{var1}, \code{var2}.
#'
#' @seealso \code{\link{default_param_grid}}
#' @examples
#' n_pairs <- matrix(c(10, 10, 20, 20), ncol = 2, byrow = TRUE)
#' generate_param_grid(var1 = c(1, 2), var2 = 1, mu2 = 0,
#'                      J = c(0.5, 0.7), n_pairs = n_pairs)
generate_param_grid <- function(var1, var2, mu2, J, n_pairs) {
  base_grid <- expand.grid(
    var1 = var1, var2 = var2, mu2 = mu2, J = J,
    KEEP.OUT.ATTRS = FALSE, stringsAsFactors = FALSE
  )

  df_params <- do.call(
    rbind,
    lapply(seq_len(nrow(n_pairs)), function(i) {
      cbind(base_grid, n1 = n_pairs[i, 1], n2 = n_pairs[i, 2])
    })
  )

  df_params$sigma1 <- sqrt(df_params$var1)
  df_params$sigma2 <- sqrt(df_params$var2)
  df_params <- df_params[, c("mu2", "sigma1", "sigma2", "J", "n1", "n2", "var1", "var2")]
  rownames(df_params) <- NULL
  df_params
}

#' Run the Monte Carlo simulation for a single parameter scenario
#'
#' Repeatedly draws paired samples under the normal or contaminated-mixture
#' model, computes GPQ-based confidence intervals for the Youden index and
#' its optimal cut-point via \code{\link{gpq_from_sample}}, and summarises
#' empirical coverage and mean interval length over \code{n_mc} replications.
#'
#' @param mu2 Numeric. Mean of the non-diseased population.
#' @param sigma1 Numeric. Standard deviation of the diseased population.
#' @param sigma2 Numeric. Standard deviation of the non-diseased population.
#' @param J Numeric. Target (nominal) Youden index.
#' @param n1 Integer. Sample size for the diseased population.
#' @param n2 Integer. Sample size for the non-diseased population.
#' @param n_mc Integer. Number of Monte Carlo replications (default: 200).
#' @param K Integer. Number of GPQ iterations per replication, passed to
#'   \code{\link{gpq_from_sample}} (default: 250).
#' @param alpha Numeric. Significance level for the confidence intervals
#'   (default: 0.05, giving 95\% intervals).
#' @param generator Character. \code{"normal"} or \code{"mixture"}.
#' @param p Numeric. Contamination probability, used only when
#'   \code{generator = "mixture"} (default: 0.1).
#'
#' @export
#' @return A list with elements \code{params} (the scenario's parameters,
#'   including the derived \code{mu1}), \code{J} and \code{c}, each a list
#'   with \code{coverage} and \code{mean_length}.
#'
#' @seealso \code{\link{run_one_row}}, \code{\link{run_parallel_grid}}
#' @examples
#' set.seed(1)
#' run_simulation(mu2 = 0, sigma1 = 1, sigma2 = 1, J = 0.6,
#'                 n1 = 20, n2 = 20, n_mc = 20, K = 50)
run_simulation <- function(mu2, sigma1, sigma2, J, n1, n2,
                           n_mc = 200, K = 250, alpha = 0.05,
                           generator = c("normal", "mixture"),
                           p = 0.1) {
  generator <- match.arg(generator)

  if (generator == "mixture") {
    mu1 <- find_mu1_from_J2(J2 = J, mu2 = mu2, sigma2 = sigma2,
                            sigma1 = sigma1, model = "norm")
  } else {
    mu1 <- find_mu1_from_J2(J2 = J, mu2 = mu2, sigma2 = sigma2,
                            sigma1 = sigma1, model = "norm")
  }

  if (!is.finite(mu1)) stop("mu1 is not finite: check parameters of scenario")
  if (mu1 <= mu2)      stop("mu1 <= mu2: invalid scenario for Youden Index")

  if (generator == "mixture") {
    J_true  <- J
    c_true  <- J_true_mixture(mu1, mu2, sigma1, sigma2, p)$c
  } else {
    J_true <- J
    a <- mu1 - mu2
    b <- sigma1 / sigma2
    c_true <- if (abs(b - 1) < 1e-10) {
      (mu1 + mu2) / 2
    } else {
      (mu2 * (b^2 - 1) - a +
         b * sqrt(a^2 + (b^2 - 1) * sigma2^2 * log(b^2))) / (b^2 - 1)
    }
  }

  cover_J <- logical(n_mc)
  len_J   <- numeric(n_mc)
  cover_c <- logical(n_mc)
  len_c   <- numeric(n_mc)

  for (m in seq_len(n_mc)) {
    if (generator == "normal") {
      y1 <- rnorm(n1, mu1, sigma1)
      y2 <- rnorm(n2, mu2, sigma2)
    } else {
      y1 <- generate_mixture_sample(n1, mu1, sigma1, p)
      y2 <- generate_mixture_sample(n2, mu2, sigma2, p)
    }

    gpq <- gpq_from_sample(y1, y2, K = K)

    ci_J <- quantile(gpq$RJ, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)
    ci_c <- quantile(gpq$Rc, probs = c(alpha / 2, 1 - alpha / 2), na.rm = TRUE)

    cover_J[m] <- J_true  >= ci_J[1] && J_true  <= ci_J[2]
    cover_c[m] <- c_true  >= ci_c[1] && c_true  <= ci_c[2]
    len_J[m]   <- diff(ci_J)
    len_c[m]   <- diff(ci_c)
  }

  list(
    params = list(mu1 = mu1, mu2 = mu2, sigma1 = sigma1, sigma2 = sigma2,
                  J = J, n1 = n1, n2 = n2),
    J = list(coverage = mean(cover_J), mean_length = mean(len_J)),
    c = list(coverage = mean(cover_c), mean_length = mean(len_c))
  )
}

#' Run the simulation for a single row of the parameter grid
#'
#' Convenience wrapper around \code{\link{run_simulation}} that accepts a
#' single-row data frame (as returned by subsetting the output of
#' \code{\link{default_param_grid}}) and returns a one-row summary data frame.
#'
#' @param row A one-row data frame with columns \code{mu2}, \code{sigma1},
#'   \code{sigma2}, \code{J}, \code{n1}, \code{n2}, \code{var1}, \code{var2}.
#' @param generator Character. \code{"normal"} or \code{"mixture"}.
#' @param n_mc Integer. Number of Monte Carlo replications (default: 200).
#' @param K Integer. Number of GPQ iterations per replication (default: 250).
#' @param p Numeric. Contamination probability. Ignored when
#'   \code{generator = "normal"} (default: 0.1).
#'
#' @export
#' @return A one-row data frame with columns \code{mu2}, \code{var1},
#'   \code{var2}, \code{J}, \code{n1}, \code{n2}, \code{coverage_J},
#'   \code{length_J}, \code{coverage_c}, \code{length_c}, \code{generator},
#'   \code{p}.
#'
#' @seealso \code{\link{run_simulation}}, \code{\link{run_parallel_grid}}
#' @examples
#' row <- default_param_grid()[1, ]
#' set.seed(1)
#' run_one_row(row, n_mc = 10, K = 30)
run_one_row <- function(row, generator = c("normal", "mixture"),
                        n_mc = 200, K = 250, p = 0.1) {
  generator <- match.arg(generator)
  res <- run_simulation(
    mu2    = row$mu2,
    sigma1 = row$sigma1,
    sigma2 = row$sigma2,
    J      = row$J,
    n1     = row$n1,
    n2     = row$n2,
    n_mc   = n_mc,
    K      = K,
    generator = generator,
    p      = p
  )

  data.frame(
    mu2        = row$mu2,
    var1       = row$var1,
    var2       = row$var2,
    J          = row$J,
    n1         = row$n1,
    n2         = row$n2,
    coverage_J = res$J$coverage,
    length_J   = res$J$mean_length,
    coverage_c = res$c$coverage,
    length_c   = res$c$mean_length,
    generator  = generator,
    p          = p
  )
}

#' Run the simulation in parallel, detecting the OS automatically
#'
#' Uses \code{mclapply} on Unix/macOS and \code{doParallel}+\code{foreach}
#' on Windows.
#'
#' @param df_params Data frame of parameter combinations (from
#'   \code{default_param_grid} or \code{generate_param_grid}).
#' @param generator Character. \code{"normal"} or \code{"mixture"}.
#' @param n_mc Integer. Number of Monte Carlo replications.
#' @param K Integer. Number of GPQ iterations per replication.
#' @param p Numeric. Mixture proportion. Ignored when
#'   \code{generator = "normal"}; passed through to \code{run_one_row} but
#'   has no effect on the normal simulation.
#' @param n_cores Integer. Number of parallel workers to use. Defaults to
#'   \code{getOption("mc.cores", 2L)}, capped at the number of cores actually
#'   detected; pass a higher value explicitly to use more (CRAN policy caps
#'   the *default* number of cores a package may use at 2, but the caller may
#'   always opt into more).
#'
#' @export
#' @return A data frame with one row per parameter combination.
#' @examples
#' grid <- default_param_grid()[1:2, ]
#' run_parallel_grid(grid, n_mc = 10, K = 30)
run_parallel_grid <- function(df_params,
                              generator = c("normal", "mixture"),
                              n_mc = 200, K = 250, p = 0.1,
                              n_cores = getOption("mc.cores", 2L)) {
  generator <- match.arg(generator)
  n_cores   <- max(1L, min(n_cores, parallel::detectCores()))

  if (.Platform$OS.type == "windows") {
    if (!requireNamespace("doParallel", quietly = TRUE) ||
        !requireNamespace("foreach", quietly = TRUE)) {
      stop(
        "Packages 'doParallel' and 'foreach' are required on Windows. ",
        "Install them with: install.packages(c('doParallel', 'foreach'))"
      )
    }
    cl <- parallel::makeCluster(n_cores)
    doParallel::registerDoParallel(cl)
    parallel::clusterEvalQ(cl, library(gpqyouden))
    on.exit(parallel::stopCluster(cl), add = TRUE)

    # foreach/doParallel are optional (Suggests, Windows-only branch), so
    # %dopar% cannot be imported at the NAMESPACE level without forcing a
    # hard dependency; bind it locally instead once we know it is available.
    `%dopar%` <- foreach::`%dopar%`

    results_df <- foreach::foreach(
      i        = seq_len(nrow(df_params)),
      .combine = rbind
    ) %dopar% {
      run_one_row(df_params[i, ], generator = generator, n_mc = n_mc,
                  K = K, p = p)
    }
  } else {
    results_list <- parallel::mclapply(
      seq_len(nrow(df_params)),
      function(i) run_one_row(df_params[i, ], generator = generator,
                              n_mc = n_mc, K = K, p = p),
      mc.cores = n_cores
    )
    results_df <- do.call(rbind, results_list)
  }

  results_df
}
