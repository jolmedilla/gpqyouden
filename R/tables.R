# Inject a short caption into a kableExtra LaTeX string.
# Replaces the first \caption{ with \caption[short]{ so the List of Tables
# entry is concise.  Only the first occurrence is replaced; continuation-page
# captions (already \caption[]{...}) are left untouched.
inject_short_caption <- function(latex_str, short) {
  # sub() replacement strings treat \\ as one literal backslash, so each
  # LaTeX backslash in `short` must be doubled before it enters the replacement.
  short_esc <- gsub("\\\\", "\\\\\\\\", short)
  sub("\\\\caption\\{",
      paste0("\\\\caption[", short_esc, "]{"),
      latex_str)
}

auc_from_params <- function(mu2, var1, var2, J) {
  mu1 <- find_mu1_from_J2(J2 = J, mu2 = mu2, sigma2 = sqrt(var2),
                           sigma1 = sqrt(var1), model = "norm")
  pnorm((mu1 - mu2) / sqrt(var1 + var2))
}

#' Generate LaTeX tables for the normal model simulation
#'
#' Produces two long-format tables (coverage probability and average length,
#' one for the Youden index and one for the cut-point), with \eqn{\sigma_1^2}
#' as the first grouping column.  These correspond to Tables 1 and 2 of
#' Lai et al. (2012).
#'
#' @param results_df Data frame as returned by \code{run_parallel_grid} with
#'   \code{generator = "normal"}.
#' @param out_dir Character. Directory where the \code{.tex} files are written.
#'
#' @export
#' @return Invisibly \code{NULL}. Writes \code{table_J_normal.tex} and
#'   \code{table_c_normal.tex} to \code{out_dir}.
#' @examples
#' \donttest{
#' grid <- default_param_grid()[1:2, ]
#' results <- run_parallel_grid(grid, n_mc = 10, K = 30)
#' generate_normal_tables(results, out_dir = tempdir())
#' }
generate_normal_tables <- function(results_df, out_dir) {

  # Table 1 — J under the normal model
  table_J <- results_df |>
    dplyr::arrange(var1, n1, n2, J) |>
    dplyr::mutate(
      AUC        = round(mapply(auc_from_params, mu2, var1, var2, J), 3),
      coverage_J = round(coverage_J, 4),
      length_J   = round(length_J,   4)
    ) |>
    dplyr::select(var1, n1, n2, J, AUC, coverage_J, length_J)

  colnames(table_J) <- c(
    "$\\sigma_1^2$", "$n_1$", "$n_2$", "$J$", "AUC", "Coverage", "Length"
  )

  latex_J <- kableExtra::kable(
    table_J,
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    escape    = FALSE,
    align     = "c",
    label     = "j_normal",
    caption   = "Coverage probabilities and mean lengths of the 95\\% confidence interval for the Youden index $J$, with the theoretical AUC for each parameter scenario."
  ) |>
    kableExtra::collapse_rows(
      columns                  = 1:3,
      valign                   = "top",
      row_group_label_position = "first"
    ) |>
    kableExtra::kable_styling(latex_options = "repeat_header") |>
    inject_short_caption("95\\% CI coverage and length for $J$")

  writeLines(latex_J, con = file.path(out_dir, "table_J_normal.tex"))

  # Table 2 — c under the normal model
  table_c <- results_df |>
    dplyr::arrange(var1, n1, n2, J) |>
    dplyr::mutate(
      AUC        = round(mapply(auc_from_params, mu2, var1, var2, J), 3),
      coverage_c = round(coverage_c, 4),
      length_c   = round(length_c,   4)
    ) |>
    dplyr::select(var1, n1, n2, J, AUC, coverage_c, length_c)

  colnames(table_c) <- c(
    "$\\sigma_1^2$", "$n_1$", "$n_2$", "$J$", "AUC", "Coverage", "Length"
  )

  latex_c <- kableExtra::kable(
    table_c,
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    escape    = FALSE,
    align     = "c",
    label     = "c_normal",
    caption   = "Coverage probabilities and mean lengths of the 95\\% confidence interval for the optimal cut-point $c^*$, with the theoretical AUC for each parameter scenario."
  ) |>
    kableExtra::collapse_rows(
      columns                  = 1:3,
      valign                   = "top",
      row_group_label_position = "first"
    ) |>
    kableExtra::kable_styling(latex_options = "repeat_header") |>
    inject_short_caption("Id.\\ for $c^*$")

  writeLines(latex_c, con = file.path(out_dir, "table_c_normal.tex"))

  invisible(NULL)
}

#' Generate LaTeX tables for the mixture model simulation
#'
#' Produces two wide-format tables (coverage probability only, with
#' \eqn{\sigma_1^2} as column headers), one for the Youden index and one for
#' the cut-point.  These correspond to Tables 3 and 4 of Lai et al. (2012)
#' for \eqn{p = 0.1}, and to the analogous tables for other values of
#' \eqn{p}.
#'
#' @param results_df Data frame as returned by \code{run_parallel_grid} with
#'   \code{generator = "mixture"}.
#' @param p Numeric. Contamination probability used in the simulation
#'   (included in captions and file names).
#' @param out_dir Character. Directory where the \code{.tex} files are written.
#'
#' @export
#' @return Invisibly \code{NULL}. Writes
#'   \code{table_J_mixture_p=<p>.tex} and
#'   \code{table_c_mixture_p=<p>.tex} to \code{out_dir}.
#' @examples
#' \donttest{
#' grid <- default_param_grid()[1:2, ]
#' results <- run_parallel_grid(grid, generator = "mixture",
#'                               n_mc = 10, K = 30, p = 0.1)
#' generate_mixture_tables(results, p = 0.1, out_dir = tempdir())
#' }
generate_mixture_tables <- function(results_df, p, out_dir) {

  # Table 3 (or equivalent) — J under the mixture model
  table_J_mix <- results_df |>
    dplyr::arrange(n1, n2, J, var1) |>
    dplyr::select(n1, n2, J, var1, coverage_J) |>
    dplyr::mutate(coverage_J = round(coverage_J, 4)) |>
    tidyr::pivot_wider(names_from = var1, values_from = coverage_J)

  colnames(table_J_mix) <- c("$n_1$", "$n_2$", "$J$", "0.5", "1", "3", "5")

  latex_J_mix <- kableExtra::kable(
    table_J_mix,
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    escape    = FALSE,
    align     = "c",
    label     = paste0("j_mixture_", sprintf("%02d", as.integer(round(p * 10)))),
    caption   = paste0(
      "Coverage probabilities of the 95\\% confidence interval for $J$ ",
      "under the mixture model ($p = ", p, "$)"
    )
  ) |>
    kableExtra::add_header_above(
      c(" " = 3, "$\\\\sigma_{1}^{2}$" = 4),
      escape = FALSE, align = "l"
    ) |>
    kableExtra::row_spec(0, extra_latex_after = "\\cmidrule(lr){4-7}") |>
    kableExtra::collapse_rows(
      columns                  = 1:3,
      valign                   = "top",
      row_group_label_position = "first"
    ) |>
    kableExtra::kable_styling(latex_options = "repeat_header") |>
    inject_short_caption(paste0("Id.\\ for $J$ under mixture model ($p = ", p, "$)"))

  writeLines(
    latex_J_mix,
    con = file.path(out_dir, paste0("table_J_mixture_p=", p, ".tex"))
  )

  # Table 4 (or equivalent) — c under the mixture model
  table_c_mix <- results_df |>
    dplyr::arrange(n1, n2, J, var1) |>
    dplyr::select(n1, n2, J, var1, coverage_c) |>
    dplyr::mutate(coverage_c = round(coverage_c, 4)) |>
    tidyr::pivot_wider(names_from = var1, values_from = coverage_c)

  colnames(table_c_mix) <- c("$n_1$", "$n_2$", "$J$", "0.5", "1", "3", "5")

  latex_c_mix <- kableExtra::kable(
    table_c_mix,
    format    = "latex",
    booktabs  = TRUE,
    longtable = TRUE,
    escape    = FALSE,
    align     = "c",
    label     = paste0("c_mixture_", sprintf("%02d", as.integer(round(p * 10)))),
    caption   = paste0(
      "Coverage probabilities of the 95\\% confidence interval for cut-point $c^*$ ",
      "under the mixture model ($p = ", p, "$)"
    )
  ) |>
    kableExtra::add_header_above(
      c(" " = 3, "$\\\\sigma_{1}^{2}$" = 4),
      escape = FALSE, align = "l"
    ) |>
    kableExtra::row_spec(0, extra_latex_after = "\\cmidrule(lr){4-7}") |>
    kableExtra::collapse_rows(
      columns                  = 1:3,
      valign                   = "top",
      row_group_label_position = "first"
    ) |>
    kableExtra::kable_styling(latex_options = "repeat_header") |>
    inject_short_caption(paste0("Id.\\ for $c^*$ under mixture model ($p = ", p, "$)"))

  writeLines(
    latex_c_mix,
    con = file.path(out_dir, paste0("table_c_mixture_p=", p, ".tex"))
  )

  invisible(NULL)
}
