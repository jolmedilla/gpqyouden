test_that("inject_short_caption replaces only the first \\caption{", {
  latex <- "\\caption{Long one}...\\caption{Long two}"
  out <- inject_short_caption(latex, "Short")
  expect_equal(out, "\\caption[Short]{Long one}...\\caption{Long two}")
})

test_that("inject_short_caption escapes backslashes in the short caption", {
  out <- inject_short_caption("\\caption{x}", "Id.\\ for $J$")
  expect_equal(out, "\\caption[Id.\\ for $J$]{x}")
})

test_that("auc_from_params matches the theoretical AUC formula", {
  auc <- auc_from_params(mu2 = 0, var1 = 1, var2 = 1, J = 0.6)
  mu1 <- find_mu1_from_J2(J2 = 0.6, mu2 = 0, sigma2 = 1, sigma1 = 1, model = "norm")
  expect_equal(unname(auc), pnorm((mu1 - 0) / sqrt(1 + 1)))
})

test_that("generate_normal_tables writes two non-empty .tex files", {
  grid <- default_param_grid()[1:2, ]
  results <- run_parallel_grid(grid, n_mc = 10, K = 30, n_cores = 1)

  out_dir <- tempfile("gpqyouden-normal-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  generate_normal_tables(results, out_dir = out_dir)

  j_file <- file.path(out_dir, "table_J_normal.tex")
  c_file <- file.path(out_dir, "table_c_normal.tex")
  expect_true(file.exists(j_file))
  expect_true(file.exists(c_file))

  j_lines <- readLines(j_file)
  expect_true(any(grepl("\\\\caption", j_lines)))
  expect_true(any(grepl("\\\\begin\\{longtable\\}", j_lines)))
})

test_that("generate_mixture_tables writes two non-empty, p-labelled .tex files", {
  grid <- default_param_grid()[1:2, ]
  results <- run_parallel_grid(grid, generator = "mixture",
                                n_mc = 10, K = 30, p = 0.1, n_cores = 1)

  out_dir <- tempfile("gpqyouden-mixture-")
  dir.create(out_dir)
  on.exit(unlink(out_dir, recursive = TRUE), add = TRUE)
  generate_mixture_tables(results, p = 0.1, out_dir = out_dir)

  j_file <- file.path(out_dir, "table_J_mixture_p=0.1.tex")
  c_file <- file.path(out_dir, "table_c_mixture_p=0.1.tex")
  expect_true(file.exists(j_file))
  expect_true(file.exists(c_file))
  expect_gt(length(readLines(j_file)), 0)
})
