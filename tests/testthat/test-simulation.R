test_that("default_param_grid reproduces the 160-scenario design of Lai et al.", {
  grid <- default_param_grid()

  expect_equal(nrow(grid), 160)
  expect_named(grid, c("mu2", "sigma1", "sigma2", "J", "n1", "n2", "var1", "var2"))
  expect_equal(sort(unique(grid$var1)), c(0.5, 1, 3, 5))
  expect_equal(sort(unique(grid$J)), c(0.2, 0.4, 0.6, 0.8, 0.9))
  expect_equal(nrow(unique(grid[, c("n1", "n2")])), 8)
})

test_that("generate_param_grid builds the full factorial of its inputs", {
  n_pairs <- matrix(c(10, 10, 20, 20), ncol = 2, byrow = TRUE)
  grid <- generate_param_grid(var1 = c(1, 2), var2 = 1, mu2 = 0,
                               J = c(0.5, 0.7), n_pairs = n_pairs)

  # 2 var1 * 1 var2 * 1 mu2 * 2 J * 2 n_pairs = 8 rows
  expect_equal(nrow(grid), 8)
  expect_equal(grid$sigma1, sqrt(grid$var1))
  expect_equal(grid$sigma2, sqrt(grid$var2))
})

test_that("run_simulation returns coverage/length summaries in sensible ranges", {
  set.seed(1)
  res <- run_simulation(mu2 = 0, sigma1 = 1, sigma2 = 1, J = 0.6,
                         n1 = 20, n2 = 20, n_mc = 20, K = 50)

  expect_named(res, c("params", "J", "c"))
  expect_true(res$params$mu1 > res$params$mu2)
  expect_true(res$J$coverage >= 0 && res$J$coverage <= 1)
  expect_true(res$c$coverage >= 0 && res$c$coverage <= 1)
  expect_gt(res$J$mean_length, 0)
  expect_gt(res$c$mean_length, 0)
})

test_that("run_simulation works for both the normal and mixture generators", {
  set.seed(1)
  res_mix <- run_simulation(mu2 = 0, sigma1 = 1, sigma2 = 1, J = 0.6,
                             n1 = 20, n2 = 20, n_mc = 10, K = 30,
                             generator = "mixture", p = 0.1)
  expect_true(res_mix$J$coverage >= 0 && res_mix$J$coverage <= 1)
})

test_that("run_simulation rejects an unsolvable scenario (non-finite mu1)", {
  expect_error(
    run_simulation(mu2 = 0, sigma1 = -1, sigma2 = 1, J = 0.6,
                    n1 = 10, n2 = 10, n_mc = 5, K = 20),
    "mu1"
  )
})

test_that("run_one_row wraps run_simulation into a one-row summary data frame", {
  row <- default_param_grid()[1, ]
  set.seed(1)
  out <- run_one_row(row, n_mc = 10, K = 30)

  expect_equal(nrow(out), 1)
  expect_named(out, c("mu2", "var1", "var2", "J", "n1", "n2",
                       "coverage_J", "length_J", "coverage_c", "length_c",
                       "generator", "p"))
  expect_equal(out$J, row$J)
  expect_equal(out$n1, row$n1)
})

test_that("run_parallel_grid returns one row per scenario and respects n_cores", {
  grid <- default_param_grid()[1:3, ]
  out <- run_parallel_grid(grid, n_mc = 10, K = 30, n_cores = 1)

  expect_equal(nrow(out), 3)
  expect_equal(out$J, grid$J)
  expect_true(all(out$coverage_J >= 0 & out$coverage_J <= 1))
})

test_that("run_parallel_grid honours an explicit, CRAN-safe n_cores", {
  grid <- default_param_grid()[1:2, ]
  out <- run_parallel_grid(grid, n_mc = 5, K = 20, n_cores = 2L)
  expect_equal(nrow(out), 2)
})
