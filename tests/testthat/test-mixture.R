test_that("mix_sd returns the two component standard deviations", {
  sds <- mix_sd(sigma = 1)
  expect_equal(unname(sds["sd_a"]), sqrt((1 / 1.1) * 1^2))
  expect_equal(unname(sds["sd_b"]), sqrt((2 / 1.1) * 1^2))
  expect_true(sds["sd_b"] > sds["sd_a"])
})

test_that("dmixnorm is a (1-p)/p mixture of the two component densities", {
  sds <- mix_sd(sigma = 1)
  expected <- 0.9 * dnorm(0, 0, sds["sd_a"]) + 0.1 * dnorm(0, 0, sds["sd_b"])
  expect_equal(unname(dmixnorm(x = 0, mu = 0, sigma = 1, p = 0.1)),
               unname(expected))
  expect_true(all(dmixnorm(seq(-3, 3), mu = 0, sigma = 1, p = 0.1) >= 0))
})

test_that("pmixnorm is a valid CDF: monotone, in [0, 1], tends to 0/1", {
  xs <- seq(-5, 5, by = 0.5)
  ps <- pmixnorm(xs, mu = 0, sigma = 1, p = 0.3)

  expect_true(all(diff(ps) >= 0))
  expect_true(all(ps >= 0 & ps <= 1))
  expect_lt(pmixnorm(-20, mu = 0, sigma = 1, p = 0.1), 1e-6)
  expect_gt(pmixnorm(20, mu = 0, sigma = 1, p = 0.1), 1 - 1e-6)
})

test_that("generate_mixture_sample returns n finite draws", {
  set.seed(1)
  y <- generate_mixture_sample(n = 200, mu = 0, sigma = 1, p = 0.1)

  expect_length(y, 200)
  expect_true(all(is.finite(y)))
  # Loose sanity check on location, not a strict statistical test.
  expect_lt(abs(mean(y)), 1)
})

test_that("generate_mixture_sample with p = 0 behaves like a plain normal draw", {
  set.seed(1)
  y <- generate_mixture_sample(n = 500, mu = 2, sigma = 1, p = 0)
  expect_equal(mean(y), 2, tolerance = 0.2)
})

test_that("J_true_mixture returns a Youden index in (0, 1) and a finite cut-point", {
  res <- J_true_mixture(mu1 = 1.68, mu2 = 0, sigma1 = 1, sigma2 = 1, p = 0.1)

  expect_named(res, c("J", "c"))
  expect_true(res$J > 0 && res$J < 1)
  expect_true(is.finite(res$c))
})

test_that("J_true_mixture matches pmixnorm at its own optimum", {
  res <- J_true_mixture(mu1 = 1.68, mu2 = 0, sigma1 = 1, sigma2 = 1, p = 0.1)
  expected_J <- pmixnorm(res$c, 0, 1, 0.1) - pmixnorm(res$c, 1.68, 1, 0.1)
  expect_equal(res$J, expected_J, tolerance = 1e-6)
})
