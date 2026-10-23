test_that("J2.theoretical matches the closed-form pnorm formula", {
  # (mu2, mu1, sigma2, sigma1), equal sigmas -> cut-point is the midpoint.
  J <- J2.theoretical(c(0, 1.68, 1, 1), model = "norm")
  c_val <- (1.68 + 0) / 2
  expected <- pnorm((1.68 - c_val) / 1) + pnorm((c_val - 0) / 1) - 1
  expect_equal(unname(J), expected)
})

test_that("J2.theoretical vectorises over rows of a parameter matrix", {
  params <- rbind(
    c(0, 1.68, 1, 1),
    c(0, 1.0,  1, 1)
  )
  J <- J2.theoretical(params, model = "norm")
  expect_length(J, 2)
  expect_true(J[1] > J[2]) # larger mu1 separation -> larger Youden index
})

test_that("J2.theoretical rejects unsupported models", {
  expect_error(J2.theoretical(c(0, 1, 1, 1), model = "lognorm"))
})

test_that("eq.J2 is J2.theoretical(mu1) minus the target", {
  val <- eq.J2(x = 1.68, parameters = c(0, 1, 1), J2 = 0.6, model = "norm")
  expected <- J2.theoretical(matrix(c(0, 1.68, 1, 1), nrow = 1), "norm") - 0.6
  expect_equal(val, unname(expected))
})

test_that("find_mu1_from_J2 round-trips through J2.theoretical", {
  mu1 <- find_mu1_from_J2(J2 = 0.6, mu2 = 0, sigma2 = 1, sigma1 = 1, model = "norm")
  J_back <- J2.theoretical(c(0, mu1, 1, 1), model = "norm")
  expect_equal(unname(J_back), 0.6, tolerance = 1e-5)
})

test_that("find_mu1_from_J2 returns NA for invalid inputs instead of erroring", {
  expect_true(is.na(find_mu1_from_J2(J2 = NA, mu2 = 0, sigma2 = 1, sigma1 = 1, model = "norm")))
  expect_true(is.na(find_mu1_from_J2(J2 = 0.6, mu2 = 0, sigma2 = -1, sigma1 = 1, model = "norm")))
})

test_that("find_mu1_from_J_mixture round-trips through J_true_mixture", {
  mu1 <- find_mu1_from_J_mixture(J_target = 0.6, mu2 = 0, sigma1 = 1, sigma2 = 1, p = 0.1)
  expect_true(is.finite(mu1))

  J_back <- J_true_mixture(mu1 = mu1, mu2 = 0, sigma1 = 1, sigma2 = 1, p = 0.1)$J
  expect_equal(J_back, 0.6, tolerance = 1e-4)
})

test_that("find_mu1_from_J_mixture returns NA for an out-of-range target", {
  expect_true(is.na(find_mu1_from_J_mixture(J_target = 1.5, mu2 = 0, sigma1 = 1, sigma2 = 1)))
  expect_true(is.na(find_mu1_from_J_mixture(J_target = 0, mu2 = 0, sigma1 = 1, sigma2 = 1)))
})

test_that("calculate_mu1_for_dataframe recovers the requested J for every row", {
  df <- data.frame(mu2 = 0, var1 = c(1, 1), var2 = c(1, 1), J = c(0.4, 0.6))
  out <- calculate_mu1_for_dataframe(df)

  expect_named(out, c("mu1", "original_J", "theoretical_J"))
  expect_equal(nrow(out), 2)
  expect_equal(out$original_J, df$J)
  expect_equal(out$theoretical_J, df$J, tolerance = 1e-5)
})
