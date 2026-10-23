test_that("gpq_mu matches the closed-form pivot", {
  expect_equal(gpq_mu(ybar = 5, s = 1.2, n = 30, t_value = 0.5),
               5 - 0.5 * 1.2 / sqrt(30))
  expect_equal(gpq_mu(ybar = 2, s = 1, n = 10, t_value = 0), 2)
})

test_that("gpq_sigma matches the closed-form pivot", {
  expect_equal(gpq_sigma(s2 = 1.5, v = 28.3, n = 30),
               sqrt(29 * 1.5 / 28.3))
})

test_that("gpq_c reduces to the midpoint when the two scales coincide", {
  expect_equal(gpq_c(R_mu1 = 2, R_mu2 = 0, R_sigma1 = 1, R_sigma2 = 1), 1)
  expect_equal(gpq_c(R_mu1 = 5, R_mu2 = -3, R_sigma1 = 2, R_sigma2 = 2), 1)
})

test_that("gpq_c matches the quadratic formula when scales differ", {
  R_mu1 <- 2; R_mu2 <- 0; R_sigma1 <- 1; R_sigma2 <- 1.5
  Ra <- R_mu1 - R_mu2
  Rb <- R_sigma1 / R_sigma2
  expected <- (R_mu2 * (Rb^2 - 1) - Ra +
                 Rb * sqrt(Ra^2 + (Rb^2 - 1) * R_sigma2^2 * log(Rb^2))) /
    (Rb^2 - 1)
  expect_equal(gpq_c(R_mu1, R_mu2, R_sigma1, R_sigma2), expected)
})

test_that("gpq_J matches the sensitivity + specificity - 1 formula", {
  expect_equal(
    gpq_J(R_mu1 = 2, R_mu2 = 0, R_sigma1 = 1, R_sigma2 = 1.2, Rc = 1),
    pnorm((2 - 1) / 1) + pnorm((1 - 0) / 1.2) - 1
  )
})

test_that("gpq_J always lies in the open interval (-1, 1)", {
  set.seed(42)
  Rc_vals <- rnorm(200, 1, 1)
  RJ <- mapply(function(rc) gpq_J(2, 0, 1, 1.2, rc), Rc_vals)
  expect_true(all(RJ > -1 & RJ < 1))
})

test_that("gpq_from_sample returns K realisations within (-1, 1)", {
  set.seed(1)
  y1 <- rnorm(30, mean = 2, sd = 1)
  y2 <- rnorm(30, mean = 0, sd = 1)
  gpq <- gpq_from_sample(y1, y2, K = 100)

  expect_named(gpq, c("RJ", "Rc"))
  expect_length(gpq$RJ, 100)
  expect_length(gpq$Rc, 100)
  expect_true(all(is.finite(gpq$RJ)))
  expect_true(all(is.finite(gpq$Rc)))
  expect_true(all(gpq$RJ > -1 & gpq$RJ < 1))
})

test_that("gpq_from_sample is consistent with gpq_mu/gpq_sigma/gpq_c/gpq_J", {
  # With K = 1 and a fixed seed, gpq_from_sample's single draw must match
  # calling the four building blocks by hand with the same random numbers.
  set.seed(7)
  y1 <- rnorm(20, mean = 3, sd = 1.5)
  y2 <- rnorm(20, mean = 0, sd = 1)

  set.seed(99)
  gpq <- gpq_from_sample(y1, y2, K = 1)

  set.seed(99)
  t1 <- rt(1, df = length(y1) - 1)
  t2 <- rt(1, df = length(y2) - 1)
  v1 <- rchisq(1, df = length(y1) - 1)
  v2 <- rchisq(1, df = length(y2) - 1)

  R_mu1 <- gpq_mu(mean(y1), sd(y1), length(y1), t1)
  R_mu2 <- gpq_mu(mean(y2), sd(y2), length(y2), t2)
  R_sigma1 <- gpq_sigma(var(y1), v1, length(y1))
  R_sigma2 <- gpq_sigma(var(y2), v2, length(y2))
  Rc <- gpq_c(R_mu1, R_mu2, R_sigma1, R_sigma2)
  RJ <- gpq_J(R_mu1, R_mu2, R_sigma1, R_sigma2, Rc)

  expect_equal(gpq$Rc[1], Rc)
  expect_equal(gpq$RJ[1], RJ)
})
