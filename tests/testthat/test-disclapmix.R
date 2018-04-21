context("disclapmix")

ESTIMATION_TOL_DUE_TO_SAMPLING <- 1e-3

set.seed(1)
p <- 0.3 # Dispersion parameter 
y <- 13 # Location parameter 
x <- disclap::rdisclap(1000, p) + y # Generate a sample using the rdisclap function

y_hat <- median(x)
mu_hat <- mean(abs(x - y_hat))
p_hat <- mu_hat^(-1) * (sqrt(mu_hat^2 + 1) - 1)

test_that("rdisclap", {
  expect_equal(y_hat, y)
  expect_equal(p_hat, p, tol = ESTIMATION_TOL_DUE_TO_SAMPLING)
})

init_y <- 12L
fit <- disclapmix(x = as.matrix(as.integer(x)), 
                  clusters = 1L, 
                  init_y = matrix(init_y), 
                  init_y_method = NULL)

test_that("disclapmix", {
  expect_equal(as.integer(fit$init_y), init_y)
  expect_equal(as.integer(fit$y), y)
  expect_equal(c(fit$disclap_parameters), p, tol = ESTIMATION_TOL_DUE_TO_SAMPLING)
  expect_equal(c(fit$disclap_parameters), p_hat)
})
