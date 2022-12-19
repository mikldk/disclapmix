library(disclapmix)
data(danes)
db <- as.matrix(danes[rep(1:nrow(danes), danes$n), 1:(ncol(danes)-1)])

get_best_BIC <- function(fits) {
  fits[[which.min(sapply(fits, function(x) x$BIC_marginal))]]
}

get_best_AICc <- function(fits) {
  fits[[which.min(sapply(fits, function(x) x$AICc_marginal))]]
}

get_c <- function(fit) {
  nrow(fit$y)
}

test_that("adaptive", {
  expect_error(disclapmix_adaptive(db, margin = 5L, criteria = 'X'))
  
  marg <- 5L
  
  fits <- disclapmix_adaptive(db, margin = marg, criteria = 'BIC')
  expect_equal(fits |> length(), 4L + marg)
  expect_equal(fits |> get_best_BIC() |> get_c(), 4L)
  
  
  fits <- disclapmix_adaptive(db, margin = marg, criteria = c('AICc', 'BIC'))
  expect_equal(fits |> length(), 9L + marg)
  expect_equal(fits |> get_best_AICc() |> get_c(), 9L)
})
