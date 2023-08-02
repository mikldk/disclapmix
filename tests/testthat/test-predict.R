test_that("predict", {
  set.seed(1)
  ps_loc <- c(0.3, 0.2, 0.4)
  ps_ys <- c(0.1, 0)
  ys <- matrix(c(13L, 19L, 20L,
                 13L, 20L, 22L), nrow = 2L, byrow = TRUE) # Location parameter 
  tau <- c(0.2, 0.8)
  N <- 1000L
  z <- sample(seq_along(tau), N, prob = tau, replace = TRUE)
  x <- do.call(cbind, lapply(seq_along(ps_loc), function(k) {
    p <- ps_loc[k] + ps_ys[z]
    return(disclap::rdisclap(N, p) + ys[z, k])
  }))
  
  fit <- disclapmix(x = x, clusters = 2L)
  newx <- simulate(fit, nsim = 5L)
  #newx
  #predict(fit, newx)
  
  newx_dropout <- newx
  drop_idx <- sample(seq_len(ncol(newx)), nrow(newx), replace = TRUE)
  for (i in seq_len(nrow(newx_dropout))) {
    newx_dropout[i, drop_idx[i]] <- NA_integer_
  }
  #newx_dropout
  
  expect_error(predict(fit, newx_dropout))

  # Manual marginalisation: One locus
  for (i in seq_len(nrow(newx_dropout))) {
    h <- newx_dropout[i, ]
    ido <- drop_idx[i]
    dh <- do.call(rbind, lapply((-100):100, function(a) {
      hdo <- h
      hdo[ido] <- a
      hdo
    }))
    hm <- newx_dropout[i, , drop = FALSE]
    expect_error(predict(fit, hm))
    ph_marg <- sum(predict(fit, dh))
    ph_do <- predict(fit, hm, marginalise = TRUE)
    stopifnot(isTRUE(all.equal(ph_marg, ph_do)))
  }
  
  # Manual marginalisation: Two locus
  #for (i in seq_len(nrow(x))) {
  for (i in 1:10) {
    dh <- expand.grid(L1 = 15+(-30):30, 
                      L2 = x[i, 2], 
                      L3 = 15+(-30):30) |> 
      as.matrix()
    
    hm <- x[i, , drop = FALSE]
    hm[c(1, 3)] <- NA_integer_
    
    ph_marg <- sum(predict(fit, dh))
    ph_do <- predict(fit, hm, marginalise = TRUE)
    stopifnot(isTRUE(all.equal(ph_marg, ph_do, tolerance = 1e-6)))
  }
})

