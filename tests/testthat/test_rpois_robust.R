context("rpois for large numbers")

test_that("rpois coincides with rpois_robust", {
  
  lambda <- c(1000, 100, 1e11)
  
  set.seed(20200317)
  rp <- rpois(3, lambda)
  
  set.seed(20200317)
  rprob <- rpois_robust(3, lambda)
  
  expect_equal(rp[1:2], rprob[1:2], tolerance = 0.01)
  
})

test_that("rpois for large numbers runs and has correct moments", {
  
  n <- 10000
  lambda <- rep(1.23e10, n)
  
  set.seed(20200317)
  rprob <- rpois_robust(n, lambda)
  
  expect_equal(mean(rprob), lambda[1], tolerance = 0.01)
  expect_equal(sd(rprob), sqrt(lambda[1]), tolerance = 0.01)
  
})