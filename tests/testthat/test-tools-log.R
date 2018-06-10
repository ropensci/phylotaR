# LIBS
library(testthat)

# DATA
ps <- parameters()

# RUNNING
phylotaR:::cleanup()
context('Testing \'log-tools\'')
test_that('info() works', {
  msg <- 'test message'
  phylotaR:::info(ps = ps, lvl = 1, msg)
  phylotaR:::info(ps = ps, lvl = 2, msg)
  res <- scan(file = 'log.txt', what = character())
  expect_true(length(res) == 5)
  phylotaR:::cleanup()
})
test_that('error() works', {
  msg <- 'test error'
  expect_error(phylotaR:::error(ps = ps, msg))
  res <- scan(file = 'log.txt', what = character())
  expect_true(length(res) == 2)
  phylotaR:::cleanup()
})
test_that('warn() works', {
  msg <- 'test warn'
  expect_warning(phylotaR:::warn(ps = ps, msg))
  res <- scan(file = 'log.txt', what = character())
  expect_true(length(res) == 3)
  phylotaR:::cleanup()
})
phylotaR:::cleanup()