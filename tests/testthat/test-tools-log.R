# LIBS
library(testthat)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'log-tools\'')
test_that('info() works', {
  msg <- 'test message'
  phylotaR:::info(ps = ps, lvl = 1, msg)
  phylotaR:::info(ps = ps, lvl = 2, msg)
  res <- scan(file = file.path(wd, 'log.txt'), what = character())
  expect_true(length(res) == 5)
  phylotaR:::cleanup(wd)
})
test_that('error() works', {
  msg <- 'test error'
  expect_error(phylotaR:::error(ps = ps, msg))
  res <- scan(file = file.path(wd, 'log.txt'), what = character())
  expect_true(length(res) == 2)
  phylotaR:::cleanup(wd)
})
test_that('warn() works', {
  msg <- 'test warn'
  expect_warning(phylotaR:::warn(ps = ps, msg))
  res <- scan(file = file.path(wd, 'log.txt'), what = character())
  expect_true(length(res) == 3)
  phylotaR:::cleanup(wd)
})
phylotaR:::cleanup(wd)