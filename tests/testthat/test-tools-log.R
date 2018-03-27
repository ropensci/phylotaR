# LIBS
library(testthat)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('log.txt')) {
    file.remove('log.txt')
  }
}

# DUMMIES
ps <- list('wd'='.', 'v'=FALSE)

# RUNNING
cleanUp()
context('Testing \'log-tools\'')
test_that('info() works', {
  msg <- 'test message'
  phylotaR:::info(ps=ps, lvl=1, msg)
  phylotaR:::info(ps=ps, lvl=2, msg)
  res <- scan(file='log.txt', what=character())
  expect_true(length(res) == 5)
  cleanUp()
})
test_that('error() works', {
  msg <- 'test error'
  expect_error(phylotaR:::error(ps=ps, msg))
  res <- scan(file='log.txt', what=character())
  expect_true(length(res) == 3)
  cleanUp()
})
test_that('warn() works', {
  msg <- 'test warn'
  expect_warning(phylotaR:::warn(ps=ps, msg))
  res <- scan(file='log.txt', what=character())
  expect_true(length(res) == 3)
  cleanUp()
})
cleanUp()