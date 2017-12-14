# LIBS
library(phylotaR)
library(testthat)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('log.txt')) {
    file.remove('log.txt')
  }
}

# RUNNING
cleanUp()
context('Testing \'log-tools\'')
test_that('info() works', {
  msg <- 'test message'
  info(wd='.', lvl=1, v=FALSE, msg)
  info(wd='.', lvl=2, v=FALSE, msg)
  res <- scan(file='log.txt', what=character())
  expect_true(length(res) == 5)
  cleanUp()
})
test_that('error() works', {
  msg <- 'test error'
  expect_error(error(wd='.', msg))
  res <- scan(file='log.txt', what=character())
  expect_true(length(res) == 3)
  cleanUp()
})
test_that('warn() works', {
  msg <- 'test warn'
  expect_warning(warn(wd='.', msg))
  res <- scan(file='log.txt', what=character())
  expect_true(length(res) == 3)
  cleanUp()
})
cleanUp()