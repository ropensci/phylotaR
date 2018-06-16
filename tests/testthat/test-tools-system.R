# LIBS
library(phylotaR)
library(testthat)

# TEST
context('Testing \'tools-system\'')
phylotaR:::cleanup()
test_that('cmdln() works', {
  res <- phylotaR:::cmdln(cmd = 'notacommand', args = list())
  expect_true(res[['status']] == 1)
  res <- phylotaR:::cmdln(cmd = 'notacommand', args = list(), lgfl = 'log')
  expect_true(file.exists('log.log'))
  expect_true(res == 1)
})
phylotaR:::cleanup()
