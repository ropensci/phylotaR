# LIBS
library(phylotaR)
library(testthat)

# VARS
wd <- tempdir()
ps <- parameters(wd = wd)

# TEST
context('Testing \'tools-system\'')
phylotaR:::cleanup(wd)
test_that('cmdln() works', {
  res <- phylotaR:::cmdln(cmd = 'notacommand', args = list(), ps = ps)
  expect_true(res[['status']] == 1)
  res <- phylotaR:::cmdln(cmd = 'notacommand', args = list(), lgfl =
                            file.path(wd, 'log'), ps = ps)
  expect_true(file.exists(file.path(wd, 'log.log')))
  expect_true(res == 1)
})
phylotaR:::cleanup(wd)
