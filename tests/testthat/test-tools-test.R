# LIBS
library(phylotaR)
library(testthat)

# RUNNING
context('Testing \'test-tools\'')
test_that('datadir_get() works', {
  res <- phylotaR:::datadir_get()
  expect_true(dir.exists(res))
})
test_that('random_phylota() works', {
  res <- phylotaR:::random_phylota()
  expect_true(inherits(res, 'Phylota'))
})
test_that('cleanup() works', {
  wd <- tempdir()
  res <- phylotaR:::cleanup(wd)
  expect_null(res)
  expect_false(dir.exists(file.path(wd, 'cache')))
})
test_that('cmdln_blastcheck() works', {
  res <- phylotaR:::cmdln_blastcheck(cmd = 'blastn', args = NULL)
  expect_true(inherits(res, 'list'))
})
test_that('testsqs_gen() works', {
  res <- phylotaR:::testsqs_gen(n = 10, l = 10)
  expect_true(inherits(res, 'SeqArc'))
  expect_true(all(res@nncltds == 10))
})
