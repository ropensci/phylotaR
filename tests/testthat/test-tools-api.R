# LIBS
library(testthat)

# DATA
load(file.path(phylotaR:::datadir_get(subdir  =  'api'), 'xml_timeout.rda'))
wd <- tempdir()
ps <- parameters(wd = wd)
ps[['mxrtry']] <- 2
ps[['wt_tms']] <- c(0, 0)

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'tools-api\'')
test_that('download_obj_check() works', {
  expect_false(phylotaR:::download_obj_check(xml_tmout))
})
test_that('safely_connect(fnm = search) works', {
  phylotaR:::cache_setup(ps  =  ps)
  args <- list('term'  =  'ncbi search term', 'db'  =  'nucleotide', 'x')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'search',
                                   ps = ps)
  expect_true(res == 1)
  myfunc <- function(...) {
    print(...)
    stop()
  }
  args <- list('term' = 'another ncbi search term', 'db' = 'nucleotide', 'x')
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'search',
                                   ps = ps)
  expect_null(res)
})
phylotaR:::cleanup(wd)
test_that('safely_connect(fnm = fetch) works', {
  phylotaR:::cache_setup(ps = ps)
  args <- list('id' = c(1, 2), 'db' = 'nucleotide', 'x')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'fetch',
                                   ps = ps)
  expect_true(res  ==  1)
  myfunc <- function(...) {
    print(...)
    stop()
  }
  args <- list('id' = c(2, 3), 'db' = 'nucleotide', 'x')
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'fetch',
                                   ps = ps)
  expect_null(res)
})
phylotaR:::cleanup(wd)
test_that('batcher() works', {
  phylotaR:::cache_setup(ps = ps)
  res <- phylotaR:::batcher(ids = 1:10, func = function(ids, ps) {
    return(ids)
  }, ps, lvl = 1)
  expect_true(all(res == 1:10))
  res <- phylotaR:::batcher(ids = 1:1000, func = function(ids, ps) {
    return(ids)
  }, ps, lvl = 1)
  expect_true(all(res == 1:1000))
})
phylotaR:::cleanup(wd)