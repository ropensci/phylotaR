# LIBS
library(phylotaR)
library(testthat)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}

# RUNNING
context('Testing \'cache-tools\'')
cleanUp()
test_that('setUpCch() works', {
  ps <- list('alovelyparameter'=1, 'wd'='.')
  setUpCch(ps=ps)
  expect_true(file.exists(file.path('cache',
                                    'prmtrs.RData')))
  expect_error(setUpCch(ps=ps))
  cleanUp()
})
test_that('ldPrmtrs() works', {
  setUpCch(ps=list('anotherlovelyparameter'=1, 'wd'='.'))
  ps <- ldPrmtrs(wd='.')
  expect_true('anotherlovelyparameter' %in% names(ps))
  cleanUp()
})
test_that('rmCch() works', {
  setUpCch(ps=list('andfinally'=1, 'wd'='.'))
  rmCch(wd='.')
  expect_false(file.exists(file.path('cache',
                                     'prmtrs.RData')))
  cleanUp()
})
test_that('chkObj() works', {
  expect_false(chkObj(wd='.', nm='mylovelyobj'))
  svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  expect_true(chkObj(wd='.', nm='mylovelyobj'))
  cleanUp()
})
test_that('svObj() works', {
  svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  expect_true(file.exists(file.path('cache',
                                    'mylovelyobj.RData')))
  cleanUp()
})
test_that('ldObj() works', {
  expect_error(ldObj(wd='.', nm='mylovelyobj'))
  svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  obj <- ldObj(wd='.', nm='mylovelyobj')
  expect_true(obj[[1]] == 'interestingthings')
  cleanUp()
})
cleanUp()