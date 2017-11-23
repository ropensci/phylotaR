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
  prmtrs <- list('alovelyparameter'=1)
  setUpCch(wd='.', prmtrs=prmtrs)
  expect_true(file.exists(file.path('cache',
                                    'prmtrs.RData')))
  expect_error(setUpCch(wd='.', prmtrs=prmtrs))
  cleanUp()
})
test_that('ldPrmtrs() works', {
  setUpCch(wd='.', prmtrs=list('anotherlovelyparameter'=1))
  prmtrs <- ldPrmtrs(wd='.')
  expect_true(names(prmtrs) == 'anotherlovelyparameter')
  cleanUp()
})
test_that('rmCch() works', {
  setUpCch(wd='.', prmtrs=list('andfinally'=1))
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