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
cleanUp()