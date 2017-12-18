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
test_that('intPrgrss() works', {
  setUpCch(ps=list('wd'='.'))
  intPrgrss(wd='.')
  expect_true(file.exists(file.path('cache', 'progress.RData')))
  cleanUp()
})
test_that('rdPrgrss() works', {
  setUpCch(ps=list('wd'='.'))
  intPrgrss(wd='.')
  stg <- rdPrgrss(wd='.')
  expect_true(stg == 'taxise')
  cleanUp()
})
test_that('svPrgrss() works', {
  setUpCch(ps=list('wd'='.'))
  intPrgrss(wd='.')
  svPrgrss(wd='.', stg='taxise')
  expect_true(rdPrgrss(wd='.') == 'download')
  svPrgrss(wd='.', stg='download')
  expect_true(rdPrgrss(wd='.') == 'cluster')
  svPrgrss(wd='.', stg='cluster')
  expect_true(rdPrgrss(wd='.') == 'align')
  svPrgrss(wd='.', stg='align')
  expect_true(is.na(rdPrgrss(wd='.')))
  cleanUp()
})
test_that('rstPrgrss() works', {
  # create prgrss
  setUpCch(ps=list('wd'='.'))
  intPrgrss(wd='.')
  svPrgrss(wd='.', stg='taxise')
  svPrgrss(wd='.', stg='download')
  svPrgrss(wd='.', stg='cluster')
  svPrgrss(wd='.', stg='align')
  # reset
  rstPrgrss(wd='.', stg='download')
  expect_true(rdPrgrss(wd='.') == 'download')
  cleanUp()
})
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