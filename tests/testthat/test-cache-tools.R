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
})
cleanUp()
test_that('rdPrgrss() works', {
  setUpCch(ps=list('wd'='.'))
  intPrgrss(wd='.')
  stg <- rdPrgrss(wd='.')
  expect_true(stg == 'taxise')
})
cleanUp()
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
})
cleanUp()
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
})
cleanUp()
test_that('setUpCch() works', {
  ps <- list('alovelyparameter'=1, 'wd'='.')
  setUpCch(ps=ps)
  expect_true(file.exists(file.path('cache',
                                    'prmtrs.RData')))
  expect_error(setUpCch(ps=ps))
})
cleanUp()
test_that('ldPrmtrs() works', {
  setUpCch(ps=list('anotherlovelyparameter'=1, 'wd'='.'))
  ps <- ldPrmtrs(wd='.')
  expect_true('anotherlovelyparameter' %in% names(ps))
})
cleanUp()
test_that('rmCch() works', {
  setUpCch(ps=list('andfinally'=1, 'wd'='.'))
  rmCch(wd='.')
  expect_false(file.exists(file.path('cache',
                                     'prmtrs.RData')))
})
cleanUp()
test_that('chkObj() works', {
  expect_false(chkObj(wd='.', nm='mylovelyobj'))
  svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  expect_true(chkObj(wd='.', nm='mylovelyobj'))
})
cleanUp()
test_that('svObj() works', {
  svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  expect_true(file.exists(file.path('cache',
                                    'mylovelyobj.RData')))
})
cleanUp()
test_that('ldObj() works', {
  expect_error(ldObj(wd='.', nm='mylovelyobj'))
  svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  obj <- ldObj(wd='.', nm='mylovelyobj')
  expect_true(obj[[1]] == 'interestingthings')
})
cleanUp()
test_that('svNcbiCch() works', {
  setUpCch(ps=list('wd'='.'))
  args <- list('term'='unique_search_term',
               'db'='nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  exptd <- file.path('cache', 'ncbi', 'search',
                     'nucleotide', '1.RData')
  svNcbiCch(fnm=fnm, args=args, wd='.', obj=obj)
  expect_true(file.exists(exptd))
})
cleanUp()
test_that('ldNcbiCch() works', {
  setUpCch(ps=list('wd'='.'))
  args <- list('term'='unique_search_term',
               'db'='nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  svNcbiCch(fnm=fnm, args=args, wd='.', obj=obj)
  res <- ldNcbiCch(fnm=fnm, args=args, wd='.')
  expect_equal(obj, res)
})
cleanUp()
test_that('svBlstCch() works', {
  setUpCch(ps=list('wd'='.'))
  sqs_1 <- list(list('gi'='1'), list('gi'='2'), list('gi'='3'))
  sqs_2 <- list(list('gi'='1'), list('gi'='2'), list('gi'='4'))
  svBlstCch(sqs=sqs_1, wd='.', obj='res_1')
  svBlstCch(sqs=sqs_2, wd='.', obj='res_2')
  expect_true(file.exists(file.path('cache', 'blast', '1.RData')))
  expect_true(file.exists(file.path('cache', 'blast', '2.RData')))
  expect_false(file.exists(file.path('cache', 'blast', '3.RData')))
})
cleanUp()
test_that('ldBlstCch() works', {
  setUpCch(ps=list('wd'='.'))
  sqs_1 <- list(list('gi'='1'), list('gi'='2'), list('gi'='3'))
  sqs_2 <- list(list('gi'='1'), list('gi'='2'), list('gi'='4'))
  sqs_3 <- list(list('gi'='1'), list('gi'='2'), list('gi'='5'))
  svBlstCch(sqs=sqs_1, wd='.', obj='res_1')
  svBlstCch(sqs=sqs_2, wd='.', obj='res_2')
  expect_true(ldBlstCch(sqs=sqs_1, wd='.') == 'res_1')
  expect_true(ldBlstCch(sqs=sqs_2, wd='.') == 'res_2')
  expect_null(ldBlstCch(sqs=sqs_3, wd='.'))
})
cleanUp()