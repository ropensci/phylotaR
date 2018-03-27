# LIBS
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
  phylotaR:::setUpCch(ps=list('wd'='.'))
  phylotaR:::intPrgrss(wd='.')
  expect_true(file.exists(file.path('cache', 'progress.RData')))
})
cleanUp()
test_that('rdPrgrss() works', {
  phylotaR:::setUpCch(ps=list('wd'='.'))
  phylotaR:::intPrgrss(wd='.')
  stg <- phylotaR:::rdPrgrss(wd='.')
  expect_true(stg == 'taxise')
})
cleanUp()
test_that('svPrgrss() works', {
  phylotaR:::setUpCch(ps=list('wd'='.'))
  phylotaR:::intPrgrss(wd='.')
  phylotaR:::svPrgrss(wd='.', stg='taxise')
  expect_true(phylotaR:::rdPrgrss(wd='.') == 'download')
  phylotaR:::svPrgrss(wd='.', stg='download')
  expect_true(phylotaR:::rdPrgrss(wd='.') == 'cluster')
  phylotaR:::svPrgrss(wd='.', stg='cluster')
  expect_true(phylotaR:::rdPrgrss(wd='.') == 'align')
  phylotaR:::svPrgrss(wd='.', stg='align')
  expect_true(is.na(phylotaR:::rdPrgrss(wd='.')))
})
cleanUp()
test_that('rstPrgrss() works', {
  # create prgrss
  phylotaR:::setUpCch(ps=list('wd'='.'))
  phylotaR:::intPrgrss(wd='.')
  phylotaR:::svPrgrss(wd='.', stg='taxise')
  phylotaR:::svPrgrss(wd='.', stg='download')
  phylotaR:::svPrgrss(wd='.', stg='cluster')
  phylotaR:::svPrgrss(wd='.', stg='align')
  # reset
  phylotaR:::rstPrgrss(wd='.', stg='download')
  expect_true(phylotaR:::rdPrgrss(wd='.') == 'download')
})
cleanUp()
test_that('setUpCch() works', {
  ps <- list('alovelyparameter'=1, 'wd'='.')
  phylotaR:::setUpCch(ps=ps)
  expect_true(file.exists(file.path('cache',
                                    'prmtrs.RData')))
  expect_error(phylotaR:::setUpCch(ps=ps))
})
cleanUp()
test_that('ldPrmtrs() works', {
  phylotaR:::setUpCch(ps=list('anotherlovelyparameter'=1, 'wd'='.'))
  ps <- ldPrmtrs(wd='.')
  expect_true('anotherlovelyparameter' %in% names(ps))
})
cleanUp()
test_that('rmCch() works', {
  phylotaR:::setUpCch(ps=list('andfinally'=1, 'wd'='.'))
  rmCch(wd='.')
  expect_false(file.exists(file.path('cache',
                                     'prmtrs.RData')))
})
cleanUp()
test_that('chkObj() works', {
  expect_false(phylotaR:::chkObj(wd='.', nm='mylovelyobj'))
  phylotaR:::svObj(wd='.', nm='mylovelyobj',
                  obj=list('interestingthings'))
  expect_true(phylotaR:::chkObj(wd='.', nm='mylovelyobj'))
})
cleanUp()
test_that('svObj() works', {
  phylotaR:::svObj(wd='.', nm='mylovelyobj',
                  obj=list('interestingthings'))
  expect_true(file.exists(file.path('cache',
                                    'mylovelyobj.RData')))
})
cleanUp()
test_that('ldObj() works', {
  expect_error(phylotaR:::ldObj(wd='.', nm='mylovelyobj'))
  phylotaR:::svObj(wd='.', nm='mylovelyobj',
        obj=list('interestingthings'))
  obj <- phylotaR:::ldObj(wd='.', nm='mylovelyobj')
  expect_true(obj[[1]] == 'interestingthings')
})
cleanUp()
test_that('svNcbiCch() works', {
  phylotaR:::setUpCch(ps=list('wd'='.'))
  args <- list('term'='unique_search_term',
               'db'='nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  exptd <- file.path('cache', 'ncbi', 'search',
                     'nucleotide', '1.RData')
  phylotaR:::svNcbiCch(fnm=fnm, args=args, wd='.', obj=obj)
  expect_true(file.exists(exptd))
})
cleanUp()
test_that('ldNcbiCch() works', {
  phylotaR:::setUpCch(ps=list('wd'='.'))
  args <- list('term'='unique_search_term',
               'db'='nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  phylotaR:::svNcbiCch(fnm=fnm, args=args, wd='.', obj=obj)
  res <- phylotaR:::ldNcbiCch(fnm=fnm, args=args, wd='.')
  expect_equal(obj, res)
})
cleanUp()
test_that('svBlstCch() works', {
  phylotaR:::setUpCch(ps=list('wd'='.'))
  sids_1 <- c('1', '2', '3')
  sids_2 <- c('1', '2', '4')
  phylotaR:::svBlstCch(sids=sids_1, wd='.', obj='res_1')
  phylotaR:::svBlstCch(sids=sids_2, wd='.', obj='res_2')
  expect_true(file.exists(file.path('cache', 'blast', '1.RData')))
  expect_true(file.exists(file.path('cache', 'blast', '2.RData')))
  expect_false(file.exists(file.path('cache', 'blast', '3.RData')))
})
cleanUp()
test_that('ldBlstCch() works', {
  phylotaR:::setUpCch(ps=list('wd'='.'))
  blst_rs <- data.frame('query.id'=as.character(c(1:5)),
                        'subject.id'='1')
  sids_1 <- c('1', '2', '3')
  sids_2 <- c('1', '2', '4')
  sids_3 <- c('1', '2', '5')
  sids_4 <- c('1', '2')
  phylotaR:::svBlstCch(sids=sids_1, wd='.', obj=blst_rs[1:3, ])
  phylotaR:::svBlstCch(sids=sids_2, wd='.', obj=blst_rs[c(1, 2, 4), ])
  expect_true(nrow(phylotaR:::ldBlstCch(sids=sids_1, wd='.')) == 3)
  expect_true(nrow(phylotaR:::ldBlstCch(sids=sids_2, wd='.')) == 3)
  expect_null(phylotaR:::ldBlstCch(sids=sids_3, wd='.'))
  # sqs 4 blast results are contained within sqs 1 and/or 2
  expect_true(nrow(phylotaR:::ldBlstCch(sids=sids_4, wd='.')) == 2)
})
cleanUp()