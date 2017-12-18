# LIBS
library(phylotaR)
library(testthat)

# VARS

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
  if(file.exists('log.txt')) {
    file.remove('log.txt')
  }
}
mckSetupNcbiTools <- function(...) {
  list('mkblstdb'='.',
       'blstn'='.')
}
mckFun <- function(...) {
  NULL
}

# RUNNING
context('Testing \'stages\'')
cleanUp()
test_that('runTaxise() works', {
  with_mock(
    `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
    `phylotaR:::dwnldTD`=mckFun,
    `phylotaR:::genTDObj`=mckFun,
    `phylotaR:::getMngblIds`=mckFun,
    `phylotaR:::genPhylotaNds`=mckFun,
    `phylotaR:::writeTax`=mckFun,
    setUp(wd='.', txid=9606),
    runTaxise(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
  cleanUp()
})
test_that('runDownload() works', {
  with_mock(
    `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
    `phylotaR:::ldObj`=mckFun,
    `phylotaR:::fltr`=mckFun,
    `phylotaR:::dwnld`=mckFun,
    setUp(wd='.', txid=9606),
    runDownload(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
  cleanUp()
})
test_that('runClusters() works', {
  with_mock(
    `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
    `phylotaR:::ldObj`=mckFun,
    `phylotaR:::calcClstrs`=mckFun,
    setUp(wd='.', txid=9606),
    runClusters(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
  cleanUp()
})
test_that('runAlign() works', {
  
})
cleanUp()