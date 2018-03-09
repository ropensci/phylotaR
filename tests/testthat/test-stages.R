# LIBS
library(testthat)

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
    `phylotaR:::getTxids`=mckFun,
    `phylotaR:::dwnldTxRcrds`=mckFun,
    `phylotaR:::genTxDct`=mckFun,
    `phylotaR:::svObj`=mckFun,
    phylotaR::setUp(wd='.', txid=9606),
    phylotaR::runTaxise(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage',
                    lglns[length(lglns) - 1]))
})
cleanUp()
test_that('runDownload() works', {
  with_mock(
    `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
    `phylotaR:::ldObj`=mckFun,
    `phylotaR:::cldIdntfy`=mckFun,
    `phylotaR:::dwnldSqRcrds`=mckFun,
    phylotaR::setUp(wd='.', txid=9606),
    phylotaR::runDownload(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage',
                    lglns[length(lglns) - 1]))
})
cleanUp()
test_that('runClusters() works', {
  with_mock(
    `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
    `phylotaR:::ldObj`=mckFun,
    `phylotaR:::calcClstrs`=mckFun,
    phylotaR::setUp(wd='.', txid=9606),
    phylotaR::runClusters(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage',
                    lglns[length(lglns) - 1]))
})
cleanUp()
test_that('runClusters2() works', {
  with_mock(
    `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
    `phylotaR:::ldObj`=mckFun,
    `phylotaR:::clstrClstrs`=mckFun,
    phylotaR::setUp(wd='.', txid=9606),
    phylotaR::runClusters2(wd='.')
  )
  lglns <- readLines('log.txt')
  expect_true(grepl('Completed stage',
                    lglns[length(lglns) - 1]))
})
cleanUp()