# LIBS
library(testthat)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
sqs <- readRDS(phylotaR:::datadir_get('sqrecs.rda'))
txdct <- readRDS(phylotaR:::datadir_get('txdct.rda'))

# RUNNING
context('Testing \'stage2\'')
phylotaR:::cleanup(wd)
test_that('download_run() works', {
  with_mock(
    `phylotaR:::parameters_load` = function(...) ps,
    `phylotaR:::clade_select` = function(...) NULL,
    `phylotaR:::seq_download` = function(...) NULL,
    `phylotaR:::obj_load` = function(...) NULL,
    download_run(wd = wd)
  )
  lglns <- readLines(file.path(wd, 'log.txt'))
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
})
phylotaR:::cleanup(wd)
test_that('clade_select() works', {
  ps[['txid']] <- '9606'
  ps[['mxsqs']] <- 10
  ps[['mxnds']] <- 10
  res <- with_mock(
    `phylotaR:::sqs_count` = function(...) 1,
    `phylotaR:::txnds_count` = function(...) 1,
    `phylotaR:::descendants_get` = function(...) NULL,
    phylotaR:::clade_select(txdct = NULL, ps = ps)
  )
  expect_true(res == ps[['txid']])
  res <- with_mock(
    `phylotaR:::sqs_count` = function(...) 11,
    `phylotaR:::txnds_count` = function(...) 1,
    `phylotaR:::descendants_get` = function(...) NULL,
    phylotaR:::clade_select(txdct = NULL, ps = ps)
  )
  expect_true(length(res) == 0)
  res <- with_mock(
    `phylotaR:::sqs_count` = function(...) 1,
    `phylotaR:::txnds_count` = function(...) 11,
    `phylotaR:::descendants_get` = function(...) NULL,
    phylotaR:::clade_select(txdct = NULL, ps = ps)
  )
  expect_true(length(res) == 0)
})
phylotaR:::cleanup(wd)
test_that('seq_download() works', {
  phylotaR:::cache_setup(ps)
  with_mock(
    `phylotaR:::hierarchic_download` = function(...) sqs,
    phylotaR:::seq_download(txids = '1', txdct = txdct, ps = ps)
  )
  expect_true(file.exists(file.path(wd, 'cache', 'sqs', '1.RData')))
})
phylotaR:::cleanup(wd)
