# LIBS
library(phylotaR)
library(testthat)

# DATA
recs <- readRDS(phylotaR:::datadir_get('txrecs.rda'))
wd <- tempdir()
ps <- parameters(wd = wd)

# RUNNING
context('Testing \'test-stage1\'')
phylotaR:::cleanup(wd)
test_that('taxise_run() works', {
  with_mock(
    `phylotaR:::blast_setup` = function(...) {
      list('mkblstdb' = '.', 'blstn' = '.')}
    ,
    `phylotaR:::txids_get` = function(...) NULL,
    `phylotaR:::batcher` = function(...) NULL,
    `phylotaR:::taxdict_gen` = function(...) NULL,
    `phylotaR:::obj_save` = function(...) NULL,
    phylotaR::setup(wd = wd, txid = 9606),
    taxise_run(wd = wd)
  )
  lglns <- readLines(file.path(wd, 'log.txt'))
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
})
phylotaR:::cleanup(wd)
test_that('txids_get() works', {
  mock_search <- function(...) {
    list('count' = 100, 'ids' = as.character(1:100))
  }
  phylotaR:::cache_setup(ps)
  res <- with_mock(
    `rentrez::entrez_search` = mock_search,
    phylotaR:::txids_get(ps = ps, retmax = 150)
  )
  expect_true(length(res) == 100)
  phylotaR:::cleanup(wd)
})
phylotaR:::cleanup(wd)
test_that('taxdict_gen() works', {
  phylotaR:::cache_setup(ps)
  txids <- vapply(X = recs, FUN = function(x) x@id, character(1))
  res <- phylotaR:::taxdict_gen(txids, recs, ps)
  expect_true(inherits(res, 'TaxDict'))
})
phylotaR:::cleanup(wd)

# Example recs
# devtools::load_all('~/Coding/phylotaR')
# wd <- file.path(getwd(), 'anisoptera')
# dir.create(wd)
# ps <- parameters(wd = wd, txid = 6962)
# cache_setup(ps)
# txids <- txids_get(ps = ps)
# recs <- batcher(ids = txids, func = tax_download, ps = ps, lvl = 2)
# saveRDS(object = recs, file = phylotaR:::datadir_get('txrecs.rda'))
# unlink(x = wd, recursive = TRUE)
