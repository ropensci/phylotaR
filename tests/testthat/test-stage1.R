# LIBS
library(phylotaR)
library(testthat)

# DATA
recs <- readRDS(phylotaR:::datadir_get("txrecs.rda"))
wd <- tempdir()
ps <- parameters(wd = wd)

# RUNNING
context("Testing 'test-stage1'")
phylotaR:::cleanup(wd)
test_that("taxise_run() works", {
  with_mocked_bindings(
    {
      phylotaR::setup(wd = wd, txid = 9606)
      taxise_run(wd = wd)
    },
    outfmt_get = function(...) "",
    blast_setup = function(...) {
      list("mkblstdb" = ".", "blstn" = ".")
    },
    txids_get = function(...) NULL,
    batcher = function(...) NULL,
    taxdict_gen = function(...) NULL,
    obj_save = function(...) NULL,
    .package = "phylotaR"
  )
  lglns <- readLines(file.path(wd, "log.txt"))
  expect_true(grepl("Completed stage", lglns[length(lglns) - 1]))
})
phylotaR:::cleanup(wd)
test_that("txids_get() works", {
  mock_search <- function(...) {
    res <- list("count" = 100, "ids" = as.character(1:100))
    class(res) <- "esearch"
    res
  }
  phylotaR:::cache_setup(ps)
  res <- with_mocked_bindings(
    phylotaR:::txids_get(ps = ps, retmax = 150),
    entrez_search = mock_search,
    .package = "rentrez"
  )
  expect_true(length(res) == 100)
  phylotaR:::cleanup(wd)
})
phylotaR:::cleanup(wd)
test_that("taxdict_gen() works", {
  phylotaR:::cache_setup(ps)
  txids <- vapply(X = recs, FUN = function(x) x@id, character(1))
  res <- phylotaR:::taxdict_gen(txids, recs, ps)
  expect_true(inherits(res, "TaxDict"))
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
