# LIBS
library(testthat)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
sqs <- readRDS(phylotaR:::datadir_get("sqrecs.rda"))
txdct <- readRDS(phylotaR:::datadir_get("txdct.rda"))

# RUNNING
context("Testing 'stage2'")
phylotaR:::cleanup(wd)
test_that("download_run() works", {
  with_mocked_bindings(
    download_run(wd = wd),
    parameters_load = function(...) ps,
    clade_select = function(...) NULL,
    seq_download = function(...) NULL,
    obj_load = function(...) NULL,
    .package = "phylotaR"
  )
  lglns <- readLines(file.path(wd, "log.txt"))
  expect_true(grepl("Completed stage", lglns[length(lglns) - 1]))
})
phylotaR:::cleanup(wd)
test_that("clade_select() works", {
  ps[["txid"]] <- "9606"
  ps[["mxsqs"]] <- 10
  ps[["mxnds"]] <- 10
  res <- with_mocked_bindings(
    phylotaR:::clade_select(txdct = NULL, ps = ps),
    sqs_count = function(...) 1,
    txnds_count = function(...) 1,
    descendants_get = function(...) NULL,
    .package = "phylotaR"
  )
  expect_true(res == ps[["txid"]])
  res <- with_mocked_bindings(
    phylotaR:::clade_select(txdct = NULL, ps = ps),
    sqs_count = function(...) 11,
    txnds_count = function(...) 1,
    descendants_get = function(...) NULL,
    .package = "phylotaR"
  )
  expect_true(length(res) == 0)
  res <- with_mocked_bindings(
    phylotaR:::clade_select(txdct = NULL, ps = ps),
    sqs_count = function(...) 1,
    txnds_count = function(...) 11,
    descendants_get = function(...) NULL,
    .package = "phylotaR"
  )
  expect_true(length(res) == 0)
})
phylotaR:::cleanup(wd)
test_that("seq_download() works", {
  phylotaR:::cache_setup(ps)
  with_mocked_bindings(
    phylotaR:::seq_download(txids = "1", txdct = txdct, ps = ps),
    hierarchic_download = function(...) sqs,
    .package = "phylotaR"
  )
  expect_true(file.exists(file.path(wd, "cache", "sqs", "1.RData")))
})
phylotaR:::cleanup(wd)
