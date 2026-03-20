# LIBS
library(testthat)

# VARS
wd <- tempdir()

# RUNNING
context("Testing 'pipeline'")
phylotaR:::cleanup(wd)
test_that("setup() works", {
  res <- with_mocked_bindings(
    phylotaR::setup(wd = wd, txid = 9606),
    outfmt_get = function(...) "",
    cmdln = phylotaR:::cmdln_blastcheck,
    .package = "phylotaR"
  )
  expect_true(file.exists(file.path(wd, "cache", "prmtrs.RData")))
  phylotaR:::cleanup(wd)
})
test_that("run() works", {
  res <- with_mocked_bindings(
    phylotaR::run(wd = wd, nstages = 4),
    stages_run = function(...) NULL,
    .package = "phylotaR"
  )
  expect_null(res)
})
test_that("restart() works", {
  phylotaR:::cache_setup(ps = list("wd" = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = "taxise")
  res <- with_mocked_bindings(
    phylotaR::restart(wd = wd, nstages = 4),
    stages_run = function(...) NULL,
    .package = "phylotaR"
  )
  res <- with_mocked_bindings(
    expect_error(phylotaR::restart(wd = wd, nstages = 1)),
    stages_run = function(...) NULL,
    .package = "phylotaR"
  )
  phylotaR:::progress_save(wd = wd, stg = "download")
  phylotaR:::progress_save(wd = wd, stg = "cluster")
  phylotaR:::progress_save(wd = wd, stg = "cluster2")
  res <- with_mocked_bindings(
    expect_error(phylotaR::restart(wd = wd, nstages = 4)),
    stages_run = function(...) NULL,
    .package = "phylotaR"
  )
  phylotaR:::cleanup(wd)
})
test_that("reset(hard=FALSE) works", {
  phylotaR:::cache_setup(ps = parameters(wd = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = "taxise")
  phylotaR:::progress_save(wd = wd, stg = "download")
  phylotaR:::progress_save(wd = wd, stg = "cluster")
  phylotaR:::reset(wd = wd, stage = "download")
  expect_true(phylotaR:::progress_read(wd = wd) == "download")
  phylotaR:::cleanup(wd)
})
test_that("reset(hard=TRUE) works", {
  phylotaR:::cache_setup(ps = parameters(wd = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = "taxise")
  phylotaR:::progress_save(wd = wd, stg = "download")
  phylotaR:::obj_save(wd = wd, obj = NULL, nm = "txdct")
  phylotaR:::sqs_save(wd = wd, txid = "1", sqs = NULL)
  phylotaR:::reset(wd = wd, stage = "taxise", hard = TRUE)
  expect_true(phylotaR:::progress_read(wd = wd) == "taxise")
  expect_null(phylotaR:::obj_load(wd = wd, nm = "txdct"))
  expect_error(phylotaR:::sqs_load(wd = wd, txid = "1"))
  phylotaR:::cleanup(wd)
})
test_that("parameters_reset() works", {
  res <- with_mocked_bindings(
    phylotaR:::setup(wd = wd, txid = 9606),
    outfmt_get = function(...) "",
    cmdln = phylotaR:::cmdln_blastcheck,
    .package = "phylotaR"
  )
  phylotaR::parameters_reset(wd = wd, parameters = "txid", values = 0000)
  expect_true(phylotaR:::parameters_load(wd = wd)[["txid"]] == 0000)
  phylotaR:::cleanup(wd)
})
phylotaR:::cleanup(wd)
