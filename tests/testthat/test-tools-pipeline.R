# LIBS
library(testthat)

# VARS
wd <- tempdir()

# RUNNING
context("Testing 'setup-tools'")
phylotaR:::cleanup(wd)
test_that("blast_setup() works", {
  # test with fake system
  res <- with_mocked_bindings(
    phylotaR:::blast_setup(d = ".", v = FALSE, wd = wd, otsdr = FALSE),
    outfmt_get = function(...) "",
    cmdln = phylotaR:::cmdln_blastcheck,
    .package = "phylotaR"
  )
  expect_true(length(res) == 2)
  # make sure wrong versions are flagged
  res <- with_mocked_bindings(
    expect_error(phylotaR:::blast_setup(
      d = "wrngvrsn", v = FALSE,
      wd = wd, otsdr = FALSE
    )),
    outfmt_get = function(...) "",
    cmdln = phylotaR:::cmdln_blastcheck,
    .package = "phylotaR"
  )
  # make sure wrong dirs are flagged
  expect_error(phylotaR:::blast_setup(
    d = ".", v = FALSE, wd = NULL,
    otsdr = FALSE
  ))
})
test_that("parameters_setup() works", {
  expect_error(phylotaR:::parameters_setup(
    wd = wd, txid = 9606,
    ncbi_execs = c("", "")
  ))
  ncbi_execs <- list("mkblstdb" = NA, "blstn" = NA)
  with_mocked_bindings(
    phylotaR:::parameters_setup(wd = wd, txid = 9606, ncbi_execs = ncbi_execs),
    outfmt_get = function(...) "",
    .package = "phylotaR"
  )
  ps <- phylotaR:::parameters_load(wd = wd)
  expect_true(length(ps) == 22)
})
phylotaR:::cleanup(wd)
test_that("stage_args_check() works", {
  expect_error(phylotaR:::stage_args_check(frm = -1, to = -1))
  expect_error(phylotaR:::stage_args_check(frm = 5, to = 5))
  expect_error(phylotaR:::stage_args_check(frm = 2, to = 1))
  res <- phylotaR:::stage_args_check(frm = 1, to = 4)
  expect_true(res == "Running stages: taxise, download, cluster, cluster2")
  res <- phylotaR:::stage_args_check(frm = 4, to = 4)
  expect_true(res == "Running stages: cluster2")
})
test_that("stages_run() works", {
  res <- with_mocked_bindings(
    {
      phylotaR::setup(wd = wd, txid = 9606)
      phylotaR:::stages_run(
        wd = wd, to = 4, frm = 1, stgs_msg = "",
        rstrt = FALSE
      )
    },
    outfmt_get = function(...) "",
    cmdln = phylotaR:::cmdln_blastcheck,
    taxise_run = function(...) {
      NULL
    },
    download_run = function(...) {
      NULL
    },
    clusters_run = function(...) {
      NULL
    },
    clusters2_run = function(...) {
      NULL
    },
    .package = "phylotaR"
  )
  expect_null(res)
  res <- with_mocked_bindings(
    phylotaR:::stages_run(
      wd = wd, to = 4, frm = 1, stgs_msg = "",
      rstrt = TRUE
    ),
    taxise_run = function(...) {
      NULL
    },
    download_run = function(...) {
      NULL
    },
    clusters_run = function(...) {
      NULL
    },
    clusters2_run = function(...) {
      NULL
    },
    .package = "phylotaR"
  )
  expect_null(res)
  phylotaR:::cleanup(wd)
})
phylotaR:::cleanup(wd)
