# LIBS
library(testthat)

# VARS
blast_res_dir <- phylotaR:::datadir_get(file.path("blast", "blast_res.rda"))
blast_res <- readRDS(file = blast_res_dir)
ps <- parameters()
ps[["wd"]] <- phylotaR:::datadir_get("")

# RUNNING
context("Testing 'blast-tools'")
test_that("blastdb_gen() works", {
  sqs <- phylotaR:::testsqs_gen(n = 100)
  res <- with_mocked_bindings(
    phylotaR:::blastdb_gen(sqs = sqs, dbfl = "testdb", ps = ps),
    cmdln = function(...) 0,
    .package = "phylotaR"
  )
  expect_null(res)
  res <- with_mocked_bindings(
    expect_error(phylotaR:::blastdb_gen(sqs = sqs, dbfl = "testdb", ps = ps)),
    cmdln = function(...) 1,
    .package = "phylotaR"
  )
})
test_that("blastn_run() works", {
  res <- with_mocked_bindings(
    phylotaR:::blastn_run(dbfl = "testdb", outfl = "testblstn", ps = ps),
    cmdln = function(...) 0,
    .package = "phylotaR"
  )
  nms <- colnames(res)
  expect_true(all(nms %in% names(blast_res)))
})
test_that("outfmt_get() works", {
  res <- with_mocked_bindings(
    phylotaR:::outfmt_get(ps = ps),
    blastdb_gen = function(...) 0,
    cmdln = function(...) 0,
    .package = "phylotaR"
  )
  expect_true(inherits(res, "character"))
})
test_that("blast_filter() works", {
  blast_res <- blast_res[-1, ]
  rnd <- sample(2:nrow(blast_res), 2)
  blast_res[rnd, "qcovs"] <- 48
  res <- phylotaR:::blast_filter(blast_res = blast_res, ps = ps)
  expect_true(nrow(res) < nrow(blast_res))
  # make sure rnd and its sq/qs pair is not present
  frst <- blast_res[rnd[1], "query.id"] == res[["query.id"]] &
    blast_res[rnd[1], "subject.id"] == res[["subject.id"]]
  expect_false(any(frst))
  scnd <- blast_res[rnd[1], "subject.id"] == res[["query.id"]] &
    blast_res[rnd[1], "query.id"] == res[["subject.id"]]
  expect_false(any(scnd))
})
