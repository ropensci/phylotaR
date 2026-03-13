# LIBS
library(phylotaR)
library(testthat)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
exclstrarc <- phylotaR:::clstrarc_gen(list())
sqs <- readRDS(phylotaR:::datadir_get("sqrecs.rda"))
sqs <- phylotaR:::seqarc_gen(sqs)
blastres_flpth <- phylotaR:::datadir_get(file.path("blast", "blast_res.rda"))
blast_res <- readRDS(blastres_flpth)
# adapt sids for blast res
sids <- unique(c(blast_res$subject.id, blast_res$query.id))
blast_res_sqs <- sample(sqs@sqs, length(sids))
blast_res_sqs <- phylotaR:::seqarc_gen(blast_res_sqs)
blast_res_sqs@ids <- sids

# RUNNING
phylotaR:::cleanup(wd)
context("Testing 'stage3-tools'")
test_that("blast_sqs() works", {
  phylotaR:::cache_setup(ps)
  res <- with_mocked_bindings(
    phylotaR:::blast_sqs(
      txid = "1", typ = "direct", sqs = sqs, ps = ps,
      lvl = 1
    ),
    blastn_run = function(...) blast_res,
    blastdb_gen = function(...) NULL,
    .package = "phylotaR"
  )
  expect_true("data.frame" %in% is(res))
})
phylotaR:::cleanup(wd)
test_that("clstr_sqs() works", {
  phylotaR:::cache_setup(ps)
  ps[["v"]] <- TRUE
  res <- with_mocked_bindings(
    phylotaR:::clstr_sqs(
      txid = "1", sqs = blast_res_sqs, ps = ps, lvl = 0,
      typ = "subtree"
    ),
    blast_sqs = function(...) blast_res,
    .package = "phylotaR"
  )
  expect_true(inherits(res, "ClstrArc"))
})
phylotaR:::cleanup(wd)
# no cache tests
test_that("clstr_all() works", {
  mock_dscdnts_get <- function(id, ...) {
    if (id == 1) {
      return(NULL)
    }
    return(1)
  }
  res <- with_mocked_bindings(
    phylotaR:::clstr_all(txid = 0, txdct = NULL, sqs = sqs, ps = ps, lvl = 0),
    clstr_subtree = function(...) exclstrarc,
    descendants_get = mock_dscdnts_get,
    .package = "phylotaR"
  )
  expect_true(inherits(res, "ClstrArc"))
})
test_that("clstr_subtree() works", {
  res <- with_mocked_bindings(
    phylotaR:::clstr_subtree(
      txid = 9479, sqs = sqs, dds = 1, ps = ps,
      lvl = 0, txdct = NULL
    ),
    clstr_sqs = function(...) exclstrarc,
    clstr_direct = function(...) exclstrarc,
    rank_get = function(...) "species",
    descendants_get = function(...) "",
    .package = "phylotaR"
  )
  expect_true(inherits(res, "ClstrArc"))
})
test_that("clstr_direct() works", {
  res <- with_mocked_bindings(
    phylotaR:::clstr_direct(
      txid = 9479, sqs = sqs, ps = ps, lvl = 0,
      txdct = NULL
    ),
    clstr_sqs = function(...) exclstrarc,
    rank_get = function(...) "species",
    .package = "phylotaR"
  )
  expect_true(inherits(res, "ClstrArc"))
})
test_that("blast_clstr() works", {
  res <- phylotaR:::blast_clstr(blast_res = blast_res)
  expect_true(length(res) == 1)
  nrws <- blast_res[1:2, ]
  nrws[["query.id"]] <- 1:2
  nrws[["subject.id"]] <- 2:1
  res <- phylotaR:::blast_clstr(blast_res = rbind(blast_res, nrws))
  expect_true(length(res) == 2)
})
test_that("clstrarc_gen() works", {
  # TODO: use more realistic cluster records
  clstrrecs <- vector(mode = "list", length = 10)
  names(clstrrecs) <- paste0(seq_along(clstrrecs))
  clstrarc <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
  expect_true(inherits(x = clstrarc, what = "ClstrArc"))
})
test_that("clstrarc_join() works", {
  clstrrecs <- vector(mode = "list", length = 10)
  names(clstrrecs) <- paste0(seq_along(clstrrecs))
  clstrarc1 <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
  names(clstrrecs) <- paste0(seq_along(clstrrecs) + length(clstrrecs))
  clstrarc2 <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
  joined <- phylotaR:::clstrarc_join(
    clstrarc_1 = clstrarc1,
    clstrarc_2 = clstrarc2
  )
  expect_true(inherits(x = joined, what = "ClstrArc"))
})
