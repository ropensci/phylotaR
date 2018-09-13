# LIBS
library(phylotaR)
library(testthat)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
raw_recs <- readRDS(phylotaR:::datadir_get('raw_seqrecs.rda'))
sqs <- readRDS(phylotaR:::datadir_get('sqrecs.rda'))
txdct <- readRDS(phylotaR:::datadir_get('txdct.rda'))

# RUNNING
context('Testing \'stage2-tools\'')
phylotaR:::cleanup(wd)
test_that('hierarchic_download() works', {
  res <- with_mock(
    `phylotaR:::descendants_get` = function(...) sample(c(1, 10), 1),
    `phylotaR:::sqs_count` = function(...) sample(c(1, 100000000000), 1),
    `phylotaR:::seqrec_get` = function(...) 1,
    phylotaR:::hierarchic_download(txid = 1, txdct = NULL, ps = ps)
  )
  expect_true(1 %in% res)
})
test_that('seqrec_augment() works', {
  seqarc <- phylotaR:::seqrec_augment(sqs = sqs, txdct = txdct)
  expect_true(inherits(seqarc, 'SeqArc'))
})
phylotaR:::cleanup(wd)
test_that('seqrec_get() works', {
  ps[['mdlthrs']] <- 100
  phylotaR:::cache_setup(ps)
  ex_sid_list <- list('none' = NULL, 'model' = 1:101, 'normal' = 1:50)
  ex_sids <- ex_sid_list[[sample(seq_along(ex_sid_list), 1)]]
  res <- with_mock(
    `rentrez::entrez_fetch` = function(...) {
      raw_recs[[sample(seq_along(raw_recs), 1)]]
    },
    `phylotaR:::sids_get` = function(...) ex_sids,
    phylotaR:::seqrec_get(txid = 1, ps = ps, direct = FALSE, lvl = 1)
  )
  res <- vapply(X = res, FUN = function(x) inherits(x, 'SeqRec'), logical(1))
  expect_true(all(res))
})
phylotaR:::cleanup(wd)

# Example recs
# devtools::load_all('~/Coding/phylotaR')
# wd <- file.path(getwd(), 'aotus')
# dir.create(wd)
# ps <- parameters(wd = wd, txid = 9504, v = TRUE)
# cache_setup(ps)
# taxise_run(wd = wd)
# sqs <- seqrec_get(txid = 9504, ps = ps, direct = FALSE, lvl = 1)
# saveRDS(object = sqs[1:100],
#         file = phylotaR:::datadir_get('sqrecs.rda'))
# txdct <- obj_load(wd = wd, nm = 'txdct')
# saveRDS(object = txdct, file = phylotaR:::datadir_get('txdct.rda'))
# unlink(x = wd, recursive = TRUE)
