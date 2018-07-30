# LIBS
library(phylotaR)
library(testthat)

# DATA
wd <- tempdir()
phylota <- phylotaR:::random_phylota()

# RUNNING
context('Testing \'user-get\'')
test_that('get_ntaxa() works', {
  ntxs <- get_clstr_slot(phylota = phylota, cid = phylota@cids, slt_nm = 'ntx')
  pssbls <- names(ntxs)[ntxs > 3]
  cid <- sample(pssbls, 1)
  res1 <- get_ntaxa(phylota = phylota, cid = cid, sid = NULL, rnk = NULL,
                    keep_higher = FALSE)
  expect_true(ntxs[[cid]] == res1)
  res2 <- get_ntaxa(phylota = phylota, cid = cid, sid = NULL, rnk = 'family',
                    keep_higher = FALSE)
  expect_true(res1 >= res2)
  res3 <- get_ntaxa(phylota = phylota, cid = cid, sid = NULL,
                    rnk = 'subspecies', keep_higher = TRUE)
  expect_true(res1 >= res3)
  res4 <- get_ntaxa(phylota = phylota, cid = NULL, sid = phylota[[cid]]@sids,
                    rnk = NULL, keep_higher = TRUE)
  expect_true(res1 == res4)
})
test_that('get_txids() works', {
  ntxs <- get_clstr_slot(phylota = phylota, cid = phylota@cids, slt_nm = 'ntx')
  pssbls <- names(ntxs)[ntxs > 3]
  cid <- sample(pssbls, 1)
  res1 <- get_txids(phylota = phylota, cid = cid, sid = NULL, rnk = NULL,
                    keep_higher = FALSE)
  expect_true(ntxs[[cid]] == length(unique(res1)))
  res2 <- get_txids(phylota = phylota, cid = cid, sid = NULL, rnk = 'family',
                    keep_higher = FALSE)
  expect_true(length(unique(res1)) >= length(unique(res2)))
  res3 <- get_txids(phylota = phylota, cid = cid, sid = NULL,
                    rnk = 'subspecies', keep_higher = TRUE)
  expect_true(length(unique(res1)) >= length(unique(res3)))
  res4 <- get_txids(phylota = phylota, cid = NULL, sid = phylota[[cid]]@sids,
                    rnk = NULL, keep_higher = TRUE)
  expect_true(all(res1 %in% res4))
})
test_that('get_nsqs() works', {
  nsq <- get_nsqs(phylota = phylota, cid = sample(phylota@cids, 1))
  expect_true(inherits(nsq[[1]], 'integer'))
})
test_that('get_sq_slot() works', {
  sltnm <- sample(list_seqrec_slots(), 1)
  sq <- sample(phylota@sqs@sqs, 1)[[1]]
  res <- get_sq_slot(phylota = phylota, sid = sq@id, slt_nm = sltnm)
  expect_equal(slot(sq, sltnm), res[[1]])
  cid <- sample(phylota@cids, 1)
  res <- get_sq_slot(phylota = phylota, cid = cid, slt_nm = sltnm)
  expect_true(length(res) > 1)
})
test_that('get_clstr_slot() works', {
  sltnm <- sample(list_clstrrec_slots(), 1)
  clstr <- sample(phylota@clstrs@clstrs, 1)[[1]]
  res <- get_clstr_slot(phylota = phylota, cid = clstr@id, slt_nm = sltnm)
  expect_equal(slot(clstr, sltnm), res[[1]])
})
test_that('get_tx_slot() works', {
  sltnm <- sample(list_taxrec_slots(), 1)
  txid <- sample(phylota@txids, 1)[[1]]
  res <- get_tx_slot(phylota = phylota, txid = txid, slt_nm = sltnm)
  rec <- phylota@txdct@recs[[txid]]
  expect_equal(slot(rec, sltnm), res[[1]])
})
phylotaR:::cleanup(wd)
test_that('get_stage_times() works', {
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
  timings <- get_stage_times(wd = wd)
  expect_true(inherits(timings, 'numeric'))
})
phylotaR:::cleanup(wd)
