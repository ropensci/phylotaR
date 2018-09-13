# LIBS
library(testthat)
library(phylotaR)

# DATA
flpth <- phylotaR:::datadir_get(file.path('taxonomy', 'txids_records.RData'))
txids_recs <- readRDS(file = flpth)
rndmn <- sample(x = seq_along(txids_recs), size = 1)
txids <- txids_recs[[rndmn]][['txids']]
recs <- txids_recs[[rndmn]][['records']]
ps <- parameters(wd = tempdir())
txdct <- phylotaR:::taxdict_gen(txids = txids, recs = recs, ps = ps)

# RUNNING
test_that('taxtree_gen() works', {
  rndm_tree <- treeman::randTree(10)
  prinds <- rndm_tree@prinds
  ids <- rndm_tree@all
  root <- rndm_tree@root
  txtree <- phylotaR:::taxtree_gen(prinds = prinds, ids = ids,
                                   root = root, ps = ps)
  expect_true(inherits(txtree, 'TreeMan'))
  expect_true(treeman::checkNdlst(txtree@ndlst, root))
})
test_that('rank_get() works', {
  txid <- sample(x = txdct@txids, size = 1)
  rank <- phylotaR:::rank_get(txid = txid, txdct = txdct)
  expect_true(is.character(rank))
})
test_that('descendants_get(direct=TRUE) works', {
  txid <- sample(x = txdct@txids, size = 1)
  dds <- phylotaR:::descendants_get(id = txid, txdct = txdct, direct = TRUE)
  expect_true(is.character(dds))
})
test_that('descendants_get(direct=FALSE) works', {
  txid <- sample(x = txdct@txids, size = 1)
  ads <- phylotaR:::descendants_get(id = txid, txdct = txdct, direct = FALSE)
  expect_true(is.character(ads))
})
test_that('parent_get() works', {
  txid <- sample(x = txdct@txids, size = 1)
  prnt <- phylotaR:::parent_get(id = txid, txdct = txdct)
  expect_true(is.numeric(as.numeric(prnt)))
})
