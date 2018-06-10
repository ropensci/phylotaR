# # LIBS
# library(testthat)
# library(phylotaR)
# 
# # DATA
# test_data_path <- file.path('tests', 'testthat', 'data',
#                             'taxonomy', 'txids_records.RData')
# txids_records <- readRDS(file = test_data_path)
# rndmn <- sample(x = seq_along(txids_records), size = 1)
# txids <- txids_records[[rndmn]][['txids']]
# rcrds <- txids_records[[rndmn]][['records']]
# dictionary <- phylotaR:::genTxDct(txids = txids, rcrds = rcrds)
# 
# # RUNNING
# test_that('taxonomic_tree_generate() works', {
#   rndm_tree <- treeman::randTree(10)
#   prinds <- rndm_tree@prinds
#   trids <- rndm_tree@all
#   root <- rndm_tree@root
#   txtree <- phylotaR:::taxonomic_tree_generate(prinds = prinds,
#                                                trids = trids,
#                                                root = root)
#   expect_true(inherits(txtree, 'TreeMan'))
#   expect_true(treeman::checkNdlst(txtree@ndlst, root))
# })
# test_that('trid_to_txid() works', {
#   trids <- sample(x = dictionary@indx, size = 5)
#   txids <- phylotaR:::trid_to_txid(trid = trids, dictionary = dictionary)
#   res <- phylotaR:::txid_to_trid(txid = txids, dictionary = dictionary)
#   expect_true(all(trids == res))
# })
# test_that('txid_to_trid() works', {
#   txids <- sample(x = dictionary@txids, size = 5)
#   trids <- phylotaR:::txid_to_trid(txid = txids, dictionary = dictionary)
#   res <- phylotaR:::trid_to_txid(trid = trids, dictionary = dictionary)
#   expect_true(all(txids == res))
# })
# test_that('txid_singletons_get() works', {
#   txid <- sample(x = dictionary@txids, size = 1)
#   sngltons <- phylotaR:::txid_singletons_get(txid = txid,
#                                              dictionary = dictionary)
#   expect_true(is.character(sngltons))
# })
# test_that('txid_rank_get() works', {
#   txid <- sample(x = dictionary@txids, size = 1)
#   rank <- phylotaR:::txid_rank_get(txid = txid,
#                                        dictionary = dictionary)
#   expect_true(is.character(rank))
# })
# test_that('txid_descendants_get(direct=TRUE) works', {
#   txid <- sample(x = dictionary@txids, size = 1)
#   dds <- phylotaR:::txid_descendants_get(txid = txid,
#                                          dictionary = dictionary,
#                                          direct = TRUE)
#   expect_true(is.character(dds))
# })
# test_that('txid_descendants_get(direct=FALSE) works', {
#   txid <- sample(x = dictionary@txids, size = 1)
#   ads <- phylotaR:::txid_descendants_get(txid = txid,
#                                          dictionary = dictionary,
#                                          direct = FALSE)
#   expect_true(is.character(dds))
# })
# test_that('txid_parent_get() works', {
#   txid <- sample(x = dictionary@txids, size = 1)
#   prnt <- phylotaR:::txid_parent_get(txid = txid,
#                                      dictionary = dictionary)
#   expect_true(is.numeric(as.numeric(dds)))
# })
