# LIBS
library(testthat)
library(phylotaR)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
sqs <- readRDS(phylotaR:::datadir_get('sqrecs.rda'))
sqs <- phylotaR:::seqarc_gen(sqs)
blastres_flpth <- phylotaR:::datadir_get(file.path('blast', 'blast_res.rda'))
blast_res <- readRDS(blastres_flpth)
exclstrarc <- readRDS(phylotaR:::datadir_get('clstrarc.rda'))
# adapt sids for blast res
sids <- unique(c(blast_res$subject.id, blast_res$query.id))
blast_res_sqs <- sample(sqs@sqs, length(sids))
blast_res_sqs <- phylotaR:::seqarc_gen(blast_res_sqs)
blast_res_sqs@ids <- sids

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'stage4-tools\'')
test_that('clstrs_join() works', {
  clstrs_jnd <- phylotaR:::clstrs_join(ps = ps, seed_ids = sids,
                                       blast_res = blast_res,
                                       all_clstrs = exclstrarc@clstrs)
  expect_true(clstrs_jnd[[1]][['typ']] == 'merged')
})
test_that('clstrs_merge() works', {
  clstrs_jnd <- phylotaR:::clstrs_join(ps = ps, seed_ids = sids,
                                       blast_res = blast_res,
                                       all_clstrs = exclstrarc@clstrs)
  res <- with_mock(
    `phylotaR:::parent_get` = function(...) '1',
    phylotaR:::clstrs_merge(jnd_clstrs = clstrs_jnd, txdct = NULL)
  )
  expect_true(inherits(res[[1]], 'ClstrRec'))
})
test_that('clstrs_renumber() works', {
  clstrarc <- phylotaR:::clstrs_renumber(exclstrarc@clstrs)
  topclstr <- clstrarc[['0']]
  bttmclstr <- clstrarc[['9']]
  expect_true(topclstr@nsqs > bttmclstr@nsqs)
})
test_that('seeds_blast() works', {
  res <- with_mock(
    `phylotaR:::blastn_run` = function(...) blast_res,
    `phylotaR:::blastdb_gen` = function(...) NULL,
    phylotaR:::seeds_blast(sqs = sqs, ps = ps)
  )
  expect_true(inherits(res, 'data.frame'))
})
phylotaR:::cleanup(wd)
