# LIBS
library(phylotaR)
library(testthat)

# DATA
ps <- parameters()
sqs <- readRDS(phylotaR:::datadir_get('sqrecs.rda'))
sqs <- phylotaR:::seqarc_gen(sqs)
blast_res <- readRDS(phylotaR:::datadir_get(file.path('blast',
                                                      'blast_res.rda')))

# RUNNING
phylotaR:::cleanup()
context('Testing \'stage3-tools\'')
test_that('blast_sqs() works', {
  phylotaR:::cache_setup(ps)
  res <- with_mock(
    `phylotaR:::blastn_run` = function(...) blast_res,
    `phylotaR:::blastdb_gen` = function(...) NULL,
    phylotaR:::blast_sqs(txid = '1', typ = 'direct', sqs = sqs, ps = ps,
                         lvl = 0)
  )
  expect_true('data.frame' %in% is(res))
})
phylotaR:::cleanup()
test_that('clstrSqs() works', {
  # setUpCch(ps=ps)
  # txid <- sample(phylt_nds[['ti']], 1)
  # res <- with_mock(
  #   `phylotaR::blstN`=mckBlstN,
  #   `phylotaR::mkBlstDB`=mckMkBlstDB,
  #   clstrSqs(txid=txid, sqs=sqs, phylt_nds=phylt_nds,
  #            ps=ps, lvl=0)
  # )
  # expect_true('list' %in% is(res))
  # res <- with_mock(
  #   `phylotaR::blstN`=mckBlstN,
  #   `phylotaR::mkBlstDB`=mckMkBlstDB,
  #   clstrSqs(txid=9479, sqs=sqs, phylt_nds=phylt_nds,
  #            ps=ps, lvl=0)
  # )
  # # using the platyrrhini txid we should get more res
  # expect_true(length(res) > 0)
})
phylotaR:::cleanup()
# no cache tests
test_that('clstrAll() works', {
  # rndn <- round(runif(min=5, max=10, n=1))
  # res <- with_mock(
  #   `phylotaR::getDDFrmPhyltNds`=function(txid, ...) {
  #     if(txid == 0) return(rep(1, rndn - 1))
  #     NULL
  #     },
  #   `phylotaR::clstrSbtr`=function(...) list(clstr),
  #   clstrAll(txid=0, sqs=sqs, phylt_nds=phylt_nds,
  #            ps=ps, lvl=0)
  # )
  # expect_true(length(res) == rndn)
})
test_that('clstrSbtr() works', {
  # # the parent should represent all sequences
  # res <- with_mock(
  #   `phylotaR::clstrSqs`=function(sqs, ...) {
  #     sapply(sqs, function(x) x[['gi']])
  #   },
  #   clstrSbtr(txid=9479, sqs=sqs,
  #             phylt_nds=phylt_nds,
  #             dds=NULL, ps=ps, lvl=0)
  # )
  # expect_true(length(res) == length(sqs))
  # # a random should repr. equal or fewer seqs.
  # txid <- sample(phylt_nds[['ti']], 1)
  # res <- with_mock(
  #   `phylotaR::clstrSqs`=function(sqs, ...) {
  #     sapply(sqs, function(x) x[['gi']])
  #   },
  #   clstrSbtr(txid=txid, sqs=sqs,
  #             phylt_nds=phylt_nds,
  #             dds=NULL, ps=ps, lvl=0)
  # )
  # expect_true(length(res) <= length(sqs))
})
test_that('clstrDrct() works', {
  # # expect empty list because 0 is not in phylt_nds
  # res <- with_mock(
  #   `phylotaR::clstrSqs`=function(...) NA,
  #   clstrDrct(txid=0, sqs=sqs,
  #             phylt_nds=phylt_nds, ps=ps, lvl=0)
  # )
  # expect_true(length(res) == 0)
  # # choose a txid from sqs
  # # should run clstrSqs -- returning TRUE
  # rnd <- sample(1:length(sqs), 1)
  # txid <- sqs[[rnd]][['ti']]
  # res <- with_mock(
  #   `phylotaR::clstrSqs`=function(...) TRUE,
  #   clstrDrct(txid=txid, sqs=sqs,
  #             phylt_nds=phylt_nds,
  #             ps=ps, lvl=0)
  # )
  # expect_true(res)
})
test_that('blast_clstr() works', {
  res <- phylotaR:::blast_clstr(blast_res = blast_res)
  expect_true(length(res) == 1)
  nrws <- blast_res[1:2, ]
  nrws[['query.id']] <- 1:2
  nrws[['subject.id']] <- 2:1
  res <- phylotaR:::blast_clstr(blast_res = rbind(blast_res, nrws))
  expect_true(length(res) == 2)
})
test_that('clstrarc_gen() works', {
  # TODO: use more realistic cluster records
  clstrrecs <- vector(mode = 'list', length = 10)
  names(clstrrecs) <- paste0(seq_along(clstrrecs))
  clstrarc <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
  expect_true(inherits(x = clstrarc, what = 'ClstrArc'))
})
test_that('clstrarc_join() works', {
  clstrrecs <- vector(mode = 'list', length = 10)
  names(clstrrecs) <- paste0(seq_along(clstrrecs))
  clstrarc1 <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
  names(clstrrecs) <- paste0(seq_along(clstrrecs) + length(clstrrecs))
  clstrarc2 <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
  joined <- phylotaR:::clstrarc_join(clstrarc_1 = clstrarc1,
                                     clstrarc_2 = clstrarc2)
  expect_true(inherits(x = joined, what = 'ClstrArc'))
})
