# LIBS
library(testthat)

# DATA
data('blst_rs')

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}
# stubs
mckBlstN <- function(dbfl, outfl, ps) {
  blst_rs
}
mckMkBlstDB <- function(sqs, dbfl, ps) {
  NULL
}

# RUNNING
cleanUp()
context('Testing \'stage-tools-cluster\'')
test_that('blstSqs() works', {
  # setUpCch(ps=ps)
  # res <- with_mock(
  #   `phylotaR::blstN`=mckBlstN,
  #   `phylotaR::mkBlstDB`=mckMkBlstDB,
  #   blstSqs(txid=1, typ='direct', sqs=sqs, ps=ps, lvl=0)
  # )
  # expect_true('data.frame' %in% is(res))
})
cleanUp()
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
cleanUp()
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
test_that('clstrBlstRs() works', {
  res <- phylotaR:::clstrBlstRs(blst_rs)
  expect_true(length(res) == 1)
  nrws <- blst_rs[1:2, ]
  nrws[['query.id']] <- 1:2
  nrws[['subject.id']] <- 2:1
  res <- phylotaR:::clstrBlstRs(rbind(blst_rs, nrws))
  expect_true(length(res) == 2)
})