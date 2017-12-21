# LIBS
library(phylotaR)
library(testthat)

# VARS
data("phylt_nds")
data('sqs')
data('blst_rs')
ps <- list('wd'='.', 'v'=FALSE)
clstr <- list("gis"=NA, "seed_gi"=NA, "ti_root"=NA, "ci"=NA,
              "cl_type"=NA, "n_gi"=1, "tis"=NA, "n_ti"=NA,
              "MinLength"=NA, "MaxLength"=NA, "n_gen"=NA,
              "n_child"=NA, "ci_anc"=NA, "unique_id"=NA)
clstrs <- list(clstr, clstr, clstr)

# FUNCTIONS
# stubs
mckBlstN <- function(dbfl, outfl, ps) {
  blst_rs
}
mckMkBlstDB <- function(sqs, dbfl, ps) {
  NULL
}

# RUNNING
context('Testing \'cluster-tools\'')
test_that('clstrPhylt() works', {
  phylt_nds <- clstrPhylt(clstrs=clstrs)
  expect_true(nrow(phylt_nds) == 3)
})
test_that('clstrCiGi() works', {
  cigi <- clstrCiGi(clstrs=clstrs)
  expect_true(nrow(cigi) == 3)
})
test_that('clstrSqs() works', {
  ps <- list('v'=FALSE)
  txid <- sample(phylt_nds[['ti']], 1)
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    clstrSqs(txid=txid, sqs=sqs, phylt_nds=phylt_nds,
             ps=ps)
  )
  expect_true('list' %in% is(res))
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    clstrSqs(txid=9479, sqs=sqs, phylt_nds=phylt_nds,
             ps=ps)
  )
  # using the platyrrhini txid we should get more res
  expect_true(length(res) > 0)
})
test_that('calcClstrs() works', {
  
})
test_that('clstrBlstRs() works', {
  # TODO: not sure how to test this well
  clstrBlstRs(blst_rs, infrmtv=TRUE)
  clstrBlstRs(blst_rs, infrmtv=FALSE)
})
test_that('addClstrInf() works', {
  # TODO: not sure how to test this well, either
  clstrs <- clstrBlstRs(blst_rs, infrmtv=TRUE)
  res <- addClstrInf(clstrs=clstrs, phylt_nds=phylt_nds,
              txid=9479, sqs=sqs, drct=FALSE)
  expect_true(length(clstrs) == length(res))
})
test_that('getADs() works', {
  # all nodes should descend from Platyrrhini
  txids <- getADs(txid=9479, phylt_nds=phylt_nds)
  expect_true(length(txids) == (nrow(phylt_nds) - 1))
})
test_that('getGnsFrmPhyltNds() works', {
  rnd <- sample(phylt_nds[['ti_anc']], 1)
  res <- getGnsFrmPhyltNds(txid=rnd, phylt_nds=phylt_nds)
  expect_true(res %in% phylt_nds[['ti_genus']])
})
test_that('blstSqs() works', {
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    blstSqs(txid=1, typ='direct', sqs=sqs, ps=ps)
  )
  expect_true('data.frame' %in% is(res))
})
