# LIBS
library(phylotaR)
library(testthat)

# VARS
data("phylt_nds")
data('sqs')
data('blst_rs')
ps <- list('wd'='.', 'v'=FALSE)

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
  
})
test_that('clstrCiGi() works', {
  
})
test_that('clstrSqs() works', {
  
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
