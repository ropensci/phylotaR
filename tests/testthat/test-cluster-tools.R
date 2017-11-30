# LIBS
library(phylotaR)
library(testthat)

# VARS
data("phylt_nds")
data('sqs')
data('blst_rs')


# FUNCTIONS
# stubs
mckBlstN <- function(dbfl, outfl, wd,
                     eval_ctoff=1.0e-10,
                     verbose=FALSE) {
  blst_rs
}
mckMkBlstDB <- function(sqs, dbfl, wd,
                        verbose=FALSE) {
  NULL
}

# RUNNING
context('Testing \'cluster-tools\'')
test_that('getSeqs() works', {
  
})
test_that('clstrIds() works', {
  
})
test_that('fltrIds() works', {
  
})
test_that('getADs() works', {
  # all nodes should descend from Platyrrhini
  txids <- getADs(txid=9479, phylt_nds=phylt_nds)
  expect_true(length(txids) == (nrow(phylt_nds) - 1))
})
test_that('writeClstr() works', {
  
})
test_that('blstSqs() works', {
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    blstSqs(txid=1, typ='direct',
            sqs=sqs, wd='.', verbose=FALSE)
  )
})
