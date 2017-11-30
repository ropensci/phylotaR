# LIBS
library(phylotaR)
library(testthat)

# VARS
data("phylt_nds")

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
