# LIBS
library(testthat)
library(phylotaR)

# RUNNING
context('Testing \'phylotaR\'')
test_that('parameters() works', {
  ps <- parameters()
  expect_true(inherits(ps, 'list'))
})
test_that('list_ncbi_ranks() works', {
  ranks <- list_ncbi_ranks()
  expect_true('species' %in% ranks)
})
test_that('list_seqrec_slots() works', {
  slts <- list_seqrec_slots()
  expect_true('id' %in% slts)
})
test_that('list_clstrrec_slots() works', {
  slts <- list_clstrrec_slots()
  expect_true('id' %in% slts)
})
test_that('list_taxrec_slots() works', {
  slts <- list_taxrec_slots()
  expect_true('id' %in% slts)
})
