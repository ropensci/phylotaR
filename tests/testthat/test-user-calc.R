# LIBS
library(phylotaR)
library(testthat)

# DATA
phylota <- phylotaR:::random_phylota()

# RUNNING
context('Testing \'user-calc\'')
test_that('calc_mad() works', {
  cid <- sample(phylota@cids, 5)
  mad_scrs <- calc_mad(phylota = phylota, cid = cid)
  expect_true(inherits(mad_scrs, 'numeric'))
  expect_true(all(mad_scrs >= 0 & mad_scrs <= 1))
})
test_that('calc_wrdfrq() works', {
  cid <- sample(phylota@cids, 5)
  wrd_frqs <- calc_wrdfrq(phylota = phylota, cid = cid, type = 'dfln')
  expect_true(inherits(unlist(wrd_frqs), 'numeric'))
  wrd_frqs <- calc_wrdfrq(phylota = phylota, cid = cid, type = 'nm')
  expect_true(inherits(unlist(wrd_frqs), 'numeric'))
})
