# LIBS
library(phylotaR)
library(testthat)

# DATA
phylota <- phylotaR:::random_phylota()

# RUNNING
context('Testing \'user-drop\'')
test_that('drop_sqs() works', {
  nsqs <- get_clstr_slot(phylota = phylota, cid = phylota@cids, slt_nm = 'nsq')
  pssbls <- names(nsqs)[nsqs > 3]
  cid <- sample(pssbls, 1)
  sids <- phylota@clstrs@clstrs[[cid]]@sids
  sid <- sample(sids, 1)
  res <- drop_sqs(phylota = phylota, cid = cid, sid = sid)
  expect_true(length(res@clstrs@clstrs[[cid]]@sids) == 1)
  # total seqs won't necessarily be affected, seqs could be represented in
  # other clusters
  expect_true(length(phylota@sids) >= length(res@sids))
})
test_that('drop_clstrs() works', {
  cid <- sample(phylota@cids, 1)
  res <- drop_clstrs(phylota = phylota, cid = cid)
  expect_true(length(phylota@sids) > length(res@sids))
  expect_true(length(res@cids) == 1)
})
test_that('drop_by_rank() works', {
  ncids <- ifelse(length(phylota@cids) >= 5, 5, length(phylota@cids))
  cids <- sample(phylota@cids, ncids)
  small_phylota <- drop_clstrs(phylota = phylota, cid = cids)
  res <- drop_by_rank(phylota = small_phylota, rnk = 'genus', n = 2,
                      keep_higher = FALSE,
                      choose_by = c('pambgs', 'age', 'nncltds'),
                      greatest = c(FALSE, FALSE, TRUE))
  expect_true(length(small_phylota@sids) > length(res@sids))
  # # most likely to be multiples of two
  # nsqs <- vapply(X = res@clstrs@clstrs, FUN = function(x) length(x@sids),
  #                FUN.VALUE = integer(1))
  # nrms <- nsqs %% 2
  # uniqs <- unique(nrms)
  # nrm_mode <- uniqs[which.max(tabulate(match(nrms, uniqs)))]
  # expect_true(nrm_mode == 0)
})
