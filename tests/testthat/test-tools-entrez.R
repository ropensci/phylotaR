# LIBS
library(testthat)

# DATA ----
wd <- tempdir()
ps <- parameters(wd = wd)
n <- integer(1)
#9607

# STUBS ----
esearch_mock <- function(db, term, retmax, retstart, ...) {
  if (n == 0) {
    ids <- NULL
  } else {
    ids <- as.character(1:n)
  }
  list('ids' = ids, 'count' = n)
}
efetch_mock <- function(retmax, ...) {
  # when rettype == 'acc'
  paste0(rep('seqid1', retmax), collapse = '\n')
}

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'entrez-tools\'')
test_that('searchterm_gen() works', {
  trm <- phylotaR:::searchterm_gen(txid = '9606', ps = ps, direct = FALSE)
  expect_true(is.character(trm))
  trm <- phylotaR:::searchterm_gen(txid = '9606', ps = ps, direct = TRUE)
  expect_true(is.character(trm))
})
test_that('seqs_count() works', {
  phylotaR:::cache_setup(ps = ps)
  res <- with_mock(
    `phylotaR:::safely_connect` = function(func, args, fnm, ps){
      list('count' = args)
    },
    phylotaR:::sqs_count(txid = '9606', direct = FALSE, ps = ps)
  )
  expect_true(grepl(':exp', res[['term']]))
  res <- with_mock(
    `phylotaR:::safely_connect` = function(func, args, fnm, ps){
      list('count' = args)
    },
    phylotaR:::sqs_count(txid = '9606', direct = TRUE, ps = ps)
  )
  expect_true(grepl(':noexp', res[['term']]))
})
phylotaR:::cleanup(wd)
test_that('txnds_count() works', {
  phylotaR:::cache_setup(ps = ps)
  res <- with_mock(
    `phylotaR:::safely_connect` = function(func, args, fnm, ps){
      list('count' = args)
    },
    phylotaR:::txnds_count(txid = '9606', ps = ps)
  )
  expect_true(res[['term']] == "txid9606[Subtree]")
})
phylotaR:::cleanup(wd)
test_that('sids_get() works', {
  phylotaR:::cleanup(wd)
  phylotaR:::cache_setup(ps = ps)
  n <<- 100
  res <- with_mock(
    `rentrez::entrez_search` = esearch_mock,
    `rentrez::entrez_fetch` = efetch_mock,
    phylotaR:::sids_get(txid = '9606', direct = FALSE, ps = ps)
  )
  expect_true(length(res) == n)
  phylotaR:::cleanup(wd)
  phylotaR:::cache_setup(ps = ps)
  n <<- 100
  res <- with_mock(
    `rentrez::entrez_search` = esearch_mock,
    `rentrez::entrez_fetch` = efetch_mock,
    phylotaR:::sids_get(txid = '9606', direct = FALSE, ps = ps, hrdmx = 20,
                        retmax = 10)
  )
  expect_true(length(res) != n)
})
phylotaR:::cleanup(wd)