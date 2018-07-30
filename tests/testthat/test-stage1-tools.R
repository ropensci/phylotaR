# LIBS
library(phylotaR)
library(testthat)

# DATA
raw_recs <- readRDS(phylotaR:::datadir_get('raw_txrecs.rda'))
wd <- tempdir()
ps <- parameters(wd = wd)

# RUNNING
context('Testing \'stage1-tools\'')
phylotaR:::cleanup(wd)
test_that('tax_download() works', {
  example_res <- raw_recs[[sample(seq_along(raw_recs), 1)]]
  res <- with_mock(
    `phylotaR:::search_and_cache` = function(...) example_res,
    phylotaR:::tax_download(ids = NULL, ps = ps)
  )
  res <- vapply(X = res, FUN = function(x) inherits(x, 'TaxRec'),
                FUN.VALUE = logical(1))
  expect_true(all(res))
})
phylotaR:::cleanup(wd)


# # Example raw_txrecs ----
# search_obj <- rentrez::entrez_search(db = 'taxonomy', retmax = 0,
#                                      term = 'txid2759[Subtree]',
#                                      use_history = TRUE)
# # random 100 taxrecords
# retmaxes <- round(runif(n = 10, min = 1,
#                         max = search_obj[['count']]))
# raw_recs <- NULL
# for (i in retmaxes) {
#   whobj <- search_obj[['web_history']]
#   raw_rec <- rentrez::entrez_fetch(db = 'taxonomy',
#                                    retstart = i, rettype = 'xml',
#                                    web_history = whobj,
#                                    retmax = sample(1:4, 1))
#   raw_recs <- c(raw_recs, raw_rec)
# }
# otfl <- phylotaR:::datadir_get('raw_txrecs.rda')
# saveRDS(object = raw_recs, file = otfl)
