# LIBS
library(phylotaR)
library(testthat)

# DATA
raw_recs <- readRDS(phylotaR:::datadir_get('raw_seqrecs.rda'))
wd <- tempdir()
ps <- parameters(wd = wd)

# RUNNING
context('Testing \'tools-sqs\'')
test_that('seqrec_convert() works', {
  raw_rec <- raw_recs[[sample(seq_along(raw_recs), 1)]]
  sqs <- phylotaR:::seqrec_convert(raw_recs = raw_rec, ps = ps)
  expect_true(inherits(sqs[[1]], 'SeqRec'))
})
test_that('seqrec_gen() works', {
  seqrec <- phylotaR:::seqrec_gen(accssn = '1', nm = '1', txid = '1',
                                  sq = 'atcg', dfln = 'defline',
                                  orgnsm = 'G. sp', ml_typ = 'DNA',
                                  rec_typ = 'Full', vrsn = '1.1',
                                  age = 30L)
  expect_true(inherits(seqrec, 'SeqRec'))
})
test_that('seqarc_gen() works', {
  seqrec1 <- phylotaR:::seqrec_gen(accssn = '1', nm = '1', txid = '1',
                                   sq = 'atcg', dfln = 'defline',
                                   orgnsm = 'G. sp', ml_typ = 'DNA',
                                   rec_typ = 'Full', vrsn = '1.1',
                                   age = 30L)
  seqrec2 <- phylotaR:::seqrec_gen(accssn = '2', nm = '1', txid = '1',
                                   sq = 'atcg', dfln = 'defline',
                                   orgnsm = 'G. sp', ml_typ = 'DNA',
                                   rec_typ = 'Full', vrsn = '2.1',
                                   age = 30L)
  seqarc <- phylotaR:::seqarc_gen(seqrecs = list(seqrec1, seqrec2))
  expect_true(inherits(seqarc, 'SeqArc'))
})

# # Example raw_recs ----
# term <- phylotaR:::searchterm_gen(txid = '2759', ps = parameters(),
#                                   direct = FALSE)
# search_obj <- rentrez::entrez_search(db = 'nucleotide', term = term,
#                                      retmax = 0, use_history = TRUE)
# # random 10 sequences
# retmaxes <- round(runif(n = 10, min = 1,
#                         max = search_obj[['count']]))
# raw_recs <- NULL
# for (i in retmaxes) {
#   whobj <- search_obj[['web_history']]
#   raw_rec <- rentrez::entrez_fetch(db = 'nucleotide',
#                                    retstart = i, retmode = 'xml',
#                                    rettype = 'gbwithparts',
#                                    web_history = whobj,
#                                    retmax = sample(1:4, 1))
#   raw_recs <- c(raw_recs, raw_rec)
# }
# otfl <- phylotaR:::datadir_get('raw_recs.rda')
# saveRDS(object = raw_recs, file = otfl)
