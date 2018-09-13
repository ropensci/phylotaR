# LIBS
library(phylotaR)
library(testthat)

# DATA
# from the wild
raw_recs <- readRDS(phylotaR:::datadir_get('raw_seqrecs.rda'))
# selected to represent record diversity
test_raw_recs <- readRDS(phylotaR:::datadir_get('test_raw_seqrecs.rda'))
wd <- tempdir()
ps <- parameters(wd = wd)

# RUNNING
context('Testing \'tools-sqs\'')
test_that('gb_extract() works', {
  for (test_raw_rec in test_raw_recs) {
    record_parts <- phylotaR:::gb_extract(record = test_raw_rec)
    test <- names(record_parts) %in% c('accession', 'version', 'definition',
                                       'sequence', 'moltype', 'features',
                                       'date', 'organism', 'keywords')
    expect_true(all(test))
  }
})
test_that('rawseqrec_breakdown() works', {
  for (test_raw_rec in test_raw_recs) {
    record_parts <- phylotaR:::gb_extract(record = test_raw_rec)
    seqrecs <- phylotaR:::rawseqrec_breakdown(record_parts = record_parts,
                                              ps = ps)
    expect_true(inherits(seqrecs[[1]], 'SeqRec'))
  }
})
test_that('seqrec_convert() works', {
  for (raw_rec in raw_recs) {
    sqs <- phylotaR:::seqrec_convert(raw_recs = raw_rec, ps = ps)
    expect_true(inherits(sqs[[1]], 'SeqRec'))
  }
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
# # random 100 sequences
# retmaxes <- round(runif(n = 100, min = 1, max = search_obj[['count']]))
# raw_recs <- NULL
# for (i in retmaxes) {
#   print(i)
#   whobj <- search_obj[['web_history']]
#   raw_rec <- rentrez::entrez_fetch(db = 'nucleotide', retstart = i,
#                                    retmode = 'text', rettype = 'gb',
#                                    web_history = whobj,
#                                    retmax = sample(1:4, 1))
#   raw_recs <- c(raw_recs, raw_rec)
# }
# otfl <- phylotaR:::datadir_get('raw_seqrecs.rda')
# saveRDS(object = raw_recs, file = otfl)
