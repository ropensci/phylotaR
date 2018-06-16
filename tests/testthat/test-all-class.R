# LIBS
library(phylotaR)
library(testthat)

# DATA
phylota <- phylotaR:::random_phylota()

# RUNNING
context('Testing \'all-class\'')
test_that('SeqRec() works', {
  seqrec <- phylota@sqs@sqs[[1]]
  show(seqrec)
  print(seqrec)
  str(seqrec)
  summary(seqrec)
  seqstrng <- as.character(seqrec)
  expect_true(validObject(seqrec))
})
test_that('SeqArc() works', {
  seqarc <- phylota@sqs
  show(seqarc)
  print(seqarc)
  str(seqarc)
  summary(seqarc)
  sqstrng <- as.character(seqarc)
  expect_true(validObject(seqarc))
  # [ and [[
  seqrec <- seqarc[[seqarc@ids[[1]]]]
  expect_true(inherits(seqrec, 'SeqRec'))
  seqarc2 <- seqarc[seqarc@ids[1:10]]
  expect_true(inherits(seqarc2, 'SeqArc'))
})
test_that('ClstrRec() works', {
  clstrrec <- phylota@clstrs@clstrs[[1]]
  show(clstrrec)
  print(clstrrec)
  str(clstrrec)
  summary(clstrrec)
  clstrng <- as.character(clstrrec)
  expect_true(validObject(clstrrec))
  # [ and [[
  expect_true(validObject(clstrrec))
})
test_that('ClstrArc() works', {
  clstrarc <- phylota@clstrs
  show(clstrarc)
  print(clstrarc)
  str(clstrarc)
  summary(clstrarc)
  clstrng <- as.character(clstrarc)
  expect_true(validObject(clstrarc))
  # [ and [[
  clstrrec <- clstrarc[[clstrarc@ids[[1]]]]
  expect_true(inherits(clstrrec, 'ClstrRec'))
  clstrarc2 <- clstrarc[clstrarc@ids[1:10]]
  expect_true(inherits(clstrarc2, 'ClstrArc'))
})
test_that('TaxRec() works', {
  txid <- phylota@txids[[1]]
  taxrec <- phylota@txdct@recs[[txid]]
  show(taxrec)
  print(taxrec)
  str(taxrec)
  summary(taxrec)
  txstrng <- as.character(taxrec)
  expect_true(validObject(taxrec))
})
test_that('TaxDct() works', {
  txdct <- phylota@txdct
  show(txdct)
  print(txdct)
  str(txdct)
  summary(txdct)
  txstrng <- as.character(txdct)
  expect_true(validObject(txdct))
})
test_that('Phylota() works', {
  show(phylota)
  print(phylota)
  str(phylota)
  cids <- sample(phylota@cids, 2)
  summary(drop_clstrs(phylota=phylota, cid=cids))
  txstrng <- as.character(phylota)
  expect_true(validObject(phylota))
  # [[
  sid <- sample(phylota@sids, 1)
  expect_true(inherits(phylota[[sid]], 'SeqRec'))
  cid <- sample(phylota@cids, 1)
  expect_true(inherits(phylota[[cid]], 'ClstrRec'))
})