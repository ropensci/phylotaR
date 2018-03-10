# LIBS
library(phylotaR)
library(testthat)

# DATA
# randomly choose one of the example phylota objects
pssbls <- c("aotus", "bromeliads", "cycads", "dragonflies",
            "sturgeons", "tardigrades", "tinamous")
rndm <- sample(pssbls, 1)
list.files('data')
data_env <- new.env()
do.call(what=data, args=list(rndm, envir=data_env))
assign(x='phylota', value=data_env[[rndm]])
rm(data_env)

# RUNNING
context('Testing \'all-class\'')
test_that('SqRcrd() works', {
  sqrcrd <- phylota@sqs@sqs[[1]]
  show(sqrcrd)
  print(sqrcrd)
  str(sqrcrd)
  summary(sqrcrd)
  sqstrng <- as.character(sqrcrd)
  expect_true(validObject(sqrcrd))
})
test_that('SqRcrdBx() works', {
  sqrcrdbx <- phylota@sqs
  show(sqrcrdbx)
  print(sqrcrdbx)
  str(sqrcrdbx)
  summary(sqrcrdbx)
  sqstrng <- as.character(sqrcrdbx)
  expect_true(validObject(sqrcrdbx))
  # [ and [[
  sqrcrd <- sqrcrdbx[[sqrcrdbx@ids[[1]]]]
  expect_true(inherits(sqrcrd, 'SqRcrd'))
  sqrcrdbx2 <- sqrcrdbx[sqrcrdbx@ids[1:10]]
  expect_true(inherits(sqrcrdbx2, 'SqRcrdBx'))
})
test_that('ClRcrd() works', {
  clrcrd <- phylota@cls@cls[[1]]
  show(clrcrd)
  print(clrcrd)
  str(clrcrd)
  summary(clrcrd)
  clstrng <- as.character(clrcrd)
  expect_true(validObject(clrcrd))
  # [ and [[
  expect_true(validObject(clrcrd))
})
test_that('ClRcrdBx() works', {
  clrcrdbx <- phylota@cls
  show(clrcrdbx)
  print(clrcrdbx)
  str(clrcrdbx)
  summary(clrcrdbx)
  clstrng <- as.character(clrcrdbx)
  expect_true(validObject(clrcrdbx))
  # [ and [[
  clrcrd <- clrcrdbx[[clrcrdbx@ids[[1]]]]
  expect_true(inherits(clrcrd, 'ClRcrd'))
  clrcrdbx2 <- clrcrdbx[clrcrdbx@ids[1:10]]
  expect_true(inherits(clrcrdbx2, 'ClRcrdBx'))
})
test_that('TxRcrd() works', {
  txid <- phylota@txids[[1]]
  txrcrd <- phylota@txdct@rcrds[[txid]]
  show(txrcrd)
  print(txrcrd)
  str(txrcrd)
  summary(txrcrd)
  txstrng <- as.character(txrcrd)
  expect_true(validObject(txrcrd))
})
test_that('TxDct() works', {
  txdct <- phylota@txdct
  show(txdct)
  print(txdct)
  str(txdct)
  summary(txdct)
  txstrng <- as.character(txdct)
  expect_true(validObject(txdct))
})
test_that('PhyLoTa() works', {
  show(phylota)
  print(phylota)
  str(phylota)
  cids <- sample(phylota@cids, 2)
  summary(drop_cls(phylota=phylota, cid=cids))
  txstrng <- as.character(phylota)
  expect_true(validObject(phylota))
  # [[
  sid <- sample(phylota@sids, 1)
  expect_true(inherits(phylota[[sid]], 'SqRcrd'))
  cid <- sample(phylota@cids, 1)
  expect_true(inherits(phylota[[cid]], 'ClRcrd'))
})