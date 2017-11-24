# TODO:
# create pseudo-nodes.dmp for testing purposes
# this would provide better testing

# LIBS
library(phylotaR)
library(testthat)

# VARS
wd <- getwd()
if(grepl('testthat', wd)) {
  data_dr <- 'data'
} else {
  # for running test at package level
  data_dr <- file.path('tests', 'testthat',
                       'data')
}
cch_dr <- file.path(data_dr, 'cache')

# FUNCTIONS
cleanUp <- function() {
  if(file.exists(cch_dr)) {
    unlink(cch_dr, recursive=TRUE)
  }
}

# RUNNING
context('Testing \'taxise\'')
cleanUp()
test_that('genTDObj() works', {
  tdobj <- genTDObj(data_dr)
  tdobj <- genTDObj(data_dr)
  two <- sum(names(tdobj) %in% c('nds', 'nms'))
  expect_true(two == 2)
  cleanUp()
})
test_that('getKids() works', {
  tdobj <- genTDObj(data_dr)
  td_nds <- tdobj[['nds']]
  # nodes.dmp is only a sample
  # can only use id with known parent
  pull <- td_nds[ ,'id'] %in% td_nds[ ,'parent']
  rid <- sample(td_nds[pull, 'id'], 1)
  kids <- getKids(rid, td_nds)
  expect_true(length(kids) >= 1)
  cleanUp()
})
test_that('nDscndnts() works', {
  tdobj <- genTDObj(data_dr)
  td_nds <- tdobj[['nds']]
  # nodes.dmp is only a sample
  # can only use id with known parent
  pull <- td_nds[ ,'id'] %in% td_nds[ ,'parent']
  rid <- sample(td_nds[pull, 'id'], 1)
  kids <- getKids(rid, td_nds)
  n <- nDscndnts(rid, td_nds)
  expect_true(n >= length(kids))
  # expect no kids
  rid <- sample(td_nds[!pull, 'id'], 1)
  n <- nDscndnts(rid, td_nds)
  expect_true(n == 0)
  cleanUp()
})
test_that('getMngblIds() works', {
  tdobj <- genTDObj(data_dr)
  td_nds <- tdobj[['nds']]
  td_nds <- td_nds[-1, ] # remove top, it is SELF-REFERENTIAL!
  pull <- td_nds[ ,'id'] %in% td_nds[ ,'parent']
  wdsndnts <- sample(td_nds[pull, 'id'], 10)
  wodsndnts <- sample(td_nds[!pull, 'id'], 10)
  res <- getMngblIds(txid=c(wdsndnts, wodsndnts),
                     td_nds=td_nds,
                     mx_dscndnts=10000,
                     tmout=10, verbose=FALSE)
  ns <- res[['ndscndnts']]
  expect_equal(ns > 0, rep(c(TRUE, FALSE), each=10))
  cleanUp()
})
cleanUp()