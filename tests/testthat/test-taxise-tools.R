# TODO:
# create pseudo-nodes.dmp for testing purposes
# this would provide better testing
# e.g. identify all mammals in nodes.dmp

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

# DATA
data('tdobj')
td_nds <- tdobj[['nds']]
td_nms <- tdobj[['nms']]
rm(tdobj)

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
  pull <- td_nds[ ,'id'] %in% td_nds[ ,'parent']
  rid <- sample(td_nds[pull, 'id'], 1)
  kids <- getKids(rid, td_nds)
  expect_true(length(kids) >= 1)
})
test_that('nDscndnts() works', {
  pull <- td_nds[ ,'id'] %in% td_nds[ ,'parent']
  rid <- sample(td_nds[pull, 'id'], 1)
  kids <- getKids(rid, td_nds)
  n <- nDscndnts(rid, td_nds)
  expect_true(n >= length(kids))
  # expect no kids
  rid <- sample(td_nds[!pull, 'id'], 1)
  n <- nDscndnts(rid, td_nds)
  expect_true(n == 0)
})
test_that('getGenus() works', {
  rid <- sample(td_nds[['id']][-1], 1)
  res <- getGenus(txid=rid, td_nds=td_nds)
  if(!is.na(res)) {
    expect_true(CHNOSZ::getrank(res, nodes=td_nds) == 'genus')
  }
})
test_that('getMngblIds() works', {
  pull <- td_nds[ ,'id'] %in% td_nds[ ,'parent']
  wdsndnts <- sample(td_nds[pull, 'id'], 10)
  wodsndnts <- sample(td_nds[!pull, 'id'], 10)
  res <- getMngblIds(txid=c(wdsndnts, wodsndnts),
                     td_nds=td_nds,
                     mx_dscndnts=10000,
                     tmout=10, verbose=FALSE)
  ns <- res[['ndscndnts']]
  expect_equal(ns > 0, rep(c(TRUE, FALSE), each=10))
})
test_that('getStats() works', {
  phylt_nds <- data.frame('ti'=NA,
                          'ti_anc'=NA,
                          'rank'=NA,
                          'n_gi_node'=NA,
                          'n_gi_sub_nonmodel'=NA,
                          'n_gi_sub_model'=NA,
                          'n_sp_desc'=NA,
                          'n_sp_model'=NA,
                          'n_leaf_desc'=NA,
                          'n_otu_desc'=NA,
                          'ti_genus'=NA,
                          'n_genera'=NA)
  rid <- sample(td_nds[['id']][-1], 1)
  res <- with_mock(
    `phylotaR::nSqs`=function(txid,
                              direct,
                              mx_len,
                              verbose){
      500},
    getStats(txid=rid, phylt_nds=phylt_nds,
             mx_sq_lngth=2000, mdl_thrshld=2000,
             td_nds=td_nds, td_nms=td_nms,
             verbose=TRUE, recursive=TRUE)
    )
  expect_true(nrow(res) > 1)
})
cleanUp()