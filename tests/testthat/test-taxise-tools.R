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
txdmpfl <- file.path(data_dr, 'taxonomy', 'taxdump.tar.gz')

# FUNCTIONS
cleanUp <- function() {
  if(file.exists(cch_dr)) {
    unlink(cch_dr, recursive=TRUE)
  }
  if(file.exists('test_writeTax.tsv')) {
    file.remove('test_writeTax.tsv')
  }
  if(file.exists('log.txt')) {
    file.remove('log.txt')
  }
  if(file.exists(file.path(data_dr, 'log.txt'))) {
    file.remove(file.path(data_dr, 'log.txt'))
  }
  if(file.exists(txdmpfl)) {
    file.remove(txdmpfl)
  }
}
mckSafeSrch <- function(func, fnm, args, ps) {
  list()
}
mckCurlFail <- function(url, path) {
  list('content'='non-existent-file',
       'status_code'='FAILED')
}
mckCurl <- function(url, path) {
  list('content'=file.path(data_dr, 'taxonomy',
                           'taxdump.tar.gz'),
       'status_code'=226)
}
mckUntar <- function(tarfile, files, exdir) {
  NULL
}
mckXmlToList <- function(rcrd) {
  list('Taxon'=list('Lineage'='A;B;C;', 'LineageEx'=list()))
}

# DATA
data('tdobj')
td_nds <- tdobj[['nds']]
td_nms <- tdobj[['nms']]
ps <- list('mxd'=10000, 'tmout'=10, 'v'=FALSE,
           'wd'=data_dr, 'mxsql'=2000, 'mdlt'=2000)
rm(tdobj)

# RUNNING
context('Testing \'taxise\'')
cleanUp()
test_that('dwnldTD() works', {
  file.create(txdmpfl)
  # pass download
  with_mock(
    `phylotaR:::.untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    dwnldTD(ps=list('wd'=data_dr, 'tdpth'=NULL,
                    'v'=FALSE))
  )
  # pass reading from file
  with_mock(
    `phylotaR:::.untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    dwnldTD(ps=list('wd'=data_dr, 'tdpth'=txdmpfl,
                    'v'=FALSE))
  )
  # fail to read from file
  with_mock(
    `phylotaR:::.untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    expect_error(dwnldTD(ps=list('wd'=data_dr,
                                 'tdpth'='non-existent-file',
                                 'v'=FALSE)))
  )
  # fail to unpack
  with_mock(
    `phylotaR:::.untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    expect_warning(
      expect_error(
        dwnldTD(ps=list('wd'='non-existent-dir',
                        'tdpth'=NULL,
                        'v'=FALSE))))
  )
  # fail to download
  with_mock(
    `phylotaR:::.untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurlFail,
    expect_error(dwnldTD(ps=list('wd'=data_dr,
                                 'tdpth'=NULL,
                                 'v'=FALSE)))
  )
  cleanUp()
})
test_that('genTDObj() works', {
  tdobj <- genTDObj(ps=ps)
  tdobj <- genTDObj(ps=ps)
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
                     td_nds=td_nds, ps=ps)
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
                              ps){
      500},
    getStats(txid=rid, phylt_nds=phylt_nds,
             td_nds=td_nds, td_nms=td_nms,
             ps=ps, rcrsv=TRUE)
    )
  expect_true(nrow(res) > 1)
})
# test_that('writeTax() works', {
#   phylt_nds <- data.frame('ti'=sample(td_nds[['id']], 1),
#                           'ti_anc'=NA,
#                           'rank'=NA,
#                           'n_gi_node'=NA,
#                           'n_gi_sub_nonmodel'=NA,
#                           'n_gi_sub_model'=NA,
#                           'n_sp_desc'=NA,
#                           'n_sp_model'=NA,
#                           'n_leaf_desc'=NA,
#                           'n_otu_desc'=NA,
#                           'ti_genus'=NA,
#                           'n_genera'=NA)
#   writeTax(phylt_nds=phylt_nds, td_nms=td_nms,
#            fl='test_writeTax.tsv', ps=ps)
#   expect_true(file.exists('test_writeTax.tsv'))
#   cleanUp()
# })
test_that('genPhylotaNds() works', {
  # ps[['v']] <- TRUE
  setUpCch(ps=ps)
  txids <- td_nds[['id']]
  nid_sets <- getMngblIds(txid=txids[1:10], td_nds=td_nds,
                          ps=ps)
  res <- with_mock(
    `phylotaR::nSqs`=function(txid,
                              direct,
                              ps){
      500},
    genPhylotaNds(nid_sets=nid_sets,
                  td_nds=td_nds,
                  td_nms=td_nms,
                  ps=ps)
  )
  expect_true('data.frame' %in% is(res))
  cleanUp()
})
test_that('genTxdct() works', {
  setUpCch(ps=ps)
  phylt_nds <- data.frame(ti=rep(NA, 10))
  res <- with_mock(
    `phylotaR::safeSrch`=mckSafeSrch,
    `XML::xmlToList`=mckXmlToList,
    genTxdct(phylt_nds=phylt_nds, ps=ps)
  )
  expect_true(length(res) == nrow(phylt_nds))
  cleanUp()
})
cleanUp()