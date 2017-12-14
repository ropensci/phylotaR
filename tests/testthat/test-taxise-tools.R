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
  if(file.exists(txdmpfl)) {
    file.remove(txdmpfl)
  }
}

# DATA
data('tdobj')
td_nds <- tdobj[['nds']]
td_nms <- tdobj[['nms']]
rm(tdobj)

# FUNCTIONS
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

# RUNNING
context('Testing \'taxise\'')
cleanUp()
test_that('dwnldTD() works', {
  file.create(txdmpfl)
  # pass download
  with_mock(
    `utils::untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    dwnldTD(wd=data_dr, tdpth=NULL)
  )
  # pass reading from file
  with_mock(
    `utils::untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    dwnldTD(wd=data_dr, tdpth=txdmpfl)
  )
  # fail to read from file
  with_mock(
    `utils::untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    expect_error(dwnldTD(wd=data_dr,
                         tdpth='non-existent-file'))
  )
  # fail to unpack
  with_mock(
    `utils::untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurl,
    expect_warning(
      expect_error(
        dwnldTD(wd='non-existent-dir',
                tdpth=NULL)))
  )
  # fail to download
  with_mock(
    `utils::untar`=mckUntar,
    `curl::curl_fetch_disk`=mckCurlFail,
    expect_error(dwnldTD(wd=data_dr, tdpth=NULL))
  )
  cleanUp()
})
test_that('genTDObj() works', {
  tdobj <- genTDObj(wd=data_dr)
  tdobj <- genTDObj(wd=data_dr)
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
                     tmout=10, verbose=FALSE,
                     wd=NULL)
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
    getStats(wd=NULL, txid=rid, phylt_nds=phylt_nds,
             mx_sq_lngth=2000, mdl_thrshld=2000,
             td_nds=td_nds, td_nms=td_nms,
             verbose=FALSE, recursive=TRUE)
    )
  expect_true(nrow(res) > 1)
})
test_that('writeTax() works', {
  phylt_nds <- data.frame('ti'=sample(td_nds[['id']], 1),
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
  writeTax(wd=NULL, phylt_nds=phylt_nds, td_nms=td_nms,
           fl='test_writeTax.tsv', verbose=FALSE)
  expect_true(file.exists('test_writeTax.tsv'))
  cleanUp()
})
test_that('genPhylotaNds() works', {
  txids <- td_nds[['id']]
  nid_sets <- getMngblIds(wd=NULL, txid=txids,
                          td_nds=td_nds,
                          mx_dscndnts=100,
                          tmout=10, verbose=FALSE)
  res <- with_mock(
    `phylotaR::nSqs`=function(txid,
                              direct,
                              mx_len,
                              verbose){
      500},
    genPhylotaNds(wd=NULL, nid_sets=nid_sets,
                  mx_sq_lngth=2000,
                  mdl_thrshld=2000,
                  td_nds=td_nds,
                  td_nms=td_nms,
                  verbose=FALSE)
  )
})
cleanUp()