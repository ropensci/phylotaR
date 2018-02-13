# LIBS
library(phylotaR)
library(testthat)

# DATA
data("phylt_nds")
ps <- parameters()

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}
# stub
mckNSqs <- function(txid, ps, direct=FALSE) {
  # randomly return either > or < mdl_thrshld
  sample(c(10, ps[['mdlthrs']] + 1), 1)
}
mckDwnldFrmNCBI <- function(txid=txid, ps, direct=FALSE) {
  nsqs <- runif(n=1, min=1, max=10)
  sqs <- list()
  for(i in 1:nsqs) {
    sq <- list("gi"=NA, "ti"=NA, "acc"=NA,
               "acc_vers"=NA, "length"=NA,
               "division"=NA, "acc_date"=NA,
               "gbrel"=NA, "def"=NA, "seq"=NA)
    sqs <- c(sqs, list(sq))
  }
  sqs
}

# RUNNING
context('Testing \'download-tools\'')
cleanUp()
test_that('getSqsByTxid() works', {
  txids <- fltr(txid=9479, phylt_nds=phylt_nds,
                ps=ps)
  txid <- sample(txids, 1)
  res <- with_mock(
    `phylotaR::nSqs`=mckNSqs,
    `phylotaR::dwnldFrmNCBI`=mckDwnldFrmNCBI,
    getSqsByTxid(txid=txid, phylt_nds=phylt_nds,
                 ps=ps)
  )
  expect_true(class(res) == 'list')
})
test_that('dwnld() works', {
  dir.create('cache')
  txids <- sample(phylt_nds[['ti']], 10)
  ps[['wd']] <- '.'
  res <- with_mock(
    `phylotaR::nSqs`=mckNSqs,
    `phylotaR::dwnldFrmNCBI`=mckDwnldFrmNCBI,
    dwnld(txids=txids, phylt_nds=phylt_nds,
          ps=ps)
  )
  expect_null(res)
  expect_true(file.exists(file.path('cache', 'sqs',
                                    paste0(txids[1], '.RData'))))
  cleanUp()
})
test_that('fltr() works', {
  # TODO: not sure what the goal of filtr is, hard to test
  txid <- 9479
  ps[['mdlthrs']] <- 3000
  ps[['mxsqs']] <- 10000
  res_1 <- fltr(txid=txid, phylt_nds=phylt_nds,
                ps=ps)
  ps[['mdlthrs']] <- 30
  ps[['mxsqs']] <- 100
  res_2 <- fltr(txid=txid, phylt_nds=phylt_nds,
                ps=ps)
  expect_true(length(res_1) < length(res_2))
})
test_that('getDDFrmPhyltNds() works', {
  rnd <- sample(1:nrow(phylt_nds), 1)
  txid <- phylt_nds[['ti_anc']][rnd]
  expctd <- phylt_nds[['ti']][rnd]
  res <- getDDFrmPhyltNds(txid=txid, phylt_nds=phylt_nds)
  expect_true(expctd %in% res)
})
cleanUp()
