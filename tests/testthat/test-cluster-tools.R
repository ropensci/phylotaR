# LIBS
library(phylotaR)
library(testthat)

# VARS
data("phylt_nds")
data('sqs')
data('blst_rs')
ps <- list('wd'='.', 'v'=FALSE)
clstr <- list("gis"=NA, "seed_gi"=NA, "ti_root"=NA, "ci"=NA,
              "cl_type"=NA, "n_gi"=1, "tis"=NA, "n_ti"=NA,
              "MinLength"=NA, "MaxLength"=NA, "n_gen"=NA,
              "n_child"=NA, "ci_anc"=NA, "unique_id"=NA)
clstrs <- list(clstr, clstr, clstr)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}
# stubs
mckBlstN <- function(dbfl, outfl, ps) {
  blst_rs
}
mckMkBlstDB <- function(sqs, dbfl, ps) {
  NULL
}

# RUNNING
cleanUp()
context('Testing \'cluster-tools\'')
test_that('blstSqs() works', {
  setUpCch(ps=ps)
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    blstSqs(txid=1, typ='direct', sqs=sqs, ps=ps, lvl=0)
  )
  expect_true('data.frame' %in% is(res))
})
cleanUp()
test_that('calcClstrs() works', {
  setUpCch(ps=ps)
  res <- with_mock(
    `phylotaR::clstrAll`=function(...) list(clstr),
    `list.files`=function(...) paste0(1:5, '.RData'),
    `readRDS`=function(...) sqs,
    calcClstrs(txid=0, phylt_nds=phylt_nds, ps=ps)
  )
  res <- list.files(file.path(ps[['wd']], 'cache',
                              'clstrs'))
  expect_true(all(res == paste0(1:5, '.RData')))
})
cleanUp()
test_that('clstrSqs() works', {
  setUpCch(ps=ps)
  txid <- sample(phylt_nds[['ti']], 1)
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    clstrSqs(txid=txid, sqs=sqs, phylt_nds=phylt_nds,
             infrmtv=FALSE, ps=ps, lvl=0)
  )
  expect_true('list' %in% is(res))
  res <- with_mock(
    `phylotaR::blstN`=mckBlstN,
    `phylotaR::mkBlstDB`=mckMkBlstDB,
    clstrSqs(txid=9479, sqs=sqs, phylt_nds=phylt_nds,
             infrmtv=FALSE, ps=ps, lvl=0)
  )
  # using the platyrrhini txid we should get more res
  expect_true(length(res) > 0)
})
cleanUp()
# no cache tests
test_that('clstrAll() works', {
  rndn <- round(runif(min=5, max=10, n=1))
  res <- with_mock(
    `phylotaR::getDDFrmPhyltNds`=function(txid, ...) {
      if(txid == 0) return(rep(1, rndn - 1))
      NULL
      },
    `phylotaR::clstrSbtr`=function(...) list(clstr),
    clstrAll(txid=0, sqs=sqs, phylt_nds=phylt_nds,
             ps=ps, lvl=0)
  )
  expect_true(length(res) == rndn)
})
test_that('clstrSbtr() works', {
  # the parent should represent all sequences
  res <- with_mock(
    `phylotaR::clstrSqs`=function(sqs, ...) {
      sapply(sqs, function(x) x[['gi']])
    },
    clstrSbtr(txid=9479, sqs=sqs,
              phylt_nds=phylt_nds,
              dds=NULL, ps=ps, lvl=0)
  )
  expect_true(length(res) == length(sqs))
  # a random should repr. equal or fewer seqs.
  txid <- sample(phylt_nds[['ti']], 1)
  res <- with_mock(
    `phylotaR::clstrSqs`=function(sqs, ...) {
      sapply(sqs, function(x) x[['gi']])
    },
    clstrSbtr(txid=txid, sqs=sqs,
              phylt_nds=phylt_nds,
              dds=NULL, ps=ps, lvl=0)
  )
  expect_true(length(res) <= length(sqs))
})
test_that('clstrDrct() works', {
  # expect empty list because 0 is not in phylt_nds
  res <- with_mock(
    `phylotaR::clstrSqs`=function(...) NA,
    clstrDrct(txid=0, sqs=sqs,
              phylt_nds=phylt_nds, ps=ps, lvl=0)
  )
  expect_true(length(res) == 0)
  # choose a txid from sqs
  # should run clstrSqs -- returning TRUE
  rnd <- sample(1:length(sqs), 1)
  txid <- sqs[[rnd]][['ti']]
  res <- with_mock(
    `phylotaR::clstrSqs`=function(...) TRUE,
    clstrDrct(txid=txid, sqs=sqs,
              phylt_nds=phylt_nds,
              ps=ps, lvl=0)
  )
  expect_true(res)
})
test_that('clstrPhylt() works', {
  phylt_nds <- clstrPhylt(clstrs=clstrs)
  expect_true(nrow(phylt_nds) == 3)
})
test_that('clstrCiGi() works', {
  cigi <- clstrCiGi(clstrs=clstrs)
  expect_true(nrow(cigi) == 3)
})
test_that('clstrBlstRs() works', {
  res <- clstrBlstRs(blst_rs)
  expect_true(length(res) == 1)
  nrws <- blst_rs[1:2, ]
  nrws[['query.id']] <- 1:2
  nrws[['subject.id']] <- 2:1
  res <- clstrBlstRs(rbind(blst_rs, nrws))
  expect_true(length(res) == 2)
})
test_that('addClstrInf() works', {
  clstrs <- clstrBlstRs(blst_rs)
  res <- addClstrInf(clstrs=clstrs, phylt_nds=phylt_nds,
              txid=9479, sqs=sqs, typ='direct')
  expect_true(length(clstrs) == length(res))
})
test_that('getADs() works', {
  # all nodes should descend from Platyrrhini
  txids <- getADs(txid=9479, phylt_nds=phylt_nds)
  expect_true(length(txids) == (nrow(phylt_nds) - 1))
})
test_that('getGnsFrmPhyltNds() works', {
  rnd <- sample(phylt_nds[['ti_anc']], 1)
  res <- getGnsFrmPhyltNds(txid=rnd, phylt_nds=phylt_nds)
  expect_true(res %in% phylt_nds[['ti_genus']])
})