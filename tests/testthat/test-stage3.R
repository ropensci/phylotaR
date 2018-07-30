# LIBS
library(testthat)
library(phylotaR)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
sqs <- readRDS(phylotaR:::datadir_get('sqrecs.rda'))
sqs <- phylotaR:::seqarc_gen(sqs)
exclstrarc <- readRDS(phylotaR:::datadir_get('clstrarc.rda'))

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'stage3\'')
test_that('clusters_run() works', {
  phylotaR:::cache_setup(ps)
  with_mock(
    `phylotaR:::obj_load` = function(...) NULL,
    `phylotaR:::clstrs_calc` = function(...) NULL,
    clusters_run(wd = wd)
  )
  lglns <- readLines(file.path(wd, 'log.txt'))
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
})
phylotaR:::cleanup(wd)
test_that('clstrs_calc() works', {
  phylotaR:::cache_setup(ps)
  saveRDS(object = sqs, file = file.path(ps[['wd']], 'cache', 'sqs', '1.RData'))
  with_mock(
    `phylotaR:::clstr_all` = function(...) exclstrarc,
    phylotaR:::clstrs_calc(txdct = NULL, ps = ps)
  )
  expect_true(list.files(file.path(wd, 'cache', 'clstrs')) == '1.RData')
})
phylotaR:::cleanup(wd)

# # Examples clusters archive
# flpth <- file.path('demos', 'aotus', 'cache', 'clstrs', '9504.RData')
# clstrrecs <- readRDS(file = flpth)
# clstrrecs <- sample(x = clstrrecs@clstrs, size = 10)
# clstrarc <- phylotaR:::clstrarc_gen(clstrrecs = clstrrecs)
# saveRDS(object = clstrarc, file = phylotaR:::datadir_get('clstrarc.rda'))
