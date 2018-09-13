# LIBS
library(testthat)
library(phylotaR)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
sqs <- readRDS(phylotaR:::datadir_get('sqrecs.rda'))
sqs <- phylotaR:::seqarc_gen(sqs)
clstrs <- readRDS(phylotaR:::datadir_get('clstrarc.rda'))
sids <- unique(unlist(lapply(clstrs@clstrs, function(x) x@sids)))
# modify sqs
sqs <- sample(sqs@sqs, length(sids))
for (i in seq_along(sqs)) {
  sqs[[i]]@id <- sids[[i]]
}
sqs <- phylotaR:::seqarc_gen(sqs)

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'cluster^2 tools\'')
test_that('clusters2_run() works', {
  phylotaR:::cache_setup(ps)
  with_mock(
    `phylotaR:::clstr2_calc` = function(...) NULL,
    phylotaR::clusters2_run(wd = ps[['wd']])
  )
  lglns <- readLines(file.path(wd, 'log.txt'))
  expect_true(grepl('Completed stage', lglns[length(lglns) - 1]))
})
phylotaR:::cleanup(wd)
test_that('clstr2_calc() works', {
  phylotaR:::cache_setup(ps)
  # skip cluster^2
  res <- with_mock(
    `phylotaR:::seeds_blast` = function(...) NA,
    `phylotaR:::clstrs_join` = function(...) NA,
    `phylotaR:::clstrs_merge` = function(...) NA,
    `phylotaR:::clstrs_renumber` = function(...) NA,
    `phylotaR:::obj_save` = function(...) NULL,
    phylotaR:::clstr2_calc(ps = ps)
  )
  expect_null(res)
  saveRDS(clstrs, file = file.path(wd, 'cache', 'clstrs', 'id1.RData'))
  saveRDS(clstrs, file = file.path(wd, 'cache', 'clstrs', 'id2.RData'))
  saveRDS(sqs, file = file.path(wd, 'cache', 'sqs', 'id1.RData'))
  saveRDS(sqs, file = file.path(wd, 'cache', 'sqs', 'id2.RData'))
  # don't skip cluster^2
  res <- with_mock(
    `phylotaR:::seeds_blast` = function(...) NA,
    `phylotaR:::clstrs_join` = function(...) NA,
    `phylotaR:::clstrs_merge` = function(...) clstrs@clstrs,
    `phylotaR:::clstrs_renumber` = function(...) clstrs,
    `phylotaR:::obj_load` = function(...) NULL,
    phylotaR:::clstr2_calc(ps = ps)
  )
  resfl <- file.path(ps[['wd']], 'cache', 'clstrs_sqs.RData')
  expect_true(file.exists(resfl))
})
phylotaR:::cleanup(wd)
