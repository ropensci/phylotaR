# LIBS
library(testthat)

# VARS
wd <- tempdir()

# RUNNING
context('Testing \'cache-tools\'')
phylotaR:::cleanup(wd)
test_that('progress_init() works', {
  phylotaR:::cache_setup(ps  =  list('wd' = wd))
  phylotaR:::progress_init(wd = wd)
  expect_true(file.exists(file.path(wd, 'cache', 'progress.RData')))
})
phylotaR:::cleanup(wd)
test_that('progress_read() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  phylotaR:::progress_init(wd = wd)
  stg <- phylotaR:::progress_read(wd = wd)
  expect_true(stg == 'taxise')
})
phylotaR:::cleanup(wd)
test_that('progress_save() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = 'taxise')
  expect_true(phylotaR:::progress_read(wd = wd) == 'download')
  phylotaR:::progress_save(wd = wd, stg = 'download')
  expect_true(phylotaR:::progress_read(wd = wd) == 'cluster')
  phylotaR:::progress_save(wd = wd, stg = 'cluster')
  expect_true(phylotaR:::progress_read(wd = wd) == 'cluster2')
  phylotaR:::progress_save(wd = wd, stg = 'cluster2')
  expect_true(is.na(phylotaR:::progress_read(wd = wd)))
})
phylotaR:::cleanup(wd)
test_that('progress_reset() works', {
  # create prgrss
  phylotaR:::cache_setup(ps = list('wd' = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = 'taxise')
  phylotaR:::progress_save(wd = wd, stg = 'download')
  phylotaR:::progress_save(wd = wd, stg = 'cluster')
  phylotaR:::progress_save(wd = wd, stg = 'cluster2')
  expect_true(is.na(phylotaR:::progress_read(wd = wd)))
  # reset
  phylotaR:::progress_reset(wd = wd, stg = 'download')
  expect_true(phylotaR:::progress_read(wd = wd) == 'download')
})
phylotaR:::cleanup(wd)
test_that('cache_setup() works', {
  ps <- list('alovelyparameter' = 1, 'wd' = wd)
  phylotaR:::cache_setup(ps = ps)
  expect_true(file.exists(file.path(wd, 'cache', 'prmtrs.RData')))
  expect_error(phylotaR:::cache_setup(ps = ps))
})
phylotaR:::cleanup(wd)
test_that('parameters_load() works', {
  phylotaR:::cache_setup(ps = list('anotherlovelyparameter' = 1, 'wd' = wd))
  ps <- phylotaR:::parameters_load(wd = wd)
  expect_true('anotherlovelyparameter' %in% names(ps))
})
phylotaR:::cleanup(wd)
test_that('cache_rm() works', {
  phylotaR:::cache_setup(ps = list('andfinally' = 1, 'wd' = wd))
  phylotaR:::cache_rm(wd = wd)
  expect_false(file.exists(file.path(wd, 'cache', 'prmtrs.RData')))
})
phylotaR:::cleanup(wd)
test_that('obj_check() works', {
  expect_false(phylotaR:::obj_check(wd = wd, nm = 'mylovelyobj'))
  phylotaR:::obj_save(wd = wd, nm = 'mylovelyobj',
                      obj = list('interestingthings'))
  expect_true(phylotaR:::obj_check(wd = wd, nm = 'mylovelyobj'))
})
phylotaR:::cleanup(wd)
test_that('obj_save() works', {
  phylotaR:::obj_save(wd = wd, nm = 'mylovelyobj',
                      obj = list('interestingthings'))
  expect_true(file.exists(file.path(wd, 'cache', 'mylovelyobj.RData')))
})
phylotaR:::cleanup(wd)
test_that('obj_load() works', {
  expect_error(phylotaR:::obj_load(wd = wd, nm = 'mylovelyobj'))
  phylotaR:::obj_save(wd = wd, nm = 'mylovelyobj',
                      obj = list('interestingthings'))
  obj <- phylotaR:::obj_load(wd = wd, nm = 'mylovelyobj')
  expect_true(obj[[1]] == 'interestingthings')
})
phylotaR:::cleanup(wd)
test_that('ncbicache_save() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  args <- list('term' = 'unique_search_term', 'db' = 'nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  exptd <- file.path(wd, 'cache', 'ncbi', 'search', 'nucleotide', '1.RData')
  phylotaR:::ncbicache_save(fnm = fnm, args = args, wd = wd, obj = obj)
  expect_true(file.exists(exptd))
})
phylotaR:::cleanup(wd)
test_that('ncbicache_load() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  args <- list('term' = 'unique_search_term', 'db' = 'nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  phylotaR:::ncbicache_save(fnm = fnm, args = args, wd = wd, obj = obj)
  res <- phylotaR:::ncbicache_load(fnm = fnm, args = args, wd = wd)
  expect_equal(obj, res)
})
phylotaR:::cleanup(wd)
test_that('blastcache_save() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  sids_1 <- c('1', '2', '3')
  sids_2 <- c('1', '2', '4')
  phylotaR:::blastcache_save(sids = sids_1, wd = wd, obj = 'res_1')
  phylotaR:::blastcache_save(sids = sids_2, wd = wd, obj = 'res_2')
  expect_true(file.exists(file.path(wd, 'cache', 'blast', '1.RData')))
  expect_true(file.exists(file.path(wd, 'cache', 'blast', '2.RData')))
  expect_false(file.exists(file.path(wd, 'cache', 'blast', '3.RData')))
})
phylotaR:::cleanup(wd)
test_that('blastcache_load() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  blst_rs <- data.frame('query.id' = as.character(c(1:5)), 'subject.id' = '1')
  sids_1 <- c('1', '2', '3')
  sids_2 <- c('1', '2', '4')
  sids_3 <- c('1', '2', '5')
  sids_4 <- c('1', '2')
  phylotaR:::blastcache_save(sids = sids_1, wd = wd, obj = blst_rs[1:3, ])
  phylotaR:::blastcache_save(sids = sids_2, wd = wd,
                             obj = blst_rs[c(1, 2, 4), ])
  expect_true(nrow(phylotaR:::blastcache_load(sids = sids_1, wd = wd)) == 3)
  expect_true(nrow(phylotaR:::blastcache_load(sids = sids_2, wd = wd)) == 3)
  expect_null(phylotaR:::blastcache_load(sids = sids_3, wd = wd))
  # sqs 4 blast results are contained within sqs 1 and/or 2
  expect_true(nrow(phylotaR:::blastcache_load(sids = sids_4, wd = wd)) == 2)
})
phylotaR:::cleanup(wd)