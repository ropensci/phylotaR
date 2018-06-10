# LIBS
library(testthat)

# RUNNING
context('Testing \'cache-tools\'')
phylotaR:::cleanup()
test_that('progress_init() works', {
  phylotaR:::cache_setup(ps  =  list('wd' = '.'))
  phylotaR:::progress_init(wd = '.')
  expect_true(file.exists(file.path('cache', 'progress.RData')))
})
phylotaR:::cleanup()
test_that('progress_read() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  phylotaR:::progress_init(wd = '.')
  stg <- phylotaR:::progress_read(wd = '.')
  expect_true(stg == 'taxise')
})
phylotaR:::cleanup()
test_that('progress_save() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  phylotaR:::progress_init(wd = '.')
  phylotaR:::progress_save(wd = '.', stg = 'taxise')
  expect_true(phylotaR:::progress_read(wd = '.') == 'download')
  phylotaR:::progress_save(wd = '.', stg = 'download')
  expect_true(phylotaR:::progress_read(wd = '.') == 'cluster')
  phylotaR:::progress_save(wd = '.', stg = 'cluster')
  expect_true(phylotaR:::progress_read(wd = '.') == 'cluster2')
  phylotaR:::progress_save(wd = '.', stg = 'cluster2')
  expect_true(is.na(phylotaR:::progress_read(wd = '.')))
})
phylotaR:::cleanup()
test_that('progress_reset() works', {
  # create prgrss
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  phylotaR:::progress_init(wd = '.')
  phylotaR:::progress_save(wd = '.', stg = 'taxise')
  phylotaR:::progress_save(wd = '.', stg = 'download')
  phylotaR:::progress_save(wd = '.', stg = 'cluster')
  phylotaR:::progress_save(wd = '.', stg = 'cluster2')
  expect_true(is.na(phylotaR:::progress_read(wd = '.')))
  # reset
  phylotaR:::progress_reset(wd = '.', stg = 'download')
  expect_true(phylotaR:::progress_read(wd = '.') == 'download')
})
phylotaR:::cleanup()
test_that('cache_setup() works', {
  ps <- list('alovelyparameter' = 1, 'wd' = '.')
  phylotaR:::cache_setup(ps = ps)
  expect_true(file.exists(file.path('cache',
                                    'prmtrs.RData')))
  expect_error(phylotaR:::cache_setup(ps = ps))
})
phylotaR:::cleanup()
test_that('parameters_load() works', {
  phylotaR:::cache_setup(ps = list('anotherlovelyparameter' = 1, 'wd' = '.'))
  ps <- phylotaR:::parameters_load(wd = '.')
  expect_true('anotherlovelyparameter' %in% names(ps))
})
phylotaR:::cleanup()
test_that('cache_rm() works', {
  phylotaR:::cache_setup(ps = list('andfinally' = 1, 'wd' = '.'))
  phylotaR:::cache_rm(wd = '.')
  expect_false(file.exists(file.path('cache', 'prmtrs.RData')))
})
phylotaR:::cleanup()
test_that('obj_check() works', {
  expect_false(phylotaR:::obj_check(wd = '.', nm = 'mylovelyobj'))
  phylotaR:::obj_save(wd = '.', nm = 'mylovelyobj',
                      obj = list('interestingthings'))
  expect_true(phylotaR:::obj_check(wd = '.', nm = 'mylovelyobj'))
})
phylotaR:::cleanup()
test_that('obj_save() works', {
  phylotaR:::obj_save(wd = '.', nm = 'mylovelyobj',
                      obj = list('interestingthings'))
  expect_true(file.exists(file.path('cache', 'mylovelyobj.RData')))
})
phylotaR:::cleanup()
test_that('obj_load() works', {
  expect_error(phylotaR:::obj_load(wd = '.', nm = 'mylovelyobj'))
  phylotaR:::obj_save(wd = '.', nm = 'mylovelyobj',
                      obj = list('interestingthings'))
  obj <- phylotaR:::obj_load(wd = '.', nm = 'mylovelyobj')
  expect_true(obj[[1]] == 'interestingthings')
})
phylotaR:::cleanup()
test_that('ncbicache_save() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  args <- list('term' = 'unique_search_term', 'db' = 'nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  exptd <- file.path('cache', 'ncbi', 'search', 'nucleotide', '1.RData')
  phylotaR:::ncbicache_save(fnm = fnm, args = args, wd = '.', obj = obj)
  expect_true(file.exists(exptd))
})
phylotaR:::cleanup()
test_that('ncbicache_load() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  args <- list('term' = 'unique_search_term', 'db' = 'nucleotide')
  fnm <- 'search'
  obj <- 'important_result'
  phylotaR:::ncbicache_save(fnm = fnm, args = args, wd = '.', obj = obj)
  res <- phylotaR:::ncbicache_load(fnm = fnm, args = args, wd = '.')
  expect_equal(obj, res)
})
phylotaR:::cleanup()
test_that('blastcache_save() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  sids_1 <- c('1', '2', '3')
  sids_2 <- c('1', '2', '4')
  phylotaR:::blastcache_save(sids = sids_1, wd = '.', obj = 'res_1')
  phylotaR:::blastcache_save(sids = sids_2, wd = '.', obj = 'res_2')
  expect_true(file.exists(file.path('cache', 'blast', '1.RData')))
  expect_true(file.exists(file.path('cache', 'blast', '2.RData')))
  expect_false(file.exists(file.path('cache', 'blast', '3.RData')))
})
phylotaR:::cleanup()
test_that('blastcache_load() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  blst_rs <- data.frame('query.id' = as.character(c(1:5)), 'subject.id' = '1')
  sids_1 <- c('1', '2', '3')
  sids_2 <- c('1', '2', '4')
  sids_3 <- c('1', '2', '5')
  sids_4 <- c('1', '2')
  phylotaR:::blastcache_save(sids = sids_1, wd = '.', obj = blst_rs[1:3, ])
  phylotaR:::blastcache_save(sids = sids_2, wd = '.',
                             obj = blst_rs[c(1, 2, 4), ])
  expect_true(nrow(phylotaR:::blastcache_load(sids = sids_1, wd = '.')) == 3)
  expect_true(nrow(phylotaR:::blastcache_load(sids = sids_2, wd = '.')) == 3)
  expect_null(phylotaR:::blastcache_load(sids = sids_3, wd = '.'))
  # sqs 4 blast results are contained within sqs 1 and/or 2
  expect_true(nrow(phylotaR:::blastcache_load(sids = sids_4, wd = '.')) == 2)
})
phylotaR:::cleanup()