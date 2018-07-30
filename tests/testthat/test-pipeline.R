# LIBS
library(testthat)

# VARS
wd <- tempdir()

# RUNNING
context('Testing \'pipeline\'')
phylotaR:::cleanup(wd)
test_that('setup() works', {
  res <- with_mock(
    `phylotaR::cmdln` = phylotaR:::cmdln_blastcheck,
    phylotaR::setup(wd = wd, txid = 9606)
  )
  expect_true(file.exists(file.path(wd, 'cache', 'prmtrs.RData')))
  phylotaR:::cleanup(wd)
})
test_that('run() works', {
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    phylotaR::run(wd = wd, nstages = 4)
  )
  expect_null(res)
})
test_that('restart() works', {
  phylotaR:::cache_setup(ps = list('wd' = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = 'taxise')
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    phylotaR::restart(wd = wd, nstages = 4)
  )
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    expect_error(phylotaR::restart(wd = wd, nstages = 1))
  )
  phylotaR:::progress_save(wd = wd, stg = 'download')
  phylotaR:::progress_save(wd = wd, stg = 'cluster')
  phylotaR:::progress_save(wd = wd, stg = 'cluster2')
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    expect_error(phylotaR::restart(wd = wd, nstages = 4))
  )
  phylotaR:::cleanup(wd)
})
test_that('reset(hard=FALSE) works', {
  phylotaR:::cache_setup(ps = parameters(wd = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = 'taxise')
  phylotaR:::progress_save(wd = wd, stg = 'download')
  phylotaR:::progress_save(wd = wd, stg = 'cluster')
  phylotaR:::reset(wd = wd, stage = 'download')
  expect_true(phylotaR:::progress_read(wd = wd) == 'download')
  phylotaR:::cleanup(wd)
})
test_that('reset(hard=TRUE) works', {
  phylotaR:::cache_setup(ps = parameters(wd = wd))
  phylotaR:::progress_init(wd = wd)
  phylotaR:::progress_save(wd = wd, stg = 'taxise')
  phylotaR:::progress_save(wd = wd, stg = 'download')
  phylotaR:::obj_save(wd = wd, obj = NULL, nm = 'txdct')
  phylotaR:::sqs_save(wd = wd, txid = '1', sqs = NULL)
  phylotaR:::reset(wd = wd, stage = 'taxise', hard = TRUE)
  expect_true(phylotaR:::progress_read(wd = wd) == 'taxise')
  expect_null(phylotaR:::obj_load(wd = wd, nm = 'txdct'))
  expect_error(phylotaR:::sqs_load(wd = wd, txid = '1'))
  phylotaR:::cleanup(wd)
})
test_that('parameters_reset() works', {
  res <- with_mock(
    `phylotaR:::cmdln` = phylotaR:::cmdln_blastcheck,
    phylotaR:::setup(wd = wd, txid = 9606)
  )
  phylotaR::parameters_reset(wd = wd, parameters = 'txid', values = 0000)
  expect_true(phylotaR:::parameters_load(wd = wd)[['txid']] == 0000)
  phylotaR:::cleanup(wd)
})
phylotaR:::cleanup(wd)