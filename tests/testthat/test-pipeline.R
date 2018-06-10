# LIBS
library(testthat)

# RUNNING
context('Testing \'pipeline\'')
phylotaR:::cleanup()
test_that('setup() works', {
  res <- with_mock(
    `phylotaR::cmdln` = phylotaR:::cmdln_blastcheck,
    phylotaR::setup(wd = '.', txid = 9606)
  )
  expect_true(file.exists(file.path('cache', 'prmtrs.RData')))
  phylotaR:::cleanup()
})
test_that('run() works', {
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    phylotaR::run(wd = '.', nstages = 4)
  )
  expect_null(res)
})
test_that('restart() works', {
  phylotaR:::cache_setup(ps = list('wd' = '.'))
  phylotaR:::progress_init(wd = '.')
  phylotaR:::progress_save(wd = '.', stg = 'taxise')
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    phylotaR::restart(wd = '.', nstages = 4)
  )
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    expect_error(phylotaR::restart(wd = '.', nstages = 1))
  )
  phylotaR:::progress_save(wd = '.', stg = 'download')
  phylotaR:::progress_save(wd = '.', stg = 'cluster')
  phylotaR:::progress_save(wd = '.', stg = 'cluster2')
  res <- with_mock(
    `phylotaR:::stages_run` = function(...) NULL,
    expect_error(phylotaR::restart(wd = '.', nstages = 4))
  )
  phylotaR:::cleanup()
})
test_that('reset() works', {
  # TODO
})
test_that('parameters_reset() works', {
  res <- with_mock(
    `phylotaR:::cmdln` = phylotaR:::cmdln_blastcheck,
    phylotaR:::setup(wd = '.', txid = 9606)
  )
  phylotaR::parameters_reset(wd = '.', parameters = 'txid', values = 0000)
  expect_true(phylotaR:::parameters_load(wd = '.')[['txid']] == 0000)
  phylotaR:::cleanup()
})
phylotaR:::cleanup()