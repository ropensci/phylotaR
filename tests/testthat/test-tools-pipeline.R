# LIBS
library(testthat)

# STUBS ----
cmdln_mock <- function(cmd, args, lgfl = NULL) {
  if (grepl('makeblastdb', cmd)) {
    out <- c("makeblastdb: 2.7.1+\nPackage: blast 2.7.1, build Oct 18 2017 19:57:24")
  }
  if (grepl('blastn', cmd)) {
    out <- c("blastn: 2.7.1+\nPackage: blast 2.7.1, build Oct 18 2017 19:57:24")
  }
  if (grepl('wrngvrsn', cmd)) {
    out <- c("blastn: 1.6.1+\nPackage: blast 1.6.1, build Oct 18 2017 19:57:24")
  }
  list(status = 0, stdout = charToRaw(out), stderr = charToRaw(''))
}

# RUNNING
context('Testing \'setup-tools\'')
phylotaR:::cleanup()
test_that('blast_setup() works', {
  # test with fake system
  res <- with_mock(
    `phylotaR:::cmdln` = cmdln_mock,
    phylotaR:::blast_setup(d  =  '.', v  =  FALSE, wd  =  NULL)
  )
  expect_true(length(res) == 2)
  # make sure wrong versions are flagged
  res <- with_mock(
    `phylotaR:::cmdln`  =  cmdln_mock,
    expect_error(phylotaR:::blast_setup(d  =  'wrngvrsn', v  =  FALSE,
                                        wd  =  NULL))
  )
  # make sure wrong dirs are flagged
  expect_error(phylotaR:::blast_setup(d  =  '.', v  =  FALSE, wd  =  NULL))
})
test_that('parameters_setup() works', {
  expect_error(phylotaR:::parameters_setup(wd = '.', txid = 9606,
                                           ncbi_execs = c('', '')))
  ncbi_execs <- list('mkblstdb' = NA, 'blstn' = NA)
  phylotaR:::parameters_setup(wd = '.', txid = 9606, ncbi_execs = ncbi_execs)
  ps <- phylotaR:::parameters_load(wd = '.')
  expect_true(length(ps) == 17)
})
phylotaR:::cleanup()
test_that('stage_args_check() works', {
  expect_error(phylotaR:::stage_args_check(frm = -1, to = -1))
  expect_error(phylotaR:::stage_args_check(frm = 5, to = 5))
  expect_error(phylotaR:::stage_args_check(frm = 2, to = 1))
  res <- phylotaR:::stage_args_check(frm = 1, to = 4)
  expect_true(res == 'Running stages: taxise, download, cluster, cluster2')
  res <- phylotaR:::stage_args_check(frm = 4, to = 4)
  expect_true(res == 'Running stages: cluster2')
})
test_that('stages_run() works', {
  res <- with_mock(
    `phylotaR::cmdln` = cmdln_mock,
    `phylotaR:::taxise_run` = function(...){NULL},
    `phylotaR:::download_run` = function(...){NULL},
    `phylotaR:::clusters_run` = function(...){NULL},
    `phylotaR:::clusters2_run` = function(...){NULL},
    phylotaR::setup(wd = '.', txid = 9606),
    phylotaR:::stages_run(wd = '.', to = 4, frm = 1, stgs_msg = '',
                          rstrt = FALSE)
  )
  expect_null(res)
  res <- with_mock(
    `phylotaR:::taxise_run` = function(...){NULL},
    `phylotaR:::download_run` = function(...){NULL},
    `phylotaR:::clusters_run` = function(...){NULL},
    `phylotaR:::clusters2_run` = function(...){NULL},
    phylotaR:::stages_run(wd = '.', to = 4, frm = 1, stgs_msg = '',
                          rstrt = TRUE)
  )
  expect_null(res)
  phylotaR:::cleanup()
})
phylotaR:::cleanup()