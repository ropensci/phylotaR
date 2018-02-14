# LIBS
library(phylotaR)
library(testthat)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}

# STUBS
mckSystem <- function(command, args, fl=NULL) {
  if(grepl('makeblastdb', command)) {
    res <- c("makeblastdb: 2.7.1+",
             " Package: blast 2.7.1, build Oct 18 2017 19:57:24")
  }
  if(grepl('blastn', command)) {
    res <- c("blastn: 2.7.1+", 
             " Package: blast 2.7.1, build Oct 18 2017 19:57:24")
  }
  res
}
mckRunStgs <- function(wd, to, frm, stgs_msg, rstrt=FALSE) {
  NULL
}
mckRunTaxise <- function(wd) {
  NULL
}
mckRunDownload <- function(wd) {
  NULL
}
mckRunClusters <- function(wd) {
  NULL
}
mckRunClusters2 <- function(wd) {
  NULL
}

# RUNNING
context('Testing \'pipeline\'')
cleanUp()
test_that('setUp() works', {
  res <- with_mock(
    `phylotaR:::.system`=mckSystem,
    setUp(wd='.', txid=9606)
  )
  expect_true(file.exists(file.path('cache',
                                    'prmtrs.RData')))
  cleanUp()
})
test_that('run() works', {
  res <- with_mock(
    `phylotaR:::runStgs`=mckRunStgs,
    run(wd='.', nstages=4)
  )
  expect_null(res)
})
test_that('runStgs() works', {
  res <- with_mock(
    `phylotaR:::.system`=mckSystem,
    `phylotaR:::runTaxise`=mckRunTaxise,
    `phylotaR:::runDownload`=mckRunDownload,
    `phylotaR:::runClusters`=mckRunClusters,
    `phylotaR:::runClusters`=mckRunClusters2,
    setUp(wd='.', txid=9606),
    runStgs(wd='.', to=4, frm=1, stgs_msg='',
            rstrt=FALSE)
  )
  expect_null(res)
  res <- with_mock(
    `phylotaR:::runTaxise`=mckRunTaxise,
    `phylotaR:::runDownload`=mckRunDownload,
    `phylotaR:::runClusters`=mckRunClusters,
    `phylotaR:::runClusters`=mckRunClusters2,
    runStgs(wd='.', to=4, frm=1, stgs_msg='',
            rstrt=TRUE)
  )
  expect_null(res)
  cleanUp()
})
test_that('restart() works', {
  setUpCch(ps=list('wd'='.'))
  intPrgrss(wd='.')
  svPrgrss(wd='.', stg='taxise')
  res <- with_mock(
    `phylotaR:::runStgs`=mckRunStgs,
    restart(wd='.', nstages=4)
  )
  res <- with_mock(
    `phylotaR:::runStgs`=mckRunStgs,
    expect_error(restart(wd='.', nstages=1))
  )
  svPrgrss(wd='.', stg='download')
  svPrgrss(wd='.', stg='cluster')
  svPrgrss(wd='.', stg='align')
  res <- with_mock(
    `phylotaR:::runStgs`=mckRunStgs,
    expect_error(restart(wd='.', nstages=4))
  )
  cleanUp()
})
test_that('chckStgs() works', {
  expect_error(chckStgs(frm=-1, to=-1))
  expect_error(chckStgs(frm=5, to=5))
  expect_error(chckStgs(frm=2, to=1))
  res <- chckStgs(frm=1, to=4)
  expect_true(res == 'Running stages: taxise, download, cluster, cluster2')
  res <- chckStgs(frm=4, to=4)
  expect_true(res == 'Running stages: cluster2')
})
test_that('reset() works', {
  # TODO
})
test_that('setParameters() works', {
  res <- with_mock(
    `phylotaR:::.system`=mckSystem,
    setUp(wd='.', txid=9606)
  )
  setParameters(wd='.', parameters='txid', values=0000)
  expect_true(ldPrmtrs(wd='.')[['txid']] == 0000)
  cleanUp()
})
cleanUp()