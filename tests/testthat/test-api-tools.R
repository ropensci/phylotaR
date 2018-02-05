# LIBS
library(phylotaR)
library(testthat)

# DATA
wd <- getwd()
if(grepl('testthat', wd)) {
  data_d <- file.path('data', 'api')
} else {
  # for running test at package level
  data_d <- file.path('tests', 'testthat',
                      'data', 'api')
}
load(file.path(data_d, 'xml_timeout.rda'))
ps <- list(wd='.', txid=9607,
           tdpth=NULL, mxd=10000,
           tmout=100, mdlt=3000,
           mxsqs=10000, mxsql=25000,
           mxretry=2, v=FALSE, ncps=1,
           wt_tms=c(0, 0))

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}

# RUNNING
cleanUp()
context('Testing \'api-tools\'')
test_that('chckSrchObj() works', {
  expect_false(chckSrchObj(xml_tmout))
})
test_that('safeSrch(fnm=search) works', {
  setUpCch(ps=ps)
  args <- list('term'='ncbi search term',
               'db'='nucleotide', 'x')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- safeSrch(func=myfunc,
                  args=args,
                  fnm='search',
                  ps=ps)
  expect_true(res == 1)
  myfunc <- function(...) {
    print(...)
    stop()
  }
  args <- list('term'='another ncbi search term',
               'db'='nucleotide',
               'x')
  res <- safeSrch(func=myfunc,
                  args=args,
                  fnm='search',
                  ps=ps)
  expect_null(res)
})
cleanUp()
test_that('safeSrch(fnm=fetch) works', {
  setUpCch(ps=ps)
  args <- list('id'=c(1, 2),
               'db'='nucleotide',
               'x')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- safeSrch(func=myfunc,
                  args=args,
                  fnm='fetch',
                  ps=ps)
  expect_true(res == 1)
  myfunc <- function(...) {
    print(...)
    stop()
  }
  args <- list('id'=c(2, 3),
               'db'='nucleotide',
               'x')
  res <- safeSrch(func=myfunc,
                  args=args,
                  fnm='fetch',
                  ps=ps)
  expect_null(res)
})
cleanUp()