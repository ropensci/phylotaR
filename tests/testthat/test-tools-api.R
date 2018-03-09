# LIBS
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
ps <- parameters()
ps[['mxrtry']] <- 2
ps[['wt_tms']] <- c(0, 0)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}

# RUNNING
cleanUp()
context('Testing \'tools-api\'')
test_that('phylotaR::chckSrchObj() works', {
  expect_false(phylotaR:::chckSrchObj(xml_tmout))
})
test_that('safeSrch(fnm=search) works', {
  phylotaR:::setUpCch(ps=ps)
  args <- list('term'='ncbi search term',
               'db'='nucleotide', 'x')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- phylotaR:::safeSrch(func=myfunc,
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
  res <- phylotaR:::safeSrch(func=myfunc,
                            args=args,
                            fnm='search',
                            ps=ps)
  expect_null(res)
})
cleanUp()
test_that('safeSrch(fnm=fetch) works', {
  phylotaR:::setUpCch(ps=ps)
  args <- list('id'=c(1, 2),
               'db'='nucleotide',
               'x')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- phylotaR:::safeSrch(func=myfunc,
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
  res <- phylotaR:::safeSrch(func=myfunc,
                            args=args,
                            fnm='fetch',
                            ps=ps)
  expect_null(res)
})
cleanUp()