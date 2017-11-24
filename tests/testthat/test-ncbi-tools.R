# LIBS
library(phylotaR)
library(testthat)

# RUNNING
context('Testing \'ncbi-tools\'')
test_that('safeSrch() works', {
  args <- list('this and that')
  myfunc <- function(...) {
    print(...)
    return(1)
  }
  res <- safeSrch(func=myfunc,
                  args=args,
                  fnm='print()',
                  verbose=TRUE,
                  mx_retry=2)
  expect_true(res == 1)
  myfunc <- function(...) {
    print(...)
    stop()
  }
  res <- safeSrch(func=myfunc,
                  args=args,
                  fnm='print()',
                  verbose=TRUE,
                  mx_retry=2)
  expect_null(res)
})
test_that('nSqs', {
  res <- with_mock(
    `phylotaR::safeSrch`=function(func,
                                  args,
                                  fnm,
                                  verbose,
                                  mx_retry){
      list('count'=args)
    },
    nSqs(txid=9606, direct=FALSE)
  )
  expect_true(grepl(':exp', res[['term']]))
  res <- with_mock(
    `phylotaR::safeSrch`=function(func,
                                  args,
                                  fnm,
                                  verbose,
                                  mx_retry){
      list('count'=args)
    },
    nSqs(txid=9606, direct=TRUE)
  )
  expect_true(grepl(':noexp', res[['term']]))
})