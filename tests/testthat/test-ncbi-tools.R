# LIBS
library(phylotaR)
library(testthat)

# FUNCTIONS
randFasta <- function() {
  seq <- sample(c('A', 'T', 'C', 'G'), size=1000,
                replace=TRUE)
  seq <- paste(seq, collapse='')
  paste0('> rand0m.name.--jg82-5969\n', seq, '\n')
}
randSummary <- function() {
  list(uid=NA, txid=NA, caption=NA,
       accessionversion=NA,
       slen=NA, createdate=NA, title=NA)
}
# stubs
mckEntrezSearch <- function(db, term,
                            use_history,
                            retmax) {
  if(n == 0) {
    ids <- NULL
  } else {
    ids <- as.character(1:n)
  }
  list('ids'=ids)
}
mckEntrezFetch <- function(db, rettype,
                           web_history,
                           retmax, retstart) {
  mx <- ifelse(retmax<n, retmax, n)
  seqs <- ''
  for(i in (retstart+1):(mx+retstart)) {
    seqs <- paste0(seqs, randFasta())
  }
  seqs
}
mckEntrezSummary <- function(db, web_history,
                             retmax, retstart) {
  if(n == 1) {
    return(randSummary())
  }
  mx <- ifelse(retmax<n, retmax, n)
  summ <- vector('list', length=n)
  for(i in (retstart+1):(mx+retstart)) {
    summ[[i]] <- randSummary()
  }
  summ
}

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
test_that('dwnldSqs', {
  # n determines the number of available seqs.
  n <- 0
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldSqs(txid=1, direct=FALSE, mx_lngth=25000,
             mx_sqs=100000, verbose=FALSE)
  )
  expect_true(class(res) == 'list')
  expect_true(length(res) == 0)
  n <<- 1
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldSqs(txid=1, direct=FALSE, mx_lngth=25000,
             mx_sqs=100000, verbose=FALSE)
  )
  expect_true(length(res) == 1)
  n <<- 100
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldSqs(txid=1, direct=FALSE, mx_lngth=25000,
             mx_sqs=100000, verbose=FALSE)
  )
  expect_true(length(res) == 100)
  n <<- 1000
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldSqs(txid=1, direct=FALSE, mx_lngth=25000,
             mx_sqs=100000, verbose=FALSE)
  )
  expect_true(length(res) == 1000)
})