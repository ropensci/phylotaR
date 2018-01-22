# LIBS
library(phylotaR)
library(testthat)

# DATA
ps <- list(wd='.', txid=9607,
           tdpth=NULL, mxd=10000,
           tmout=100, mdlt=3000,
           mxsqs=10000, mxsql=25000,
           mxretry=2, v=FALSE, ncps=1)

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
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}
# stubs
mckEntrezSearch <- function(db, term, retmax, retstart) {
  if(n == 0) {
    ids <- NULL
  } else {
    ids <- as.character(1:n)
  }
  list('ids'=ids, 'count'=n)
}
mckEntrezFetch <- function(db, rettype, id) {
  seqs <- ''
  for(i in 1:length(id)) {
    seqs <- paste0(seqs, randFasta())
  }
  seqs
}
mckEntrezSummary <- function(db, id) {
  if(n == 1) {
    return(randSummary())
  }
  summ <- vector('list', length=n)
  for(i in 1:length(id)) {
    summ[[i]] <- randSummary()
  }
  summ
}

# RUNNING
cleanUp()
context('Testing \'ncbi-tools\'')
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
test_that('nSqs() works', {
  setUpCch(ps=ps)
  res <- with_mock(
    `phylotaR::safeSrch`=function(func,
                                  args,
                                  fnm,
                                  ps){
      list('count'=args)
    },
    nSqs(txid=9606, direct=FALSE, ps=ps)
  )
  expect_true(grepl(':exp', res[['term']]))
  res <- with_mock(
    `phylotaR::safeSrch`=function(func,
                                  args,
                                  fnm,
                                  ps){
      list('count'=args)
    },
    nSqs(txid=9606, direct=TRUE, ps=ps)
  )
  expect_true(grepl(':noexp', res[['term']]))
})
cleanUp()
test_that('getGIs() works', {
  cleanUp()
  setUpCch(ps=ps)
  n <<- 100
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    getGIs(txid=9606, direct=FALSE, sqcnt=100, ps=ps)
  )
  expect_true(length(res) == n)
  cleanUp()
  setUpCch(ps=ps)
  n <<- 100
  # TODO: develop mock entrez search to improve
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    getGIs(txid=9606, direct=FALSE, sqcnt=100, ps=ps,
           hrdmx=20)
  )
  expect_true(length(res) != n)
})
test_that('dwnldFrmNCBI() works', {
  # n determines the number of available seqs.
  cleanUp()
  n <<- 0
  setUpCch(ps=ps)
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldFrmNCBI(txid=1, direct=FALSE, ps=ps)
  )
  expect_true(class(res) == 'list')
  expect_true(length(res) == 0)
  cleanUp()
  setUpCch(ps=ps)
  n <<- 1
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldFrmNCBI(txid=1, direct=FALSE, ps=ps)
  )
  expect_true(length(res) == 1)
  cleanUp()
  setUpCch(ps=ps)
  n <<- 100
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldFrmNCBI(txid=1, direct=FALSE, ps=ps)
  )
  expect_true(length(res) == 100)
  cleanUp()
  setUpCch(ps=ps)
  n <<- 1000
  res <- with_mock(
    `rentrez::entrez_search`=mckEntrezSearch,
    `rentrez::entrez_fetch`=mckEntrezFetch,
    `rentrez::entrez_summary`=mckEntrezSummary,
    dwnldFrmNCBI(txid=1, direct=FALSE, ps=ps)
  )
  expect_true(length(res) == 1000)
  cleanUp()
})
cleanUp()