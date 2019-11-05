# LIBS
library(testthat)

# DATA
# error examples
failed_search <- 'Error in entrez_check(...'
load(file.path(phylotaR:::datadir_get(subdir  =  'api'), 'xml_timeout.rda'))
# real example: search tax
flpth <- file.path(phylotaR:::datadir_get(subdir  =  'api'),
                   'entrez_search_tax.rda')
# trm <- 'txid10151[Subtree]'
# entrez_search_tax <- rentrez::entrez_search(db = 'taxonomy', retmax = 10,
#                                             term = trm)
# save(entrez_search_tax, file = flpth)
load(flpth)
# real example: search sids
flpth <- file.path(phylotaR:::datadir_get(subdir  =  'api'),
                   'entrez_search_sid.rda')
# trm <- searchterm_gen(txid = '10151', ps = parameters(), direct = FALSE)
# entrez_search_sid <- rentrez::entrez_search(db = 'nucleotide', retmax = 0,
#                                             term = trm)
# save(entrez_search_sid, file = flpth)
load(flpth)
# real example: fetch sids
flpth <- file.path(phylotaR:::datadir_get(subdir  =  'api'),
                   'entrez_fetch_sid.rda')
# trm <- phylotaR:::searchterm_gen(txid = '10151', ps = parameters(),
#                                  direct = FALSE)
# srch <- rentrez::entrez_search(db = 'nucleotide', retmax = 0, term = trm,
#                                use_history = TRUE)
# entrez_fetch_sid <- rentrez::entrez_fetch(db = 'nucleotide',
#                                           web_history = srch[['web_history']],
#                                           retmax = 1000L, rettype = 'acc',
#                                           retstart = 0L)
# save(entrez_fetch_sid, file = flpth)
load(flpth)
# real example: fetch seqs
flpth <- file.path(phylotaR:::datadir_get(subdir  =  'api'),
                   'entrez_fetch_seqs.rda')
# ids <- sample(strsplit(entrez_fetch_sid, '\n')[[1]], 100)
# entrez_fetch_seqs <- rentrez::entrez_fetch(db = 'nucleotide', id = ids,
#                                            rettype = 'gb', retmode = 'text')
# save(entrez_fetch_seqs, file = flpth)
load(flpth)
wd <- tempdir()
ps <- parameters(wd = wd)
ps[['mxrtry']] <- 2
ps[['wt_tms']] <- c(0, 0)

# RUNNING
phylotaR:::cleanup(wd)
context('Testing \'tools-api\'')
test_that('download_obj_check() works', {
  expect_true(phylotaR:::download_obj_check(entrez_search_tax))
  expect_true(phylotaR:::download_obj_check(entrez_search_sid))
  expect_true(phylotaR:::download_obj_check(entrez_fetch_sid))
  expect_true(phylotaR:::download_obj_check(entrez_fetch_seqs))
  expect_false(phylotaR:::download_obj_check(xml_tmout))
  expect_false(phylotaR:::download_obj_check(failed_search))
})
test_that('safely_connect(fnm = search) works', {
  phylotaR:::cache_setup(ps  =  ps)
  args <- list('term'  =  'ncbi search term', 'db'  =  'nucleotide', 'x')
  myfunc <- function(...) {
    print(...)
    return('NCBI record')
  }
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'search',
                                   ps = ps)
  expect_true(res == 'NCBI record')
  myfunc <- function(...) {
    print(...)
    stop()
  }
  args <- list('term' = 'another ncbi search term', 'db' = 'nucleotide', 'x')
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'search',
                                   ps = ps)
  expect_null(res)
})
phylotaR:::cleanup(wd)
test_that('safely_connect(fnm = fetch) works', {
  phylotaR:::cache_setup(ps = ps)
  args <- list('id' = c(1, 2), 'db' = 'nucleotide', 'x')
  myfunc <- function(...) {
    print(...)
    return('NCBI record')
  }
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'fetch',
                                   ps = ps)
  expect_true(res  ==  'NCBI record')
  myfunc <- function(...) {
    print(...)
    stop()
  }
  args <- list('id' = c(2, 3), 'db' = 'nucleotide', 'x')
  res <- phylotaR:::safely_connect(func = myfunc, args = args, fnm = 'fetch',
                                   ps = ps)
  expect_null(res)
})
phylotaR:::cleanup(wd)
test_that('batcher() works', {
  phylotaR:::cache_setup(ps = ps)
  res <- phylotaR:::batcher(ids = 1:10, func = function(ids, ps) {
    return(ids)
  }, ps, lvl = 1)
  expect_true(all(res == 1:10))
  res <- phylotaR:::batcher(ids = 1:1000, func = function(ids, ps) {
    return(ids)
  }, ps, lvl = 1)
  expect_true(all(res == 1:1000))
})
phylotaR:::cleanup(wd)