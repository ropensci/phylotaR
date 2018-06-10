# # LIBS
# library(testthat)
# 
# # DATA
# ps <- parameters()
# #9607
# 
# # FUNCTIONS
# randFasta <- function() {
#   seq <- sample(c('A', 'T', 'C', 'G'), size=1000,
#                 replace=TRUE)
#   seq <- paste(seq, collapse='')
#   paste0('> rand0m.name.--jg82-5969\n', seq, '\n')
# }
# randSummary <- function() {
#   list(uid=NA, txid=NA, caption=NA,
#        accessionversion=NA,
#        slen=NA, createdate=NA, title=NA)
# }
# cleanUp <- function() {
#   if(file.exists('cache')) {
#     unlink('cache', recursive=TRUE)
#   }
# }
# # stubs
# mckGetGIs <- function(...) {
#   if(n == 0) {
#     ids <- NULL
#   } else {
#     ids <- as.character(1:n)
#   }
#   ids
# }
# mckEntrezSearch <- function(db, term, retmax, retstart, ...) {
#   if(n == 0) {
#     ids <- NULL
#   } else {
#     ids <- as.character(1:n)
#   }
#   list('ids'=ids, 'count'=n)
# }
# mckEntrezFetch <- function(db, rettype, id) {
#   seqs <- ''
#   for(i in 1:length(id)) {
#     seqs <- paste0(seqs, randFasta())
#   }
#   seqs
# }
# mckEntrezFetchIDs <- function(...) {
#   '11111\n11111\n11111\n11111\n11111'
# }
# mckEntrezSummary <- function(db, id) {
#   if(n == 1) {
#     return(randSummary())
#   }
#   summ <- vector('list', length=n)
#   for(i in 1:length(id)) {
#     summ[[i]] <- randSummary()
#   }
#   summ
# }
# 
# # RUNNING
# cleanUp()
# context('Testing \'entrez-tools\'')
# test_that('nSqs() works', {
#   setUpCch(ps=ps)
#   res <- with_mock(
#     `phylotaR:::safeSrch`=function(func,
#                                   args,
#                                   fnm,
#                                   ps){
#       list('count'=args)
#     },
#     phylotaR:::nSqs(txid=9606, drct=FALSE, ps=ps)
#   )
#   expect_true(grepl(':exp', res[['term']]))
#   res <- with_mock(
#     `phylotaR:::safeSrch`=function(func,
#                                   args,
#                                   fnm,
#                                   ps){
#       list('count'=args)
#     },
#     phylotaR:::nSqs(txid=9606, drct=TRUE, ps=ps)
#   )
#   expect_true(grepl(':noexp', res[['term']]))
# })
# cleanUp()
# test_that('nNds() works', {
#   setUpCch(ps=ps)
#   res <- with_mock(
#     `phylotaR:::safeSrch`=function(func,
#                                   args,
#                                   fnm,
#                                   ps){
#       list('count'=args)
#     },
#     phylotaR:::nNds(txid=9606, ps=ps)
#   )
#   expect_true(res[['term']] == "txid9606[Subtree]")
# })
# cleanUp()
# test_that('getGIs() works', {
#   cleanUp()
#   setUpCch(ps=ps)
#   n <<- 100
#   res <- with_mock(
#     `rentrez::entrez_search`=mckEntrezSearch,
#     `rentrez::entrez_fetch`=mckEntrezFetchIDs,
#     phylotaR:::getGIs(txid=9606, drct=FALSE,
#                      sqcnt=100, ps=ps)
#   )
#   expect_true(length(res) == n)
#   cleanUp()
#   setUpCch(ps=ps)
#   n <<- 100
#   res <- with_mock(
#     `rentrez::entrez_search`=mckEntrezSearch,
#     `rentrez::entrez_fetch`=mckEntrezFetchIDs,
#     phylotaR:::getGIs(txid=9606, drct=FALSE,
#                      sqcnt=100, ps=ps, hrdmx=20,
#                      retmax=10)
#   )
#   expect_true(length(res) != n)
# })
# cleanUp()