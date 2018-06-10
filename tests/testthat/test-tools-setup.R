# # LIBS
# library(testthat)
# 
# # FUNCTIONS
# cleanUp <- function() {
#   if(file.exists('cache')) {
#     unlink('cache', recursive=TRUE)
#   }
#   if(file.exists('log.txt')) {
#     file.remove('log.txt')
#   }
# }
# # stub
# mckCmdLn <- function(cmd, args, lgfl=NULL) {
#   if(grepl('makeblastdb', cmd)) {
#     out <- c("makeblastdb: 2.7.1+\nPackage: blast 2.7.1, build Oct 18 2017 19:57:24")
#   }
#   if(grepl('blastn', cmd)) {
#     out <- c("blastn: 2.7.1+\nPackage: blast 2.7.1, build Oct 18 2017 19:57:24")
#   }
#   if(grepl('wrngvrsn', cmd)) {
#     out <- c("blastn: 1.6.1+\nPackage: blast 1.6.1, build Oct 18 2017 19:57:24")
#   }
#   list(status=0, stdout=charToRaw(out), stderr=charToRaw(''))
# }
# 
# # RUNNING
# context('Testing \'setup-tools\'')
# cleanUp()
# test_that('setUpNcbiTools() works', {
#   # test with fake system
#   res <- with_mock(
#     `phylotaR:::cmdLn`=mckCmdLn,
#     phylotaR:::setUpNcbiTools(d='.', v=FALSE, wd=NULL)
#   )
#   expect_true(length(res) == 2)
#   # make sure wrong versions are flagged
#   res <- with_mock(
#     `phylotaR:::cmdLn`=mckCmdLn,
#     expect_error(phylotaR:::setUpNcbiTools(d='wrngvrsn',
#                                           v=FALSE,
#                                           wd=NULL))
#   )
#   # make sure wrong dirs are flagged
#   expect_error(phylotaR:::setUpNcbiTools(d='.',
#                                         v=FALSE,
#                                         wd=NULL))
# })
# test_that('setUpPrmtrs() works', {
#   expect_error(phylotaR:::setUpPrmtrs(wd='.',
#                                      txid=9606,
#                                      ncbi_execs=c('', '')))
#   ncbi_execs <- list('mkblstdb'=NA,
#                      'blstn'=NA)
#   phylotaR:::setUpPrmtrs(wd='.', txid=9606,
#               ncbi_execs=ncbi_execs)
#   ps <- ldPrmtrs(wd='.')
#   expect_true(length(ps) == 17)
# })
# cleanUp()
