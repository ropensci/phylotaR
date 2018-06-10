# # LIBS
# library(testthat)
# 
# # FUNCTIONS
# cleanUp <- function() {
#   if(file.exists('cache')) {
#     unlink('cache', recursive=TRUE)
#   }
# }
# 
# # STUBS
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
# mckRunStgs <- function(wd, to, frm, stgs_msg, rstrt=FALSE) {
#   NULL
# }
# mckRunTaxise <- function(wd) {
#   NULL
# }
# mckRunDownload <- function(wd) {
#   NULL
# }
# mckRunClusters <- function(wd) {
#   NULL
# }
# mckRunClusters2 <- function(wd) {
#   NULL
# }
# 
# # RUNNING
# context('Testing \'pipeline\'')
# cleanUp()
# test_that('setUp() works', {
#   res <- with_mock(
#     `phylotaR::cmdLn`=mckCmdLn,
#     phylotaR::setUp(wd='.', txid=9606)
#   )
#   expect_true(file.exists(file.path('cache',
#                                     'prmtrs.RData')))
#   cleanUp()
# })
# test_that('run() works', {
#   res <- with_mock(
#     `phylotaR:::runStgs`=mckRunStgs,
#     phylotaR::run(wd='.', nstages=4)
#   )
#   expect_null(res)
# })
# test_that('runStgs() works', {
#   res <- with_mock(
#     `phylotaR::cmdLn`=mckCmdLn,
#     `phylotaR:::runTaxise`=mckRunTaxise,
#     `phylotaR:::runDownload`=mckRunDownload,
#     `phylotaR:::runClusters`=mckRunClusters,
#     `phylotaR:::runClusters`=mckRunClusters2,
#     phylotaR::setUp(wd='.', txid=9606),
#     phylotaR::runStgs(wd='.', to=4, frm=1,
#                       stgs_msg='', rstrt=FALSE)
#   )
#   expect_null(res)
#   res <- with_mock(
#     `phylotaR:::runTaxise`=mckRunTaxise,
#     `phylotaR:::runDownload`=mckRunDownload,
#     `phylotaR:::runClusters`=mckRunClusters,
#     `phylotaR:::runClusters`=mckRunClusters2,
#     phylotaR::runStgs(wd='.', to=4, frm=1,
#                       stgs_msg='', rstrt=TRUE)
#   )
#   expect_null(res)
#   cleanUp()
# })
# test_that('restart() works', {
#   phylotaR:::setUpCch(ps=list('wd'='.'))
#   phylotaR:::intPrgrss(wd='.')
#   phylotaR:::svPrgrss(wd='.', stg='taxise')
#   res <- with_mock(
#     `phylotaR:::runStgs`=mckRunStgs,
#     phylotaR::restart(wd='.', nstages=4)
#   )
#   res <- with_mock(
#     `phylotaR:::runStgs`=mckRunStgs,
#     expect_error(phylotaR::restart(wd='.',
#                                    nstages=1))
#   )
#   phylotaR:::svPrgrss(wd='.', stg='download')
#   phylotaR:::svPrgrss(wd='.', stg='cluster')
#   phylotaR:::svPrgrss(wd='.', stg='align')
#   res <- with_mock(
#     `phylotaR:::runStgs`=mckRunStgs,
#     expect_error(phylotaR::restart(wd='.',
#                                    nstages=4))
#   )
#   cleanUp()
# })
# test_that('chckStgs() works', {
#   expect_error(phylotaR:::chckStgs(frm=-1, to=-1))
#   expect_error(phylotaR:::chckStgs(frm=5, to=5))
#   expect_error(phylotaR:::chckStgs(frm=2, to=1))
#   res <- phylotaR:::chckStgs(frm=1, to=4)
#   expect_true(res == 'Running stages: taxise, download, cluster, cluster2')
#   res <- phylotaR:::chckStgs(frm=4, to=4)
#   expect_true(res == 'Running stages: cluster2')
# })
# test_that('reset() works', {
#   # TODO
# })
# test_that('setParameters() works', {
#   res <- with_mock(
#     `phylotaR::cmdLn`=mckCmdLn,
#     phylotaR::setUp(wd='.', txid=9606)
#   )
#   phylotaR::setParameters(wd='.', parameters='txid',
#                           values=0000)
#   expect_true(phylotaR:::ldPrmtrs(wd='.')[['txid']] == 0000)
#   cleanUp()
# })
# cleanUp()