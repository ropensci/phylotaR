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
# # DATA
# # ps <- parameters()
# # # mock clstrs 
# # clstr <- list("gis"='101', "seed_gi"=NA, "ti_root"='101',
# #               "ci"=NA, "cl_type"=NA, "n_gi"=1,
# #               "tis"='101', "n_ti"=1, "MinLength"=0,
# #               "MaxLength"=0, "n_gen"=1, "n_child"=1,
# #               "ci_anc"=NA, "unique_id"=NA)
# # all_clstrs <- lapply(1:100, function(x) {
# #   clstr[['seed_gi']] <- as.character(x)
# #   clstr[['ci']] <- as.character(x+100)
# #   clstr[['unique_id']] <- as.character(x)
# #   clstr
# #   })
# # seed_ids <- as.character(1:100)
# # sqs <- vector('list', length=200)
# # names(sqs) <- as.character(1:length(sqs))
# # sqs <- sample(sqs)
# # # mock blast results from seed blast
# # blst_rs <- data.frame(query.id=NA, subject.id=NA,
# #                       qcovs=runif(min=1, max=100, n=100))
# # blst_rs[['query.id']] <- as.character(sample(1:100, 50,
# #                                              replace=FALSE))
# # blst_rs[['subject.id']] <- as.character(sample(1:100, 50,
# #                                                replace=FALSE))
# # pull <- blst_rs[['query.id']] == blst_rs[['subject.id']]
# # blst_rs[pull, 'qcovs'] <- 100
# 
# # RUNNING
# cleanUp()
# context('Testing \'cluster^2 tools\'')
# test_that('clstrClstrs() works', {
#   # setUpCch(ps)
#   # # skip cluster^2
#   # res <- with_mock(
#   #   `phylotaR:::getSeedSqs`=function(...) NA,
#   #   `phylotaR:::blstSeeds`=function(...) NA,
#   #   `phylotaR:::jnClstrs`=function(...) NA,
#   #   `phylotaR:::mrgClstrs`=function(...) NA,
#   #   `phylotaR:::rnmbrClstrs`=function(...) NA,
#   #   `phylotaR:::svObj`=function(...) NULL,
#   #   clstrClstrs(ps=ps)
#   # )
#   # expect_null(res)
#   # saveRDS(NA, file=file.path('cache', 'clstrs',
#   #                            'id1.RData'))
#   # saveRDS(NA, file=file.path('cache', 'clstrs',
#   #                            'id2.RData'))
#   # saveRDS(NA, file=file.path('cache', 'sqs',
#   #                            'id1.RData'))
#   # saveRDS(NA, file=file.path('cache', 'sqs',
#   #                            'id2.RData'))
#   # # don't skip cluster^2
#   # res <- with_mock(
#   #   `phylotaR:::getSeedSqs`=function(...) NA,
#   #   `phylotaR:::blstSeeds`=function(...) NA,
#   #   `phylotaR:::jnClstrs`=function(...) NA,
#   #   `phylotaR:::mrgClstrs`=function(...) NA,
#   #   `phylotaR:::rnmbrClstrs`=function(...) NA,
#   #   `phylotaR:::svObj`=function(...) NULL,
#   #   clstrClstrs(ps=ps)
#   # )
#   # expect_null(res)
# })
# cleanUp()

#   with_mock(
#     `phylotaR:::setUpNcbiTools`=mckSetupNcbiTools,
#     `phylotaR:::ldObj`=mckFun,
#     `phylotaR:::clstrClstrs`=mckFun,
#     phylotaR::setUp(wd='.', txid=9606),
#     phylotaR::runClusters2(wd='.')
#   )
#   lglns <- readLines('log.txt')
#   expect_true(grepl('Completed stage',
#                     lglns[length(lglns) - 1]))
# })
# cleanUp()