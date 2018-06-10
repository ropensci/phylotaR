# # LIBS
# library(testthat)
# 
# # DATA
# ps <- parameters()
# 
# # FUNCTIONS
# cleanUp <- function() {
#   if(file.exists('cache')) {
#     unlink('cache', recursive=TRUE)
#   }
# }
# # stub
# mckNSqs <- function(txid, ps, direct=FALSE) {
#   # randomly return either > or < mdl_thrshld
#   sample(c(10, ps[['mdlthrs']] + 1), 1)
# }
# 
# # RUNNING
# context('Testing \'stage-tools-download\'')
# cleanUp()
# test_that('hrrchcDwnld() works', {
#   
# })
# test_that('prtDwnldSqRcrds() works', {
#   
# })
# test_that('btchDwnldSqRcrds() works', {
#   
# })
# test_that('agmntDwnldSqRcrds() works', {
#   
# })
# cleanUp()