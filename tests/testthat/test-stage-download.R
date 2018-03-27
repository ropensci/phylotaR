# LIBS
library(testthat)

# DATA
ps <- parameters()

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}

# RUNNING
context('Testing \'stage-download\'')
cleanUp()
test_that('cldIdntfy() works', {
  
})
test_that('dwnldSqRcrds() works', {
  # dir.create('cache')
  # txids <- sample(phylt_nds[['ti']], 10)
  # ps[['wd']] <- '.'
  # res <- with_mock(
  #   `phylotaR::nSqs`=mckNSqs,
  #   `phylotaR::dwnldFrmNCBI`=mckDwnldFrmNCBI,
  #   dwnld(txids=txids, phylt_nds=phylt_nds,
  #         ps=ps)
  # )
  # expect_null(res)
  # expect_true(file.exists(file.path('cache', 'sqs',
  #                                   paste0(txids[1], '.RData'))))
  # cleanUp()
})
