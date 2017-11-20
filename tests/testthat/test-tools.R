# LIBS
library(phylotaR)
library(testthat)

# RUNNING
# pretty lame tests.... any better ideas would be appreciated
context('Testing \'tools\'')
test_that('setUpNcbiTools() works', {
  # can't really test if BLAST tools are installed
  # instead just test whether the tools fail
  expect_error(setUpNcbiTools(dr='.'))
})
test_that('mkPrmtrs() works', {
  # simple simple test
  prmtrs <- mkPrmtrs(ncbi_execs=c('', ''))
  expect_true(length(prmtrs) == 8)
})