# LIBS
library(phylotaR)
library(testthat)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
}

# RUNNING
# pretty lame tests.... any better ideas would be appreciated
context('Testing \'tools\'')
cleanUp()
test_that('setUpNcbiTools() works', {
  # can't really test if BLAST tools are installed
  # instead just test whether the tools fail
  expect_error(setUpNcbiTools(dr='.'))
})
test_that('setUpPrmtrs() works', {
  expect_error(setUpPrmtrs(wd='.', ncbi_execs=c('', '')))
  ncbi_execs <- list('mkblstdb'=NA,
                     'blstn'=NA)
  setUpPrmtrs(wd='.', ncbi_execs=ncbi_execs)
  prmtrs <- ldPrmtrs(wd='.')
  expect_true(length(prmtrs) == 7)
})
cleanUp()