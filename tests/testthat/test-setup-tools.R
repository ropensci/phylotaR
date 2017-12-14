# LIBS
library(phylotaR)
library(testthat)

# FUNCTIONS
cleanUp <- function() {
  if(file.exists('cache')) {
    unlink('cache', recursive=TRUE)
  }
  if(file.exists('log.txt')) {
    file.remove('log.txt')
  }
}
# stub
mckSystem <- function(cmd, intern, ignore.stderr) {
  if(grepl('makeblastdb', cmd)) {
    res <- c("makeblastdb: 2.7.1+",
             " Package: blast 2.7.1, build Oct 18 2017 19:57:24")
  }
  if(grepl('blastn', cmd)) {
    res <- c("blastn: 2.7.1+", 
             " Package: blast 2.7.1, build Oct 18 2017 19:57:24")
  }
  if(grepl('wrngvrsn', cmd)) {
    res <- res <- c("blastn: 2.6.1+", 
                    " Package: blast 2.6.1, build Oct 18 2017 19:57:24")
  }
  res
}

# RUNNING
# pretty lame tests.... any better ideas would be appreciated
context('Testing \'setup-tools\'')
cleanUp()
test_that('setUpNcbiTools() works', {
  # test with fake system
  res <- with_mock(
    `base::system`=mckSystem,
    setUpNcbiTools(d='.',
                   verbose=FALSE, wd=NULL)
  )
  expect_true(length(res) == 2)
  # make sure wrong versions are flagged
  res <- with_mock(
    `base::system`=mckSystem,
    expect_error(setUpNcbiTools(d='wrngvrsn',
                                verbose=FALSE,
                                wd=NULL))
  )
  # make sure wrong dirs are flagged
  expect_error(setUpNcbiTools(d='.',
                              verbose=FALSE,
                              wd=NULL))
})
test_that('setUpPrmtrs() works', {
  expect_error(setUpPrmtrs(wd='.', txid=9606,
                           ncbi_execs=c('', '')))
  ncbi_execs <- list('mkblstdb'=NA,
                     'blstn'=NA)
  setUpPrmtrs(wd='.', txid=9606,
              ncbi_execs=ncbi_execs)
  prmtrs <- ldPrmtrs(wd='.')
  expect_true(length(prmtrs) == 12)
})
cleanUp()