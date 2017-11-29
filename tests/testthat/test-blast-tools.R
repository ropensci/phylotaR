# LIBS
library(phylotaR)
library(testthat)

# VARS
data('sqs')  # TODO: make sure sqs match blast results
wd <- getwd()
if(grepl('testthat', wd)) {
  data_d <- file.path('data')
} else {
  # for running test at package level
  data_d <- file.path('tests', 'testthat',
                       'data')
}

# FUNCTIONS
# stubs
mckLdPrmtrs <- function(wd) {
  list('mkblstdb'='',
       'blstn'='')
}
mckSystem <- function(cmd) {
  0
}

# RUNNING
context('Testing \'blast-tools\'')
test_that('mkBlstDB() works', {
  res <- with_mock(
    `phylotaR::ldPrmtrs`=mckLdPrmtrs,
    `base::system`=mckSystem,
    mkBlstDB(sqs=sqs, dbfl='testdb', wd=data_d,
             verbose=FALSE)
  )
  expect_true(res == file.path(data_d, 'blast',
                               'testdb'))
})
test_that('blstN() works', {
  res <- with_mock(
    `phylotaR::ldPrmtrs`=mckLdPrmtrs,
    `base::system`=mckSystem,
    blstN(dbfl='testdb', outfl='testblstn', wd=data_d)
  )
  nms <- colnames(res)
  exp_nms <- c('query.id', 'subject.id', 'identity',
               'alignment.length', 'mismatches',
               'gap.opens', 'q.start', 'q.end',
               's.start', 's.end', 'evalue',
               'bit.score', 'qcovs', 'qcovhsp')
  expect_true(all(nms %in% exp_nms))
})