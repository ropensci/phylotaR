# LIBS
library(phylotaR)
library(testthat)

# VARS
data('sqs')  # TODO: make sure sqs match blast results
data('blst_rs')
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
  expect_true(all(nms %in% names(blst_rs)))
})
test_that('fltrBlstRs() works', {
  rnd <- sample(1:nrow(blst_rs), 1)
  blst_rs[rnd, 'qcovs'] <- 0.1
  res <- fltrBlstRs(blst_rs=blst_rs, min_cvrg=0.49)
  expect_true(nrow(res) < nrow(blst_rs))
  # make sure rnd and its sq/qs pair is not present
  frst <- blst_rs[rnd, 'query.id'] == res[['query.id']] &
    blst_rs[rnd, 'subject.id'] == res[['subject.id']]
  expect_false(any(frst))
  scnd <- blst_rs[rnd, 'subject.id'] == res[['query.id']] &
    blst_rs[rnd, 'query.id'] == res[['subject.id']]
  expect_false(any(scnd))
})