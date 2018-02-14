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
ps <- parameters()
ps[['wd']] <- data_d

# FUNCTIONS
# stubs
mckCmdLn <- function(cmd, args, lgfl=NULL) {
  0
}

# RUNNING
context('Testing \'blast-tools\'')
test_that('mkBlstDB() works', {
  res <- with_mock(
    `phylotaR::cmdLn`=mckCmdLn,
    mkBlstDB(sqs=sqs, dbfl='testdb', ps=ps)
  )
  expect_null(res)
})
test_that('blstN() works', {
  res <- with_mock(
    `phylotaR::cmdLn`=mckCmdLn,
    blstN(dbfl='testdb', outfl='testblstn', ps=ps)
  )
  nms <- colnames(res)
  expect_true(all(nms %in% names(blst_rs)))
})
test_that('fltrBlstRs() works', {
  blst_rs <- blst_rs[-1, ]
  rnd <- sample(2:nrow(blst_rs), 2)
  blst_rs[rnd, 'qcovs'] <- 48
  res <- fltrBlstRs(blst_rs=blst_rs, ps=ps)
  expect_true(nrow(res) < nrow(blst_rs))
  # make sure rnd and its sq/qs pair is not present
  frst <- blst_rs[rnd[1], 'query.id'] == res[['query.id']] &
    blst_rs[rnd[1], 'subject.id'] == res[['subject.id']]
  expect_false(any(frst))
  scnd <- blst_rs[rnd[1], 'subject.id'] == res[['query.id']] &
    blst_rs[rnd[1], 'query.id'] == res[['subject.id']]
  expect_false(any(scnd))
})