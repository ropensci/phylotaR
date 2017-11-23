# LIBS
library(phylotaR)
library(testthat)

# VARS
wd <- getwd()
if(grepl('testthat', wd)) {
  data_dr <- 'data'
} else {
  # for running test at package level
  data_dr <- file.path('tests', 'testthat',
                       'data')
}
cch_dr <- file.path(data_dr, 'cache')

# FUNCTIONS
cleanUp <- function() {
  if(file.exists(cch_dr)) {
    unlink(cch_dr, recursive=TRUE)
  }
}

# RUNNING
context('Testing \'taxise\'')
cleanUp()
test_that('genTDObj() works', {
  tdobj <- genTDObj(data_dr)
  tdobj <- genTDObj(data_dr)
  two <- sum(names(tdobj) %in% c('nds', 'nms'))
  expect_true(two == 2)
  cleanUp()
})
test_that('nDscndnts() works', {
  # stub
})
test_that('getMngblIds() works', {
  # stub
})
cleanUp()