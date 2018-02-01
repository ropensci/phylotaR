# LIBS
library(phylotaR)
library(testthat)

# DATA
ps <- list('mncvrg'=51, 'v'=FALSE)
# mock clstrs 
clstr <- list("gis"='101', "seed_gi"=NA, "ti_root"='101',
              "ci"=NA, "cl_type"=NA, "n_gi"=1,
              "tis"='101', "n_ti"=1, "MinLength"=0,
              "MaxLength"=0, "n_gen"=1, "n_child"=1,
              "ci_anc"=NA, "unique_id"=NA)
all_clstrs <- lapply(1:100, function(x) {
  clstr[['seed_gi']] <- as.character(x)
  clstr[['ci']] <- as.character(x+100)
  clstr[['unique_id']] <- as.character(x)
  clstr
  })
seed_ids <- as.character(1:100)
sqs <- vector('list', length=200)
names(sqs) <- as.character(1:length(sqs))
sqs <- sample(sqs)
# mock blast results from seed blast
blst_rs <- data.frame(query.id=NA, subject.id=NA,
                      qcovs=runif(min=1, max=100, n=100))
blst_rs[['query.id']] <- as.character(sample(1:100, 50,
                                             replace=FALSE))
blst_rs[['subject.id']] <- as.character(sample(1:100, 50,
                                               replace=FALSE))
pull <- blst_rs[['query.id']] == blst_rs[['subject.id']]
blst_rs[pull, 'qcovs'] <- 100

# RUNNING
context('Testing \'cluster^2 tools\'')
test_that('clstrClstrs() works', {
  
})
# no cache needed
test_that('jnClstrs() works', {
  jnd_clstrs <- jnClstrs(blst_rs=blst_rs, ps=ps,
                         seed_ids=seed_ids,
                         all_clstrs=all_clstrs)
  ids <- unlist(lapply(jnd_clstrs, function(x)
    x[['unique_id']]))
  expect_false(any(duplicated(ids)))
  nids <- vapply(jnd_clstrs, function(x)
    length(x[['unique_id']]), numeric(1))
  ncids <- vapply(jnd_clstrs, function(x)
    length(x[['ci']]), numeric(1))
  expect_true(all(nids == ncids))
  # unique cluster IDs cannot be greater than the
  #  max n IDs in blst_rs
  pssbl_mx <- length(unique(c(blst_rs[['subject.id']],
                              blst_rs[['query.id']])))
  expect_true(length(ids) <= pssbl_mx)
})
test_that('mrgClstrs() works', {
  jnd_clstrs <- jnClstrs(blst_rs=blst_rs, ps=ps,
                         seed_ids=seed_ids,
                         all_clstrs=all_clstrs)
  mrg_clstrs <- mrgClsts(jnd_clstrs=jnd_clstrs)
  expect_true(length(jnd_clstrs) == length(mrg_clstrs))
  rnd <- sample(1:length(jnd_clstrs), 1)
  exmpl <- mrg_clstrs[[rnd]]
  lngths <- vapply(exmpl, length, numeric(1))
  rdclbls <- c('seed_gi', 'ti_root', 'ci', 'cl_type', 'n_gi',
               'n_ti', 'MinLength', 'MaxLength', 'n_gen',
               'n_child', 'ci_anc', 'unique_id')
  tst <- vapply(rdclbls, function(x) lngths[[x]] == 1,
                logical(1))
  expect_true(all(tst))
  expect_true(length(exmpl[['tis']]) == exmpl[['n_ti']])
})
test_that('getSeedSqs() works', {
  sqs <- getSeedSqs(clstrs=all_clstrs, sqs=sqs)
  expect_true(length(sqs) == length(all_clstrs))
})
test_that('rnmbrClstrs() works', {
  clstrs <- all_clstrs
  rndsqs <- sample(letters, size=100, replace=TRUE)
  clstrs[[1]][['gis']] <- rndsqs
  clstrs[[1]][['n_gi']] <- 100
  clstrs <- rnmbrClstrs(clstrs=sample(clstrs))
  expect_true(all(clstrs[[1]][['gis']] == rndsqs))
  cis <- vapply(clstrs, function(x) x[['ci']], numeric(1))
  expect_true(all(cis == 0:99))
})
test_that('blstSeeds() works', {
  res <- with_mock(
    `phylotaR:::mkBlstDB`=function(...) NA,
    `phylotaR:::blstN`=function(...) data.frame(NA),
    blstSeeds(sqs=sqs, ps=ps)
  )
  expect_true('data.frame' %in% is(res))
})