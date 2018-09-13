# LIBS
library(phylotaR)
library(testthat)

# DATA
wd <- tempdir()
ps <- parameters(wd = wd)
phylota <- phylotaR:::random_phylota()

# RUNNING
context('Testing \'user-special\'')
phylotaR:::cleanup(wd)
test_that('read_phylota() works', {
  # setup clstrs_sqs
  phylotaR:::cache_setup(ps)
  phylotaR:::obj_save(wd = ps[['wd']], nm = 'clstrs_sqs',
                      obj = list('clstrs' = phylota@clstrs,
                                 'sqs' = phylota@sqs))
  phylotaR:::obj_save(wd = ps[['wd']], obj = phylota@txdct, nm = 'txdct')
  res <- read_phylota(wd = ps[['wd']])
  expect_true(inherits(res, 'Phylota'))
})
phylotaR:::cleanup(wd)
test_that('write_sqs() works', {
  cid <- sample(phylota@cids, 1)
  sids <- phylota@clstrs@clstrs[[cid]]@sids
  n <- ifelse(length(sids) > 50, 50, length(sids))
  sids <- sample(sids, n)
  write_sqs(phylota = phylota, outfile = file.path(wd, 'test.fasta'),
            sid = sids, sq_nm = 'myfavouritegene')
  expect_true(file.exists(file.path(wd, 'test.fasta')))
})
phylotaR:::cleanup(wd)
test_that('plot_phylota_treemap() works', {
  n <- ifelse(length(phylota@cids) > 5, 5, length(phylota@cids))
  cids <- phylota@cids[1:n]
  area <- sample(c('ntx', 'nsq'), 1)
  fill <- sample(c('NULL', 'typ', 'ntx', 'nsq'), 1)
  p1 <- plot_phylota_treemap(phylota = phylota, cids = cids, area = area,
                             fill = fill)
  p1
  expect_true(inherits(p1, 'gg'))
  n <- ifelse(length(phylota@txids) > 5, 5, length(phylota@txids))
  txids <- phylota@txids[1:n]
  area <- sample(c('ncl', 'nsq'), 1)
  fill <- sample(c('NULL', 'nsq', 'ncl'), 1)
  p2 <- plot_phylota_treemap(phylota = phylota, txids = txids, area = area,
                             fill = fill)
  p2
  expect_true(inherits(p2, 'gg'))
})
test_that('plot_phylota_pa() works', {
  n <- ifelse(length(phylota@cids) > 10, 10, length(phylota@cids))
  cids <- phylota@cids[1:n]
  spp <- unique(get_txids(phylota = phylota, txids = phylota@txids,
                          rnk = 'species'))
  spp <- spp[spp != '']
  p <- plot_phylota_pa(phylota = phylota, cids = cids, txids = spp)
  expect_true(inherits(p, 'gg'))
})
test_that('mk_txid_in_sq_mtrx() works', {
  n <- ifelse(length(phylota@sids) > 50, 50, length(phylota@sids))
  sids <- sample(phylota@sids, n)
  n <- ifelse(length(phylota@txids) > 10, 10, length(phylota@txids))
  txids <- sample(phylota@txids, n)
  res <- phylotaR:::mk_txid_in_sq_mtrx(phylota = phylota, txids = txids,
                                       sids = sids)
  expect_true(inherits(res, 'matrix'))
  sid <- sample(colnames(res), 1)
  txid <- sample(rownames(res), 1)
  ads <- phylotaR:::descendants_get(id = txid, txdct = phylota@txdct,
                                    direct = FALSE)
  all_ids <- c(ads, txid)
  sq <- phylota@sqs@sqs[[which(phylota@sqs@ids == sid)]]
  expect_true(res[txid, sid] == (sq@txid %in% all_ids))
})
test_that('is_txid_in_sq() works', {
  txid <- sample(phylota@txids, 1)
  sid <- sample(phylota@sids, 1)
  res <- is_txid_in_sq(phylota = phylota, txid = txid, sid = sid)
  expect_true(inherits(res, 'logical'))
})
test_that('is_txid_in_clstr() works', {
  txid <- sample(phylota@txids, 1)
  cid <- sample(phylota@cids, 1)
  res <- is_txid_in_clstr(phylota = phylota, txid = txid, cid = cid)
  expect_true(inherits(res, 'logical'))
})
test_that('summary_phylota() works', {
  test_phylota <- phylota
  n <- ifelse(length(phylota@cids) > 3, 3, length(phylota@cids))
  cids <- sample(phylota@cids, n)
  test_phylota <- drop_clstrs(phylota = test_phylota, cid = cids)
  res <- phylotaR:::summary_phylota(test_phylota)
  expect_true(inherits(res, 'data.frame'))
})
test_that('update_phylota() works', {
  test_phylota <- phylota
  n <- ifelse(length(phylota@cids) > 10, 10, length(phylota@cids))
  cltrrecs <- sample(test_phylota@clstrs@clstrs, n)
  test_phylota@clstrs <- phylotaR:::clstrarc_gen(clstrrecs = cltrrecs)
  expect_true(length(test_phylota@cids) == length(phylota@cids))
  test_phylota <- phylotaR:::update_phylota(test_phylota)
  expect_false(length(test_phylota@cids) == length(phylota@cids))
  expect_true(length(test_phylota@cids) == 10)
})
