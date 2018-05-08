# LIBS
library(testthat)
library(phylotaR)

# RUNNING
test_that('archive_cluster_generate() works', {
  # TODO: use more realistic cluster records
  clstr_rcrds <- vector(mode = 'list', length = 10)
  names(clstr_rcrds) <- paste0(seq_along(clstr_rcrds))
  archive_cluster <- phylotaR:::archive_cluster_generate(clstr_rcrds =
                                                           clstr_rcrds)
  expect_true(inherits(x = archive_cluster, what = 'ArchiveCluster'))
})
test_that('archive_cluster_join() works', {
  clstr_rcrds <- vector(mode = 'list', length = 10)
  names(clstr_rcrds) <- as.character(seq_along(clstr_rcrds))
  ac_1 <- phylotaR:::archive_cluster_generate(clstr_rcrds =  clstr_rcrds)
  names(clstr_rcrds) <- paste0(seq_along(clstr_rcrds) + length(clstr_rcrds))
  ac_2 <- phylotaR:::archive_cluster_generate(clstr_rcrds =  clstr_rcrds)
  joined <- phylotaR:::archive_cluster_join(ac_1 = ac_1, ac_2 = ac_2)
  expect_true(inherits(x = joined, what = 'ArchiveCluster'))
})