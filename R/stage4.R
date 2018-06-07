#' @name clusters2_run
#' @title Cluster2
#' @description Cluster the clusters
#' @param wd Working directory
#' @export
clusters2_run <- function(wd) {
  ps <- ldPrmtrs(wd)
  # stage print
  msg <- paste0('Starting stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  # generate clusters
  cluster2_calc(ps = ps)
  # stage print
  msg <- paste0('Completed stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name clusterclusters
#' @title Cluster sets of clusters identified in cluster stage
#' @description Loads cluster sets from cache. Extracts seed sequences
#' and runs all-v-all BLAST of seeds to identify sister clusters.
#' Sisters are then merged. An object of all sequences and clusters
#' is then saved in cache.
#' @param ps Parameters
cluster2_calc <- function(ps) {
  # TODO: break up to manageable blasting sizes
  info(lvl = 1, ps = ps, 'Loading clusters ...')
  clusterpth <- file.path(ps[['wd']], 'cache', 'clusters')
  sqpth <- file.path(ps[['wd']], 'cache', 'sqs')
  clusterfls <- list.files(clusterpth)
  seeds <- NULL
  all_sqs <- all_clusters <- list()
  for (i in seq_along(clusterfls)) {
    clusters <- readRDS(file.path(clusterpth, clusterfls[i]))
    sqs <- readRDS(file.path(sqpth, clusterfls[i]))
    all_sqs <- c(all_sqs, sqs@sqs)
    all_clusters <- c(all_clusters, clusters@cls)
  }
  all_sqs <- seqarc_gen(all_sqs)
  if (length(clusterfls) > 1) {
    info(lvl = 1, ps = ps, 'Cluster-clustering ...')
    seed_ids <- vapply(all_clusters, function(x) x@seed, '')
    non_dups <- seed_ids[!duplicated(seed_ids)]
    seeds <- all_sqs[non_dups]
    blast_res <- seeds_blast(sqs = seeds, ps = ps)
    info(lvl = 1, ps = ps, 'Merging ...')
    jnd_clusters <- clusters_join(blast_res = blast_res, ps = ps,
                                seed_ids = seed_ids,
                                all_clusters = all_clusters)
    txdct <- ldObj(wd = ps[['wd']], nm = 'txdct')
    mrg_clusters <- clusters_merge(jnd_clusters = jnd_clusters,
                                 txdct = txdct)
    all_clusters <- c(all_clusters, mrg_clusters)
  } else {
    info(lvl = 1, ps = ps,
         'Done. Only one cluster set -- skipping cluster^2')
  }
  info(lvl = 1, ps = ps, 'Dropping all clusters of < 3 sqs ...')
  nsqs <- vapply(all_clusters, function(x) length(x@sids), 1L)
  all_clusters <- all_clusters[nsqs >= 3]
  info(lvl = 1, ps = ps, 'Renumbering clusters ...')
  # returns cluster record box
  all_clusters <- clusters_renumber(cluster_rcrds = all_clusters)
  info(lvl = 1, ps = ps, 'Saving ...')
  svObj(wd = ps[['wd']], obj = list('clusters' = all_clusters,
                                'sqs' = all_sqs),
        nm = 'clusters_sqs')
}
