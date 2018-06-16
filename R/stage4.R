#' @name clusters2_run
#' @title Run the cluster2 stage
#' @description Run the fourth stage of the phylotaR pipeline,
#' cluster2. Identify clusters at higher taxonomic levels by
#' merging sister clusters.
#' @param wd Working directory
#' @export
#' @family run-public
#' @example examples/clusters2_run.R
#' @return NULL
clusters2_run <- function(wd) {
  ps <- parameters_load(wd)
  # stage print
  msg <- paste0('Starting stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  # generate clusters
  clstr2_calc(ps = ps)
  # stage print
  msg <- paste0('Completed stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name clstr2_calc
#' @title Cluster sets of clusters identified in cluster stage
#' @description Loads cluster sets from cache. Extracts seed sequences
#' and runs all-v-all BLAST of seeds to identify sister clusters.
#' Sisters are then merged. An object of all sequences and clusters
#' is then saved in cache.
#' @template ps
#' @return NULL
#' @family run-private
clstr2_calc <- function(ps) {
  # TODO: break up to manageable blasting sizes
  info(lvl = 1, ps = ps, 'Loading clusters ...')
  clstrpth <- file.path(ps[['wd']], 'cache', 'clstrs')
  sqpth <- file.path(ps[['wd']], 'cache', 'sqs')
  clstrfls <- list.files(clstrpth)
  seeds <- NULL
  all_sqs <- all_clstrs <- list()
  for (i in seq_along(clstrfls)) {
    clstrs <- readRDS(file.path(clstrpth, clstrfls[i]))
    sqs <- readRDS(file.path(sqpth, clstrfls[i]))
    all_sqs <- c(all_sqs, sqs@sqs)
    all_clstrs <- c(all_clstrs, clstrs@clstrs)
  }
  all_sqs <- seqarc_gen(all_sqs)
  if (length(clstrfls) > 1) {
    info(lvl = 1, ps = ps, 'Cluster-clustering ...')
    seed_ids <- vapply(all_clstrs, function(x) x@seed, '')
    non_dups <- seed_ids[!duplicated(seed_ids)]
    seeds <- all_sqs[non_dups]
    blast_res <- seeds_blast(sqs = seeds, ps = ps)
    info(lvl = 1, ps = ps, 'Merging ...')
    jnd_clstrs <- clstrs_join(blast_res = blast_res, ps = ps,
                              seed_ids = seed_ids, all_clstrs = all_clstrs)
    txdct <- obj_load(wd = ps[['wd']], nm = 'txdct')
    mrg_clstrs <- clstrs_merge(jnd_clstrs = jnd_clstrs, txdct = txdct)
    all_clstrs <- c(all_clstrs, mrg_clstrs)
  } else {
    info(lvl = 1, ps = ps, 'Done. Only one cluster set',
         ' -- skipping cluster^2')
  }
  info(lvl = 1, ps = ps, 'Dropping all clusters of < 3 sqs ...')
  nsqs <- vapply(all_clstrs, function(x) length(x@sids), 1L)
  all_clstrs <- all_clstrs[nsqs >= 3]
  info(lvl = 1, ps = ps, 'Renumbering clusters ...')
  # returns cluster record box
  all_clstrs <- clstrs_renumber(clstrrecs = all_clstrs)
  info(lvl = 1, ps = ps, 'Saving ...')
  obj_save(wd = ps[['wd']], obj = list('clstrs' = all_clstrs,
                                       'sqs' = all_sqs), nm = 'clstrs_sqs')
}
