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
  clstrClstrs(ps = ps)
  # stage print
  msg <- paste0('Completed stage CLUSTER^2: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name clstrClstrs
#' @title Cluster sets of clusters identified in cluster stage
#' @description Loads cluster sets from cache. Extracts seed sequences
#' and runs all-v-all BLAST of seeds to identify sister clusters.
#' Sisters are then merged. An object of all sequences and clusters
#' is then saved in cache.
#' @param ps Parameters
clstrClstrs <- function(ps) {
  # TODO: break up to manageable blasting sizes
  info(lvl=1, ps=ps, 'Loading clusters ...')
  clstrpth <- file.path(ps[['wd']], 'cache', 'clstrs')
  sqpth <- file.path(ps[['wd']], 'cache', 'sqs')
  clstrfls <- list.files(clstrpth)
  seeds <- NULL
  all_sqs <- all_clstrs <- list()
  for(i in seq_along(clstrfls)) {
    clstrs <- readRDS(file.path(clstrpth, clstrfls[i]))
    sqs <- readRDS(file.path(sqpth, clstrfls[i]))
    all_sqs <- c(all_sqs, sqs@sqs)
    all_clstrs <- c(all_clstrs, clstrs@cls)
  }
  all_sqs <- genSqRcrdBx(all_sqs)
  if(length(clstrfls) > 1) {
    info(lvl=1, ps=ps, 'Cluster-clustering ...')
    seed_ids <- vapply(all_clstrs, function(x) x@seed, '')
    non_dups <- seed_ids[!duplicated(seed_ids)]
    seeds <- all_sqs[non_dups]
    blst_rs <- blstSeeds(sqs=seeds, ps=ps)
    info(lvl=1, ps=ps, 'Merging ...')
    jnd_clstrs <- jnClstrs(blst_rs=blst_rs, ps=ps,
                           seed_ids=seed_ids,
                           all_clstrs=all_clstrs)
    txdct <- ldObj(wd=ps[['wd']], nm='txdct')
    mrg_clstrs <- mrgClstrs(jnd_clstrs=jnd_clstrs,
                            txdct=txdct)
    all_clstrs <- c(all_clstrs, mrg_clstrs)
  } else {
    info(lvl=1, ps=ps,
         'Done. Only one cluster set -- skipping cluster^2')
  }
  info(lvl=1, ps=ps, 'Dropping all clusters of < 3 sqs ...')
  nsqs <- vapply(all_clstrs, function(x) length(x@sids), 1L)
  all_clstrs <- all_clstrs[nsqs >= 3]
  info(lvl=1, ps=ps, 'Renumbering clusters ...')
  # returns cluster record box
  all_clstrs <- rnmbrClstrs(clstr_rcrds=all_clstrs)
  info(lvl=1, ps=ps, 'Saving ...')
  svObj(wd=ps[['wd']], obj=list('clstrs'=all_clstrs,
                                'sqs'=all_sqs),
        nm='clstrs_sqs')
}
