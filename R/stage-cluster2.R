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
  sqfls <- list.files(sqpth)
  seeds <- NULL
  all_sqs <- all_clstrs <- list()
  for(i in seq_along(clstrfls)) {
    clstrs <- readRDS(file.path(clstrpth, clstrfls[i]))
    sqs <- readRDS(file.path(sqpth, sqfls[i]))
    all_sqs <- c(all_sqs, sqs@sqs)
    all_clstrs <- c(all_clstrs, clstrs)
  }
  all_sqs <- genSqsRcrd(all_sqs)
  if(length(clstrfls) > 1) {
    info(lvl=1, ps=ps, 'Cluster-clustering ...')
    seed_ids <- vapply(clstrs, function(x) x@seed, '')
    seeds <- all_sqs[seed_ids[!duplicated(seed_ids)]]
    blst_rs <- blstSeeds(sqs=seeds, ps=ps)
    info(lvl=1, ps=ps, 'Merging ...')
    jnd_clstrs <- jnClstrs(blst_rs=blst_rs, ps=ps,
                           seed_ids=seed_ids,
                           all_clstrs=all_clstrs)
    mrg_clstrs <- mrgClstrs(jnd_clstrs=jnd_clstrs)
    info(lvl=2, ps=ps, "Generated [", length(mrg_clstrs),
         "] merged clusters")
    all_clstrs <- c(all_clstrs, mrg_clstrs)
  } else {
    info(lvl=1, ps=ps,
         'Done. Only one cluster set -- skipping cluster^2')
  }
  info(lvl=1, ps=ps, 'Renumbering clusters ...')
  all_clstrs <- rnmbrClstrs(clstrs=all_clstrs)
  info(lvl=1, ps=ps, 'Saving ...')
  svObj(wd=ps[['wd']], obj=list('clstrs'=all_clstrs,
                                'sqs'=all_sqs),
        nm='clstrs_sqs')
}
