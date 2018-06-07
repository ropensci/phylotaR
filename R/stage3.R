#' @name clusters_run
#' @title Cluster
#' @description Identify sequence clusters
#' for an initial PhyLoTa cluster table
#' generated
#' @param wd Working directory
#' @export
clusters_run <- function(wd) {
  ps <- parameters_load(wd)
  msg <- paste0('Starting stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  txdct <- obj_load(wd = wd, nm = 'txdct')
  clusters_calc(ps = ps, txdct = txdct)
  msg <- paste0('Completed stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name calcclusters
#' @title Calculate clusters for all sequences in WD
#' @description Loop through downloaded sequences
#' for each clade and hierarchically find clusters
#' using BLAST.
#' @param txdct Taxonomic dictionary
#' @param ps Parameters
clusters_calc <- function(txdct, ps) {
  # load sequences
  sq_fls <- list.files(file.path(ps[['wd']], 'cache', 'sqs'))
  sq_fls <- sq_fls[!grepl('paraphyly', sq_fls)]
  fld <- NULL
  for (i in seq_along(sq_fls)) {
    sq_fl <- sq_fls[i]
    # TODO: use the cache tool
    clfl <- file.path(file.path(file.path(ps[['wd']], 'cache',
                                          'clusters', sq_fl)))
    if (file.exists(clfl)) {
      next
    }
    sqs <- readRDS(file = file.path(file.path(ps[['wd']], 'cache',
                                            'sqs', sq_fl)))
    txid <- as.character(sub('\\.RData', '', sq_fl))
    info(lvl = 1, ps = ps, "Working on [id ", txid, "]")
    clusters <- cluster_all(txid = txid, sqs = sqs, txdct = txdct,
                           ps = ps)
    if (length(clusters@ids) > 0) {
      clusters_save(wd = ps[['wd']], txid = txid, clusters = clusters)
    } else {
      fld <- c(fld, i)
    }
    info(lvl = 1, ps = ps, "[", i, "/", length(sq_fls), "]")
  }
  if (length(fld) > 1) {
    info(lvl = 1, ps = ps, "Paraphyletic retry with unsuccessful clades ...")
    all_sqs <- NULL
    for (i in fld) {
      sq_fl <- sq_fls[i]
      # TODO: use the cache tool
      sqs <- readRDS(file = file.path(file.path(ps[['wd']], 'cache',
                                              'sqs', sq_fl)))
      all_sqs <- c(all_sqs, sqs@sqs)
    }
    all_sqs <- seqarc_gen(all_sqs)
    clusters <- cluster_sqs(txid = '', sqs = all_sqs, ps = ps,
                           lvl = 1, typ = 'paraphyly')
    if (length(clusters@ids) > 0) {
      clusters_save(wd = ps[['wd']], txid = 'paraphyly', clusters = clusters)
      sqs_save(wd = ps[['wd']], txid = 'paraphyly', sqs = all_sqs)
    }
  }
}
