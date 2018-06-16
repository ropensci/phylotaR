#' @name clusters_run
#' @title Run the cluster stage
#' @description Run the third stage of the phylotaR pipeline, cluster.
#' This stage hierarchically traverses the taxonomy identifying all
#' direct and subtree clusters from downloaded sequences. Any
#' taxonomic nodes too small for cluster identification are placed
#' into paraphyletic clusters.
#' @param wd Working directory
#' @family run-public
#' @example examples/clusters_run.R
#' @export
clusters_run <- function(wd) {
  ps <- parameters_load(wd)
  msg <- paste0('Starting stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
  txdct <- obj_load(wd = wd, nm = 'txdct')
  clstrs_calc(ps = ps, txdct = txdct)
  msg <- paste0('Completed stage CLUSTER: [', Sys.time(), ']')
  .stgMsg(ps = ps, msg = msg)
}

#' @name clstrs_calc
#' @title Calculate clusters for all sequences in wd
#' @description Loop through downloaded sequences for each clade and
#' hierarchically find clusters using BLAST.
#' @param txdct Taxonomic dictionary
#' @template ps
#' @family run-private
#' @return NULL
clstrs_calc <- function(txdct, ps) {
  # load sequences
  sq_fls <- list.files(file.path(ps[['wd']], 'cache', 'sqs'))
  sq_fls <- sq_fls[!grepl('paraphyly', sq_fls)]
  fld <- NULL
  for (i in seq_along(sq_fls)) {
    sq_fl <- sq_fls[i]
    # TODO: use the cache tool
    clfl <- file.path(ps[['wd']], 'cache', 'clusters', sq_fl)
    if (file.exists(clfl)) {
      next
    }
    sqs <- readRDS(file = file.path(ps[['wd']], 'cache', 'sqs', sq_fl))
    txid <- as.character(sub('\\.RData', '', sq_fl))
    info(lvl = 1, ps = ps, "Working on [id ", txid, "]")
    clstrs <- clstr_all(txid = txid, sqs = sqs, txdct = txdct, ps = ps)
    if (length(clstrs@ids) > 0) {
      clstrs_save(wd = ps[['wd']], txid = txid, clstrs = clstrs)
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
      sqs <- readRDS(file = file.path(ps[['wd']], 'cache', 'sqs', sq_fl))
      all_sqs <- c(all_sqs, sqs@sqs)
    }
    all_sqs <- seqarc_gen(all_sqs)
    clstrs <- clstr_sqs(txid = '', sqs = all_sqs, ps = ps, lvl = 1,
                        typ = 'paraphyly')
    if (length(clstrs@ids) > 0) {
      clstrs_save(wd = ps[['wd']], txid = 'paraphyly', clstrs = clstrs)
      sqs_save(wd = ps[['wd']], txid = 'paraphyly', sqs = all_sqs)
    }
  }
}
